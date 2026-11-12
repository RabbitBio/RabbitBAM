// Compare canonical SAM content, using SAM stdin or an independent HTSlib reader.
#ifdef SWBAM_DIGEST_PORTABLE_SHA256
#include "vendor/picosha2/picosha2.h"
#else
#include <openssl/evp.h>
#endif
#ifdef SWBAM_DIGEST_WITH_HTSLIB
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#endif
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Hash = std::array<unsigned char, 32>;

struct Sha256 {
#ifdef SWBAM_DIGEST_PORTABLE_SHA256
    picosha2::hash256_one_by_one context;
    void Reset() { context.init(); }
    void Add(const std::string &text) { context.process(text.begin(), text.end()); }
    Hash Finish() {
        Hash hash{};
        context.finish();
        context.get_hash_bytes(hash.begin(), hash.end());
        return hash;
    }
#else
    EVP_MD_CTX *context = EVP_MD_CTX_new();
    Sha256() { Reset(); }
    ~Sha256() { EVP_MD_CTX_free(context); }
    void Reset() {
        if (!context || EVP_DigestInit_ex(context, EVP_sha256(), nullptr) != 1)
            throw std::runtime_error("SHA-256 initialization failed");
    }
    void Add(const std::string &text) {
        if (EVP_DigestUpdate(context, text.data(), text.size()) != 1)
            throw std::runtime_error("SHA-256 update failed");
    }
    Hash Finish() {
        Hash hash{};
        unsigned length = 0;
        if (EVP_DigestFinal_ex(context, hash.data(), &length) != 1 || length != hash.size())
            throw std::runtime_error("SHA-256 finalization failed");
        return hash;
    }
#endif
    Hash One(const std::string &text) { Reset(); Add(text); return Finish(); }
};

#ifdef SWBAM_DIGEST_WITH_HTSLIB
class BamSamReader {
public:
    BamSamReader() : input_(nullptr), header_(nullptr), record_(nullptr), text_{0, 0, nullptr} {}
    ~BamSamReader() {
        free(text_.s);
        bam_destroy1(record_);
        sam_hdr_destroy(header_);
        if (input_) sam_close(input_);
    }
    void Open(const std::string &path, int threads) {
        input_ = sam_open(path.c_str(), "r");
        if (!input_) throw std::runtime_error("cannot open BAM/SAM input");
        if (hts_set_threads(input_, threads) < 0)
            throw std::runtime_error("cannot initialize HTSlib threads");
        const int eof = hts_check_EOF(input_);
        if (eof == 0 || eof == -1) throw std::runtime_error("missing BGZF EOF or input error");
        header_ = sam_hdr_read(input_);
        record_ = bam_init1();
        if (!header_ || !record_) throw std::runtime_error("cannot read header or allocate record");
        const char *header_text = sam_hdr_str(header_);
        if (!header_text) throw std::runtime_error("cannot format SAM header");
        header_lines_.str(std::string(header_text, sam_hdr_length(header_)));
    }
    void SetRange(uint64_t begin, uint64_t end) {
        if (!input_ || input_->format.format != bam || end <= begin || end >= (uint64_t(1) << 47))
            throw std::runtime_error("invalid compressed BAM range");
        if (begin && bgzf_seek(input_->fp.bgzf, static_cast<int64_t>(begin << 16), SEEK_SET) < 0)
            throw std::runtime_error("cannot seek to BAM range");
        end_ = end << 16;
    }
    bool Next(std::string *line) {
        if (std::getline(header_lines_, *line)) return true;
        if (end_ && static_cast<uint64_t>(bgzf_tell(input_->fp.bgzf)) >= end_) return false;
        const int status = sam_read1(input_, header_, record_);
        if (status == -1) return false;
        if (status < -1) throw std::runtime_error("malformed or truncated BAM record");
        if (sam_format1(header_, record_, &text_) < 0)
            throw std::runtime_error("cannot format SAM record");
        line->assign(text_.s, text_.l);
        return true;
    }
    void Close() {
        samFile *input = input_;
        input_ = nullptr;
        if (sam_close(input) != 0) throw std::runtime_error("BAM close failed");
    }
private:
    samFile *input_;
    sam_hdr_t *header_;
    bam1_t *record_;
    kstring_t text_;
    std::istringstream header_lines_;
    uint64_t end_ = 0;
};
#endif

static void Sum(Hash *sum, const Hash &hash) {
    unsigned carry = 0;
    for (size_t i = 0; i < hash.size(); ++i) {
        carry += (*sum)[i] + hash[i];
        (*sum)[i] = static_cast<unsigned char>(carry);
        carry >>= 8;
    }
}

static std::string Hex(const Hash &hash) {
    std::ostringstream out;
    for (unsigned char byte : hash)
        out << std::hex << std::setfill('0') << std::setw(2) << static_cast<unsigned>(byte);
    return out.str();
}

static std::string NameHex(const std::string &name) {
    static const char digits[] = "0123456789abcdef";
    std::string result;
    for (unsigned char byte : name) {
        result += digits[byte >> 4]; result += digits[byte & 15];
    }
    return result;
}

static std::vector<std::string> Fields(const std::string &line) {
    std::vector<std::string> fields;
    size_t start = 0;
    for (;;) {
        size_t end = line.find('\t', start);
        fields.push_back(line.substr(start, end == std::string::npos ? end : end - start));
        if (end == std::string::npos) return fields;
        start = end + 1;
    }
}

static uint64_t Number(const std::string &text) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("invalid numeric SAM field");
    size_t end = 0;
    uint64_t result = std::stoull(text, &end);
    if (end != text.size()) throw std::runtime_error("invalid numeric SAM field");
    return result;
}

static std::string Join(const std::vector<std::string> &fields) {
    std::string result;
    for (size_t i = 0; i < fields.size(); ++i) {
        if (i) result += '\t';
        result += fields[i];
    }
    return result + '\n';
}

int main(int argc, char **argv) {
    try {
        std::ofstream groups;
        std::string bam_path;
        int threads = 2;
        uint64_t progress_every = 0;
        uint64_t range_begin = 0, range_end = 0;
        for (int i = 1; i < argc; ++i) {
            const std::string option(argv[i]);
            if (i + 1 == argc) throw std::runtime_error("option requires a value: " + option);
            const std::string value(argv[++i]);
            if (option == "--groups") {
                groups.open(value, std::ios::out | std::ios::trunc);
                if (!groups) throw std::runtime_error("cannot open QNAME group file");
            } else if (option == "--bam") {
                bam_path = value;
            } else if (option == "--threads") {
                const uint64_t number = Number(value);
                if (number < 1 || number > 64) throw std::runtime_error("threads must be 1..64");
                threads = static_cast<int>(number);
            } else if (option == "--progress-every") {
                progress_every = Number(value);
            } else if (option == "--range-begin") {
                range_begin = Number(value);
            } else if (option == "--range-end") {
                range_end = Number(value);
            } else {
                throw std::runtime_error("unknown option: " + option);
            }
        }
        if ((range_begin || range_end) && (bam_path.empty() || !range_end))
            throw std::runtime_error("ranges require --bam and --range-end");
#ifdef SWBAM_DIGEST_WITH_HTSLIB
        BamSamReader reader;
        if (!bam_path.empty()) reader.Open(bam_path, threads);
        if (range_end) reader.SetRange(range_begin, range_end);
#else
        if (!bam_path.empty()) throw std::runtime_error("this build has no HTSlib input support");
        if (range_begin || range_end) throw std::runtime_error("ranges require HTSlib input");
        (void)threads;
#endif
        Sha256 header, ordered, record;
        Hash multiset{}, neutral{}, duplicate_set{};
        uint64_t count = 0, duplicates = 0, group_count = 0, coordinate_errors = 0, read_order_errors = 0;
        std::map<std::string, int> references;
        std::pair<int, uint64_t> previous(-1, 0);
        std::string line, last_qname;
        unsigned previous_read_order = 0;
        std::pair<int, uint64_t> first_coordinate(-1, 0);
        std::string first_qname;
        unsigned first_read_order = 0;
        const auto start = std::chrono::steady_clock::now();
        auto next = [&]() {
#ifdef SWBAM_DIGEST_WITH_HTSLIB
            if (!bam_path.empty()) return reader.Next(&line);
#endif
            return static_cast<bool>(std::getline(std::cin, line));
        };
        while (next()) {
            if (!count && !line.empty() && line[0] == '@') {
                if (line.compare(0, 4, "@PG\t") != 0) header.Add(line + '\n');
                if (line.compare(0, 4, "@SQ\t") == 0) {
                    for (const auto &field : Fields(line)) {
                        if (field.compare(0, 3, "SN:") == 0) {
                            const int tid = static_cast<int>(references.size());
                            if (!references.emplace(field.substr(3), tid).second)
                                throw std::runtime_error("duplicate reference name");
                        }
                    }
                }
                continue;
            }
            auto fields = Fields(line);
            if (fields.size() < 11) throw std::runtime_error("SAM record has fewer than 11 fields");
            const uint64_t flag = Number(fields[1]);
            if (flag > 65535) throw std::runtime_error("SAM FLAG exceeds 16 bits");
            int tid = std::numeric_limits<int>::max();
            if (fields[2] != "*") {
                auto ref = references.find(fields[2]);
                if (ref == references.end()) throw std::runtime_error("unknown reference name");
                tid = ref->second;
            }
            const auto coordinate = std::make_pair(tid, Number(fields[3]));
            if (count && coordinate < previous) ++coordinate_errors;
            previous = coordinate;
            const unsigned read_order = (flag >> 6) & 3;
            if (!count) {
                first_coordinate = coordinate;
                first_qname = fields[0];
                first_read_order = read_order;
            }
            if (!count || fields[0] != last_qname) {
                ++group_count;
                last_qname = fields[0];
                previous_read_order = 0;
                if (groups.is_open()) groups << last_qname << '\n';
            } else if (read_order < previous_read_order) ++read_order_errors;
            previous_read_order = read_order;
            std::sort(fields.begin() + 11, fields.end());
            const std::string canonical = Join(fields);
            const Hash full = record.One(canonical);
            Sum(&multiset, full);
            ordered.Add(canonical);
            if (flag & 1024) {
                ++duplicates;
                fields[1] = std::to_string(flag & ~uint64_t(1024));
                const Hash identity = record.One(Join(fields));
                Sum(&neutral, identity);
                Sum(&duplicate_set, identity);
            } else Sum(&neutral, full);
            ++count;
            if (progress_every && count % progress_every == 0)
                std::cerr << "[digest] records=" << count << " seconds="
                          << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                          << '\n';
        }
        if (bam_path.empty() && std::cin.bad()) throw std::runtime_error("SAM input read failed");
#ifdef SWBAM_DIGEST_WITH_HTSLIB
        if (!bam_path.empty()) reader.Close();
#endif
        if (groups.is_open()) {
            groups.close();
            if (!groups) throw std::runtime_error("QNAME group write failed");
        }
        std::cout << "{\n  \"records\": " << count << ",\n  \"duplicates\": " << duplicates
                  << ",\n  \"qname_groups\": " << group_count
                  << ",\n  \"coordinate_errors\": " << coordinate_errors
                  << ",\n  \"read_order_errors\": " << read_order_errors
                  << ",\n  \"header_without_pg_sha256\": \"" << Hex(header.Finish())
                  << "\",\n  \"ordered_sha256\": \"" << Hex(ordered.Finish())
                  << "\",\n  \"record_multiset_sha256_sum\": \"" << Hex(multiset)
                  << "\",\n  \"without_dup_multiset_sha256_sum\": \"" << Hex(neutral)
                  << "\",\n  \"duplicate_members_sha256_sum\": \"" << Hex(duplicate_set)
                  << "\",\n  \"scope\": \"11 SAM fields + all AUX; AUX order normalized; sums modulo 2^256, little endian\"";
        if (range_end) {
            std::cout << ",\n  \"first_coordinate\": [" << first_coordinate.first << ',' << first_coordinate.second << ']'
                      << ",\n  \"last_coordinate\": [" << previous.first << ',' << previous.second << ']'
                      << ",\n  \"first_qname_hex\": \"" << NameHex(first_qname) << '"'
                      << ",\n  \"last_qname_hex\": \"" << NameHex(last_qname) << '"'
                      << ",\n  \"first_read_order\": " << first_read_order
                      << ",\n  \"last_read_order\": " << previous_read_order;
        }
        std::cout << "\n}\n";
        if (!std::cout) throw std::runtime_error("JSON output write failed");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "sam_correctness_digest: " << error.what() << '\n';
        return 1;
    }
}
