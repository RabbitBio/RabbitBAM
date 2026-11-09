// Input is a complete SAM stream emitted by an independent SAMtools reader.
#include <openssl/evp.h>
#include <algorithm>
#include <array>
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
    Hash One(const std::string &text) { Reset(); Add(text); return Finish(); }
};

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
        if (argc == 3 && std::string(argv[1]) == "--groups") {
            groups.open(argv[2], std::ios::out | std::ios::trunc);
            if (!groups) throw std::runtime_error("cannot open QNAME group file");
        } else if (argc != 1) {
            throw std::runtime_error("usage: samtools view -h input.bam | sam_correctness_digest [--groups FILE]");
        }
        Sha256 header, ordered, record;
        Hash multiset{}, neutral{}, duplicate_set{};
        uint64_t count = 0, duplicates = 0, group_count = 0, coordinate_errors = 0, read_order_errors = 0;
        std::map<std::string, int> references;
        std::pair<int, uint64_t> previous(-1, 0);
        std::string line, last_qname;
        unsigned previous_read_order = 0;
        while (std::getline(std::cin, line)) {
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
        }
        if (std::cin.bad()) throw std::runtime_error("SAM input read failed");
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
                  << "\",\n  \"scope\": \"11 SAM fields + all AUX; AUX order normalized; sums modulo 2^256, little endian\"\n}\n";
        if (!std::cout) throw std::runtime_error("JSON output write failed");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "sam_correctness_digest: " << error.what() << '\n';
        return 1;
    }
}
