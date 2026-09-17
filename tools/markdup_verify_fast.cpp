#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

const uint64_t kFnvOffset = 1469598103934665603ull;
const uint64_t kFnvPrime = 1099511628211ull;

struct Field {
    const char *data;
    size_t size;
};

struct UnorderedFingerprint {
    uint64_t count;
    uint64_t sum1;
    uint64_t sum2;
    uint64_t xors;

    UnorderedFingerprint()
        : count(0), sum1(0), sum2(0), xors(0) {}
};

struct OrderedFingerprint {
    uint64_t count;
    uint64_t h1;
    uint64_t h2;

    OrderedFingerprint()
        : count(0), h1(0x123456789abcdef0ull),
          h2(0xfedcba9876543210ull) {}
};

struct Summary {
    uint64_t records;
    uint64_t primary;
    uint64_t secondary;
    uint64_t supplementary;
    uint64_t duplicates;
    uint64_t primary_duplicates;
    uint64_t mapped;
    uint64_t primary_mapped;
    uint64_t paired;
    uint64_t read1;
    uint64_t read2;
    uint64_t properly_paired;
    uint64_t with_mate_mapped;
    uint64_t singletons;
    uint64_t mate_diff_chr;
    uint64_t mate_diff_chr_mapq5;
    uint64_t qc_failed;
    OrderedFingerprint header_without_pg;
    OrderedFingerprint ordered_body;
    UnorderedFingerprint record_multiset;
    UnorderedFingerprint mandatory_fields;
    UnorderedFingerprint duplicate_key;

    Summary()
        : records(0), primary(0), secondary(0),
          supplementary(0), duplicates(0),
          primary_duplicates(0), mapped(0),
          primary_mapped(0), paired(0), read1(0), read2(0),
          properly_paired(0), with_mate_mapped(0),
          singletons(0), mate_diff_chr(0),
          mate_diff_chr_mapq5(0), qc_failed(0) {}
};

static uint64_t SplitMix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

static uint64_t HashBytes(const char *data, size_t size,
                          uint64_t seed = kFnvOffset) {
    uint64_t h = seed;
    for (size_t i = 0; i < size; ++i) {
        h ^= (unsigned char)data[i];
        h *= kFnvPrime;
    }
    return h;
}

struct HashBuilder {
    uint64_t h;

    HashBuilder() : h(kFnvOffset) {}

    void AddByte(unsigned char c) {
        h ^= c;
        h *= kFnvPrime;
    }

    void AddField(const Field &field) {
        h = HashBytes(field.data, field.size, h);
        AddByte('\t');
    }

    void AddUInt(uint64_t value) {
        char buf[32];
        int n = 0;
        do {
            buf[n++] = (char)('0' + value % 10);
            value /= 10;
        } while (value);
        while (n-- > 0) AddByte((unsigned char)buf[n]);
        AddByte('\t');
    }
};

static void AddUnordered(UnorderedFingerprint *fp, uint64_t hash) {
    uint64_t a = SplitMix64(hash);
    uint64_t b = SplitMix64(hash ^ 0xd6e8feb86659fd93ull);
    fp->count++;
    fp->sum1 += a;
    fp->sum2 += b;
    fp->xors ^= SplitMix64(a ^ b);
}

static void AddOrdered(OrderedFingerprint *fp, uint64_t hash) {
    uint64_t a = SplitMix64(hash + fp->count);
    fp->count++;
    fp->h1 = SplitMix64(fp->h1 ^ a ^
                        (0x9e3779b97f4a7c15ull + fp->count));
    fp->h2 = fp->h2 * 1099511628211ull ^
             SplitMix64(hash ^ (fp->count << 1));
}

static std::string Hex64(uint64_t value) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16)
        << value;
    return out.str();
}

static std::string FingerprintString(
        const UnorderedFingerprint &fp) {
    std::ostringstream out;
    out << fp.count << ":" << Hex64(fp.sum1) << ":"
        << Hex64(fp.sum2) << ":" << Hex64(fp.xors);
    return out.str();
}

static std::string FingerprintString(
        const OrderedFingerprint &fp) {
    std::ostringstream out;
    out << fp.count << ":" << Hex64(fp.h1) << ":"
        << Hex64(fp.h2);
    return out.str();
}

static bool FieldEquals(const Field &field, const char *literal) {
    size_t n = std::strlen(literal);
    return field.size == n &&
           std::memcmp(field.data, literal, n) == 0;
}

static uint64_t ParseUInt(const Field &field) {
    uint64_t value = 0;
    for (size_t i = 0; i < field.size; ++i) {
        unsigned char c = (unsigned char)field.data[i];
        if (c < '0' || c > '9') break;
        value = value * 10u + (uint64_t)(c - '0');
    }
    return value;
}

static int SplitSamFields(const std::string &line,
                          Field *fields, int max_fields) {
    const char *base = line.data();
    size_t start = 0;
    int count = 0;
    for (size_t i = 0; i <= line.size(); ++i) {
        if (i == line.size() || line[i] == '\t') {
            if (count < max_fields) {
                fields[count].data = base + start;
                fields[count].size = i - start;
            }
            count++;
            start = i + 1;
        }
    }
    return count;
}

static uint64_t DuplicateKeyHash(const Field *fields,
                                 uint64_t flag) {
    uint64_t dup = (flag & 1024u) ? 1u : 0u;
    uint64_t flag_no_dup = flag & ~1024ull;
    HashBuilder hb;
    hb.AddField(fields[0]);
    hb.AddUInt(flag_no_dup);
    hb.AddField(fields[2]);
    hb.AddField(fields[3]);
    hb.AddField(fields[4]);
    hb.AddField(fields[5]);
    hb.AddField(fields[6]);
    hb.AddField(fields[7]);
    hb.AddField(fields[8]);
    hb.AddField(fields[9]);
    hb.AddUInt(dup);
    return hb.h;
}

static uint64_t MandatoryFieldsHash(const Field *fields) {
    HashBuilder hb;
    for (int i = 0; i < 11; ++i) {
        hb.AddField(fields[i]);
    }
    return hb.h;
}

static bool ScanSam(const std::string &path, Summary *summary) {
    std::ifstream in(path.c_str());
    if (!in) {
        std::cerr << "ERROR: cannot read SAM: " << path << "\n";
        return false;
    }

    std::string line;
    uint64_t line_no = 0;
    while (std::getline(in, line)) {
        line_no++;
        if (!line.empty() && line[line.size() - 1] == '\r') {
            line.erase(line.size() - 1);
        }
        if (line.empty()) continue;
        uint64_t line_hash = HashBytes(line.data(), line.size());
        if (line[0] == '@') {
            if (line.compare(0, 3, "@PG") != 0) {
                AddOrdered(&summary->header_without_pg,
                           line_hash);
            }
            continue;
        }

        Field fields[11];
        int n_fields = SplitSamFields(line, fields, 11);
        if (n_fields < 10) {
            std::cerr << "ERROR: malformed SAM body line "
                      << line_no << " in " << path << "\n";
            return false;
        }

        uint64_t flag = ParseUInt(fields[1]);
        bool is_paired = (flag & 1u) != 0;
        bool is_proper = (flag & 2u) != 0;
        bool is_unmapped = (flag & 4u) != 0;
        bool mate_unmapped = (flag & 8u) != 0;
        bool is_read1 = (flag & 64u) != 0;
        bool is_read2 = (flag & 128u) != 0;
        bool is_secondary = (flag & 256u) != 0;
        bool is_qcfail = (flag & 512u) != 0;
        bool is_dup = (flag & 1024u) != 0;
        bool is_supp = (flag & 2048u) != 0;

        summary->records++;
        if (is_secondary) summary->secondary++;
        if (is_supp) summary->supplementary++;
        if (!is_secondary && !is_supp) summary->primary++;
        if (is_dup) summary->duplicates++;
        if (is_dup && !is_secondary && !is_supp) {
            summary->primary_duplicates++;
        }
        if (!is_unmapped) summary->mapped++;
        if (!is_unmapped && !is_secondary && !is_supp) {
            summary->primary_mapped++;
        }
        if (is_paired) summary->paired++;
        if (is_read1) summary->read1++;
        if (is_read2) summary->read2++;
        if (is_proper) summary->properly_paired++;
        if (is_paired && !is_unmapped && !mate_unmapped) {
            summary->with_mate_mapped++;
        }
        if (is_paired && !is_unmapped && mate_unmapped) {
            summary->singletons++;
        }
        if (is_paired && !is_unmapped && !mate_unmapped &&
            !FieldEquals(fields[6], "=") &&
            !FieldEquals(fields[6], "*") &&
            !FieldEquals(fields[2], "*") &&
            !(fields[6].size == fields[2].size &&
              std::memcmp(fields[6].data, fields[2].data,
                          fields[2].size) == 0)) {
            summary->mate_diff_chr++;
            if (ParseUInt(fields[4]) >= 5u) {
                summary->mate_diff_chr_mapq5++;
            }
        }
        if (is_qcfail) summary->qc_failed++;

        AddOrdered(&summary->ordered_body, line_hash);
        AddUnordered(&summary->record_multiset, line_hash);
        AddUnordered(&summary->mandatory_fields,
                     MandatoryFieldsHash(fields));
        AddUnordered(&summary->duplicate_key,
                     DuplicateKeyHash(fields, flag));
    }

    return true;
}

static void PrintSummary(const char *label,
                         const Summary &summary) {
    std::cout << "[" << label << "]\n";
    std::cout << "records " << summary.records << "\n";
    std::cout << "primary " << summary.primary << "\n";
    std::cout << "secondary " << summary.secondary << "\n";
    std::cout << "supplementary " << summary.supplementary << "\n";
    std::cout << "duplicates " << summary.duplicates << "\n";
    std::cout << "primary_duplicates "
              << summary.primary_duplicates << "\n";
    std::cout << "mapped " << summary.mapped << "\n";
    std::cout << "primary_mapped "
              << summary.primary_mapped << "\n";
    std::cout << "paired " << summary.paired << "\n";
    std::cout << "read1 " << summary.read1 << "\n";
    std::cout << "read2 " << summary.read2 << "\n";
    std::cout << "properly_paired "
              << summary.properly_paired << "\n";
    std::cout << "with_mate_mapped "
              << summary.with_mate_mapped << "\n";
    std::cout << "singletons " << summary.singletons << "\n";
    std::cout << "mate_diff_chr "
              << summary.mate_diff_chr << "\n";
    std::cout << "mate_diff_chr_mapq5 "
              << summary.mate_diff_chr_mapq5 << "\n";
    std::cout << "qc_failed " << summary.qc_failed << "\n";
    std::cout << "header_without_pg_fp "
              << FingerprintString(summary.header_without_pg)
              << "\n";
    std::cout << "ordered_body_fp "
              << FingerprintString(summary.ordered_body) << "\n";
    std::cout << "record_multiset_fp "
              << FingerprintString(summary.record_multiset)
              << "\n";
    std::cout << "mandatory_fields_fp "
              << FingerprintString(summary.mandatory_fields)
              << "\n";
    std::cout << "duplicate_key_fp "
              << FingerprintString(summary.duplicate_key) << "\n";
}

template <typename T>
static bool CompareValue(const char *key, const T &base,
                         const T &test, bool required,
                         bool *failed) {
    bool ok = base == test;
    if (!ok && required) *failed = true;
    std::cout << "compare " << std::left << std::setw(24)
              << key << (ok ? "OK" : "DIFF")
              << " base=" << base << " test=" << test << "\n";
    return ok;
}

static bool CompareFp(const char *key, const std::string &base,
                      const std::string &test, bool required,
                      bool *failed) {
    bool ok = base == test;
    if (!ok && required) *failed = true;
    std::cout << "compare " << std::left << std::setw(24)
              << key << (ok ? "OK" : "DIFF")
              << " base=" << base << " test=" << test << "\n";
    return ok;
}

} // namespace

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0]
                  << " BASE.sam TEST.sam\n";
        std::cerr << "Set STRICT_ORDER=1 to require body order "
                  << "fingerprint equality.\n";
        return 2;
    }

    Summary base;
    Summary test;
    if (!ScanSam(argv[1], &base) ||
        !ScanSam(argv[2], &test)) {
        return 1;
    }

    PrintSummary("base", base);
    PrintSummary("test", test);

    bool failed = false;
    bool strict_order = false;
    bool strict_tags = false;
    const char *strict_env = std::getenv("STRICT_ORDER");
    if (strict_env && std::strcmp(strict_env, "0") != 0) {
        strict_order = true;
    }
    const char *strict_tags_env = std::getenv("STRICT_TAGS");
    if (strict_tags_env && std::strcmp(strict_tags_env, "0") != 0) {
        strict_tags = true;
    }

    std::cout << "[comparison]\n";
    CompareValue("records", base.records, test.records, true,
                 &failed);
    CompareValue("primary", base.primary, test.primary, true,
                 &failed);
    CompareValue("secondary", base.secondary, test.secondary,
                 true, &failed);
    CompareValue("supplementary", base.supplementary,
                 test.supplementary, true, &failed);
    CompareValue("duplicates", base.duplicates, test.duplicates,
                 true, &failed);
    CompareValue("primary_duplicates",
                 base.primary_duplicates,
                 test.primary_duplicates, true, &failed);
    CompareValue("mapped", base.mapped, test.mapped, true,
                 &failed);
    CompareValue("primary_mapped", base.primary_mapped,
                 test.primary_mapped, true, &failed);
    CompareValue("paired", base.paired, test.paired, true,
                 &failed);
    CompareValue("read1", base.read1, test.read1, true,
                 &failed);
    CompareValue("read2", base.read2, test.read2, true,
                 &failed);
    CompareValue("properly_paired", base.properly_paired,
                 test.properly_paired, true, &failed);
    CompareValue("with_mate_mapped", base.with_mate_mapped,
                 test.with_mate_mapped, true, &failed);
    CompareValue("singletons", base.singletons,
                 test.singletons, true, &failed);
    CompareValue("mate_diff_chr", base.mate_diff_chr,
                 test.mate_diff_chr, true, &failed);
    CompareValue("mate_diff_chr_mapq5",
                 base.mate_diff_chr_mapq5,
                 test.mate_diff_chr_mapq5, true, &failed);
    CompareValue("qc_failed", base.qc_failed, test.qc_failed,
                 true, &failed);
    CompareFp("header_without_pg_fp",
              FingerprintString(base.header_without_pg),
              FingerprintString(test.header_without_pg), true,
              &failed);
    CompareFp("record_multiset_fp",
              FingerprintString(base.record_multiset),
              FingerprintString(test.record_multiset),
              strict_tags,
              &failed);
    CompareFp("mandatory_fields_fp",
              FingerprintString(base.mandatory_fields),
              FingerprintString(test.mandatory_fields), true,
              &failed);
    CompareFp("duplicate_key_fp",
              FingerprintString(base.duplicate_key),
              FingerprintString(test.duplicate_key), true,
              &failed);
    CompareFp("ordered_body_fp",
              FingerprintString(base.ordered_body),
              FingerprintString(test.ordered_body),
              strict_order, &failed);

    if (failed) {
        std::cerr << "markdup_verify_fast: FAILED\n";
        return 1;
    }
    std::cerr << "markdup_verify_fast: PASS\n";
    return 0;
}
