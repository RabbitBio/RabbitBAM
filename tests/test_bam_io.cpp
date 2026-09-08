#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "BamReader.h"
#include "rb_bam_sort.h"

static void require(bool ok, const char *message) {
    if (!ok) {
        std::fprintf(stderr, "%s\n", message);
        // Stop immediately on failure: the imported reader requires draining
        // its queues before destruction, which a failed comparison cannot do.
        std::exit(1);
    }
}

static std::vector<std::string> reference_records(const char *path, int count) {
    samFile *input = sam_open(path, "rb");
    require(input != nullptr, "Cannot open the written BAM");
    sam_hdr_t *header = sam_hdr_read(input);
    require(header != nullptr, "Cannot read the written header");
    if (count > 0) {
        require(input->fp.bgzf->block_offset < input->fp.bgzf->block_length,
                "Regression fixture must pack the header and records together");
    }
    std::vector<std::string> records;
    bam1_t *record = bam_init1();
    kstring_t line = {0, 0, nullptr};
    int result;
    while ((result = sam_read1(input, header, record)) >= 0) {
        line.l = 0;
        require(sam_format1(header, record, &line) >= 0, "Cannot format a record");
        records.emplace_back(line.s, line.l);
    }
    require(result == -1 && records.size() == static_cast<size_t>(count),
            "HTSlib did not recover every written record");
    std::free(line.s);
    bam_destroy1(record);
    sam_hdr_destroy(header);
    require(sam_close(input) == 0, "Cannot close the reference reader");
    return records;
}

static void check_reader(BamReader &reader,
                         const std::vector<std::string> &expected) {
    bam1_t *record = bam_init1();
    kstring_t line = {0, 0, nullptr};
    size_t count = 0;
    while (reader.getBam1_t(record)) {
        require(count < expected.size(), "RabbitBAM returned extra records");
        line.l = 0;
        require(sam_format1(reader.getHeader(), record, &line) >= 0,
                "Cannot format a RabbitBAM record");
        require(std::string(line.s, line.l) == expected[count],
                "RabbitBAM changed record content or order");
        ++count;
    }
    require(count == expected.size(), "RabbitBAM lost buffered records");
    std::free(line.s);
    bam_destroy1(record);
}

static void roundtrip(const char *path, int count) {
    const char *text = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n";
    sam_hdr_t *header = sam_hdr_parse(std::strlen(text), text);
    require(header != nullptr, "Cannot create the test header");
    std::vector<bam1_t *> records;
    for (int i = 0; i < count; ++i) {
        const std::string text = "read" + std::to_string(i) + "\t0\tchr1\t" +
            std::to_string(count - i) + "\t60\t75M\t*\t0\t0\t" +
            std::string(75, 'A') + "\t" + std::string(75, 'I');
        kstring_t line = {0, 0, nullptr};
        require(kputs(text.c_str(), &line) >= 0, "Cannot allocate the test record");
        bam1_t *record = bam_init1();
        require(sam_parse1(&line, header, record) >= 0, "Cannot parse the test record");
        records.push_back(record);
        std::free(line.s);
    }
    omp_set_num_threads(2);
    rb_stable_coord_sort(records);
    require(rb_write_sorted_bam(header, records, path, 2, 6, true) == 0,
            "Cannot write and index the sorted BAM");
    rb_free_records(records, 2);
    sam_hdr_destroy(header);

    const auto expected = reference_records(path, count);
    {
        BamReader reader(path, 2, true);
        check_reader(reader, expected);
    }
    {
        BamReader reader(path, 8, 8, 8, 2, true);
        check_reader(reader, expected);
    }
    std::remove(path);
    std::remove((std::string(path) + ".bai").c_str());
}

int main(int argc, char *argv[]) {
    require(argc == 2, "Expected an output BAM path");
    // Cover EOF, records wholly buffered after the header, and records spanning
    // several BGZF blocks. Exercise both reader constructors for every case.
    for (int count : {0, 1, 6000}) roundtrip(argv[1], count);
    return 0;
}
