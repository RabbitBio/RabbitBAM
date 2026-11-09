#!/usr/bin/env python3
"""Exercise the real fast-reference lookup with an injected dictionary miss."""

from pathlib import Path
import subprocess
import tempfile
import unittest


SOURCE = Path(__file__).resolve().parents[1] / "swbam-sdk/slave/sam_parse.c"


def function(text, name):
    start = text.index("static inline int " + name + "(")
    brace = text.index("{", start)
    depth = 1
    end = brace + 1
    while depth:
        depth += (text[end] == "{") - (text[end] == "}")
        end += 1
    return text[start:end]


HARNESS = r'''
#include <assert.h>
#include <stdio.h>
#include <string.h>
typedef struct { int n_targets; char **target_name; } sam_hdr_t;
typedef struct {
    const sam_hdr_t *hdr[2];
    int name_len[2], tid[2];
    char name[2][64];
} MpiSamParseFastCache;
static int result = -1, calls = 0;
static int bam_name2id(sam_hdr_t *h, const char *name) {
    (void)h; (void)name; ++calls; return result;
}
'''

CASES = r'''
int main(void) {
    char *names[] = {"chr1", "chr2", "1234567890123456789012345678901234567890123456789012345678901234567890"};
    sam_hdr_t header = {3, names};
    MpiSamParseFastCache cache = {0};
    int tid = -9;
    char field[] = "chr2\tNEXT";
    assert(mpi_sam_name2id_fast(&header, field, field + 4, &tid, &cache, 0));
    assert(tid == 1 && calls == 1 && field[4] == '\t');
    assert(mpi_sam_name2id_fast(&header, field, field + 4, &tid, &cache, 0));
    assert(tid == 1 && calls == 1); /* Cached fallback result. */
    char mate[] = "chr1\t";
    assert(mpi_sam_name2id_fast(&header, mate, mate + 4, &tid, &cache, 1));
    assert(tid == 0 && calls == 2);
    char star[] = "*\t";
    assert(mpi_sam_name2id_fast(&header, star, star + 1, &tid, &cache, 0));
    assert(tid == -1 && calls == 2);
    result = 2;
    assert(mpi_sam_name2id_fast(&header, field, field + 4, &tid, NULL, 0));
    assert(tid == 2); /* Preserve a successful fast lookup. */
    result = -2;
    assert(!mpi_sam_name2id_fast(&header, field, field + 4, &tid, NULL, 0));
    result = -1;
    char unknown[] = "missing\t";
    assert(mpi_sam_name2id_fast(&header, unknown, unknown + 7, &tid, NULL, 0));
    assert(tid == -1); /* Preserve existing unknown-reference semantics. */
    char long_name[80];
    strcpy(long_name, names[2]);
    size_t length = strlen(long_name);
    long_name[length] = '\t';
    assert(mpi_sam_name2id_fast(&header, long_name, long_name + length, &tid, &cache, 0));
    assert(tid == 2 && long_name[length] == '\t');
    sam_hdr_t empty = {0, NULL};
    assert(!mpi_sam_name2id_fast(&empty, mate, mate + 4, &tid, &cache, 0));
    assert(mpi_sam_name2id_fast(&empty, star, star + 1, &tid, &cache, 0));
    sam_hdr_t no_dictionary = {1, NULL};
    assert(mpi_sam_name2id_fast(&no_dictionary, mate, mate + 4, &tid, NULL, 0));
    assert(tid == -1);
    puts("PASS: forced miss, cache reuse, mate slot, *, successful lookup, error, unknown and long names");
    return 0;
}
'''


class ReferenceLookupTest(unittest.TestCase):
    def test_actual_lookup(self):
        text = SOURCE.read_text()
        source = HARNESS + "\n" + function(text, "mpi_sam_name_equal_cached") + "\n" + \
            function(text, "mpi_sam_name2id_fast") + "\n" + CASES
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "lookup.c").write_text(source)
            subprocess.run(["cc", "-std=c11", "-Wall", "-Wextra", "-Werror",
                            "-fsanitize=address,undefined", "-fno-omit-frame-pointer",
                            str(path / "lookup.c"), "-o", str(path / "lookup")], check=True)
            subprocess.run([str(path / "lookup")], check=True)


if __name__ == "__main__":
    unittest.main()
