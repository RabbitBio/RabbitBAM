# RabbitBAM: sortedbam

This branch extracts the BAM I/O modules from RabbitBin into a standalone
library and command-line tool (source commit: `4acde498998c3f5cfe8481b2a0505f8ff69ee95d`).

- `src/bamsort`: stable coordinate sorting and parallel BAM/BGZF writing with libdeflate.
- `src/bamindex`: BAI index generation through HTSlib.
- `src/rabbitbam`: the RabbitBAM parallel BGZF reader used in RabbitBin.

## Build

Requires Linux on x86-64, GCC 8.5 or newer (GNU parallel algorithms/OpenMP),
CMake 3.18 or newer, HTSlib 1.15 or newer, libdeflate, zlib, and pthreads.
Install the libraries and their development headers first.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
(cd build && ctest --output-on-failure)
```

For dependencies outside the system paths, add `-DHTSLIB_ROOT=/path/to/htslib`
and `-DLIBDEFLATE_ROOT=/path/to/libdeflate` to the configure command.
These paths should contain `include/` and `lib/` or `lib64/`.

## Usage

```bash
# Sort and also create output.sorted.bam.bai.
./build/rabbitbam sortbam -t 8 --write-index -o output.sorted.bam input.bam

# Index an existing coordinate-sorted BAM.
./build/rabbitbam bai -t 8 -o output.bai output.sorted.bam

./build/rabbitbam --help
./build/rabbitbam sortbam --help
./build/rabbitbam bai --help
```

`sortbam` also accepts standard input (`-`) and writes BAM to standard output
when `-o` is omitted. Use a file output with `--write-index`.
Sorting keeps all records and output buffers in memory; external sorting and
a memory limit are not implemented. The imported writer uses the normal BAM
CIGAR layout (at most 65,535 operations per record).

For C++ reuse, link the CMake target `rabbitbam_io`. It exposes the sorting and
writing API in `rb_bam_sort.h` and the reader in `BamReader.h`. Construct the
sequential reader as `BamReader reader(path, threads, true)` and consume
`getBam1_t()` through EOF before destroying it. The `sortbam` command itself
reads via HTSlib; the RabbitBAM reader is available separately through the library.
This extraction also preserves records buffered with the BAM header, so the
reader can consume files produced by the parallel writer.

The RabbitBin license is retained in `LICENSE`; the notice for HTSlib-derived
code is in `licenses/htslib.txt`.
