#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?Usage: build_rabbitbam_library_benchmark.sh ROOT}
HTS=${HTS:-/home/user_home/wzs/software/src/samtools-1.20/htslib-1.20}
DEFLATE=${DEFLATE:-/home/user_home/wzs/software/libdeflate-1.20}
COMMIT=2c2bf3572fa76ae775787d324f6660f06d6f6519
[[ $(git -C "$ROOT/upstream" rev-parse HEAD) == "$COMMIT" ]]
git -C "$ROOT/upstream" diff --exit-code
mkdir -p "$ROOT/bin" "$ROOT/logs" "$ROOT/out" "$ROOT/validation"
cd "$ROOT/upstream"
sources=(BamTools.cpp Buffer.cpp BamCompleteBlock.cpp BamCompress.cpp BamBlock.cpp
         BamRead.cpp BamReader.cpp BamWriteCompress.cpp BamWriter.cpp)
flags=(-O3 -g -std=c++11 -ffast-math -fopenmp -pthread -I. -I./htslib -I"$HTS"
       -I"$DEFLATE/include")
libs=("$HTS/libhts.a" -L"$DEFLATE/lib64" -Wl,-rpath,"$DEFLATE/lib64"
      -ldeflate -lz -lbz2 -llzma -lcurl -lcrypto -lm)
command=(g++ "${flags[@]}" "$ROOT/scripts/rabbitbam_library_benchmark.cpp"
         "${sources[@]}" "${libs[@]}" -o "$ROOT/bin/rabbitbam-library")
printf '%q ' "${command[@]}" > "$ROOT/build-command.txt"
printf '\n' >> "$ROOT/build-command.txt"
"${command[@]}"
gcc -O3 -I"$HTS" "$ROOT/scripts/sam_stream.c" "$HTS/libhts.a" \
    -L"$DEFLATE/lib64" -Wl,-rpath,"$DEFLATE/lib64" \
    -ldeflate -lz -lbz2 -llzma -lcurl -lcrypto -lpthread -lm \
    -o "$ROOT/bin/sam-stream"
git rev-parse HEAD > "$ROOT/upstream-commit.txt"
git status --porcelain > "$ROOT/upstream-status.txt"
sha256sum "$ROOT/bin/rabbitbam-library" "$ROOT/bin/sam-stream" \
    "$HTS/libhts.a" "$DEFLATE/lib64/libdeflate.so.0" \
    "$ROOT"/scripts/*.cpp "$ROOT"/scripts/*.py "$ROOT"/scripts/*.sh \
    "$ROOT"/scripts/*.c > "$ROOT/frozen.sha256"
ldd "$ROOT/bin/rabbitbam-library" > "$ROOT/linkage.txt"
g++ --version > "$ROOT/compiler.txt"
