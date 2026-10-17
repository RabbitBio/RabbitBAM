#!/usr/bin/env bash
set -Eeuo pipefail
export TZ=Asia/Shanghai LC_ALL=C

ROOT=${D3_ROOT:-/home/user_home/wzs/wzs_data/D3-mouse-SRR7768905}
TOOLS=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BWA=/home/user_home/wzs/software/bwa-0.7.18-clean/bwa
SAMTOOLS=/home/user_home/wzs/software/samtools/bin/samtools
HTSLIB=/home/user_home/wzs/BAM/samtools/samtools-1.20/htslib-1.20
DEFLATE=/home/user_home/wzs/software/libdeflate_dist/lib
export LD_LIBRARY_PATH="$DEFLATE:${LD_LIBRARY_PATH:-}"
BWA_THREADS=${BWA_THREADS:-24}
SORT_THREADS=${SORT_THREADS:-8}
SORT_MEMORY=${SORT_MEMORY:-8G}
ID=D3-mouse-SRR7768905
REF_NAME=Mus_musculus.GRCm39.dna.primary_assembly.fa
REF_URL=https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna
FASTQ_URL=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR776/005/SRR7768905
stage() { printf '%s\n' "$1" > "$ROOT/stage"; printf '[%(%F %T %Z)T] %s\n' -1 "$1"; }

download() {
    local url=$1 output=$2
    if [[ -f "$output" ]]; then return; fi
    if [[ ${D3_DOWNLOAD_METHOD:-aria2} == aspera && "$url" == "$FASTQ_URL/"* ]]; then
        [[ -f ${D3_ASPERA_KEY:-} ]] || return 1
        timeout 12h bash "$TOOLS/d3_ascp.sh" -T -l "${D3_ASPERA_RATE:-200m}" \
            -P 33001 -i "$D3_ASPERA_KEY" -k 1 --mode=recv --user=era-fasp \
            -L "$ROOT/logs" --host=fasp.sra.ebi.ac.uk \
            "${url#https://ftp.sra.ebi.ac.uk/}" \
            "$output.aspera.part" || return
        mv -- "$output.aspera.part" "$output"
        return
    fi
    timeout 12h aria2c --continue=true \
        --max-connection-per-server="${D3_DOWNLOAD_CONNECTIONS:-1}" \
        --split="${D3_DOWNLOAD_CONNECTIONS:-1}" \
        --min-split-size=16M --file-allocation=none --auto-file-renaming=false \
        --allow-overwrite=false --max-tries=8 --retry-wait=20 --connect-timeout=30 \
        --timeout=60 --summary-interval=60 --console-log-level=notice --enable-color=false \
        --dir="$(dirname "$output")" --out="$(basename "$output").part" "$url" || return
    mv -- "$output.part" "$output"
}
verify_fastq() {
    local file=$1 bytes=$2 md5=$3
    [[ "$md5" =~ ^[[:xdigit:]]{32}$ ]] || return 1
    [[ $(stat -c %s "$file") == "$bytes" ]] || return 1
    printf '%s  %s\n' "$md5" "$file" | md5sum -c -
}

ensure_fastq() {
    local url=$1 output=$2 bytes=$3 md5=$4 attempt quarantine path
    for attempt in 1 2 3; do
        if download "$url" "$output"; then
            if verify_fastq "$output" "$bytes" "$md5"; then return 0; fi
            quarantine=$(mktemp -d "$ROOT/quarantine/checksum-XXXXXXXX") || return
            md5sum "$output" > "$quarantine/actual.md5" || return
            printf 'expected_bytes=%s\nexpected_md5=%s\nurl=%s\n' \
                "$bytes" "$md5" "$url" > "$quarantine/expected.txt"
            for path in "$output" "$output.part" "$output.part.aria2" \
                        "$output.aspera.part" "$output.aspera.part.aspera-ckpt"; do
                if [[ -e "$path" ]]; then mv -- "$path" "$quarantine/" || return; fi
            done
            printf 'FASTQ checksum failed; preserved file in %s\n' "$quarantine"
        else
            printf 'FASTQ download attempt %s failed; retaining partial data for resume\n' "$attempt"
        fi
        if (( attempt < 3 )); then sleep "${D3_RETRY_WAIT:-30}"; fi
    done
    printf 'ERROR: FASTQ not verified after three attempts: %s\n' "$output" >&2
    return 1
}

main() {
mkdir -p "$ROOT"/{fastq,reference,initial,tmp,logs,bin,quarantine}
exec 9>"$ROOT/.prepare.lock"
flock -n 9
exec > >(tee -a "$ROOT/logs/prepare.log") 2>&1
trap 'printf "FAILED line=%s exit=%s\n" "$LINENO" "$?" > "$ROOT/stage"' ERR
stage PREFLIGHT
[[ $(df --output=avail -B1 "$ROOT" | tail -n 1) -ge $((200 * 1024**3)) ]]

g++ -O3 -std=c++11 -Wall -Wextra -I"$HTSLIB" "$TOOLS/d3_prepare_audit.cpp" \
    "$HTSLIB/libhts.a" -L"$DEFLATE" -Wl,-rpath,"$DEFLATE" \
    -ldeflate -lz -lbz2 -llzma -lcurl -lcrypto -lm -lpthread \
    -o "$ROOT/bin/d3_prepare_audit"
{
    hostname
    date --iso-8601=seconds
    "$SAMTOOLS" --version
    "$BWA" 2>&1 || true
    printf 'reference_url=%s/%s.gz\nbwa_threads=%s\nsort_threads=%s\nsort_memory=%s\n' \
        "$REF_URL" "$REF_NAME" "$BWA_THREADS" "$SORT_THREADS" "$SORT_MEMORY"
    printf 'download_method=%s\naspera_rate_cap=%s\n' \
        "${D3_DOWNLOAD_METHOD:-aria2}" "${D3_ASPERA_RATE:-200m}"
    if [[ ${D3_DOWNLOAD_METHOD:-aria2} == aspera ]]; then
        bash "$TOOLS/d3_ascp.sh" -A
        sha256sum "$TOOLS/d3_ascp.sh"
    fi
    sha256sum "$BWA" "$SAMTOOLS" "$TOOLS/prepare_nonhuman_d3.sh" \
        "$TOOLS/d3_prepare_audit.cpp" "$ROOT/bin/d3_prepare_audit"
} > "$ROOT/initial/preparation_environment.txt"

stage DOWNLOAD_R1
R1="$ROOT/fastq/SRR7768905_1.fastq.gz"
R2="$ROOT/fastq/SRR7768905_2.fastq.gz"
ensure_fastq "$FASTQ_URL/$(basename "$R1")" "$R1" 7311488932 38b11ff91acce49be79e25c2340c115e
stage DOWNLOAD_R2
ensure_fastq "$FASTQ_URL/$(basename "$R2")" "$R2" 7295993483 21bc6f6433e3976f2bc7f63d1106dc15
stage CHECK_FASTQ
"$ROOT/bin/d3_prepare_audit" fastq "$R1" "$R2" 82318164 \
    > "$ROOT/initial/fastq_audit.json"

stage DOWNLOAD_REFERENCE
REF="$ROOT/reference/$REF_NAME"
download "$REF_URL/$REF_NAME.gz" "$REF.gz"
[[ $(stat -c %s "$REF.gz") == 806418890 ]]
read -r ref_sum ref_blocks _ < <(sum -r "$REF.gz")
[[ "$ref_sum" == 16996 && "$ref_blocks" == 787519 ]]
download "$REF_URL/CHECKSUMS" "$ROOT/reference/CHECKSUMS"
[[ $(awk -v name="$REF_NAME.gz" '$3 == name {print $1, $2}' "$ROOT/reference/CHECKSUMS") == '16996 787519' ]]
if [[ ! -f "$REF" ]]; then
    gzip -dc "$REF.gz" > "$REF.part"
    mv -- "$REF.part" "$REF"
fi
sha256sum "$REF.gz" "$REF" > "$ROOT/initial/reference.sha256"
stage INDEX_REFERENCE
if [[ ! -f "$ROOT/reference/.index.complete" ]]; then
    "$BWA" index -a bwtsw "$REF"
    "$SAMTOOLS" faidx "$REF"
    touch "$ROOT/reference/.index.complete"
fi

BAM="$ROOT/initial/$ID.coordinate.bam"
if [[ ! -f "$BAM" ]]; then
    stage ALIGN_AND_SORT
    "$BWA" mem -t "$BWA_THREADS" -K 100000000 \
        -R '@RG\tID:SRR7768905\tSM:MPR116\tLB:MWNMPRL007_116_index4\tPL:ILLUMINA' \
        "$REF" "$R1" "$R2" 2> "$ROOT/logs/bwa.stderr.log" |
        "$SAMTOOLS" sort -@ "$SORT_THREADS" -m "$SORT_MEMORY" -l 6 \
            -T "$ROOT/tmp/$ID.sort" -o "$BAM.part" - \
            2> "$ROOT/logs/sort.stderr.log"
    "$SAMTOOLS" quickcheck -v "$BAM.part"
    mv -- "$BAM.part" "$BAM"
fi
stage VALIDATE_BAM
"$SAMTOOLS" quickcheck -v "$BAM"
"$ROOT/bin/d3_prepare_audit" bam "$BAM" 164636328 > "$ROOT/initial/bam_audit.json"
"$SAMTOOLS" index -@ "$SORT_THREADS" "$BAM"
"$SAMTOOLS" view -H "$BAM" > "$ROOT/initial/$ID.header.sam"
"$SAMTOOLS" flagstat -@ "$SORT_THREADS" "$BAM" > "$ROOT/initial/flagstat.txt"
"$SAMTOOLS" idxstats "$BAM" > "$ROOT/initial/idxstats.tsv"
stage HASH_OUTPUT
(
    cd "$ROOT/initial"
    sha256sum "$ID.coordinate.bam" "$ID.coordinate.bam.bai" bam_audit.json \
        fastq_audit.json "$ID.header.sam" flagstat.txt idxstats.tsv \
        preparation_environment.txt reference.sha256 > SHA256SUMS
    sha256sum -c SHA256SUMS
)
touch "$ROOT/.ready"
stage READY
}

if [[ ${BASH_SOURCE[0]} == "$0" ]]; then main; fi
