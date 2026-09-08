// rb_sort.cpp - `rabbitbam sortbam`: coordinate-sort a BAM.
// Extracted from RabbitBin. Records are stable-sorted in memory by reference,
// position and strand, then serialized and BGZF-compressed in parallel using
// libdeflate. External sorting and a memory limit are not implemented.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <string>
#include <vector>

#include <parallel/algorithm>

#include <omp.h>
#include <libdeflate.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "rb_bam_cli.h"
#include "rb_bam_sort.h"

static void sort_usage() {
  fprintf(stderr,
      "\nrabbitbam sortbam: coordinate-sort a BAM (parallel BGZF write)\n\n"
      "Usage: rabbitbam sortbam [options] <in.bam>\n"
      "       cat in.bam | rabbitbam sortbam -o out.sorted.bam -\n\n"
      "Options:\n"
      "  -o, --out FILE     Output sorted BAM (default: stdout)\n"
      "  -t, --threads N    Sorting and BGZF I/O threads (default: 8)\n"
      "  -l, --level N      libdeflate compression level 1-12 (default: 6)\n"
      "      --write-index  Also write <out>.bai (output must be a file)\n"
      "  -h, --help         Show this help\n\n"
      "Output records are ordered identically to `samtools sort`.\n");
}

// samtools coordinate-sort key (htslib 1.x bam1_lt):
//   ((uint64_t)tid << 32) | ((pos+1) << 1) | is_reverse
// so at equal (tid,pos) forward reads precede reverse reads, and tid<0 (unmapped)
// sign-extends to a huge value and therefore sorts last. Combined with a STABLE
// sort (equal keys keep input order) this reproduces `samtools sort` exactly.
static inline uint64_t rb_sort_key(const bam1_t *b) {
  uint64_t pos1 = (uint64_t)(b->core.pos + 1);
  uint64_t rev = (b->core.flag & BAM_FREVERSE) ? 1u : 0u;
  return ((uint64_t)b->core.tid << 32) | (pos1 << 1) | rev;
}

// Stable coordinate sort of an in-memory vector of bam1_t*, exposed for reuse.
void rb_stable_coord_sort(std::vector<bam1_t *> &recs) {
  // Parallel STABLE sort: __gnu_parallel::stable_sort preserves equal-key input
  // order exactly like std::stable_sort (uses OpenMP threads; thread count is
  // omp_get_max_threads(), which the CLI sets from -t). The output order
  // is therefore identical to the serial sort -> byte-identical to samtools.
  __gnu_parallel::stable_sort(recs.begin(), recs.end(),
                              [](const bam1_t *a, const bam1_t *b) {
                                return rb_sort_key(a) < rb_sort_key(b);
                              });
}

// ── Parallel BGZF writer ────────────────────────────────────────────────────
// htslib's writer caps at ~8 cores because sam_write1 serialises every record
// on one thread. We instead (a) serialise all bam1_t records to BAM wire bytes
// in parallel, (b) cut the stream into 64 KiB BGZF blocks and libdeflate-compress
// them in parallel, then (c) write the blocks in order. The DECODED bytes are
// identical to htslib's output (same header + records, same record encoding as
// bam_write1), so `samtools view` / depth / binning see identical content; only
// the BGZF block framing (and hence the .bam md5) differs. The .bai is still
// built from the finished file with sam_index_build3 (a fast parallel re-read).

// On-disk size of one BAM record (block_size i32 + 32 core bytes + payload,
// where the in-memory l_extranul qname padding is dropped, exactly as
// bam_write1 does).
static inline size_t bam_wire_size(const bam1_t *b) {
  return 4u + 32u + (size_t)(b->l_data - b->core.l_extranul);
}

// Serialise one record to `dst` (must hold bam_wire_size(b) bytes); little-endian
// host (x86) so the int32 fields are written by memcpy, matching bam_write1.
static inline void bam_wire_write(const bam1_t *b, unsigned char *dst) {
  const bam1_core_t *c = &b->core;
  uint8_t ex = c->l_extranul;
  uint32_t qn = (uint32_t)c->l_qname - ex;           // on-disk read-name length
  int32_t block_len = 32 + (int32_t)(b->l_data - ex);
  uint32_t x[8];
  x[0] = (uint32_t)c->tid;
  x[1] = (uint32_t)c->pos;
  x[2] = ((uint32_t)c->bin << 16) | ((uint32_t)c->qual << 8) | qn;
  x[3] = ((uint32_t)c->flag << 16) | (c->n_cigar & 0xffff);
  x[4] = (uint32_t)c->l_qseq;
  x[5] = (uint32_t)c->mtid;
  x[6] = (uint32_t)c->mpos;
  x[7] = (uint32_t)c->isize;
  memcpy(dst, &block_len, 4);
  memcpy(dst + 4, x, 32);
  memcpy(dst + 36, b->data, qn);                      // qname (no extranul)
  uint32_t rest = (uint32_t)b->l_data - (uint32_t)c->l_qname;  // cigar+seq+qual+aux
  memcpy(dst + 36 + qn, b->data + c->l_qname, rest);
}

// Standard 28-byte empty BGZF block that marks end-of-file.
static const unsigned char kBgzfEof[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
    0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00};

// Build the uncompressed BAM header block bytes ("BAM\1", text, refs).
static void build_bam_header_bytes(sam_hdr_t *hdr, std::vector<unsigned char> &H) {
  const char *text = sam_hdr_str(hdr);
  size_t l_text = sam_hdr_length(hdr);
  int n_ref = sam_hdr_nref(hdr);
  if (!text) { text = ""; l_text = 0; }
  H.clear();
  H.reserve(12 + l_text + (size_t)n_ref * 24 + 64);
  auto put32 = [&](int32_t v) {
    unsigned char p[4];
    memcpy(p, &v, 4);
    H.insert(H.end(), p, p + 4);
  };
  H.push_back('B'); H.push_back('A'); H.push_back('M'); H.push_back(1);
  put32((int32_t)l_text);
  H.insert(H.end(), text, text + l_text);
  put32(n_ref);
  for (int i = 0; i < n_ref; ++i) {
    const char *nm = sam_hdr_tid2name(hdr, i);
    int l_ref = (int)sam_hdr_tid2len(hdr, i);
    size_t nl = nm ? strlen(nm) : 0;
    put32((int32_t)(nl + 1));
    if (nl) H.insert(H.end(), nm, nm + nl);
    H.push_back(0);
    put32(l_ref);
  }
}

// Parallel BGZF write of `recs` (already coordinate-sorted) to `out_path`.
// Returns 0 on success, non-zero on error. (The .bai is built separately by the
// caller via sam_index_build3, which is fast and parallel; see rb_write_sorted_bam.)
static int write_bam_parallel(sam_hdr_t *hdr, std::vector<bam1_t *> &recs,
                              const std::string &out_path, int threads,
                              int level) {
  int clevel = level;
  if (clevel < 1) clevel = 1;
  if (clevel > 12) clevel = 12;  // libdeflate range
  if (threads < 1) threads = 1;

  std::vector<unsigned char> H;
  build_bam_header_bytes(hdr, H);
  const size_t header_len = H.size();
  const size_t nrec = recs.size();

  // Per-record offsets (prefix sum of wire sizes).
  std::vector<size_t> roff(nrec + 1);
  const long long nl = (long long)nrec;
#pragma omp parallel for num_threads(threads) schedule(static)
  for (long long i = 0; i < nl; ++i) roff[i + 1] = bam_wire_size(recs[i]);
  roff[0] = 0;
  for (size_t i = 0; i < nrec; ++i) roff[i + 1] += roff[i];
  const size_t rec_total = roff[nrec];

  // Assemble the full uncompressed BAM stream U = header ++ records.
  const size_t U_size = header_len + rec_total;
  std::vector<unsigned char> U(U_size);
  if (header_len) memcpy(U.data(), H.data(), header_len);
  std::vector<unsigned char>().swap(H);
  unsigned char *Ud = U.data() + header_len;
#pragma omp parallel for num_threads(threads) schedule(static)
  for (long long i = 0; i < nl; ++i) bam_wire_write(recs[i], Ud + roff[i]);

  // Cut U into BGZF blocks (<=0xff00 uncompressed) and compress in parallel.
  const size_t UBS = 0xff00;                 // 65280 uncompressed per block
  const size_t STRIDE = 0x10000 + 64;        // max compressed block slot
  const size_t nblk = (U_size + UBS - 1) / UBS;
  std::vector<unsigned char> ob(nblk ? nblk * STRIDE : STRIDE);
  std::vector<uint32_t> blen(nblk, 0);
  std::atomic<int> bad{0};
#pragma omp parallel num_threads(threads)
  {
    struct libdeflate_compressor *comp = libdeflate_alloc_compressor(clevel);
    if (!comp) {
      bad.fetch_add(1);
    } else {
#pragma omp for schedule(dynamic, 8)
      for (long long bi = 0; bi < (long long)nblk; ++bi) {
        size_t s = (size_t)bi * UBS;
        size_t len = (U_size - s < UBS) ? (U_size - s) : UBS;
        unsigned char *blk = ob.data() + (size_t)bi * STRIDE;
        size_t csize = libdeflate_deflate_compress(comp, U.data() + s, len,
                                                   blk + 18, STRIDE - 26);
        if (csize == 0) { bad.fetch_add(1); continue; }
        uint32_t crc = libdeflate_crc32(0, U.data() + s, len);
        size_t bsize = 18 + csize + 8;  // total block bytes
        static const unsigned char hh[16] = {
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
            0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00};
        memcpy(blk, hh, 16);
        uint16_t bf = (uint16_t)(bsize - 1);
        blk[16] = (unsigned char)(bf & 0xff);
        blk[17] = (unsigned char)((bf >> 8) & 0xff);
        unsigned char *tr = blk + 18 + csize;
        memcpy(tr, &crc, 4);
        uint32_t isz = (uint32_t)len;
        memcpy(tr + 4, &isz, 4);
        blen[bi] = (uint32_t)bsize;
      }
      libdeflate_free_compressor(comp);
    }
  }
  std::vector<unsigned char>().swap(U);
  if (bad.load() != 0) {
    fprintf(stderr, "[Error!] parallel BGZF compression failed\n");
    return 1;
  }

  // Write blocks in order + EOF marker.
  const bool to_stdout = out_path.empty() || out_path == "-";
  FILE *f = to_stdout ? stdout : fopen(out_path.c_str(), "wb");
  if (!f) {
    fprintf(stderr, "[Error!] cannot open output: %s\n", out_path.c_str());
    return 1;
  }
  int rc = 0;
  for (size_t bi = 0; bi < nblk; ++bi) {
    if (fwrite(ob.data() + bi * STRIDE, 1, blen[bi], f) != blen[bi]) {
      rc = 1;
      break;
    }
  }
  if (!rc && fwrite(kBgzfEof, 1, sizeof(kBgzfEof), f) != sizeof(kBgzfEof))
    rc = 1;
  if (to_stdout) fflush(f);
  else if (fclose(f) != 0) rc = 1;
  if (rc) fprintf(stderr, "[Error!] failed writing BAM output\n");
  return rc;
}

// Free a large record vector in parallel (66M bam_destroy1 calls are ~18s
// single-threaded). glibc free() is thread-safe; the records are independent.
void rb_free_records(std::vector<bam1_t *> &recs, int threads) {
  if (threads < 1) threads = 1;
  const long long n = (long long)recs.size();
#pragma omp parallel for num_threads(threads) schedule(static)
  for (long long i = 0; i < n; ++i)
    if (recs[i]) bam_destroy1(recs[i]);
  recs.clear();
}

// Build a .bai for a finished BAM via htslib's (fast, parallel) indexer.
// Exposed so callers can run it on the main thread after a background
// write (htslib's index thread pool does not engage from a worker thread).
int rb_index_bam(const std::string &bam_path, int threads) {
  std::string idx = bam_path + ".bai";
  int ir = sam_index_build3(bam_path.c_str(), idx.c_str(), /*min_shift=*/0,
                            threads);
  if (ir < 0) {
    fprintf(stderr, "[Warn] could not build index %s (code %d)\n", idx.c_str(),
            ir);
    return ir;
  }
  return 0;
}

// Write a sorted record set to `out_path` ("-"/empty = stdout) as BAM. Updates
// the header sort-order to coordinate. Optionally builds a .bai. Returns 0 ok.
int rb_write_sorted_bam(sam_hdr_t *hdr, std::vector<bam1_t *> &recs,
                        const std::string &out_path, int threads, int level,
                        bool write_index) {
  // Mark header as coordinate-sorted (samtools sets @HD SO:coordinate).
  sam_hdr_update_hd(hdr, "SO", "coordinate");
  int rc = write_bam_parallel(hdr, recs, out_path, threads, level);
  if (rc != 0) return rc;
  const bool to_stdout = out_path.empty() || out_path == "-";
  if (write_index && !to_stdout) rc = rb_index_bam(out_path, threads);
  return rc;
}

int rb_cmd_sortbam(int ac, char *av[]) {
  std::string in_path, out_path;
  int threads = 8, level = 6;
  bool write_index = false;

  enum { OPT_WRITE_INDEX = 1000 };
  static const struct option longopts[] = {
      {"out", required_argument, 0, 'o'},
      {"threads", required_argument, 0, 't'},
      {"level", required_argument, 0, 'l'},
      {"write-index", no_argument, 0, OPT_WRITE_INDEX},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};
  optind = 1;
  int c;
  while ((c = getopt_long(ac, av, "o:t:l:h", longopts, nullptr)) != -1) {
    switch (c) {
      case 'o': out_path = optarg; break;
      case 't': threads = atoi(optarg); break;
      case 'l': level = atoi(optarg); break;
      case OPT_WRITE_INDEX: write_index = true; break;
      case 'h': sort_usage(); return 0;
      default: sort_usage(); return 1;
    }
  }
  if (optind < ac) in_path = av[optind];
  if (in_path.empty()) in_path = "-";  // read BAM from stdin
  if (threads < 1) threads = 1;
  omp_set_num_threads(threads);

  fprintf(stderr, "[rabbitbam sortbam] reading %s (%d threads)\n",
          in_path.c_str(), threads);

  samFile *in = sam_open(in_path.c_str(), "rb");
  if (!in) {
    fprintf(stderr, "[Error!] cannot open input: %s\n", in_path.c_str());
    return 1;
  }
  if (threads > 1) hts_set_threads(in, threads);
  sam_hdr_t *hdr = sam_hdr_read(in);
  if (!hdr) {
    fprintf(stderr, "[Error!] cannot read BAM header\n");
    sam_close(in);
    return 1;
  }

  std::vector<bam1_t *> recs;
  recs.reserve(1u << 20);
  int64_t n = 0;
  for (;;) {
    bam1_t *b = bam_init1();
    int r = sam_read1(in, hdr, b);
    if (r < -1) {
      fprintf(stderr, "[Error!] truncated/corrupt BAM at record %lld\n",
              (long long)n);
      bam_destroy1(b);
      // fall through to free + fail
      for (auto *p : recs) bam_destroy1(p);
      sam_hdr_destroy(hdr);
      sam_close(in);
      return 1;
    }
    if (r == -1) {  // EOF
      bam_destroy1(b);
      break;
    }
    recs.push_back(b);
    ++n;
  }
  sam_close(in);
  fprintf(stderr, "[rabbitbam sortbam] read %lld records; sorting\n",
          (long long)n);

  rb_stable_coord_sort(recs);

  fprintf(stderr, "[rabbitbam sortbam] writing %s\n",
          (out_path.empty() || out_path == "-") ? "<stdout>"
                                                 : out_path.c_str());
  int rc = rb_write_sorted_bam(hdr, recs, out_path, threads, level,
                               write_index);

  for (auto *p : recs) bam_destroy1(p);
  sam_hdr_destroy(hdr);
  if (rc == 0) fprintf(stderr, "[rabbitbam sortbam] done\n");
  return rc;
}
