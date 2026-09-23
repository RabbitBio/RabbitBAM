#include "swbam/bam1.h"

#include <cerrno>
#include <climits>
#include <cstring>
#include <vector>

namespace swbam {
namespace cpe {
namespace {

uint32_t ReadLe32(const unsigned char *data) {
    return (uint32_t)data[0] |
           ((uint32_t)data[1] << 8) |
           ((uint32_t)data[2] << 16) |
           ((uint32_t)data[3] << 24);
}

int32_t ReadLe32Signed(const unsigned char *data) {
    const uint32_t value = ReadLe32(data);
    int32_t result = 0;
    memcpy(&result, &value, sizeof(result));
    return result;
}

void WriteHost32(unsigned char *data, uint32_t value) {
    memcpy(data, &value, sizeof(value));
}

int GrowBamData(bam1_t *record, size_t desired,
                std::vector<unsigned char> *scratch) {
    if (!record || !scratch || desired > INT_MAX) return -1;
    if (desired <= record->m_data) return 0;

    const int old_length = record->l_data;
    try {
        scratch->assign(desired, 0);
    } catch (...) {
        return -1;
    }
    if (old_length > 0) {
        if (!record->data || (size_t)old_length > desired) return -1;
        memcpy(scratch->data(), record->data, (size_t)old_length);
    }

    bam1_t source = *record;
    source.data = scratch->data();
    source.l_data = (int)desired;
    source.m_data = desired > UINT32_MAX ? UINT32_MAX : (uint32_t)desired;
    if (!bam_copy1(record, &source)) return -1;
    record->l_data = old_length;
    return 0;
}

int NormalizeLongCigar(bam1_t *record,
                       std::vector<unsigned char> *scratch) {
    bam1_core_t *core = &record->core;
    if (core->n_cigar == 0 || core->tid < 0 || core->pos < 0) return 0;

    uint32_t *cigar = bam_get_cigar(record);
    if (bam_cigar_op(cigar[0]) != BAM_CSOFT_CLIP ||
        bam_cigar_oplen(cigar[0]) != core->l_qseq) {
        return 0;
    }

    const int saved_errno = errno;
    uint8_t *cg_value = bam_aux_get(record, "CG");
    if (!cg_value) {
        if (errno != ENOENT) return -1;
        errno = saved_errno;
        return 0;
    }
    if (cg_value[0] != 'B' || cg_value[1] != 'I') return 0;

    const size_t cg_value_offset =
        (size_t)(cg_value - record->data);
    if (cg_value_offset < 2 || cg_value_offset > (size_t)record->l_data ||
        (size_t)record->l_data - cg_value_offset < 6) {
        return -1;
    }
    const uint32_t cg_count = ReadLe32(cg_value + 2);
    if (cg_count < core->n_cigar || cg_count >= (1U << 29)) return 0;

    const size_t old_cigar_bytes = (size_t)core->n_cigar * 4;
    const size_t new_cigar_bytes = (size_t)cg_count * 4;
    const size_t cigar_offset =
        (size_t)((uint8_t *)cigar - record->data);
    const size_t cg_offset = cg_value_offset - 2;
    const size_t cg_end = cg_offset + 8 + new_cigar_bytes;
    const size_t old_length = (size_t)record->l_data;

    if (cigar_offset > old_length || old_cigar_bytes > old_length - cigar_offset ||
        cg_offset > old_length || cg_end > old_length ||
        new_cigar_bytes < old_cigar_bytes ||
        new_cigar_bytes - old_cigar_bytes > SIZE_MAX - old_length) {
        return -1;
    }
    if (GrowBamData(record,
                    old_length + new_cigar_bytes - old_cigar_bytes,
                    scratch) != 0) {
        return -1;
    }

    cigar = bam_get_cigar(record);
    memmove(record->data + cigar_offset + new_cigar_bytes,
            record->data + cigar_offset + old_cigar_bytes,
            old_length - cigar_offset - old_cigar_bytes);
    memcpy(record->data + cigar_offset,
           record->data + cg_offset +
               (new_cigar_bytes - old_cigar_bytes) + 8,
           new_cigar_bytes);
    if (old_length > cg_end) {
        memmove(record->data + cg_offset +
                    new_cigar_bytes - old_cigar_bytes,
                record->data + cg_end +
                    new_cigar_bytes - old_cigar_bytes,
                old_length - cg_end);
    }
    core->n_cigar = cg_count;
    record->l_data = (int)(old_length - old_cigar_bytes - 8);
    return 0;
}

int ValidateAndUpdateCigar(bam1_t *record) {
    bam1_core_t *core = &record->core;
    if (core->n_cigar == 0) return 0;

    hts_pos_t reference_length = 0;
    hts_pos_t query_length = 0;
    const uint32_t *cigar = bam_get_cigar(record);
    for (uint32_t i = 0; i < core->n_cigar; ++i) {
        const int type = bam_cigar_type(bam_cigar_op(cigar[i]));
        const hts_pos_t length = bam_cigar_oplen(cigar[i]);
        if (type & 1) query_length += length;
        if (type & 2) reference_length += length;
    }
    if ((core->flag & BAM_FUNMAP) || reference_length == 0) {
        reference_length = 1;
    }
    if (core->tid >= 0 && core->pos >= 0) {
        core->bin = hts_reg2bin(
            core->pos, core->pos + reference_length, 14, 5);
    }
    if (core->l_qseq > 0 && !(core->flag & BAM_FUNMAP) &&
        query_length != core->l_qseq) {
        return -1;
    }
    return 0;
}

int MaterializeRecord(const RawBamRecordView &view, bam1_t *output,
                      std::vector<unsigned char> *resize_scratch) {
    if (!output || !resize_scratch || !view.encoded ||
        view.encoded_size < 36 || view.block_size < 32 ||
        (uint64_t)view.block_size + 4 != view.encoded_size) {
        return -1;
    }

    const unsigned char *encoded = view.encoded;
    const unsigned char *fields = encoded + 4;
    const unsigned char *payload = encoded + 36;
    const size_t payload_size = (size_t)view.block_size - 32;
    const uint32_t bin_mq_nl = ReadLe32(fields + 8);
    const uint32_t flag_nc = ReadLe32(fields + 12);
    const uint32_t raw_qname_length = bin_mq_nl & 0xff;
    const int32_t query_length = ReadLe32Signed(fields + 16);

    if (raw_qname_length < 1 || raw_qname_length > payload_size ||
        query_length < 0) {
        return -1;
    }
    const uint32_t cigar_count = flag_nc & 0xffff;
    const uint64_t minimum_payload =
        (uint64_t)raw_qname_length + (uint64_t)cigar_count * 4 +
        (((uint64_t)(uint32_t)query_length + 1) >> 1) +
        (uint64_t)(uint32_t)query_length;
    if (minimum_payload > payload_size) return -1;

    const bool qname_terminated = payload[raw_qname_length - 1] == '\0';
    uint32_t qname_padding = (0U - raw_qname_length) & 3U;
    if (!qname_terminated && qname_padding == 0) qname_padding = 4;
    const uint32_t final_qname_length = raw_qname_length + qname_padding;
    if (payload_size > (size_t)INT_MAX - qname_padding) return -1;
    const size_t normalized_size = payload_size + qname_padding;

    if (GrowBamData(output, normalized_size, resize_scratch) != 0) {
        return -1;
    }
    unsigned char *normalized = output->data;
    memcpy(normalized, payload, raw_qname_length);
    memset(normalized + raw_qname_length, 0, qname_padding);
    memcpy(normalized + final_qname_length,
           payload + raw_qname_length,
           payload_size - raw_qname_length);

    // BAM stores CIGAR words little-endian; bam1_t accessors expect host words.
    unsigned char *normalized_cigar =
        normalized + final_qname_length;
    const unsigned char *encoded_cigar = payload + raw_qname_length;
    for (uint32_t i = 0; i < cigar_count; ++i) {
        WriteHost32(normalized_cigar + (size_t)i * 4,
                    ReadLe32(encoded_cigar + (size_t)i * 4));
    }

    memset(&output->core, 0, sizeof(output->core));
    output->core.tid = ReadLe32Signed(fields);
    output->core.pos = ReadLe32Signed(fields + 4);
    output->core.bin = bin_mq_nl >> 16;
    output->core.qual = (bin_mq_nl >> 8) & 0xff;
    output->core.l_qname = final_qname_length;
    output->core.l_extranul = qname_terminated
        ? qname_padding : qname_padding - 1;
    output->core.flag = flag_nc >> 16;
    output->core.n_cigar = cigar_count;
    output->core.l_qseq = query_length;
    output->core.mtid = ReadLe32Signed(fields + 20);
    output->core.mpos = ReadLe32Signed(fields + 24);
    output->core.isize = ReadLe32Signed(fields + 28);
    output->l_data = (int)normalized_size;
    output->id = 0;

    if (NormalizeLongCigar(output, resize_scratch) != 0 ||
        ValidateAndUpdateCigar(output) != 0) {
        return -1;
    }
    output->id = 0;
    return 0;
}

class Bam1BatchPostProcessorAdapter : public RawBamBatchPostProcessor {
public:
    Bam1BatchPostProcessorAdapter(Bam1BatchPostProcessor *post_processor,
                     ComposableBam1Metrics *metrics)
        : post_processor_(post_processor), metrics_(metrics) {}

    ~Bam1BatchPostProcessorAdapter() {
        for (size_t i = 0; i < pool_.size(); ++i) {
            bam_destroy1(pool_[i]);
        }
    }

    int PostProcessRawBatch(const RawBamRecordView *records, size_t count) {
        if (!post_processor_ || (!records && count != 0)) return -1;
        if (EnsurePool(count) != 0) return -1;

        long long encoded_bytes = 0;
        long long data_bytes = 0;
        const double materialize_t0 = GetTime();
        for (size_t i = 0; i < count; ++i) {
            if (MaterializeRecord(records[i], pool_[i],
                                  &resize_scratch_) != 0) {
                return -1;
            }
            batch_[i] = pool_[i];
            encoded_bytes += records[i].encoded_size;
            data_bytes += pool_[i]->l_data;
        }
        const double materialize = GetTime() - materialize_t0;

        const double post_process_t0 = GetTime();
        if (post_processor_->PostProcessBam1Batch(
                count == 0 ? nullptr : batch_.data(), count) != 0) {
            return -1;
        }
        const double post_process = GetTime() - post_process_t0;

        if (metrics_) {
            metrics_->records += (long long)count;
            metrics_->encoded_bytes += encoded_bytes;
            metrics_->data_bytes += data_bytes;
            if (count > metrics_->peak_batch_records) {
                metrics_->peak_batch_records = count;
            }
            metrics_->materialize += materialize;
            metrics_->post_process += post_process;
        }
        return 0;
    }

private:
    int EnsurePool(size_t count) {
        try {
            batch_.resize(count);
        } catch (...) {
            return -1;
        }
        while (pool_.size() < count) {
            bam1_t *record = bam_init1();
            if (!record) return -1;
            try {
                pool_.push_back(record);
            } catch (...) {
                bam_destroy1(record);
                return -1;
            }
        }
        return 0;
    }

    Bam1BatchPostProcessor *post_processor_;
    ComposableBam1Metrics *metrics_;
    std::vector<bam1_t *> pool_;
    std::vector<const bam1_t *> batch_;
    std::vector<unsigned char> resize_scratch_;
};

} // namespace

ComposableBam1Metrics::ComposableBam1Metrics()
    : records(0), encoded_bytes(0), data_bytes(0),
      peak_batch_records(0), materialize(0.0), post_process(0.0) {}

int RunComposableBam1Pipeline(
        const BamInputBackend &input,
        const BgzfBlockSpan *spans,
        size_t span_count,
        Bam1BatchPostProcessor *post_processor,
        CpeReadPipelineTiming *timing,
        ComposableDecodeMetrics *decode_metrics,
        ComposableBam1Metrics *bam1_metrics,
        const CpeReadPipelineOptions &options) {
    if (!post_processor) return -1;
    if (bam1_metrics) *bam1_metrics = ComposableBam1Metrics();
    Bam1BatchPostProcessorAdapter adapter(post_processor, bam1_metrics);
    ComposableRawBamMetrics raw_metrics;
    return RunComposableRawBamPipeline(
        input, spans, span_count, &adapter,
        timing, decode_metrics, &raw_metrics, options);
}

} // namespace cpe
} // namespace swbam
