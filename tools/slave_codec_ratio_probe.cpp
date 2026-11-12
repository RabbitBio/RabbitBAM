#include "codec_ratio_probe.h"
#include "swbam/cpe_codec.h"
#include "swbam/slave/libdeflate.h"
#include <slave.h>

extern "C" void slave_codec_ratio_probe(CodecRatioPara *parameters) {
    const int id = _PEN;
    CodecRatioPara &p = parameters[id];
    p.status = 0;
    if (!p.input) return;
    struct libdeflate_decompressor *decoder = swbam_cpe_get_decompressor(id, nullptr);
    struct libdeflate_compressor *encoder = swbam_cpe_get_compressor(id, 1, 1, nullptr);
    if (!decoder || !encoder ||
        swbam_cpe_decode_bgzf(p.input, p.decoded, decoder, nullptr, nullptr)) {
        p.status = -1;
        return;
    }

    // Restore the bundled compressor's original level-1 defaults (nice=32).
    libdeflate_set_level1_nice_match_length(encoder, 32);
    libdeflate_set_level1_single_probe(encoder, 0);
    size_t original_size = BGZF_MAX_BLOCK_SIZE;
    if (swbam_cpe_compress_bgzf(p.original->data, &original_size,
            p.decoded->data, p.decoded->length, 1, encoder, nullptr, nullptr)) {
        p.status = -2;
        return;
    }
    p.original->length = static_cast<unsigned int>(original_size);

    // Identical payload and BGZF envelope, with the production MPI policy.
    libdeflate_set_level1_nice_match_length(encoder, 4);
    libdeflate_set_level1_single_probe(encoder, 1);
    libdeflate_set_level1_stride2_probe(encoder, 1);
    size_t optimized_size = BGZF_MAX_BLOCK_SIZE;
    if (swbam_cpe_compress_bgzf(p.optimized->data, &optimized_size,
            p.decoded->data, p.decoded->length, 1, encoder, nullptr, nullptr)) {
        p.status = -3;
        return;
    }
    p.optimized->length = static_cast<unsigned int>(optimized_size);
}
