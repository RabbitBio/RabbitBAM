#include "swbam/operators/record_count.h"
#include "swbam/cpe_codec.h"

#include <cstring>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

extern "C" void slave_swbam_record_count(
        swbam::operators::RecordCountPara paras[64]) {
    swbam::operators::RecordCountPara &para = paras[_PEN];
    if (!para.compressed) return;
    struct libdeflate_decompressor *decoder =
        swbam_cpe_get_decompressor(_PEN, nullptr);
    if (!decoder || swbam_cpe_decode_bgzf(
            para.compressed, para.decoded, decoder, nullptr, nullptr) != 0) {
        para.status = -1;
        return;
    }

    const unsigned char *data = para.decoded->data;
    const size_t length = para.decoded->length;
    size_t offset = 0;
    long long records = 0;
    while (offset < length) {
        if (length - offset < 36) {
            para.status = -2;
            return;
        }
        int32_t record_bytes;
        memcpy(&record_bytes, data + offset, sizeof(record_bytes));
        if (record_bytes < 32 ||
            static_cast<size_t>(record_bytes) > length - offset - 4) {
            para.status = -2;
            return;
        }
        offset += static_cast<size_t>(record_bytes) + 4;
        ++records;
    }
    para.records = records;
    para.status = 0;
}
