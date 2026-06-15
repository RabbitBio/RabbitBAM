
#include "CmdInfo.h"

CmdInfo::CmdInfo() {
    in_file_name_ = "";
    out_file_name_ = "";
    verbose_ = false;
    validate_bounds_ = false;
    min_mapq_ = -1;
    max_mapq_ = -1;
    require_flag_ = 0;
    exclude_flag_ = 0;
    ref_name_ = "";
    min_read_len_ = -1;
    max_read_len_ = -1;
    compress_level_ = 1;
    sort_memory_ = "";
    sort_temp_prefix_ = "";
    collate_bins_ = 64;
    collate_memory_ = "";
    collate_temp_prefix_ = "";
    markdup_remove_dups_ = false;
    markdup_clear_ = false;
    markdup_include_fails_ = false;
    markdup_memory_ = "";
    fixmate_mate_score_ = false;
}
