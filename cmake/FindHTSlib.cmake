find_path(HTSlib_INCLUDE_DIR htslib/sam.h
  HINTS "${HTSLIB_ROOT}" ENV HTSLIB_ROOT
  PATH_SUFFIXES include)
find_library(HTSlib_LIBRARY NAMES hts
  HINTS "${HTSLIB_ROOT}" ENV HTSLIB_ROOT
  PATH_SUFFIXES lib lib64)

# Older headers lack HTS_VERSION and cannot provide the required BAM APIs.
set(HTSlib_VERSION "0.0.0")
if(HTSlib_INCLUDE_DIR AND EXISTS "${HTSlib_INCLUDE_DIR}/htslib/hts.h")
  file(STRINGS "${HTSlib_INCLUDE_DIR}/htslib/hts.h" _hts_version_line
    REGEX "^#define HTS_VERSION +[0-9]+")
  if(_hts_version_line)
    string(REGEX MATCH "[0-9]+" _hts_version_number "${_hts_version_line}")
    math(EXPR _hts_major "${_hts_version_number} / 100000")
    math(EXPR _hts_minor "(${_hts_version_number} / 100) % 1000")
    math(EXPR _hts_patch "${_hts_version_number} % 100")
    set(HTSlib_VERSION "${_hts_major}.${_hts_minor}.${_hts_patch}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HTSlib
  REQUIRED_VARS HTSlib_INCLUDE_DIR HTSlib_LIBRARY
  VERSION_VAR HTSlib_VERSION)
mark_as_advanced(HTSlib_INCLUDE_DIR HTSlib_LIBRARY)
