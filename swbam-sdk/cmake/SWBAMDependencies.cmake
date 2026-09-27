function(swbam_define_dependencies hts_include hts_library deflate_include deflate_library)
    foreach(path IN ITEMS "${hts_include}/htslib/sam.h" "${hts_library}"
                          "${deflate_include}/libdeflate.h" "${deflate_library}"
                          "${SUNWAY_MPI_ROOT}/include/mpi.h"
                          "${SUNWAY_MPI_ROOT}/lib/single_static/libmpi.a")
        if(NOT EXISTS "${path}")
            message(FATAL_ERROR "Missing SWBAM dependency: ${path}")
        endif()
    endforeach()
    if(NOT TARGET swbam::htslib)
        add_library(swbam::htslib STATIC IMPORTED GLOBAL)
        set_target_properties(swbam::htslib PROPERTIES
            IMPORTED_LOCATION "${hts_library}"
            INTERFACE_INCLUDE_DIRECTORIES "${hts_include}")
    endif()
    if(NOT TARGET swbam::deflate)
        add_library(swbam::deflate STATIC IMPORTED GLOBAL)
        set_target_properties(swbam::deflate PROPERTIES
            IMPORTED_LOCATION "${deflate_library}"
            INTERFACE_INCLUDE_DIRECTORIES "${deflate_include}")
    endif()
    if(NOT TARGET swbam::mpi_platform)
        add_library(swbam::mpi_platform INTERFACE IMPORTED GLOBAL)
        set_target_properties(swbam::mpi_platform PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${SUNWAY_MPI_ROOT}/include"
            INTERFACE_LINK_DIRECTORIES "${SUNWAY_MPI_ROOT}/lib/single_static;/usr/sw/lib"
            INTERFACE_LINK_LIBRARIES "mpicxx;mpi;swverbs_ft_perf;swutils;swtm;m;dl")
    endif()
endfunction()
