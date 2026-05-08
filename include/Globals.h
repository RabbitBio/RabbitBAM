#ifndef H_GLOBALS
#define H_GLOBALS

// 修复神威平台从核缺失 max_align_t 的 编译Bug
#if defined(__sw_slave__)
#ifndef _MAX_ALIGN_T
#define _MAX_ALIGN_T
typedef struct {
  long long __max_align_ll __attribute__((__aligned__(__alignof__(long long))));
  long double __max_align_ld __attribute__((__aligned__(__alignof__(long double))));
} max_align_t;
#endif
#endif

#include <sys/time.h>
#include <cstddef>

inline double GetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}

namespace rabbit {

    // basic types
    typedef char int8;
    typedef unsigned char uchar, byte, uint8;
    typedef short int int16;
    typedef unsigned short int uint16;
    typedef int int32;
    typedef unsigned int uint32;
    typedef long long int64;
    typedef unsigned long long uint64;


}// namespace rabbit

#endif // H_GLOBALS