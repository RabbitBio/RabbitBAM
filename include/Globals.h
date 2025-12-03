#ifndef H_GLOBALS
#define H_GLOBALS


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