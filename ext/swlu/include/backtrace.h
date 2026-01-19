#ifndef SWLU_BACKTRACE_H
#define SWLU_BACKTRACE_H

#include <stdint.h>
#include <ucontext.h>

#ifdef __cplusplus
extern "C" int swlu_backtrace(uint64_t*, int);
extern "C" int swlu_backtrace_context(ucontext_t*, uint64_t*, int);
extern "C" int swlu_backtrace_symbols(uint64_t*, int);
extern "C" void swlu_show_callstack(int);
#else
int swlu_backtrace(uint64_t*, int);
int swlu_backtrace_context(ucontext_t*, uint64_t*, int);
int swlu_backtrace_symbols(uint64_t*, int);
void swlu_show_callstack(int);
#endif

#endif
