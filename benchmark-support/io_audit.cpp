#include "storage.h"
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <limits.h>
#include <unistd.h>

namespace {
bool active = false;
int calls = 0;
bool Matches(const char *path, const std::string &name) {
    const size_t length = strlen(path);
    return !name.empty() && length >= name.size() &&
        !strcmp(path + length - name.size(), name.c_str());
}
void Observe(int fd) {
    if (!active || !swbam_usage::IsRam() || fd < 0) return;
    char link[64], path[PATH_MAX];
    snprintf(link, sizeof(link), "/proc/self/fd/%d", fd);
    const ssize_t n = readlink(link, path, sizeof(path) - 1);
    if (n < 0) return;
    path[n] = 0;
    const swbam_usage::Session *session = swbam_usage::CurrentSession();
    if (Matches(path, session->input_path()) || Matches(path, session->output_path()) ||
        strstr(path, "rabbitbam-") || strstr(path, "swbam-conversion.body") || strstr(path, "sdk-body")) {
        ++calls;
        fprintf(stderr, "ERROR: physical user-data I/O in RAM MPI interval: %s\n", path);
    }
}
}
namespace swbam_usage {
void BeginMpiInterval() { calls = 0; active = true; }
void EndMpiInterval() { active = false; }
int PhysicalDataIoCalls() { return calls; }
}
extern "C" {
int __real_open(const char *, int, ...);
int __real_open64(const char *, int, ...);
int __real_close(int);
ssize_t __real_read(int, void *, size_t);
ssize_t __real_write(int, const void *, size_t);
ssize_t __real_pread(int, void *, size_t, off_t);
ssize_t __real_pwrite(int, const void *, size_t, off_t);
FILE *__real_fopen64(const char *, const char *);
size_t __real_fread(void *, size_t, size_t, FILE *);
size_t __real_fwrite(const void *, size_t, size_t, FILE *);
int __wrap_open(const char *path, int flags, ...) {
    mode_t mode = 0;
    if (flags & O_CREAT) { va_list args; va_start(args, flags); mode = va_arg(args, int); va_end(args); }
    int fd = __real_open(path, flags, mode); Observe(fd); return fd;
}
int __wrap_open64(const char *path, int flags, ...) {
    mode_t mode = 0;
    if (flags & O_CREAT) { va_list args; va_start(args, flags); mode = va_arg(args, int); va_end(args); }
    int fd = __real_open64(path, flags, mode); Observe(fd); return fd;
}
int __wrap_close(int fd) { Observe(fd); return __real_close(fd); }
ssize_t __wrap_read(int fd, void *data, size_t n) { Observe(fd); return __real_read(fd, data, n); }
ssize_t __wrap_write(int fd, const void *data, size_t n) { Observe(fd); return __real_write(fd, data, n); }
ssize_t __wrap_pread(int fd, void *data, size_t n, off_t pos) { Observe(fd); return __real_pread(fd, data, n, pos); }
ssize_t __wrap_pwrite(int fd, const void *data, size_t n, off_t pos) { Observe(fd); return __real_pwrite(fd, data, n, pos); }
FILE *__wrap_fopen64(const char *path, const char *mode) {
    FILE *file = __real_fopen64(path, mode); if (file) Observe(fileno(file)); return file;
}
size_t __wrap_fread(void *data, size_t size, size_t n, FILE *file) {
    Observe(fileno(file)); return __real_fread(data, size, n, file);
}
size_t __wrap_fwrite(const void *data, size_t size, size_t n, FILE *file) {
    Observe(fileno(file)); return __real_fwrite(data, size, n, file);
}
}
