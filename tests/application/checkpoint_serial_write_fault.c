#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Test-only write error for the unit fixture's first checkpoint destination. */
ssize_t write(int fd, const void *buffer, size_t count)
{
    typedef ssize_t (*write_type)(int, const void *, size_t);
    write_type real_write = (write_type)dlsym(RTLD_NEXT, "write");
    char link[64], path[4096];
    if (!real_write) _exit(97);
    snprintf(link, sizeof(link), "/proc/self/fd/%d", fd);
    ssize_t n = readlink(link, path, sizeof(path)-1);
    if (n >= 0) {
        path[n] = '\0';
        const char *name = strrchr(path, '/');
        name = name ? name+1 : path;
        if ((!strcmp(name,"tree.chk") || !strcmp(name,"tree.chk.tmp")) &&
            lseek(fd, 0, SEEK_CUR) >= 8192) {
            const char marker[] = "FAULT_SHIM serial checkpoint write ENOSPC\n";
            real_write(STDERR_FILENO, marker, sizeof(marker)-1);
            errno = ENOSPC;
            return -1;
        }
    }
    return real_write(fd, buffer, count);
}
