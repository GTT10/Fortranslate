#define _GNU_SOURCE
#include <mpi.h>
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Test-only fault injection for OpenMPI's mpi_f08 -> PMPI C binding. */
static int requested_rank(void)
{
    const char *value = getenv("PELEF_TEST_FAIL_MPI_RANK");
    return value ? atoi(value) : 0;
}
static int requested_mode(const char *mode)
{
    const char *value = getenv("PELEF_TEST_FAIL_MPI_MODE");
    return !strcmp(value ? value : "write", mode);
}
int PMPI_File_write_at(MPI_File file, MPI_Offset offset, const void *buffer,
                      int count, MPI_Datatype datatype, MPI_Status *status)
{
    static int calls;
    int rank = -1;
    const char *requested = getenv("PELEF_TEST_FAIL_MPI_WRITE_CALL");
    typedef int (*write_at_type)(MPI_File, MPI_Offset, const void *, int,
                                 MPI_Datatype, MPI_Status *);
    write_at_type real_write = (write_at_type)dlsym(RTLD_NEXT, "PMPI_File_write_at");
    if (!real_write) { fprintf(stderr, "FAULT_SHIM missing PMPI_File_write_at\n"); abort(); }
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ++calls;
    if (requested_mode("write") && requested && rank == requested_rank() && calls == atoi(requested)) {
        fprintf(stderr, "FAULT_SHIM mode=write rank=%d call=%d offset=%lld bytes=%d MPI_ERR_IO\n",
                rank, calls, (long long)offset, count);
        fflush(stderr);
        return MPI_ERR_IO;
    }
    return real_write(file, offset, buffer, count, datatype, status);
}
int PMPI_File_close(MPI_File *file)
{
    static int calls;
    int rank = -1;
    typedef int (*close_type)(MPI_File *);
    close_type real_close = (close_type)dlsym(RTLD_NEXT, "PMPI_File_close");
    if (!real_close) { fprintf(stderr, "FAULT_SHIM missing PMPI_File_close\n"); abort(); }
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int result = real_close(file);
    ++calls;
    if (requested_mode("close") && rank == requested_rank() && calls == 1 && result == MPI_SUCCESS) {
        fprintf(stderr, "FAULT_SHIM mode=close rank=%d call=%d MPI_ERR_IO after real close\n", rank, calls);
        fflush(stderr);
        return MPI_ERR_IO;
    }
    return result;
}
