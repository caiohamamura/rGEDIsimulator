#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "functionWrappers.h"

int msgf(const char *format, ...) {
    va_list argptr;
    int result;
    va_start(argptr, format);
    result = vfprintf(stdout, format, argptr);
    va_end(argptr);
    return(result);
}

int errorf(const char *format, ...) {
    va_list argptr;
    int result;
    va_start(argptr, format);
    result = vfprintf(stderr, format, argptr);
    va_end(argptr);
    return(result);
}

float frand() {
    return (float)unif_rand();
}

int override_fprintf(FILE* stream, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    int res;
    if (stream == r_stderr) REvprintf(fmt, args);
    else if (stream == r_stdout) Rvprintf(fmt, args);
    else res = vfprintf(stream, fmt, args);
    va_end(args);
    return res;
}

int override_fflush(FILE* stream) {
    if (stream == r_stdout || stream == r_stderr) return 0;
    return fflush(stream);
}

void override_srand(unsigned int seed) {
    return;
}