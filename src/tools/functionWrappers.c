#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "functionWrappers.h"

FILE* r_stderr = NULL;
FILE* r_stdout = NULL;

/* Helper buffer size -- keep reasonably large but finite */
#define FW_BUF_SIZE 8192

int msgf(const char *format, ...) {
    va_list args;
    char buf[FW_BUF_SIZE];
    int n;

    va_start(args, format);
    n = vsnprintf(buf, sizeof(buf), format, args);
    va_end(args);

    /* Ensure null-termination (vsnprintf does that when size>0) */
    if (n < 0) return n;
    Rprintf("%s", buf);
    return n;
}

int errorf(const char *format, ...) {
    va_list args;
    char buf[FW_BUF_SIZE];
    int n;

    va_start(args, format);
    n = vsnprintf(buf, sizeof(buf), format, args);
    va_end(args);

    if (n < 0) return n;
    REprintf("%s", buf);
    return n;
}

float frand(void) {
    return (float)unif_rand();
}

int override_fprintf(FILE* stream, const char* fmt, ...) {
    va_list args;
    int res = 0;

    va_start(args, fmt);

    if (stream == r_stderr) {
        char buf[FW_BUF_SIZE];
        int n = vsnprintf(buf, sizeof(buf), fmt, args);
        if (n < 0) {
            res = n;
        } else {
            REprintf("%s", buf);
            /* return number of characters that would have been written */
            res = n;
        }
    } else if (stream == r_stdout) {
        char buf[FW_BUF_SIZE];
        int n = vsnprintf(buf, sizeof(buf), fmt, args);
        if (n < 0) {
            res = n;
        } else {
            Rprintf("%s", buf);
            res = n;
        }
    } else {
        /* other FILE* streams: pass through to vfprintf */
        res = vfprintf(stream, fmt, args);
    }

    va_end(args);
    return res;
}

int override_fflush(FILE* stream) {
    if (stream == r_stdout || stream == r_stderr) return 0;
    return fflush(stream);
}

void override_srand(unsigned int seed) {
    /* ignore seeding of C RNG; R's RNG should be used via set_seed/unif_rand */
    (void) seed;
    return;
}
