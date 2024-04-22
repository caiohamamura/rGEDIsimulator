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
    float max=0;
    if(RAND_MAX>0)max=(float)RAND_MAX;
    else          max=-1.0*(float)RAND_MAX;
    return((float)rand()/max);
}

int rand2() {
    return rand();
}

void srand2(int seed) {
    srand(seed);
}
