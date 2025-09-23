#ifndef FUNCTION_WRAPPERS_H
#define FUNCTION_WRAPPERS_H

#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rmath.h>

#undef exit
#define exit(status) Rf_error("Library attempted to call exit(%d)", status)

#undef rand
#define rand() ((int)(unif_rand() * RAND_MAX))

#undef srand
#define srand override_srand

#undef printf
#define printf(...) Rprintf(__VA_ARGS__)

#undef fprintf
#define fprintf override_fprintf

int override_fprintf(FILE* stream, const char* fmt, ...);

#undef stderr
#define stderr r_stderr

#undef stdout
#define stdout r_stdout

extern FILE* r_stderr;
extern FILE* r_stdout;

#undef fflush
#define fflush override_fflush

int override_fflush(FILE* stream);


int msgf(const char*, ...);
int errorf(const char*, ...);
float frand();
void override_srand(unsigned int seed);

#endif