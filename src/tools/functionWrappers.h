#ifndef HANCOCK_MSGHANDLING_H
#define HANCOCK_MSGHANDLING_H

#include <stdio.h>
#include <stdarg.h>

int msgf(const char*, ...);
int errorf(const char*, ...);
float frand();
int rand2();
void srand2(int);

#endif /* HANCOCK_MSGHANDLING_H */
