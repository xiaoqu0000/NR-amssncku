#ifndef AHFINDERDIRECT__CONFIG_H
#define AHFINDERDIRECT__CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

size_t Util_Strlcat(char* dst, const char* src, size_t dst_size);
size_t Util_Strlcpy(char* dst, const char* src, size_t dst_size);

typedef CCTK_REAL fp;

typedef CCTK_INT integer;

#endif	/* AHFINDERDIRECT__CONFIG_H */
