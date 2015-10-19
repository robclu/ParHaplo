// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo debug macros 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_DEBUG_H
#define PARHAPLO_DEBUG_H

#include <stdio.h>
#include <errno.h>
#include <string.h>

#ifdef RELEASE
    #define DEBUG(M, ...)
#else
    #define DEBUG(M, ...) printf("DEBUG %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#endif              // PARAHAPLO_DEBUG_H