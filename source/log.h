/*
@(#)File:           log.h
@(#)Purpose:        Definitions for logging primitives
@(#)Author:         Carlos Montemuiño
*/

#ifndef LOG_H
#define LOG_H

#include <stdio.h>

// Determine whether logging macros are active at compile time.
#ifdef LOG
// https://stackoverflow.com/questions/26932749/how-can-wgnu-zero-variadic-macro-arguments-warning-be-turned-off-with-clang
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#define log_info(M, ...) fprintf(stderr, "[INFO] (%s:%d) " M "\n",__FILE__, __LINE__, ##__VA_ARGS__)
#pragma clang diagnostic pop
#else
#define log_info(M, ...)
#endif // LOG

#endif // LOG_H
