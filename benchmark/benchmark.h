#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include <stdio.h> // needed for "FILE"
#include "benchmark_types.h"

// Starts a stopwatch and returns the starting time.
struct timespec stopwatch_start(void);

// Stops the stopwatch and returns the time difference in seconds given the starting time.
long stopwatch_stop(struct timespec start);

int collect_memory_stats(struct memory_stats_t *memory_stats_out);

FILE *create_benchmark_file(char *template);

void close_benchmark_file(FILE *benchmark_file);

void report_elapsed_time(long elapsed_ms, FILE *benchmark_file);

void report_memory_stats(struct memory_stats_t memory_stats, FILE *benchmark_file);

#endif /* BENCHMARK_H_ */
