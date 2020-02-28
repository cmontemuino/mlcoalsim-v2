#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include <stdio.h> // needed for "FILE"
#include "benchmark_types.h"

// Starts a stopwatch and returns the starting time.
struct timespec stopwatch_start(void);

// Stops the stopwatch and returns the time difference in seconds given the starting time.
long stopwatch_stop(struct timespec start);

// Add a memory stats metric point to the head of a list.
// This works like a stack. The passed `head` will point ot the newly added metric point.
// It returns either BENCHMARK_SUCCESS or BENCHMARK_MEMORY_ERROR in case of problems collecting the memory metric.
int *add_memory_metric_point(memory_stats_set_t **head);

FILE *create_benchmark_file(char *template);

void close_benchmark_file(FILE *benchmark_file);

void report_elapsed_time(long elapsed_ms, FILE *benchmark_file);

void report_memory_stats(memory_stats_set_t *head, FILE *benchmark_file);

#endif /* BENCHMARK_H_ */
