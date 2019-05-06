#ifndef BENCHMARK_H_
#define BENCHMARK_H_
#endif

#include "benchmark_types.h"

long now_in_seconds(void);

int collect_memory_stats(struct memory_stats_t *memory_stats_out);

FILE *create_benchmark_file(char *template);

void close_benchmark_file(FILE *benchmark_file);

void report_elapsed_time(long elapsed_ms, FILE *benchmark_file);

void report_memory_stats(struct memory_stats_t memory_stats, FILE *benchmark_file);
