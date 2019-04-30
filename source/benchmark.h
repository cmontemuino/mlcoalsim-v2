#ifndef BENCHMARK_H_
#define BENCHMARK_H_
#endif

long now_in_seconds();

FILE *create_benchmark_file(char *template);

void close_benchmark_file(FILE *benchmark_file);

void report_elapsed_time(long elapsed_ms, FILE *benchmark_file);
