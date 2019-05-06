#ifndef BENCHMARK_TYPES_H_
#define BENCHMARK_TYPES_H_
#endif

#define BENCHMARK_SUCCESS       0
#define BENCHMARK_MEMORY_ERROR  10
#define BENCHMARK_UNIMPLEMENTED 30

struct memory_stats_t
{
    long long rss; /* Resident memory size (MB)*/
    long long max_rss; /* Maximum resident memory size (MB)*/ 
};
