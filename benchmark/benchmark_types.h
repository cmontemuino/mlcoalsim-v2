#ifndef BENCHMARK_TYPES_H_
#define BENCHMARK_TYPES_H_

#define BENCHMARK_SUCCESS       0
#define BENCHMARK_MEMORY_ERROR  10
#define BENCHMARK_UNIMPLEMENTED 30

#ifdef __APPLE__
struct memory_stats_t
{
    unsigned long long rss; /* Resident memory size (MB)*/
    unsigned long long max_rss; /* Maximum resident memory size (MB)*/
};
#elif __linux__
struct memory_stats_t
{
    long rss; /* Resident memory size (MB)*/
    long max_rss; /* Maximum resident memory size (MB)*/

};
#endif

#endif /* BENCHMARK_TYPES_H_ */
