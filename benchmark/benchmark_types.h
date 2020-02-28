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

// We store a set of memory_stats value in this structure
typedef struct memory_stats_set {
    struct memory_stats_t val;
    struct memory_stats_set *next;
} memory_stats_set_t;

#elif __linux__
struct memory_stats_t
{
    long rss; /* Resident memory size (MB)*/
    long max_rss; /* Maximum resident memory size (MB)*/

};
#endif

#endif /* BENCHMARK_TYPES_H_ */
