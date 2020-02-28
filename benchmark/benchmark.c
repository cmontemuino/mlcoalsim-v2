#include <sys/time.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> // so that we can use free()
#include <string.h>
#include "zf_log.h"
#include "benchmark_types.h"

#ifdef __APPLE__
#include <mach/mach.h>
#elif __linux__
#include <sys/resource.h>
#endif

// Implementation note: we use a CLOCK_MONOTONIC although it is subject to NTP adjustments. Such adjustments will not
// cause the clock to jump.
// see https://www.systutorials.com/docs/linux/man/2-clock_gettime/
// see https://www.cs.rutgers.edu/~pxk/416/notes/c-tutorials/gettime.html
#define CLOCK_ID CLOCK_MONOTONIC

// We use this constant definition for byte -> megabyte transformations
#define MEGABYTE 1048576L

struct timespec stopwatch_start(void)
{
    struct timespec tp;
    clock_settime(CLOCK_ID, &tp);
    return tp;
}

long stopwatch_stop(struct timespec start)
{
    struct timespec end;
    uint64_t diff;

    clock_gettime(CLOCK_ID, &end);

    diff = end.tv_sec - start.tv_sec;
    return diff;
}

int collect_memory_stats(struct memory_stats_t *memory_stats_out)
{
#ifdef __APPLE__
    struct mach_task_basic_info t_info;
    mach_msg_type_number_t t_info_count = MACH_TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS == task_info(mach_task_self(),
                                  MACH_TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    {
        memory_stats_out->rss = t_info.resident_size/MEGABYTE;
        memory_stats_out->max_rss = t_info.resident_size_max/MEGABYTE;
    }
    else
    {
        ZF_LOGW("task_info function has failed. Memory cannot be collected.");
        return BENCHMARK_MEMORY_ERROR;
    }
#elif __linux__
    struct rusage usage;
    int errno = getrusage(RUSAGE_SELF, &usage);
    if (errno == 0)
    {
        memory_stats_out->max_rss = usage.ru_maxrss;
    }
    else
    {
        ZF_LOGW("getrusage function has failed. Memory cannot be collected.");
        return BENCHMARK_MEMORY_ERROR;
    }
#elif
    ZF_LOGW("OS not recognized. Memory cannot be collected");
    return BENCHMARK_UNIMPLEMENTED;
#endif
    return BENCHMARK_SUCCESS;
}

FILE *create_benchmark_file(char *template)
{
    // Note: using snprintf safely inspired from https://dzone.com/articles/detecting-string-truncation-with-gcc-8

    // First we make a copy of the template because we need to massage it and we don't want to produce a side-effect
    int bytes = snprintf(0, 0, "%s", template);
    // Implementation Note:
    // Somehow we need to copy one extra byte, otherwise the source template file name gets truncated.
    // For example, `somefile.out` gets truncated to `somefile.ou`. It could be related to how the template file
    // is managed by the original ms code.
    char *template_cpy = (char *) malloc (bytes + 2);
    snprintf(template_cpy, bytes + 1, "%s", template);

    // Usually all templates are `something.txt`, so we get `.txt` here
    char *extension = strrchr(template_cpy, '.');

    // Build our custom "extension", that is going to be something like `_benchmark.txt`
    bytes = snprintf(0, 0, "_benchmark.%s", extension);
    char *benchmark_extension = (char *) malloc (bytes + 1);
    snprintf(benchmark_extension, bytes, "_benchmark%s", extension);

    // Now we can build the file where benchmarks are going to be written
    bytes = snprintf (0, 0, "%s%s", template_cpy, benchmark_extension);
    char *benchmark_path = (char *) malloc (bytes + 1);
    *extension = '\0'; // so that we remove the extension from `template` string
    snprintf (benchmark_path, bytes, "%s%s", template_cpy, benchmark_extension);

    // Finally we create the file
    ZF_LOGD("Going to create the benchmark file |%s|", benchmark_path);
    FILE *benchmark_file = fopen(benchmark_path, "w");

    // Be a good citizen
    free(template_cpy);
    free(benchmark_extension);
    free(benchmark_path);

    if(benchmark_file)
    {
        fprintf(benchmark_file, "{");
    }
    else
    {
        ZF_LOGE("Benchmark file |%s| cannot be opened for write.", benchmark_path);
    }

    return benchmark_file;
}

void close_benchmark_file(FILE *benchmark_file)
{
    if(benchmark_file)
    {
        fprintf(benchmark_file, "}");
        fclose(benchmark_file);
    }
    else
    {
        ZF_LOGE("Benchmark file is not valid.");
    }
}

void report_elapsed_time(long elapsed_ms, FILE *benchmark_file)
{
    ZF_LOGD("Elapsed time: %ld s\n", elapsed_ms);
    if(benchmark_file)
    {
        fprintf(benchmark_file, ", \"elapsed_s\": %ld", elapsed_ms);
    }
    else
    {
        ZF_LOGE("Benchmark file cannot be opened for write.");
    }
}

void report_memory_stats(struct memory_stats_t memory_stats, FILE *benchmark_file)
{
#ifdef __APPLE__
    ZF_LOGD("max_rss = %lld MB", memory_stats.max_rss);
#elif __linux__
    ZF_LOGD("max_rss = %ld MB", memory_stats.max_rss);
#endif

    if(benchmark_file)
    {
#ifdef __APPLE__
        fprintf(benchmark_file, "\"max_rss\": %lld", memory_stats.max_rss);
#elif __linux__
        fprintf(benchmark_file, "\"max_rss\": %ld", memory_stats.max_rss);
#endif
    }
    else
    {
        fprintf(stderr, "Benchmark file cannot be opened for write.");
    }
}
