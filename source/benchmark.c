#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include "benchmark_types.h"
#include <stdarg.h>

// https://stackoverflow.com/questions/26932749/how-can-wgnu-zero-variadic-macro-arguments-warning-be-turned-off-with-clang
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#define log_info(M, ...) fprintf(stderr, "[INFO] (%s:%d) " M "\n",__FILE__, __LINE__, ##__VA_ARGS__)
#pragma clang diagnostic pop

#ifdef __APPLE__
#include <mach/mach.h>
#elif __linux__
#include <sys/resource.h>
#endif

long now_in_seconds(void)
{
    struct timeval tval;
    gettimeofday( &tval, NULL );
    return tval.tv_sec;
}

int collect_memory_stats(struct memory_stats_t *memory_stats_out)
{   
#ifdef __APPLE__
    struct mach_task_basic_info t_info;
    mach_msg_type_number_t t_info_count = MACH_TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS == task_info(mach_task_self(),
                                  MACH_TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    {
        memory_stats_out->rss = (uint64_t)t_info.resident_size/1024/1024;
        memory_stats_out->max_rss = (uint64_t)t_info.resident_size_max/1024/1024;
    }
    else
    {
        log_info("task_info function has failed.");
        return BENCHMARK_MEMORY_ERROR;
    }
#elif __linux__
    struct rusage usage;
    int errno = getrusage(RUSAGE_SELF, &usage);
    if (errno == 0)
    {
        log_info("ru_maxrss --> %ld", usage.ru_maxrss);
        memory_stats_out->max_rss = usage.ru_maxrss;
    }
    else
    {
        log_info("getrusage function has failed.");
        return BENCHMARK_MEMORY_ERROR;
    }
#elif
    return BENCHMARK_UNIMPLEMENTED;
#endif
    return BENCHMARK_SUCCESS;
}

FILE *create_benchmark_file(char *template)
{
    char benchmark_out[420] = "\0";
	strncat(benchmark_out,template,420);
	char *extension = strrchr(benchmark_out, '.');
	char benchmark_name[20] = "\0";
	strcat(benchmark_name, "_benchmark");
	strncat(benchmark_name, extension, 10);
	*extension = '\0';
	strncat(benchmark_out, benchmark_name, 420);

    FILE *out = fopen(benchmark_out, "w");

    if(out)
    {
        fprintf(out, "{");
    }
    else
    {
        fprintf(stderr, "Benchmark file '%s' cannot be opened for write.\n", benchmark_out);
    }

    return out;
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
        fprintf(stderr, "Benchmark file is not valid.\n");
    }   
}

void report_elapsed_time(long elapsed_ms, FILE *benchmark_file)
{
    log_info("Elapsed time in s: %ld\n", elapsed_ms);
	if(benchmark_file)
	{
		fprintf(benchmark_file, ", \"elapsed_s\": %ld", elapsed_ms);
	}
    else
    {
        fprintf(stderr, "Benchmark file cannot be opened for write.\n");
    }
}

void report_memory_stats(struct memory_stats_t memory_stats, FILE *benchmark_file)
{
    log_info("max_rss = %lld MB\n", memory_stats.max_rss);
	if(benchmark_file)
	{
		fprintf(benchmark_file, "\"max_rss\": %lld", memory_stats.max_rss);
	}
    else
    {
        fprintf(stderr, "Benchmark file cannot be opened for write.\n");
    }
}
