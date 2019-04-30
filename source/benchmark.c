#include <sys/time.h>
#include <stdio.h>
#include <string.h>

long now_in_seconds()
{
    struct timeval tval;
    gettimeofday( &tval, NULL );
    return tval.tv_sec;
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
    fprintf(stdout, "Elapsed time in ms: %ld\n", elapsed_ms);
	if(benchmark_file)
	{
		fprintf(benchmark_file, "\"elapsed_ms\": %ld", elapsed_ms);
	}
    else
    {
        fprintf(stderr, "Benchmark file cannot be opened for write.\n");
    }   
}
