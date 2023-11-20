#ifndef __GET_TIME__
#define __GET_TIME__

#include "debug_print.h"
#if _WIN32
#include <Windows.h>
#endif
#else 
#include <sys/time.h>
#endif

#if DEBUG 

#if _WIN32

#define TINIT \
	LARGE_INTEGER cpu_freqence; \
	LARGE_INTEGER start; \
	LARGE_INTEGER end; \
	double run_time = 0.0; \
	QueryPerformanceFrequency(&cpu_freqence);

#define TIC \
  QueryPerformanceCounter(&start);

#define TOC(x) \
	QueryPerformanceCounter(&end); \
	run_time = (((end.QuadPart - start.QuadPart) * 1000.0f) / cpu_freqence.QuadPart); \
    Debug_info("%s taskes %f ms\r\n", x, run_time);

#elif __linux__
#define TINIT \
    struct timeval start, stop; \
    double elapsed_time;

#define TIC \
    gettimeofday(&start, NULL); 

#define TOC(x) \
    gettimeofday(&stop, NULL); \
    elapsed_time = (stop.tv_sec - start.tv_sec) * 1000. + \
    (stop.tv_usec - start.tv_usec) / 1000.; \
    Debug_info("%s taskes %f ms\r\n", x, elapsed_time);

#else

#define TINIT 
#define TIC 
#define TOC(x) 

#endif

#endif
