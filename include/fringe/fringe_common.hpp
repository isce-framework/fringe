#ifndef FRINGE_COMMON_H
#define FRINGE_COMMON_H

#include <time.h>
#include <sys/time.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double getWallTime()
{
    struct timeval time;
    if (gettimeofday(&time,NULL))
    {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int numberOfThreads()
{
    int nthreads = 0;

#ifdef _OPENMP
    #pragma omp parallel
    if (omp_get_thread_num() == 1)
    {
        nthreads = omp_get_num_threads();
    }
#endif

    std::cout << "Processing with " << nthreads << " threads \n";

    if (nthreads == 0)
    {
        std::cout << "Looks like the code was not linked with openmp. \n";
        std::cout << "Recompile with the right linker flags. \n";
        throw;
    }
    
    if (nthreads == 1)
    {
        std::cout << "This code has been designed to work with multiple threads. \n";
        std::cout << "Looks like sequential execution due to compilation or environment settings. \n";
        std::cout << "Check your settings for optimal performace. Continuing ... \n";
    }

    return nthreads;
}

double daysSince1900(std::string indate)
{
    tm tm0 = {0};
    tm0.tm_year = 0;
    tm0.tm_mon = 0;
    tm0.tm_mday = 1;
    time_t time0 = mktime(&tm0);

    tm tm1 = {0};
    tm1.tm_year = stoi( indate.substr(0,4)) - 1900;
    tm1.tm_mon = stoi( indate.substr(4,2)) - 1;
    tm1.tm_mday = stoi( indate.substr(6,2));

    time_t time1 = mktime(&tm1);

    return difftime(time1, time0) / (24.0*60.0*60.0);
}

#endif //FRINGE_COMMON_H
