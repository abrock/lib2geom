#ifndef PARALLELTIME_H
#define PARALLELTIME_H

#include <chrono>
#include <string>

class ParallelTime
{
    std::chrono::high_resolution_clock::time_point start_hr, stop_hr;
    clock_t start_clock, stop_clock;
    bool isStopped = false;
public:
    ParallelTime();
    void start();
    void stop();
    std::string print();
    double realTime();
    double cpuTime();
    double parallelisation();
};

#endif // PARALLELTIME_H
