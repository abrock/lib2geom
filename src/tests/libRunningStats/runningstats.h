#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

class RunningStats {

public:

    template<class StatsVec>
    static std::vector<double> getMean(const StatsVec& vec);

    template<class StatsVec>
    static std::vector<double> getStddev(const StatsVec& vec);

    template<class StatsVec>
    static std::vector<double> getMin(const StatsVec& vec);

    template<class StatsVec>
    static std::vector<double> getMax(const StatsVec& vec);

    template<class StatsVec>
    static std::vector<size_t> getCount(const StatsVec& vec);

    template<class T, class StatsVec>
    static std::vector<T> getMedian(const StatsVec& vec);

    void clear();

    void push(const double value);

    double getMean() const;

    double getLogMean() const;

    double getVar() const;

    double getLogVar() const;

    double getStddev() const;

    double getLogStddev() const;

    double getMin() const;

    double getMax() const;

    size_t getCount() const;

    void print(std::ostream& out);

    std::string print();

    void printLog(std::ostream& out);

    std::string printLog();

    std::string printBoth();

    double sum = 0;
    double squaresum = 0;
    double min = 0;
    double max = 0;
    size_t n = 0;

    bool calcLog = true;

    double log_sum = 0;
    double log_square_sum = 0;
    double log_min = 0;
    double log_max = 0;
    size_t log_n = 0;
};

template<class T>
class QuantileStats : public RunningStats {
public:
    void push(const double value);

    T getQuantile(const double quantile) const;

    T getMedian() const;

    double getAccurateVariance() const;

    double getAccurateStddev() const;

    void reserve(const size_t size);

    std::vector<T> getData();

private:
    void sort() const;

    mutable std::vector<T> values;
    mutable bool sorted = true;
};

#endif // RUNNINGSTATS_H
