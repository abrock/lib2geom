#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "runningstats.h"

void RunningStats::clear() {
    sum = 0;
    squaresum = 0;
    n = 0;
}

void RunningStats::push(const double value) {
    if (!std::isfinite(value)) {
        return;
    }
#pragma omp critical
    {
        if (n == 0) {
            min = value;
            max = value;
        }
        else {
            min = std::min(min, (double)value);
            max = std::max(max, (double)value);
        }
        sum += value;
        squaresum += value * value;
        n++;

        if (calcLog && value > 0) {
            double log_val = std::log(static_cast<double>(value)) / std::log(10);
            log_sum += log_val;
            log_square_sum += log_val * log_val;
            if (log_n == 0) {
                log_min = log_val;
                log_max = log_val;
            }
            else {
                log_min = std::min(log_min, log_val);
                log_max = std::max(log_max, log_val);
            }
            log_n++;
        }
    }
}

double RunningStats::getMean() const {
    if (n < 1) {
        return 0;
    }
    return sum / n;
}
double RunningStats::getLogMean() const {
    if (log_n < 1) {
        return 0;
    }
    return log_sum / log_n;
}
double RunningStats::getVar() const {
    if (n < 2) {
        return 0;
    }
    return 1.0/(n-1) * (squaresum - sum*sum / n);
}
double RunningStats::getLogVar() const {
    if (log_n < 2) {
        return 0;
    }
    return 1.0/(log_n-1) * (log_square_sum - log_sum*log_sum / log_n);
}
double RunningStats::getStddev() const {
    return std::sqrt(getVar());
}
double RunningStats::getLogStddev() const {
    return std::sqrt(getLogVar());
}
void RunningStats::print(std::ostream& out) {
    out << getMean() << " +- " << getStddev() << ", " << n << " Samples, range: [" << min << ", " << max << "]";
}

std::string RunningStats::print() {
    std::stringstream out;
    print(out);
    return out.str();
}

void RunningStats::printLog(std::ostream& out) {
    out << getLogMean() << " +- " << getLogStddev() << ", " << n << " Samples, range: [" << log_min << ", " << log_max << "]";
}

std::string RunningStats::printLog() {
    std::stringstream out;
    printLog(out);
    return out.str();
}

std::string RunningStats::printBoth() {
    return print() + "\nLogarithmic: " + printLog();
}

template<class StatsVec>
std::vector<double> RunningStats::getMean(const StatsVec& vec) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getMean());
    }
    return result;
}

template<class StatsVec>
std::vector<double> RunningStats::getStddev(const StatsVec& vec) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getStddev());
    }
    return result;
}

template<class StatsVec>
std::vector<double> RunningStats::getMin(const StatsVec& vec) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getMin());
    }
    return result;
}

template<class StatsVec>
std::vector<double> RunningStats::getMax(const StatsVec& vec) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getMax());
    }
    return result;
}

template<class StatsVec>
std::vector<size_t> RunningStats::getCount(const StatsVec& vec) {
    std::vector<size_t> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getCount());
    }
    return result;
}

template<class T, class StatsVec>
std::vector<T> RunningStats::getMedian(const StatsVec& vec) {
    std::vector<T> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getMedian());
    }
    return result;
}

template std::vector<double> RunningStats::getStddev(const std::vector<RunningStats>& vec);
template std::vector<double> RunningStats::getMean(const std::vector<RunningStats>& vec);
template std::vector<double> RunningStats::getMin(const std::vector<RunningStats>& vec);
template std::vector<double> RunningStats::getMax(const std::vector<RunningStats>& vec);
template std::vector<size_t> RunningStats::getCount(const std::vector<RunningStats>& vec);

template std::vector<double> RunningStats::getStddev(const std::vector<QuantileStats<double> >& vec);
template std::vector<double> RunningStats::getMean(const std::vector<QuantileStats<double> >& vec);
template std::vector<double> RunningStats::getMin(const std::vector<QuantileStats<double> >& vec);
template std::vector<double> RunningStats::getMax(const std::vector<QuantileStats<double> >& vec);
template std::vector<double> RunningStats::getMedian(const std::vector<QuantileStats<double> >& vec);
template std::vector<size_t> RunningStats::getCount(const std::vector<QuantileStats<double> >& vec);

size_t RunningStats::getCount() const {
    return n;
}

double RunningStats::getMin() const {
    return min;
}

double RunningStats::getMax() const {
    return max;
}

template<class T>
T QuantileStats<T>::getMedian() const {
    return getQuantile(0.5);
}

template<class T>
void QuantileStats<T>::push(const double value){
    sorted = false;
    values.push_back(value);
    RunningStats::push(value);
}


template<class T>
T QuantileStats<T>::getQuantile(const double quantile) const {
    if (quantile <= 0) {
        return min;
    }
    if (quantile >= 1) {
        return max;
    }
    if (values.size() == 0) {
        return 0;
    }
    if (values.size() == 1) {
        return values[0];
    }
    sort();
    return values[static_cast<size_t>(quantile * (values.size()-1))];
}

template<class T>
void QuantileStats<T>::sort() const {
    if (!sorted) {
        std::sort(values.begin(), values.end());
        sorted = true;
    }
}

template<class T>
double QuantileStats<T>::getAccurateVariance() const {
    if (n < 2) {
        return 0;
    }
    double square_sum = 0;
    const double mean = getMean();
    for (const T val : values) {
        square_sum += (val - mean) * (val - mean);
    }
    return square_sum / (n-1);
}

template<class T>
double QuantileStats<T>::getAccurateStddev() const {
    return std::sqrt(getAccurateVariance());
}

template<class T>
void QuantileStats<T>::reserve(const size_t size) {
    values.reserve(size);
}

template<class T>
std::vector<T> QuantileStats<T>::getData() {
    sort();
    return values;
}


template class QuantileStats<double>;
template class QuantileStats<float>;
