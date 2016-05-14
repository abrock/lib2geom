#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

class RunningStats {

public:
    void clear() {
        sum = 0;
        squaresum = 0;
        n = 0;
    }
    template<class T>
    void push(T value) {
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
            double logVal = std::log(value) / std::log(10);
            logN++;
            logSum += logVal;
            logSquareSum += logVal * logVal;
        }
    }

    double getMean() {
        if (n < 1) {
            return 0;
        }
        return sum / n;
    }
    double getLogMean() {
        if (logN < 1) {
            return 0;
        }
        return logSum / logN;
    }
    double getVar() {
        if (n < 2) {
            return 0;
        }
        return 1.0/(n-1) * (squaresum - sum*sum / n);
    }
    double getLogVar() {
        if (logN < 2) {
            return 0;
        }
        return 1.0/(logN-1) * (logSquareSum - logSum*logSum / logN);
    }
    double getStddev() {
        return std::sqrt(getVar());
    }
    double getLogStddev() {
        return std::sqrt(getLogVar());
    }
    void print(std::ostream& out) {
        out << getMean() << " +- " << getStddev() << ", " << n << " Samples, range: [" << min << ", " << max << "]";
    }

    std::string print() {
        std::stringstream out;
        print(out);
        return out.str();
    }

    void printLog(std::ostream& out) {
        out << getLogMean() << " +- " << getLogStddev() << ", " << n << " Samples";
    }

    std::string printLog() {
        std::stringstream out;
        printLog(out);
        return out.str();
    }

    std::string printBoth() {
        return print() + "\nLogarithmic: " + printLog();
    }

    double sum = 0;
    double squaresum = 0;
    double min = 0;
    double max = 0;
    size_t n = 0;

    bool calcLog = true;

    double logSum = 0;
    double logSquareSum = 0;
    size_t logN = 0;
};

class QuantileStats : public RunningStats {
public:
    template<class T>
    void push(T value){
        sorted = false;
        values.push_back(value);
        RunningStats::push(value);
    }

    template<class T>
    float getQuantile(const T & quantile) {
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
        return values[(size_t)(quantile * (values.size()-1))];
    }

private:
    void sort() {
        if (!sorted) {
            std::sort(values.begin(), values.end());
            sorted = true;
        }
    }

    std::vector<float> values;
    bool sorted = true;
};

#endif // RUNNINGSTATS_H
