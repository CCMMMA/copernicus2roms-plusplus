//
// Created by Ciro De Vita on 28/03/24.
//

#ifndef MYOCEAN2ROMS_INTERPOLATORBASE_H
#define MYOCEAN2ROMS_INTERPOLATORBASE_H

#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>

class Point {
public:
    size_t j;
    size_t i;
    double dist;
    double value;
    double w;

    Point(size_t j_, size_t i_, double dist_, double value_) : j(j_), i(i_), dist(dist_), value(value_) {}
    Point(size_t j_, size_t i_, double value_) : j(j_), i(i_), dist(-1), value(value_) {}
};

class InterpolatorBase {
public:
    InterpolatorBase(std::vector<std::vector<double>> srcLAT,
                     std::vector<std::vector<double>> srcLON,
                     std::vector<std::vector<double>> dstLAT,
                     std::vector<std::vector<double>> dstLON,
                     std::vector<std::vector<int>> dstMASK);
    std::vector<std::vector<double>> interp(std::vector<std::vector<double>> values, double srcMissingValue,
                                            double dstMissingValue);

private:
    std::vector<std::vector<double>> srcLAT;
    std::vector<std::vector<double>> srcLON;
    std::vector<std::vector<double>> dstLAT;
    std::vector<std::vector<double>> dstLON;
    bool USE_IDW = true;
    const int EARTH_RADIUS = 6371;

    double distance(double, double, double, double);
    double toRadians(double);
    double haversin(double);

protected:
    int dstSNDim;
    int dstWEDim;
    std::vector<std::vector<int>> dstMASK;
};


#endif //MYOCEAN2ROMS_INTERPOLATORBASE_H
