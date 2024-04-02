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

class InterpolatorBase {
public:
    InterpolatorBase(std::vector<double> srcLAT,
                     std::vector<double> srcLON,
                     std::vector<std::vector<double>> dstLAT,
                     std::vector<std::vector<double>> dstLON,
                     std::vector<std::vector<int>> dstMASK);
    std::vector<std::vector<double>> interp(std::vector<std::vector<double>> values, double srcMissingValue,
                                            double dstMissingValue);

private:
    std::vector<double> srcLAT;
    std::vector<double> srcLON;
    std::vector<std::vector<double>> dstLAT;
    std::vector<std::vector<double>> dstLON;
    bool USE_IDW = true;

    double haversine(double, double, double, double);

protected:
    int dstSNDim;
    int dstWEDim;
    std::vector<std::vector<int>> dstMASK;
};


#endif //MYOCEAN2ROMS_INTERPOLATORBASE_H
