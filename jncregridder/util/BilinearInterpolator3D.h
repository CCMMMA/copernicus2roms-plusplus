//
// Created by Ciro De Vita on 29/03/24.
//

#ifndef MYOCEAN2ROMS_BILINEARINTERPOLATOR3D_H
#define MYOCEAN2ROMS_BILINEARINTERPOLATOR3D_H

#include "InterpolatorBase.h"
#include "../roms/ROMSGrid.h"
#include "BilinearInterpolator.h"

struct Weight3DStruct {
    double w[2];
    int KK[2];
};

class BilinearInterpolator3D : public InterpolatorBase {
public:
    BilinearInterpolator3D(std::vector<std::vector<double>> srcLAT,
                           std::vector<std::vector<double>> srcLON,
                           std::vector<std::vector<std::vector<double>>> srcZ,
                           std::vector<std::vector<double>> dstLAT,
                           std::vector<std::vector<double>> dstLON,
                           std::vector<std::vector<std::vector<double>>> dstZ,
                           std::vector<std::vector<int>> dstMASK,
                           const ROMSGrid& romsGrid);

    std::vector<std::vector<std::vector<double>>> interp(std::vector<std::vector<std::vector<double>>>, double, double);

private:
    size_t srcLevs;
    size_t dstLevs;
    int maxK = 0;
    std::vector<std::vector<std::vector<double>>> srcZ, dstZ;
    std::vector<Weight3DStruct> weight3Ds;
    BilinearInterpolator bilinearInterpolator;

    void prepare();
};


#endif //MYOCEAN2ROMS_BILINEARINTERPOLATOR3D_H
