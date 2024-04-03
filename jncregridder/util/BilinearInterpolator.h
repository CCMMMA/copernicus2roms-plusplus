//
// Created by Ciro De Vita on 28/03/24.
//

#ifndef MYOCEAN2ROMS_BILINEARINTERPOLATOR_H
#define MYOCEAN2ROMS_BILINEARINTERPOLATOR_H

#include "InterpolatorBase.h"

class BilinearInterpolator : public InterpolatorBase {
public:
    BilinearInterpolator(std::vector<std::vector<double>> srcLAT,
                         std::vector<std::vector<double>> srcLON,
                         std::vector<std::vector<double>> dstLAT,
                         std::vector<std::vector<double>> dstLON,
                         std::vector<std::vector<int>> dstMASK);
};


#endif //MYOCEAN2ROMS_BILINEARINTERPOLATOR_H
