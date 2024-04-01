//
// Created by Ciro De Vita on 28/03/24.
//

#include "BilinearInterpolator.h"

BilinearInterpolator::BilinearInterpolator(std::vector<double> srcLAT,
                                           std::vector<double> srcLON,
                                           std::vector<std::vector<double>> dstLAT,
                                           std::vector<std::vector<double>> dstLON,
                                           std::vector<std::vector<int>> dstMASK)
                                           : InterpolatorBase(srcLAT, srcLON, dstLAT, dstLON, dstMASK) {

}