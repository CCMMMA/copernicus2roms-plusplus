//
// Created by Ciro De Vita on 28/03/24.
//

#include "InterpolatorBase.h"

InterpolatorBase::InterpolatorBase(std::vector<double> srcLAT,
                                   std::vector<double> srcLON,
                                   std::vector<std::vector<double>> dstLAT,
                                   std::vector<std::vector<double>> dstLON,
                                   std::vector<std::vector<int>> dstMASK)
        : srcLAT(srcLAT), srcLON(srcLON), dstLAT(dstLAT), dstLON(dstLON), dstMASK(dstMASK),
          dstSNDim(dstLAT.size()), dstWEDim(dstLAT.empty() ? 0 : dstLAT[0].size()) {}

std::vector<std::vector<double>> InterpolatorBase::interp(std::vector<std::vector<double>> values, double srcMissingValue,
                                                          double dstMissingValue) {
    size_t dstEta = dstSNDim;
    size_t dstXi = dstWEDim;
    size_t srcEta = srcLAT.size();
    size_t srcXi = srcLON.size();

    double srcLonMin = srcLON[0];
    double srcLatMin = srcLAT[0];
    double srcLonMax = srcLON[srcXi - 1];
    double srcLatMax = srcLAT[srcEta - 1];
    double srcLatDelta = srcLatMax - srcLatMin;
    double srcLonDelta = srcLonMax - srcLonMin;
    double srcLatStep = srcEta == 0 ? 0 : srcLatDelta / srcEta;
    double srcLonStep = srcXi == 0 ? 0 : srcLonDelta / srcXi;

    std::vector<std::vector<double>> dst(dstEta, std::vector<double>(dstXi, dstMissingValue));
    for (size_t dstJ = 0; dstJ < dstEta; ++dstJ) {
        for (size_t dstI = 0; dstI < dstXi; ++dstI) {
            if (dstMASK[dstJ][dstI] == 1) {
                double dstLon = dstLON[dstJ][dstI];
                double dstLat = dstLAT[dstJ][dstI];

                double srcII = (dstLon - srcLonMin) / srcLonStep;
                double srcJJ = (dstLat - srcLatMin) / srcLatStep;

                int iR = static_cast<int>(srcII);
                int jR = static_cast<int>(srcJJ);

                std::vector<std::tuple<size_t, size_t, double>> pointsBilinear;
                for (int j : {jR - 1, jR + 1}) {
                    for (int i : {iR - 1, iR + 1}) {
                        size_t jj = std::min(std::max(j, 0), static_cast<int>(srcEta) - 1);
                        size_t ii = std::min(std::max(i, 0), static_cast<int>(srcXi) - 1);
                        // TODO: fix this workaround
                        if (std::fabs(values[jj][ii] - srcMissingValue) > 1e13) {
                            pointsBilinear.emplace_back(jj, ii, values[jj][ii]);
                        }
                    }
                }

                if (pointsBilinear.size() == 4) {
                    double lon1 = srcLON[std::get<1>(pointsBilinear[0])];
                    double lat1 = srcLAT[std::get<0>(pointsBilinear[0])];
                    double lon2 = srcLON[std::get<1>(pointsBilinear[1])];
                    double lat2 = srcLAT[std::get<0>(pointsBilinear[2])];

                    double dLon = lon2 - lon1;
                    double dLat = lat2 - lat1;

                    double FXY1 = ((lon2 - dstLon) / dLon) * std::get<2>(pointsBilinear[0]) +
                                  ((dstLon - lon1) / dLon) * std::get<2>(pointsBilinear[1]);
                    double FXY2 = ((lon2 - dstLon) / dLon) * std::get<2>(pointsBilinear[2]) +
                                  ((dstLon - lon1) / dLon) * std::get<2>(pointsBilinear[3]);

                    dst[dstJ][dstI] = ((lat2 - dstLat) / dLat) * FXY1 + ((dstLat - lat1) / dLat) * FXY2;
                } else if (USE_IDW) {
                    std::vector<std::tuple<size_t, size_t, double, double>> pointsIDW;
                    size_t size = 0;

                    while (pointsIDW.empty() && size < 4) {
                        size++;
                        for (int j = jR - size; j <= jR + size; ++j) {
                            size_t jj = std::min(std::max(j, 0), static_cast<int>(srcEta) - 1);
                            for (int i = iR - size; i <= iR + size; ++i) {
                                size_t ii = std::min(std::max(i, 0), static_cast<int>(srcXi) - 1);
                                // if (!std::isnan(values[jj][ii]) && values[jj][ii] != srcMissingValue) {
                                // TODO: fix this workaround
                                if (std::fabs(values[jj][ii] - srcMissingValue) > 1e13) {
                                    pointsIDW.emplace_back(jj,
                                                           ii,
                                                           haversine(dstLat, dstLon,
                                                                     srcLAT[jj], srcLON[ii]),
                                                                     values[jj][ii]);
                                }
                            }
                        }
                    }

                    if (!pointsIDW.empty()) {
                        double weighted_values_sum = 0.0;
                        double sum_of_weights = 0.0;
                        for (const auto& point : pointsIDW) {
                            double weight = 1 / std::get<2>(point);
                            sum_of_weights += weight;
                            weighted_values_sum += weight * std::get<3>(point);
                        }
                        dst[dstJ][dstI] = weighted_values_sum / sum_of_weights;
                    }
                }
            }
        }
    }

    return dst;
}

double InterpolatorBase::haversine(double lat1, double lon1, double lat2, double lon2) {
    // convert decimal degrees to radians
    lat1 = lat1 * M_PI / 180.0;
    lon1 = lon1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;
    lon2 = lon2 * M_PI / 180.0;

    // haversine formula
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;
    double a = std::sin(dlat / 2.0) * std::sin(dlat / 2.0) +
               std::cos(lat1) * std::cos(lat2) *
               std::sin(dlon / 2.0) * std::sin(dlon / 2.0);
    double c = 2 * std::asin(std::sqrt(a));
    double km = 6367 * c;
    return km;
}