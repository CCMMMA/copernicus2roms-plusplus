//
// Created by Ciro De Vita on 28/03/24.
//

#include "InterpolatorBase.h"

InterpolatorBase::InterpolatorBase(std::vector<std::vector<double>> srcLAT,
                                   std::vector<std::vector<double>> srcLON,
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
    size_t srcXi = srcLAT[0].size();

    double srcLonMin = srcLON[0][0];
    double srcLonMax = srcLON[srcEta - 1][srcXi - 1];
    double srcLatMin = srcLAT[0][0];
    double srcLatMax = srcLAT[srcEta - 1][srcXi - 1];
    double srcLatDelta = srcLatMax - srcLatMin;
    double srcLonDelta = srcLonMax - srcLonMin;
    double srcLatStep = srcEta == 0 ? 0 : srcLatDelta / srcEta;
    double srcLonStep = srcXi == 0 ? 0 : srcLonDelta / srcXi;

    std::vector<std::vector<double>> dst(dstEta, std::vector<double>(dstXi, dstMissingValue));

    #pragma omp parallel for collapse(2) default(none) shared(dstEta, dstXi, srcLonMin, srcLatMin, srcLonStep, srcLatStep, srcEta, srcXi, values, srcMissingValue, dst)
    for (size_t dstJ = 0; dstJ < dstEta; ++dstJ) {
        for (size_t dstI = 0; dstI < dstXi; ++dstI) {
            if (dstMASK[dstJ][dstI] == 1) {
                double dstLon = dstLON[dstJ][dstI];
                double dstLat = dstLAT[dstJ][dstI];

                double srcII = (dstLon - srcLonMin) / srcLonStep;
                double srcJJ = (dstLat - srcLatMin) / srcLatStep;

                int iR = static_cast<int>(srcII);
                int jR = static_cast<int>(srcJJ);

                std::vector<Point> pointsBilinear;

                /*
                for (int j = jR; j <= jR + 1; j++) {
                    size_t jj = j >= values.size() ? values.size() - 1 : j;
                    for (int i = iR; i <= iR + 1; i++) {
                        size_t ii = i >= values[0].size() ? values[0].size() - 1 : i;
                        if (std::fabs(values[jj][ii] - srcMissingValue) > 1e13) {
                            pointsBilinear.emplace_back(Point(jj, ii, values[jj][ii]));
                        }
                    }
                }
                */

                for (int j : {jR - 1, jR + 1}) {
                    for (int i : {iR - 1, iR + 1}) {
                        size_t jj = std::min(std::max(j, 0), static_cast<int>(srcEta) - 1);
                        size_t ii = std::min(std::max(i, 0), static_cast<int>(srcXi) - 1);
                        // TODO: fix this workaround
                        if (std::fabs(values[jj][ii] - srcMissingValue) > 1e13) {
                            pointsBilinear.emplace_back(Point(jj, ii, values[jj][ii]));
                        }
                    }
                }

                if (pointsBilinear.size() == 4) {
                    double lon1 = srcLON[pointsBilinear[0].j][pointsBilinear[0].i];
                    double lat1 = srcLAT[pointsBilinear[0].j][pointsBilinear[0].i];
                    double lon2 = srcLON[pointsBilinear[1].j][pointsBilinear[1].i];
                    double lat2 = srcLAT[pointsBilinear[2].j][pointsBilinear[2].i];

                    double dLon = lon2 - lon1;
                    double dLat = lat2 - lat1;

                    double FXY1 = ((lon2 - dstLon) / dLon) * pointsBilinear[0].value + ((dstLon - lon1) / dLon) * pointsBilinear[1].value;
                    double FXY2 = ((lon2 - dstLon) / dLon) * pointsBilinear[2].value + ((dstLon - lon1) / dLon) * pointsBilinear[3].value;

                    dst[dstJ][dstI] = ((lat2 - dstLat) / dLat) * FXY1 + ((dstLat - lat1) / dLat) * FXY2;
                } else if (USE_IDW) {
                    std::vector<Point> pointsIDW;
                    int size = 0;

                    do {
                        size++;
                        for (int j = jR - size; j <= jR + size; j++) {
                            size_t jj = j >= values.size() ? values.size() - 1 : (j < 0 ? 0 : j);
                            for (int i = iR - size; i <= iR + size; i++) {
                                size_t ii = i >= values[0].size() ? values[0].size() - 1 : (i < 0 ? 0 : i);
                                // TODO: fix this workaround
                                if (std::fabs(values[jj][ii] - srcMissingValue) > 1e13) {
                                    pointsIDW.emplace_back(Point(jj, ii, distance(dstLAT[dstJ][dstI], dstLON[dstJ][dstI], srcLAT[jj][ii], srcLON[jj][ii]), values[jj][ii]));
                                }
                            }
                        }
                    } while (pointsIDW.size() == 0 && size < 4);

                    if (pointsIDW.size() > 0) {
                        double weighted_values_sum = 0.0;
                        double sum_of_weights = 0.0;
                        double weight;
                        for (const auto &point: pointsIDW) {
                            weight = 1 / point.dist;
                            sum_of_weights += weight;
                            weighted_values_sum += weight * point.value;
                        }
                        dst[dstJ][dstI] = weighted_values_sum / sum_of_weights;
                    }

                    /*
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
                                                           distance(dstLat, dstLon,srcLAT[jj][ii],srcLON[jj][ii]),
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
                     */
                }
            }
        }
    }

    return dst;
}

double InterpolatorBase::distance(double startLat, double startLong,
        double endLat, double endLong) {

    double dLat  = toRadians((endLat - startLat));
    double dLong = toRadians((endLong - startLong));

    startLat = toRadians(startLat);
    endLat   = toRadians(endLat);

    double a = haversin(dLat) + cos(startLat) * cos(endLat) * haversin(dLong);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    return EARTH_RADIUS * c;
}

double InterpolatorBase::toRadians(double degrees) {
    return degrees * M_PI / 180.0;
}

double InterpolatorBase::haversin(double val) {
    return pow(sin(val / 2), 2);
}