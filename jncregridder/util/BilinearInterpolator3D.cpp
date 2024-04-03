//
// Created by Ciro De Vita on 29/03/24.
//

#include "BilinearInterpolator3D.h"

BilinearInterpolator3D::BilinearInterpolator3D(std::vector<std::vector<double>> srcLAT,
                                               std::vector<std::vector<double>> srcLON,
                                               std::vector<std::vector<std::vector<double>>> srcZ,
                                               std::vector<std::vector<double>> dstLAT,
                                               std::vector<std::vector<double>> dstLON,
                                               std::vector<std::vector<std::vector<double>>> dstZ,
                                               std::vector<std::vector<int>> dstMASK,
                                               const ROMSGrid& romsGrid)
        : InterpolatorBase(srcLAT, srcLON, dstLAT, dstLON, dstMASK),
        bilinearInterpolator(srcLAT, srcLON, dstLAT, dstLON, dstMASK),
        srcZ(srcZ),
        dstZ(dstZ){
    srcLevs = srcZ.size();
    dstLevs = romsGrid.getS_rho().size();
    auto H = romsGrid.getH();

    for (int k = 0; k < srcLevs; ++k) {
        double copernicusDepth = srcZ[k][0][0];

        auto maxH = *std::max_element(H.begin(), H.end());
        if (copernicusDepth > maxH) {
            maxK = k + 1;
            break;
        }
    }

    srcLevs = maxK;
    // std::cout << maxK << std::endl;

    prepare();
}

void BilinearInterpolator3D::prepare() {
    for (int dstK = 0; dstK < dstLevs; dstK++) {
        for (int dstJ = 0; dstJ < dstSNDim; dstJ++) {
            for (int dstI = 0; dstI < dstWEDim; dstI++) {
                double dstZatKJI = std::abs(dstZ[dstK][dstJ][dstI]);

                int srcK = 0;
                double srcZat00 = 0.0;
                while (srcK < srcLevs) {
                    srcZat00 = std::abs(srcZ[srcK][0][0]);
                    if (srcZat00 > std::abs(dstZatKJI)) {
                        break;
                    }
                    srcK++;
                }

                int srcKmin = srcK - 1;
                int srcKmax = srcK;

                if (srcKmax == srcLevs) {
                    srcKmax = srcKmin;
                }

                if (srcKmin < 0) {
                    srcKmin = 0;
                }

                double srcZmin = std::abs(srcZ[srcKmin][0][0]);
                double srcZmax = std::abs(srcZ[srcKmax][0][0]);
                double delta = srcZmax - srcZmin;

                Weight3DStruct weight3D{};

                if (delta != 0.0) {
                    weight3D.w[0] = (dstZatKJI - srcZmin) / delta;
                    weight3D.w[1] = (srcZmax - dstZatKJI) / delta;
                } else {
                    weight3D.w[0] = 0.0;
                    weight3D.w[1] = 1.0;
                }

                weight3D.KK[0] = srcKmin;
                weight3D.KK[1] = srcKmax;

                weight3Ds.push_back(weight3D);
            }
        }
    }
}

std::vector<std::vector<std::vector<double>>> BilinearInterpolator3D::interp(std::vector<std::vector<std::vector<double>>> values, double srcMissingValue,
                                                                             double dstMissingValue) {
    std::vector<std::vector<std::vector<double>>> tSrc(maxK, std::vector<std::vector<double>>(dstSNDim, std::vector<double>(dstWEDim)));
    for (int k = 0; k < maxK; k++) {
        std::cout << "<k=" << k << " depth:" << srcZ[k][0][0] << " " << std::endl;
        tSrc[k] = bilinearInterpolator.interp(values[k], srcMissingValue, dstMissingValue);
    }

    std::cout << "Interpolating 3d..." << std::endl;

    std::vector<std::vector<std::vector<double>>> dst(dstLevs, std::vector<std::vector<double>>(dstSNDim, std::vector<double>(dstWEDim)));
    for (int dstK = 0; dstK < dstLevs; dstK++) {
        for (int dstJ = 0; dstJ < dstSNDim; dstJ++) {
            for (int dstI = 0; dstI < dstWEDim; dstI++) {
                if (dstMASK[dstJ][dstI] == 1) {
                    int srcI = dstI;
                    int srcJ = dstJ;

                    int srcK0 = weight3Ds[dstK * dstSNDim * dstWEDim + dstJ * dstWEDim + dstI].KK[0];
                    int srcK1 = weight3Ds[dstK * dstSNDim * dstWEDim + dstJ * dstWEDim + dstI].KK[1];
                    double srcW0 = weight3Ds[dstK * dstSNDim * dstWEDim + dstJ * dstWEDim + dstI].w[0];
                    double srcW1 = weight3Ds[dstK * dstSNDim * dstWEDim + dstJ * dstWEDim + dstI].w[1];
                    double tSrcK0 = tSrc[srcK0][srcJ][srcI];
                    double tSrcK1 = tSrc[srcK1][srcJ][srcI];

                    if ((tSrcK0 != dstMissingValue) && (tSrcK1 != dstMissingValue)) {
                        dst[dstK][dstJ][dstI] = (tSrcK0 * srcW0 + tSrcK1 * srcW1);
                    } else if (tSrcK0 != dstMissingValue) {
                        dst[dstK][dstJ][dstI] = tSrcK0;
                    } else if (tSrcK1 != dstMissingValue) {
                        dst[dstK][dstJ][dstI] = tSrcK1;
                    } else {
                        dst[dstK][dstJ][dstI] = dstMissingValue;
                    }
                } else {
                    dst[dstK][dstJ][dstI] = dstMissingValue;
                }
            }
        }
    }

    return dst;
}