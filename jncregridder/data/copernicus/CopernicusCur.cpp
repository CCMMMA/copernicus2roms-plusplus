//
// Created by Ciro De Vita on 28/03/24.
//

#include "CopernicusCur.h"

CopernicusCur::CopernicusCur(const std::string &url) : CopernicusBase(url) {
    // Load UO variable
    std::vector<double> vo = loadVariable("vo");
    // Load VO variable
    std::vector<double> uo = loadVariable("uo");

    // Get time, depth, lat, and lon sizes
    size_t time_size = getTime();
    size_t depth_size = getDepth();
    size_t lat_size = getLat();
    size_t lon_size = getLon();

    // Reshape UO into a 4D vector
    UO.resize(time_size);
    for (size_t t = 0; t < time_size; ++t) {
        UO[t].resize(depth_size);
        for (size_t d = 0; d < depth_size; ++d) {
            UO[t][d].resize(lat_size);
            for (size_t i = 0; i < lat_size; ++i) {
                UO[t][d][i].resize(lon_size);
                for (size_t j = 0; j < lon_size; ++j) {
                    size_t index = (t * depth_size * lat_size * lon_size) + (d * lat_size * lon_size) + (i * lon_size) + j;
                    UO[t][d][i][j] = uo[index];
                }
            }
        }
    }

    // Reshape UO into a 4D vector
    VO.resize(time_size);
    for (size_t t = 0; t < time_size; ++t) {
        VO[t].resize(depth_size);
        for (size_t d = 0; d < depth_size; ++d) {
            VO[t][d].resize(lat_size);
            for (size_t i = 0; i < lat_size; ++i) {
                VO[t][d][i].resize(lon_size);
                for (size_t j = 0; j < lon_size; ++j) {
                    size_t index = (t * depth_size * lat_size * lon_size) + (d * lat_size * lon_size) + (i * lon_size) + j;
                    VO[t][d][i][j] = uo[index];
                }
            }
        }
    }
}

CopernicusCur::~CopernicusCur() {}

std::vector<std::vector<std::vector<std::vector<double>>>>CopernicusCur::getVO() const {
    return VO;
}

std::vector<std::vector<std::vector<std::vector<double>>>> CopernicusCur::getUO() const {
    return UO;
}