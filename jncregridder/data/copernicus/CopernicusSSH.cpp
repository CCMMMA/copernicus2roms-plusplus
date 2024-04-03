//
// Created by Ciro De Vita on 28/03/24.
//

#include "CopernicusSSH.h"

CopernicusSSH::CopernicusSSH(const std::string& url) : CopernicusBase(url) {
    // Load ZOS variable
    std::vector<double> zos = loadVariable<double>("zos");

    // Get time, lat, and lon sizes
    size_t time_size = getTime();
    size_t lat_size = getLat();
    size_t lon_size = getLon();

    // Reshape ZOS into a 3D vector
    ZOS.resize(time_size);
    for (size_t t = 0; t < time_size; ++t) {
        ZOS[t].resize(lat_size);
        for (size_t i = 0; i < lat_size; ++i) {
            ZOS[t][i].resize(lon_size);
            for (size_t j = 0; j < lon_size; ++j) {
                size_t index = t * lat_size * lon_size + i * lon_size + j;
                ZOS[t][i][j] = zos[index];
            }
        }
    }
}

CopernicusSSH::~CopernicusSSH() {}

std::vector<std::vector<std::vector<double>>> CopernicusSSH::getZOS() const {
    return ZOS;
}