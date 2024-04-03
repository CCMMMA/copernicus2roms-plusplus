//
// Created by Ciro De Vita on 28/03/24.
//

#include "CopernicusTem.h"

CopernicusTem::CopernicusTem(const std::string& url) : CopernicusBase(url) {
    // Load THETAO variable
    std::vector<double> thetao = loadVariable<double>("thetao");
    // Load BOTTOMT variable
    std::vector<double> bottomT = loadVariable<double>("bottomT");

    // Get time, depth, lat, and lon sizes
    size_t time_size = getTime();
    size_t depth_size = getDepth();
    size_t lat_size = getLat();
    size_t lon_size = getLon();

    // Reshape THETAO into a 4D vector
    THETAO.resize(time_size);
    for (size_t t = 0; t < time_size; ++t) {
        THETAO[t].resize(depth_size);
        for (size_t d = 0; d < depth_size; ++d) {
            THETAO[t][d].resize(lat_size);
            for (size_t i = 0; i < lat_size; ++i) {
                THETAO[t][d][i].resize(lon_size);
                for (size_t j = 0; j < lon_size; ++j) {
                    size_t index = (t * depth_size * lat_size * lon_size) + (d * lat_size * lon_size) + (i * lon_size) + j;
                    THETAO[t][d][i][j] = thetao[index];
                }
            }
        }
    }

    // Reshape BOTTOMT into a 3D vector
    BOTTOMT.resize(time_size);
    for (size_t t = 0; t < time_size; ++t) {
        BOTTOMT[t].resize(lat_size);
        for (size_t i = 0; i < lat_size; ++i) {
            BOTTOMT[t][i].resize(lon_size);
            for (size_t j = 0; j < lon_size; ++j) {
                size_t index = t * lat_size * lon_size + i * lon_size + j;
                BOTTOMT[t][i][j] = bottomT[index];
            }
        }
    }
}

CopernicusTem::~CopernicusTem() {}

std::vector<std::vector<std::vector<std::vector<double>>>> CopernicusTem::getTHETAO() const {
    return THETAO;
}

std::vector<std::vector<std::vector<double>>> CopernicusTem::getBOTTOMT() const {
    return BOTTOMT;
}