//
// Created by Ciro De Vita on 28/03/24.
//

#include "CopernicusSal.h"

CopernicusSal::CopernicusSal(const std::string& url) : CopernicusBase(url) {
    // Load SO variable
    std::vector<double> so = loadVariable<double>("so");

    // Get time, depth, lat, and lon sizes
    size_t time_size = getTime();
    size_t depth_size = getDepth();
    size_t lat_size = getLat();
    size_t lon_size = getLon();

    // Reshape UO into a 4D vector
    SO.resize(time_size);
    for (size_t t = 0; t < time_size; ++t) {
        SO[t].resize(depth_size);
        for (size_t d = 0; d < depth_size; ++d) {
            SO[t][d].resize(lat_size);
            for (size_t i = 0; i < lat_size; ++i) {
                SO[t][d][i].resize(lon_size);
                for (size_t j = 0; j < lon_size; ++j) {
                    size_t index = (t * depth_size * lat_size * lon_size) + (d * lat_size * lon_size) + (i * lon_size) + j;
                    SO[t][d][i][j] = so[index];
                }
            }
        }
    }
}

CopernicusSal::~CopernicusSal() {}

std::vector<std::vector<std::vector<std::vector<double>>>> CopernicusSal::getSO() const {
    return SO;
}