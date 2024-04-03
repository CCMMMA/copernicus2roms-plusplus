//
// Created by Ciro De Vita on 28/03/24.
//

#include "CopernicusBase.h"

template<>
std::vector<double> CopernicusBase::loadVariable<double>(const std::string& variable_name) {
    std::vector<double> variable_data;

    auto variable = ncFile->getVar(variable_name);
    // Get the dimensions of the variable
    auto dims = variable.getDims();
    size_t totalSize = 1;
    for (const auto& dim : dims) {
        totalSize *= dim.getSize();
    }

    // Resize the vector to hold the variable data
    variable_data.resize(totalSize);

    // Read the variable data from the file
    variable.getVar(variable_data.data());

    return variable_data;
}

template<>
std::vector<std::vector<double>> CopernicusBase::loadVariable<std::vector<double>>(const std::string& variable_name) {
    if (variable_name == "LAT2") {
        std::vector<std::vector<double>> lat;
        lat.resize(LAT.size());
        for (int j = 0; j < LAT.size(); j++) {
            lat[j].resize(LON.size());
            for (int i = 0; i < LON.size(); i++) {
                lat[j][i] = LAT[j];
            }
        }

        return lat;
    } else if (variable_name == "LON2") {
        std::vector<std::vector<double>> lon;
        lon.resize(LAT.size());
        for (int j = 0; j < LAT.size(); j++) {
            lon[j].resize(LON.size());
            for (int i = 0; i < LON.size(); i++) {
                lon[j][i] = LON[i];
            }
        }
        return lon;
    }

    return {};
}

CopernicusBase::CopernicusBase(const std::string& url) {
    try {
        ncFile = new netCDF::NcFile(url, netCDF::NcFile::read);

        // Load dimensions
        auto dimTime = ncFile->getDim("time");
        auto dimLat = ncFile->getDim("latitude");
        auto dimLon = ncFile->getDim("longitude");
        auto dimDepth = ncFile->getDim("depth");
        b3D = !dimDepth.isNull();

        // Set dimension sizes
        TIME = loadVariable<double>("time");
        LAT = loadVariable<double>("latitude");
        LON = loadVariable<double>("longitude");
        DEPTH = b3D ? loadVariable<double>("depth") : std::vector<double>();
        Z = b3D ? calculateZ() : std::vector<std::vector<std::vector<double>>>();
        LAT2 = loadVariable<std::vector<double>>("LAT2");
        LON2 = loadVariable<std::vector<double>>("LON2");
    } catch (const netCDF::exceptions::NcException& e) {
        std::cerr << "NetCDF error: " << e.what() << std::endl;
    }
}

CopernicusBase::~CopernicusBase() {}

int CopernicusBase::getTime() const {
    return TIME.size();
}

int CopernicusBase::getLat() const {
    return LAT.size();
}

int CopernicusBase::getLon() const {
    return LON.size();
}

bool CopernicusBase::hasDepth() const {
    return b3D;
}

int CopernicusBase::getDepth() const {
    return DEPTH.size();
}

std::vector<double> CopernicusBase::getLAT() const {
    return LAT;
}

std::vector<std::vector<double>> CopernicusBase::getLAT2() const {
    return LAT2;
}

std::vector<double> CopernicusBase::getLON() const {
    return LON;
}

std::vector<std::vector<double>> CopernicusBase::getLON2() const {
    return LON2;
}

std::vector<double> CopernicusBase::getTIME() const {
    return TIME;
}

std::vector<double> CopernicusBase::getDEPTH() const {
    return DEPTH;
}

std::vector<std::vector<std::vector<double>>> CopernicusBase::getZ() const {
    return Z;
}

std::vector<std::vector<std::vector<double>>> CopernicusBase::calculateZ() const {
    // Get time, depth, lat, and lon sizes
    size_t depth_size = getDepth();
    size_t lat_size = getLat();
    size_t lon_size = getLon();

    std::vector<std::vector<std::vector<double>>> z(depth_size, std::vector<std::vector<double>>(lat_size, std::vector<double>(lon_size, 0.0)));

    for (size_t d = 0; d < depth_size; ++d) {
        double depth_value = DEPTH[d];
        for (size_t i = 0; i < lat_size; ++i) {
            for (size_t j = 0; j < lon_size; ++j) {
                z[d][i][j] = depth_value;
            }
        }
    }

    return z;
}