//
// Created by Ciro De Vita on 28/03/24.
//

#ifndef MYOCEAN2ROMS_COPERNICUS_BASE_H
#define MYOCEAN2ROMS_COPERNICUS_BASE_H

#include <netcdf>
#include <string>
#include <vector>

class CopernicusBase {
public:
    CopernicusBase(const std::string& url);
    virtual ~CopernicusBase();

    int getTime() const;
    int getLat() const;
    int getLon() const;
    bool hasDepth() const;
    int getDepth() const;

    std::vector<double> getLAT() const;
    std::vector<std::vector<double>> getLAT2() const;
    std::vector<double> getLON() const;
    std::vector<std::vector<double>> getLON2() const;
    std::vector<double> getTIME() const;
    std::vector<double> getDEPTH() const;
    std::vector<std::vector<std::vector<double>>> getZ() const;

protected:
    template<typename T>
    std::vector<T> loadVariable(const std::string& variable_name);

private:
    netCDF::NcFile* ncFile;
    std::vector<double> LAT;
    std::vector<std::vector<double>> LAT2;
    std::vector<double> LON;
    std::vector<std::vector<double>> LON2;
    std::vector<double> TIME;
    std::vector<double> DEPTH;
    std::vector<std::vector<std::vector<double>>> Z;
    bool b3D;

    std::vector<std::vector<std::vector<double>>> calculateZ() const;
};


#endif //MYOCEAN2ROMS_COPERNICUS_BASE_H
