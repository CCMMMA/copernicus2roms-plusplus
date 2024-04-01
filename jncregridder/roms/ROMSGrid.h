//
// Created by Ciro De Vita on 28/03/24.
//

#ifndef MYOCEAN2ROMS_ROMSGRID_H
#define MYOCEAN2ROMS_ROMSGRID_H

#include <iostream>
#include <netcdf>


class ROMSGrid {
public:
    ROMSGrid(const std::string& url);

    // Getters
    std::string getURL() const;
    size_t getEtaRho() const;
    size_t getXiRho() const;
    size_t getEtaPsi() const;
    size_t getXiPsi() const;
    size_t getEtaU() const;
    size_t getXiU() const;
    size_t getEtaV() const;
    size_t getXiV() const;
    double getNoData() const;
    double getZeta() const;
    std::vector<double> getS_rho() const;
    std::vector<double> getCs_r() const;
    std::vector<double> getS_w() const;
    std::vector<double> getCs_w() const;
    std::vector<double> getTheta_s() const;
    std::vector<double> getTheta_b() const;
    std::vector<double> getANGLE() const;
    std::vector<double> getLATRHO() const;
    std::vector<double> getLONRHO() const;
    std::vector<double> getLATPSI() const;
    std::vector<double> getLONPSI() const;
    std::vector<double> getLATU() const;
    std::vector<double> getLONU() const;
    std::vector<double> getLATV() const;
    std::vector<double> getLONV() const;
    std::vector<double> getH() const;
    std::vector<double> getHC() const;
    std::vector<double> getMASKRHO() const;
    std::vector<double> getMASKU() const;
    std::vector<double> getMASKV() const;
    std::vector<double> getZ() const;
    std::vector<double> getTCLINE() const;

    // Setters
    void setNoData(double value);
    void setZeta(double value);

private:
    std::string url;
    netCDF::NcFile ncDataset;
    size_t etaRho, xiRho, etaPsi, xiPsi, etaU, xiU, etaV, xiV;
    double no_data = 1e37;
    double zeta = 0;
    std::vector<double> s_rho, cs_r, s_w, cs_w, theta_s, theta_b, ANGLE,
            LATRHO, LONRHO, LATPSI, LONPSI, LATU, LONU, LATV, LONV,
            H, HC, MASKRHO, MASKU, MASKV, Z, TCLINE;

    std::vector<double> loadVariable(const std::string& variable_name);
};


#endif //MYOCEAN2ROMS_ROMSGRID_H
