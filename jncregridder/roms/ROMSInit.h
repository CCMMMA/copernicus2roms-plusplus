//
// Created by Ciro De Vita on 29/03/24.
//

#ifndef MYOCEAN2ROMS_ROMSINIT_H
#define MYOCEAN2ROMS_ROMSINIT_H

#include <netcdf>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <omp.h>
#include "ROMSGrid.h"

class ROMSInit {
public:
    ROMSInit(const std::string& url, const ROMSGrid& romsGrid, size_t forcingTimeSteps);
    ~ROMSInit();

    void setOceanTime(std::vector<double> ocean_time);
    void setScrumTime(std::vector<int> scrumTime);
    void setZETA(std::vector<std::vector<double>> zeta);
    void setSALT(std::vector<std::vector<std::vector<double>>> salt);
    void setTEMP(std::vector<std::vector<std::vector<double>>> temp);
    void setU(std::vector<std::vector<std::vector<double>>> u);
    void setV(std::vector<std::vector<std::vector<double>>> v);
    void setVBAR(std::vector<std::vector<double>> vbar);
    void setUBAR(std::vector<std::vector<double>> ubar);

    void make();
    void write(size_t time);

private:
    std::string url;
    netCDF::NcFile* ncfWritable;
    netCDF::NcDim dimOceanTime, dimScrumTime, dimEtaRho, dimXiRho, dimEtaU, dimXiU, dimEtaV, dimXiV, dimSRho, dimOne;
    netCDF::NcVar h, lat_rho, lon_rho, lat_u, lon_u, lat_v, lon_v, temp, salt, ubar, vbar, u, v, zeta, theta_b, theta_s, Tcline, ocean_time, hc, scrum_time, tend, Cs_r, sc_r, s_rho;
    float missing_value = 1e37;

    std::vector<int> SCRUM_TIME;
    std::vector<double> H, LATRHO, LONRHO, LATU, LONU, LATV, LONV, THETA_S, THETA_B, TCLINE, OCEAN_TIME, SC_R, S_RHO, HC, CS_R;
    std::vector<std::vector<double>> ZETA, UBAR, VBAR;
    std::vector<std::vector<std::vector<double>>> SALT, TEMP, U, V;

    std::string getCurrentDateTime() const;
};


#endif //MYOCEAN2ROMS_ROMSINIT_H
