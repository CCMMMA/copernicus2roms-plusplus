//
// Created by Ciro De Vita on 29/03/24.
//

#ifndef MYOCEAN2ROMS_ROMSBOUNDARY_H
#define MYOCEAN2ROMS_ROMSBOUNDARY_H

#include <iostream>
#include "ROMSGrid.h"

class ROMSBoundary {
public:
    ROMSBoundary(const std::string& url, const ROMSGrid& romsGrid, size_t forcingTimeSteps);
    ~ROMSBoundary();

    void setOceanTime(std::vector<double> ocean_time);
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
    float missing_value = 1e37;
    netCDF::NcDim dimOceanTime, dimEtaRho, dimXiRho, dimEtaU, dimXiU, dimEtaV, dimXiV, dimSRho, dimSw, dimOne, dimTwo, dimFour, dimBath;
    netCDF::NcVar ocean_time, angle, theta_b, theta_s, Tcline, z_r, hc, Cs_w, Cs_r, s_w, sc_r, s_rho, h, lat_rho, lon_rho, lat_u, lon_u, lat_v, lon_v,
    temp_west, temp_east, temp_south, temp_north, salt_west, salt_east, salt_south, salt_north, zeta_west, zeta_east, zeta_south, zeta_north,
    u_west, u_east, u_south, u_north, v_west, v_east, v_south, v_north, vbar_west, vbar_east, vbar_south, vbar_north, ubar_west, ubar_east, ubar_south, ubar_north;

    std::vector<double> H, LATRHO, LONRHO, LATU, LONU, LATV, LONV, THETA_S, THETA_B, TCLINE, OCEAN_TIME, SC_R, S_RHO, HC, CS_R;
    std::vector<std::vector<double>> ZETA, UBAR, VBAR;
    std::vector<std::vector<std::vector<double>>> SALT, TEMP, U, V;


    std::string getCurrentDateTime() const;
};


#endif //MYOCEAN2ROMS_ROMSBOUNDARY_H
