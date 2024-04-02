//
// Created by Ciro De Vita on 29/03/24.
//

#include "ROMSInit.h"

ROMSInit::ROMSInit(const std::string& url, const ROMSGrid& romsGrid, size_t forcingTimeSteps) {
    H = romsGrid.getH();

    LATRHO = romsGrid.getLATRHO();
    LONRHO = romsGrid.getLONRHO();
    LATV = romsGrid.getLATV();
    LONV = romsGrid.getLONV();
    LATU = romsGrid.getLATU();
    LONU = romsGrid.getLONU();

    THETA_S = romsGrid.getTheta_s();
    THETA_B = romsGrid.getTheta_b();
    TCLINE = romsGrid.getTCLINE();

    SC_R = romsGrid.getS_rho();
    S_RHO = romsGrid.getS_rho();
    HC = romsGrid.getHC();
    CS_R = romsGrid.getCs_r();

    ncfWritable = new netCDF::NcFile(url, netCDF::NcFile::replace);
    ncfWritable->putAtt("type", "Initial file");
    ncfWritable->putAtt("title", "Initialization file (INI) used for forcing of the ROMS model");
    ncfWritable->putAtt("grd_file", romsGrid.getURL());
    ncfWritable->putAtt("source", "University of Napoli Parthenope Weather Centre http://meteo.uniparthenope.it");
    ncfWritable->putAtt("date", getCurrentDateTime());

    dimOceanTime = ncfWritable->addDim("ocean_time", forcingTimeSteps);
    dimScrumTime = ncfWritable->addDim("scrum_time", forcingTimeSteps);
    dimEtaRho = ncfWritable->addDim("eta_rho", romsGrid.getEtaRho());
    dimXiRho = ncfWritable->addDim("xi_rho", romsGrid.getXiRho());
    dimEtaU = ncfWritable->addDim("eta_u", romsGrid.getEtaU());
    dimXiU = ncfWritable->addDim("xi_u", romsGrid.getXiU());
    dimEtaV = ncfWritable->addDim("eta_v", romsGrid.getEtaV());
    dimXiV = ncfWritable->addDim("xi_v", romsGrid.getXiV());
    dimSRho = ncfWritable->addDim("s_rho", romsGrid.getS_rho().size());

    dimOne = ncfWritable->addDim("one", 1);

    h = ncfWritable->addVar("h", netCDF::ncDouble, {dimEtaRho, dimXiRho});
    h.putAtt("long_name", "Final bathymetry at RHO-points");
    h.putAtt("units", "meters");
    h.putAtt("field", "bath, scalar");

    lat_rho = ncfWritable->addVar("lat_rho", netCDF::ncDouble, {dimEtaRho, dimXiRho});
    lat_rho.putAtt("long_name", "latitude of RHO-points");
    lat_rho.putAtt("units", "degree_north");
    lat_rho.putAtt("field", "lat_rho, scalar");
    lat_rho.putAtt("standard_name", "latitude");
    lat_rho.putAtt("_CoordinateAxisType", "Lat");

    lon_rho = ncfWritable->addVar("lon_rho", netCDF::ncDouble, {dimEtaRho, dimXiRho});
    lon_rho.putAtt("long_name", "longitude of RHO-points");
    lon_rho.putAtt("units", "degree_east");
    lon_rho.putAtt("field", "lon_rho, scalar");
    lon_rho.putAtt("standard_name", "longitude");
    lon_rho.putAtt("_CoordinateAxisType", "Lon");

    lat_u = ncfWritable->addVar("lat_u", netCDF::ncDouble, {dimEtaU, dimXiU});
    lat_u.putAtt("long_name", "latitude of U-points");
    lat_u.putAtt("units", "degree_north");
    lat_u.putAtt("standard_name", "latitude");
    lat_u.putAtt("_CoordinateAxisType", "Lat");

    lon_u = ncfWritable->addVar("lon_u", netCDF::ncDouble, {dimEtaU, dimXiU});
    lon_u.putAtt("long_name", "longitude of U-points");
    lon_u.putAtt("units", "degree_east");
    lon_u.putAtt("standard_name", "longitude");
    lon_u.putAtt("_CoordinateAxisType", "Lon");

    lat_v = ncfWritable->addVar("lat_v", netCDF::ncDouble, {dimEtaV, dimXiV});
    lat_v.putAtt("long_name", "latitude of V-points");
    lat_v.putAtt("units", "degree_north");
    lat_v.putAtt("standard_name", "latitude");
    lat_v.putAtt("_CoordinateAxisType", "Lat");

    lon_v = ncfWritable->addVar("lon_v", netCDF::ncDouble, {dimEtaV, dimXiV});
    lon_v.putAtt("long_name", "longitude of V-points");
    lon_v.putAtt("units", "degree_east");
    lon_v.putAtt("standard_name", "longitude");
    lon_v.putAtt("_CoordinateAxisType", "Lon");

    temp = ncfWritable->addVar("temp", netCDF::ncFloat, {dimOceanTime, dimSRho, dimEtaRho, dimXiRho});
    temp.putAtt("long_name", "potential temperature");
    temp.putAtt("units", "Celsius");
    temp.putAtt("coordinates", "lon_rho lat_rho sc_r ocean_time");
    temp.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    temp.putAtt("time", "ocean_time");

    salt = ncfWritable->addVar("salt", netCDF::ncFloat, {dimOceanTime, dimSRho, dimEtaRho, dimXiRho});
    salt.putAtt("long_name", "salinity");
    salt.putAtt("units", "PSU");
    salt.putAtt("coordinates", "lon_rho lat_rho sc_r ocean_time");
    salt.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    salt.putAtt("time", "ocean_time");

    ubar = ncfWritable->addVar("ubar", netCDF::ncFloat, {dimOceanTime, dimEtaU, dimXiU});
    ubar.putAtt("long_name", "vertically integrated u-momentum component");
    ubar.putAtt("units", "meter second-1");
    ubar.putAtt("coordinates", "lon_u lat_u ocean_time");
    ubar.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    ubar.putAtt("time", "ocean_time");

    vbar = ncfWritable->addVar("vbar", netCDF::ncFloat, {dimOceanTime, dimEtaV, dimXiV});
    vbar.putAtt("long_name", "vertically integrated v-momentum component");
    vbar.putAtt("units", "meter second-1");
    vbar.putAtt("coordinates", "lon_v lat_v ocean_time");
    vbar.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    vbar.putAtt("time", "ocean_time");

    u = ncfWritable->addVar("u", netCDF::ncFloat, {dimOceanTime, dimSRho, dimEtaU, dimXiU});
    u.putAtt("long_name", "u-momentum component");
    u.putAtt("units", "meter second-1");
    u.putAtt("coordinates", "lon_u lat_u sc_r ocean_time");
    u.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    u.putAtt("time", "ocean_time");

    v = ncfWritable->addVar("v", netCDF::ncFloat, {dimOceanTime, dimSRho, dimEtaV, dimXiV});
    v.putAtt("long_name", "v-momentum component");
    v.putAtt("units", "meter second-1");
    v.putAtt("coordinates", "lon_v lat_v sc_r ocean_time");
    v.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    v.putAtt("time", "ocean_time");

    zeta = ncfWritable->addVar("zeta", netCDF::ncFloat, {dimOceanTime, dimEtaRho, dimXiRho});
    zeta.putAtt("long_name", "free-surface");
    zeta.putAtt("units", "meter");
    zeta.putAtt("coordinates", "lon_rho lat_rho ocean_time");
    zeta.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    zeta.putAtt("time", "ocean_time");

    theta_b = ncfWritable->addVar("theta_b", netCDF::ncFloat, std::initializer_list<netCDF::NcDim>{dimOne});
    theta_b.putAtt("long_name", "S-coordinate surface control parameter");
    theta_b.putAtt("units", "nondimensional");

    theta_s = ncfWritable->addVar("theta_s", netCDF::ncFloat, std::initializer_list<netCDF::NcDim>{dimOne});
    theta_s.putAtt("long_name", "S-coordinate bottom control parameter");
    theta_s.putAtt("units", "nondimensional");

    Tcline = ncfWritable->addVar("Tcline", netCDF::ncFloat, std::initializer_list<netCDF::NcDim>{dimOne});

    ocean_time = ncfWritable->addVar("ocean_time", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime});
    ocean_time.putAtt("long_name", "ocean forcing time");
    ocean_time.putAtt("units", "days since 1968-05-23 00:00:00 GMT");
    ocean_time.putAtt("calendar", "gregorian");

    hc = ncfWritable->addVar("hc", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOne});
    hc.putAtt("long_name", "S-coordinate parameter, critical depth");
    hc.putAtt("units", "meter");

    scrum_time = ncfWritable->addVar("scrum_time", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimScrumTime});
    scrum_time.putAtt("long_name", "time since initialization");
    scrum_time.putAtt("units", "second");

    tend = ncfWritable->addVar("tend", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOne});
    tend.putAtt("long_name", "end processing day");
    tend.putAtt("units", "day");

    Cs_r = ncfWritable->addVar("Cs_r", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimSRho});
    Cs_r.putAtt("long_name", "S-coordinate stretching curves at RHO-points");
    Cs_r.putAtt("units", "nondimensional");

    sc_r = ncfWritable->addVar("sc_r", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimSRho});
    sc_r.putAtt("long_name", "S-coordinate at RHO-points");
    sc_r.putAtt("units", "nondimensional");

    s_rho = ncfWritable->addVar("s_rho", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimSRho});
    s_rho.putAtt("long_name", "oS-coordinate at RHO-points");
    s_rho.putAtt("valid_min", netCDF::ncDouble, -1.0);
    s_rho.putAtt("valid_max", netCDF::ncDouble, 0.0);
    s_rho.putAtt("positive", "up");
    s_rho.putAtt("standard_name", "ocean_s_coordinate_g1");
    s_rho.putAtt("formula_terms", "");
    s_rho.putAtt("field", "s_rho, scalar");
    s_rho.putAtt("_CoordinateTransformType", "Vertical");
    s_rho.putAtt("_CoordinateAxes", "s_rho");
    s_rho.putAtt("_CoordinateAxisType", "GeoZ");
    s_rho.putAtt("_CoordinateZisPositive", "up");

}

ROMSInit::~ROMSInit() {
    if (ncfWritable) {
        ncfWritable->close();
        delete ncfWritable;
    }
}

void ROMSInit::make() {
    h.putVar(&H.front());

    lat_rho.putVar(&LATRHO.front());
    lon_rho.putVar(&LONRHO.front());
    lat_u.putVar(&LATU.front());
    lon_u.putVar(&LONU.front());
    lat_v.putVar(&LATV.front());
    lon_v.putVar(&LONV.front());

    theta_s.putVar(&THETA_S.front());
    theta_b.putVar(&THETA_B.front());
    Tcline.putVar(&TCLINE.front());

    ocean_time.putVar(&OCEAN_TIME.front());
    scrum_time.putVar(&SCRUM_TIME.front());

    sc_r.putVar(&SC_R.front());
    s_rho.putVar(&S_RHO.front());
    hc.putVar(&HC.front());
    Cs_r.putVar(&CS_R.front());
}

void ROMSInit::write(size_t time) {
    std::vector<double> buffer(dimSRho.getSize() * dimEtaRho.getSize() * dimXiRho.getSize() * 2 +
                              dimSRho.getSize() * dimEtaV.getSize() * dimXiV.getSize() +
                              dimSRho.getSize() * dimEtaU.getSize() * dimXiU.getSize() +
                              dimEtaRho.getSize() * dimXiRho.getSize() +
                              dimEtaV.getSize() * dimXiV.getSize() +
                              dimEtaU.getSize() * dimXiU.getSize());

    size_t idx = 0;

    for (int k = 0; k < dimSRho.getSize(); k++) {
        for (int j = 0; j < dimEtaRho.getSize(); j++) {
            for (int i = 0; i < dimXiRho.getSize(); i++) {
                buffer[idx++] = SALT[k][j][i];
            }
        }
    }

    for (int k = 0; k < dimSRho.getSize(); k++) {
        for (int j = 0; j < dimEtaRho.getSize(); j++) {
            for (int i = 0; i < dimXiRho.getSize(); i++) {
                buffer[idx++] = TEMP[k][j][i];
            }
        }
    }

    for (int k = 0; k < dimSRho.getSize(); k++) {
        for (int j = 0; j < dimEtaV.getSize(); j++) {
            for (int i = 0; i < dimXiV.getSize(); i++) {
                buffer[idx++] = V[k][j][i];
            }
        }
    }

    for (int k = 0; k < dimSRho.getSize(); k++) {
        for (int j = 0; j < dimEtaU.getSize(); j++) {
            for (int i = 0; i < dimXiU.getSize(); i++) {
                buffer[idx++] = U[k][j][i];
            }
        }
    }

    for (int j = 0; j < dimEtaRho.getSize(); j++) {
        for (int i = 0; i < dimXiRho.getSize(); i++) {
            buffer[idx++] = ZETA[j][i];
        }
    }

    for (int j = 0; j < dimEtaV.getSize(); j++) {
        for (int i = 0; i < dimXiV.getSize(); i++) {
            buffer[idx++] = VBAR[j][i];
        }
    }

    for (int j = 0; j < dimEtaU.getSize(); j++) {
        for (int i = 0; i < dimXiU.getSize(); i++) {
            buffer[idx++] = UBAR[j][i];
        }
    }

    std::vector<size_t> start_salt = {time, 0, 0, 0};
    std::vector<size_t> count_salt = {1, dimSRho.getSize(), dimEtaRho.getSize(), dimXiRho.getSize()};
    salt.putVar(start_salt, count_salt, buffer.data());

    idx = dimSRho.getSize() * dimEtaRho.getSize() * dimXiRho.getSize();

    std::vector<size_t> start_temp = {time, 0, 0, 0};
    std::vector<size_t> count_temp = {1, dimSRho.getSize(), dimEtaRho.getSize(), dimXiRho.getSize()};
    temp.putVar(start_temp, count_temp, buffer.data() + idx);

    idx += dimSRho.getSize() * dimEtaRho.getSize() * dimXiRho.getSize();

    std::vector<size_t> start_v = {time, 0, 0, 0};
    std::vector<size_t> count_v = {1, dimSRho.getSize(), dimEtaV.getSize(), dimXiV.getSize()};
    v.putVar(start_v, count_v, buffer.data() + idx);

    idx += dimSRho.getSize() * dimEtaV.getSize() * dimXiV.getSize();

    std::vector<size_t> start_u = {time, 0, 0, 0};
    std::vector<size_t> count_u = {1, dimSRho.getSize(), dimEtaU.getSize(), dimXiU.getSize()};
    u.putVar(start_u, count_u, buffer.data() + idx);

    idx += dimSRho.getSize() * dimEtaU.getSize() * dimXiU.getSize();

    std::vector<size_t> start_zeta = {time, 0, 0};
    std::vector<size_t> count_zeta = {1, dimEtaRho.getSize(), dimXiRho.getSize()};
    zeta.putVar(start_zeta, count_zeta, buffer.data() + idx);

    idx += dimEtaRho.getSize() * dimXiRho.getSize();

    std::vector<size_t> start_vbar = {time, 0, 0};
    std::vector<size_t> count_vbar = {1, dimEtaV.getSize(), dimXiV.getSize()};
    vbar.putVar(start_vbar, count_vbar, buffer.data() + idx);

    idx += dimEtaV.getSize() * dimXiV.getSize();

    std::vector<size_t> start_ubar = {time, 0, 0};
    std::vector<size_t> count_ubar = {1, dimEtaU.getSize(), dimXiU.getSize()};
    ubar.putVar(start_ubar, count_ubar, buffer.data() + idx);
}

std::string ROMSInit::getCurrentDateTime() const {
    std::time_t now = std::time(nullptr);
    char buf[sizeof "YYYY-MM-DD HH:MM:SS"];
    std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return buf;
}

void ROMSInit::setOceanTime(std::vector<double> ocean_time){
    OCEAN_TIME = ocean_time;
}

void ROMSInit::setScrumTime(std::vector<int> scrumTime){
    SCRUM_TIME = scrumTime;
}

void ROMSInit::setZETA(std::vector<std::vector<double>> zeta) {
    ZETA = zeta;
}

void ROMSInit::setSALT(std::vector<std::vector<std::vector<double>>> salt) {
    SALT = salt;
}

void ROMSInit::setTEMP(std::vector<std::vector<std::vector<double>>> temp) {
    TEMP = temp;
}

void ROMSInit::setU(std::vector<std::vector<std::vector<double>>> u) {
    U = u;
}

void ROMSInit::setV(std::vector<std::vector<std::vector<double>>> v) {
    V = v;
}

void ROMSInit::setUBAR(std::vector<std::vector<double>> ubar) {
    UBAR = ubar;
}

void ROMSInit::setVBAR(std::vector<std::vector<double>> vbar) {
    VBAR = vbar;
}