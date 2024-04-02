//
// Created by Ciro De Vita on 29/03/24.
//

#include "ROMSBoundary.h"

ROMSBoundary::ROMSBoundary(const std::string& url, const ROMSGrid& romsGrid, size_t forcingTimeSteps) {
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
    ncfWritable->putAtt("type", "Boundary forcing file");
    ncfWritable->putAtt("title", "Boundary forcing file (BRY) used for forcing of the ROMS model");
    ncfWritable->putAtt("grd_file", romsGrid.getURL());
    ncfWritable->putAtt("source", "University of Napoli Parthenope Weather Centre http://meteo.uniparthenope.it");
    ncfWritable->putAtt("date", getCurrentDateTime());

    dimEtaRho = ncfWritable->addDim("eta_rho", romsGrid.getEtaRho());
    dimXiRho = ncfWritable->addDim("xi_rho", romsGrid.getXiRho());
    dimEtaU = ncfWritable->addDim("eta_u", romsGrid.getEtaU());
    dimXiU = ncfWritable->addDim("xi_u", romsGrid.getXiU());
    dimEtaV = ncfWritable->addDim("eta_v", romsGrid.getEtaV());
    dimXiV = ncfWritable->addDim("xi_v", romsGrid.getXiV());
    dimSRho = ncfWritable->addDim("s_rho", romsGrid.getS_rho().size());
    dimSw = ncfWritable->addDim("s_w", romsGrid.getS_w().size());

    dimOceanTime = ncfWritable->addDim("ocean_time", forcingTimeSteps);

    dimOne = ncfWritable->addDim("one", 1);
    dimTwo = ncfWritable->addDim("two", 2);
    dimFour = ncfWritable->addDim("four", 4);
    dimBath = ncfWritable->addDim("bath", 1);

    ocean_time = ncfWritable->addVar("ocean_time", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime});
    ocean_time.putAtt("long_name", "ocean forcing time");
    ocean_time.putAtt("units", "days since 1968-05-23 00:00:00 GMT");
    ocean_time.putAtt("calendar", "gregorian");

    angle = ncfWritable->addVar("angle", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimEtaRho, dimXiRho});
    angle.putAtt("long_name", "angle between xu axis and east");
    angle.putAtt("units", "radian");

    theta_b = ncfWritable->addVar("theta_b", netCDF::ncFloat, std::initializer_list<netCDF::NcDim>{dimOne});
    theta_b.putAtt("long_name", "S-coordinate surface control parameter");
    theta_b.putAtt("units", "nondimensional");

    theta_s = ncfWritable->addVar("theta_s", netCDF::ncFloat, std::initializer_list<netCDF::NcDim>{dimOne});
    theta_s.putAtt("long_name", "S-coordinate bottom control parameter");
    theta_s.putAtt("units", "nondimensional");

    Tcline = ncfWritable->addVar("Tcline", netCDF::ncFloat, std::initializer_list<netCDF::NcDim>{dimOne});

    z_r = ncfWritable->addVar("z_r", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimSRho, dimEtaRho, dimXiRho});
    z_r.putAtt("long_name", "Sigma layer to depth matrix");
    z_r.putAtt("units", "meter");

    hc = ncfWritable->addVar("hc", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOne});
    hc.putAtt("long_name", "S-coordinate parameter, critical depth");
    hc.putAtt("units", "meter");

    Cs_w = ncfWritable->addVar("Cs_w", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimSw});
    Cs_w.putAtt("long_name", "S-coordinate stretching curves at W-points");
    Cs_w.putAtt("valid_min", netCDF::ncDouble, -1.0);
    Cs_w.putAtt("valid_max", netCDF::ncDouble, 0.0);
    Cs_w.putAtt("field", "s_w, scalar");

    Cs_r = ncfWritable->addVar("Cs_r", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimSRho});
    Cs_r.putAtt("long_name", "S-coordinate stretching curves at RHO-points");
    Cs_r.putAtt("units", "nondimensional");

    s_w = ncfWritable->addVar("s_w", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimSw});
    s_w.putAtt("long_name", "S-coordinate at W-points");
    s_w.putAtt("valid_min", netCDF::ncDouble, -1);
    s_w.putAtt("valid_max", netCDF::ncDouble, 0);
    s_w.putAtt("standard_name", "ocean_s_coordinate_g1");
    s_w.putAtt("formula_terms", "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc");
    s_w.putAtt("field", "s_w, scalar");

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

    h = ncfWritable->addVar("h", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimEtaRho, dimXiRho});
    h.putAtt("long_name", "Final bathymetry at RHO-points");
    h.putAtt("units", "meters");
    h.putAtt("field", "bath, scalar");

    lat_rho = ncfWritable->addVar("lat_rho", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimEtaRho, dimXiRho});
    lat_rho.putAtt("long_name", "latitude of RHO-points");
    lat_rho.putAtt("units", "degree_north");
    lat_rho.putAtt("field", "lat_rho, scalar");
    lat_rho.putAtt("standard_name", "latitude");
    lat_rho.putAtt("_CoordinateAxisType", "Lat");

    lon_rho = ncfWritable->addVar("lon_rho", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimEtaRho, dimXiRho});
    lon_rho.putAtt("long_name", "longitude of RHO-points");
    lon_rho.putAtt("units", "degree_east");
    lon_rho.putAtt("field", "lon_rho, scalar");
    lon_rho.putAtt("standard_name", "longitude");
    lon_rho.putAtt("_CoordinateAxisType", "Lon");

    lat_u = ncfWritable->addVar("lat_u", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimEtaU, dimXiU});
    lat_u.putAtt("long_name", "latitude of U-points");
    lat_u.putAtt("units", "degree_north");
    lat_u.putAtt("standard_name", "latitude");
    lat_u.putAtt("_CoordinateAxisType", "Lat");

    lon_u = ncfWritable->addVar("lon_u", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimEtaU, dimXiU});
    lon_u.putAtt("long_name", "longitude of U-points");
    lon_u.putAtt("units", "degree_east");
    lon_u.putAtt("standard_name", "longitude");
    lon_u.putAtt("_CoordinateAxisType", "Lon");

    lat_v = ncfWritable->addVar("lat_v", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimEtaV, dimXiV});
    lat_v.putAtt("long_name", "latitude of V-points");
    lat_v.putAtt("units", "degree_north");
    lat_v.putAtt("standard_name", "latitude");
    lat_v.putAtt("_CoordinateAxisType", "Lat");

    lon_v = ncfWritable->addVar("lon_v", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimEtaV, dimXiV});
    lon_v.putAtt("long_name", "longitude of V-points");
    lon_v.putAtt("units", "degree_east");
    lon_v.putAtt("standard_name", "longitude");
    lon_v.putAtt("_CoordinateAxisType", "Lon");

    temp_west = ncfWritable->addVar("temp_west", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimEtaRho});
    temp_west.putAtt("long_name", "potential temperature western boundary conditions");
    temp_west.putAtt("units", "Celsius");
    temp_west.putAtt("field", "temp_west, scalar, series");
    temp_west.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    temp_west.putAtt("time", "ocean_time");

    temp_east = ncfWritable->addVar("temp_east", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimEtaRho});
    temp_east.putAtt("long_name", "potential temperature eastern boundary conditions");
    temp_east.putAtt("units", "Celsius");
    temp_east.putAtt("field", "temp_east, scalar, series");
    temp_east.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    temp_east.putAtt("time", "ocean_time");

    temp_south = ncfWritable->addVar("temp_south", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimXiRho});
    temp_south.putAtt("long_name", "potential temperature southern boundary conditions");
    temp_south.putAtt("units", "Celsius");
    temp_south.putAtt("field", "temp_south, scalar, series");
    temp_south.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    temp_south.putAtt("time", "ocean_time");

    temp_north = ncfWritable->addVar("temp_north", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimXiRho});
    temp_north.putAtt("long_name", "potential temperature northern boundary conditions");
    temp_north.putAtt("units", "Celsius");
    temp_north.putAtt("field", "temp_north, scalar, series");
    temp_north.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    temp_north.putAtt("time", "ocean_time");

    salt_west = ncfWritable->addVar("salt_west", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimEtaRho});
    salt_west.putAtt("long_name", "salinity western boundary conditions");
    salt_west.putAtt("units", "PSU");
    salt_west.putAtt("field", "salt_west, scalar, series");
    salt_west.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    salt_west.putAtt("time", "ocean_time");

    salt_east = ncfWritable->addVar("salt_east", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimEtaRho});
    salt_east.putAtt("long_name", "salinity eastern boundary conditions");
    salt_east.putAtt("units", "PSU");
    salt_east.putAtt("field", "salt_east, scalar, series");
    salt_east.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    salt_east.putAtt("time", "ocean_time");

    salt_south = ncfWritable->addVar("salt_south", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimXiRho});
    salt_south.putAtt("long_name", "salinity southern boundary conditions");
    salt_south.putAtt("units", "PSU");
    salt_south.putAtt("field", "salt_south, scalar, series");
    salt_south.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    salt_south.putAtt("time", "ocean_time");

    salt_north = ncfWritable->addVar("salt_north", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimXiRho});
    salt_north.putAtt("long_name", "salinity northern boundary conditions");
    salt_north.putAtt("units", "PSU");
    salt_north.putAtt("field", "salt_north, scalar, series");
    salt_north.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    salt_north.putAtt("time", "ocean_time");

    zeta_west = ncfWritable->addVar("zeta_west", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimEtaRho});
    zeta_west.putAtt("long_name", "free-surface western boundary conditions");
    zeta_west.putAtt("units", "meter");
    zeta_west.putAtt("field", "zeta_west, scalar, series");
    zeta_west.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    zeta_west.putAtt("time", "ocean_time");

    zeta_east = ncfWritable->addVar("zeta_east", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimEtaRho});
    zeta_east.putAtt("long_name", "free-surface eastern boundary conditions");
    zeta_east.putAtt("units", "meter");
    zeta_east.putAtt("field", "zeta_east, scalar, series");
    zeta_east.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    zeta_east.putAtt("time", "ocean_time");

    zeta_south = ncfWritable->addVar("zeta_south", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimXiRho});
    zeta_south.putAtt("long_name", "free-surface southern boundary conditions");
    zeta_south.putAtt("units", "meter");
    zeta_south.putAtt("field", "zeta_south, scalar, series");
    zeta_south.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    zeta_south.putAtt("time", "ocean_time");

    zeta_north = ncfWritable->addVar("zeta_north", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimXiRho});
    zeta_north.putAtt("long_name", "free-surface northern boundary conditions");
    zeta_north.putAtt("units", "meter");
    zeta_north.putAtt("field", "zeta_north, scalar, series");
    zeta_north.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    zeta_north.putAtt("time", "ocean_time");

    u_west = ncfWritable->addVar("u_west", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimEtaU});
    u_west.putAtt("long_name", "3D U-momentum western boundary conditions");
    u_west.putAtt("units", "meter second-1");
    u_west.putAtt("field", "u_west, scalar, series");
    u_west.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    u_west.putAtt("time", "ocean_time");

    u_east = ncfWritable->addVar("u_east", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimEtaU});
    u_east.putAtt("long_name", "3D U-momentum eastern boundary conditions");
    u_east.putAtt("units", "meter second-1");
    u_east.putAtt("field", "u_east, scalar, series");
    u_east.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    u_east.putAtt("time", "ocean_time");

    u_south = ncfWritable->addVar("u_south", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimXiU});
    u_south.putAtt("long_name", "3D U-momentum southern boundary conditions");
    u_south.putAtt("units", "meter second-1");
    u_south.putAtt("field", "u_south, scalar, series");
    u_south.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    u_south.putAtt("time", "ocean_time");

    u_north = ncfWritable->addVar("u_north", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimXiU});
    u_north.putAtt("long_name", "3D U-momentum northern boundary conditions");
    u_north.putAtt("units", "meter second-1");
    u_north.putAtt("field", "u_north, scalar, series");
    u_north.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    u_north.putAtt("time", "ocean_time");

    v_west = ncfWritable->addVar("v_west", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimEtaV});
    v_west.putAtt("long_name", "3D V-momentum western boundary conditions");
    v_west.putAtt("units", "meter second-1");
    v_west.putAtt("field", "v_west, scalar, series");
    v_west.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    v_west.putAtt("time", "ocean_time");

    v_east = ncfWritable->addVar("v_east", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimEtaV});
    v_east.putAtt("long_name", "3D V-momentum eastern boundary conditions");
    v_east.putAtt("units", "meter second-1");
    v_east.putAtt("field", "v_east, scalar, series");
    v_east.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    v_east.putAtt("time", "ocean_time");

    v_south = ncfWritable->addVar("v_south", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimXiV});
    v_south.putAtt("long_name", "3D V-momentum southern boundary conditions");
    v_south.putAtt("units", "meter second-1");
    v_south.putAtt("field", "v_south, scalar, series");
    v_south.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    v_south.putAtt("time", "ocean_time");

    v_north = ncfWritable->addVar("v_north", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimSRho, dimXiV});
    v_north.putAtt("long_name", "3D V-momentum northern boundary conditions");
    v_north.putAtt("units", "meter second-1");
    v_north.putAtt("field", "v_north, scalar, series");
    v_north.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    v_north.putAtt("time", "ocean_time");

    vbar_west = ncfWritable->addVar("vbar_west", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimEtaV});
    vbar_west.putAtt("long_name", "2D V-momentum western boundary conditions");
    vbar_west.putAtt("units", "meter second-1");
    vbar_west.putAtt("field", "vbar_west, scalar, series");
    vbar_west.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    vbar_west.putAtt("time", "ocean_time");

    vbar_east = ncfWritable->addVar("vbar_east", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimEtaV});
    vbar_east.putAtt("long_name", "2D V-momentum eastern boundary conditions");
    vbar_east.putAtt("units", "meter second-1");
    vbar_east.putAtt("field", "vbar_east, scalar, series");
    vbar_east.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    vbar_east.putAtt("time", "ocean_time");

    vbar_south = ncfWritable->addVar("vbar_south", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimXiV});
    vbar_south.putAtt("long_name", "2D V-momentum southern boundary conditions");
    vbar_south.putAtt("units", "meter second-1");
    vbar_south.putAtt("field", "vbar_south, scalar, series");
    vbar_south.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    vbar_south.putAtt("time", "ocean_time");

    vbar_north = ncfWritable->addVar("vbar_north", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimXiV});
    vbar_north.putAtt("long_name", "2D V-momentum northern boundary conditions");
    vbar_north.putAtt("units", "meter second-1");
    vbar_north.putAtt("field", "vbar_north, scalar, series");
    vbar_north.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    vbar_north.putAtt("time", "ocean_time");

    ubar_west = ncfWritable->addVar("ubar_west", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimEtaU});
    ubar_west.putAtt("long_name", "2D U-momentum western boundary conditions");
    ubar_west.putAtt("units", "meter second-1");
    ubar_west.putAtt("field", "ubar_west, scalar, series");
    ubar_west.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    ubar_west.putAtt("time", "ocean_time");

    ubar_east = ncfWritable->addVar("ubar_east", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimEtaU});
    ubar_east.putAtt("long_name", "2D U-momentum eastern boundary conditions");
    ubar_east.putAtt("units", "meter second-1");
    ubar_east.putAtt("field", "ubar_east, scalar, series");
    ubar_east.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    ubar_east.putAtt("time", "ocean_time");

    ubar_south = ncfWritable->addVar("ubar_south", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimXiU});
    ubar_south.putAtt("long_name", "2D U-momentum southern boundary conditions");
    ubar_south.putAtt("units", "meter second-1");
    ubar_south.putAtt("field", "ubar_south, scalar, series");
    ubar_south.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    ubar_south.putAtt("time", "ocean_time");

    ubar_north = ncfWritable->addVar("ubar_north", netCDF::ncDouble, std::initializer_list<netCDF::NcDim>{dimOceanTime, dimXiU});
    ubar_north.putAtt("long_name", "2D U-momentum northern boundary conditions");
    ubar_north.putAtt("units", "meter second-1");
    ubar_north.putAtt("field", "ubar_north, scalar, series");
    ubar_north.putAtt("missing_value", netCDF::ncFloat, 1, &missing_value);
    ubar_north.putAtt("time", "ocean_time");
}

ROMSBoundary::~ROMSBoundary() {
    if (ncfWritable) {
        ncfWritable->close();
        delete ncfWritable;
    }
}

void ROMSBoundary::make() {
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

    sc_r.putVar(&SC_R.front());
    s_rho.putVar(&S_RHO.front());
    hc.putVar(&HC.front());
    Cs_r.putVar(&CS_R.front());
}

void ROMSBoundary::write(size_t time) {
    const size_t bufferSize_eta_rho = dimSRho.getSize() * dimEtaRho.getSize();
    const size_t bufferSize_xi_rho = dimSRho.getSize() * dimXiRho.getSize();
    const size_t bufferSize_eta_u = dimSRho.getSize() * dimEtaU.getSize();
    const size_t bufferSize_xi_u = dimSRho.getSize() * dimXiU.getSize();
    const size_t bufferSize_eta_v = dimSRho.getSize() * dimEtaV.getSize();
    const size_t bufferSize_xi_v = dimSRho.getSize() * dimXiV.getSize();

    std::vector<double> TempWestBuffer(bufferSize_eta_rho);
    std::vector<double> TempEastBuffer(bufferSize_eta_rho);
    std::vector<double> TempSouthBuffer(bufferSize_xi_rho);
    std::vector<double> TempNorthBuffer(bufferSize_xi_rho);

    std::vector<double> SaltWestBuffer(bufferSize_eta_rho);
    std::vector<double> SaltEastBuffer(bufferSize_eta_rho);
    std::vector<double> SaltSouthBuffer(bufferSize_xi_rho);
    std::vector<double> SaltNorthBuffer(bufferSize_xi_rho);

    std::vector<double> UWestBuffer(bufferSize_eta_u);
    std::vector<double> UEastBuffer(bufferSize_eta_u);
    std::vector<double> USouthBuffer(bufferSize_xi_u);
    std::vector<double> UNorthBuffer(bufferSize_xi_u);

    std::vector<double> VWestBuffer(bufferSize_eta_v);
    std::vector<double> VEastBuffer(bufferSize_eta_v);
    std::vector<double> VSouthBuffer(bufferSize_xi_v);
    std::vector<double> VNorthBuffer(bufferSize_xi_v);

    std::vector<double> ZetaWestBuffer(dimEtaRho.getSize());
    std::vector<double> ZetaEastBuffer(dimEtaRho.getSize());
    std::vector<double> ZetaSouthBuffer(dimXiRho.getSize());
    std::vector<double> ZetaNorthBuffer(dimXiRho.getSize());

    std::vector<double> UBarWestBuffer(dimEtaU.getSize());
    std::vector<double> UBarEastBuffer(dimEtaU.getSize());
    std::vector<double> UBarSouthBuffer(dimXiU.getSize());
    std::vector<double> UBarNorthBuffer(dimXiU.getSize());

    std::vector<double> VBarWestBuffer(dimEtaV.getSize());
    std::vector<double> VBarEastBuffer(dimEtaV.getSize());
    std::vector<double> VBarSouthBuffer(dimXiV.getSize());
    std::vector<double> VBarNorthBuffer(dimXiV.getSize());

    for (int k = 0; k < dimSRho.getSize(); ++k) {
        for (int j = 0; j < dimEtaRho.getSize(); ++j) {
            size_t index = k * dimEtaRho.getSize() + j;

            TempWestBuffer[index] = TEMP[k][j][0];
            SaltWestBuffer[index] = SALT[k][j][0];
            TempEastBuffer[index] = TEMP[k][j][dimXiRho.getSize()-1];
            SaltEastBuffer[index] = SALT[k][j][dimXiRho.getSize()-1];
        }

        for (int i = 0; i < dimXiRho.getSize(); ++i) {
            size_t index = k * dimXiRho.getSize() + i;

            TempSouthBuffer[index] = TEMP[k][0][i];
            SaltSouthBuffer[index] = SALT[k][0][i];
            TempNorthBuffer[index] = TEMP[k][dimEtaRho.getSize()-1][i];
            SaltNorthBuffer[index] = SALT[k][dimEtaRho.getSize()-1][i];
        }

        for (int j = 0; j < dimEtaU.getSize(); ++j) {
            size_t index = k * dimEtaU.getSize() + j;

            UWestBuffer[index] = U[k][j][0];
            UEastBuffer[index] = U[k][j][dimXiU.getSize()-1];
        }

        for (int i = 0; i < dimXiU.getSize(); ++i) {
            size_t index = k * dimXiU.getSize() + i;

            USouthBuffer[index] = U[k][0][i];
            UNorthBuffer[index] = U[k][dimEtaU.getSize()-1][i];
        }

        for (int j = 0; j < dimEtaV.getSize(); ++j) {
            size_t index = k * dimEtaV.getSize() + j;

            VWestBuffer[index] = V[k][j][0];
            VEastBuffer[index] = V[k][j][dimXiV.getSize()-1];
        }

        for (int i = 0; i < dimXiV.getSize(); ++i) {
            size_t index = k * dimXiV.getSize() + i;

            VSouthBuffer[index] = V[k][0][i];
            VNorthBuffer[index] = V[k][dimEtaV.getSize()-1][i];
        }
    }

    for (int j = 0; j < dimEtaRho.getSize(); j++) {
        ZetaWestBuffer[j] = ZETA[j][0];
        ZetaEastBuffer[j] = ZETA[j][dimXiRho.getSize()-1];
    }

    for (int i = 0; i < dimXiRho.getSize(); i++) {
        ZetaSouthBuffer[i] = ZETA[0][i];
        ZetaNorthBuffer[i] = ZETA[dimEtaRho.getSize()-1][i];
    }

    for (int j = 0; j < dimEtaU.getSize(); ++j) {
        UBarWestBuffer[j] = UBAR[j][0];
        UBarEastBuffer[j] = UBAR[j][dimXiU.getSize()-1];
    }

    for (int i = 0; i < dimXiU.getSize(); ++i) {
        UBarSouthBuffer[i] = UBAR[0][i];
        UBarNorthBuffer[i] = UBAR[dimEtaU.getSize()-1][i];
    }

    for (int j = 0; j < dimEtaV.getSize(); ++j) {
        VBarWestBuffer[j] = VBAR[j][0];
        VBarEastBuffer[j] = VBAR[j][dimXiV.getSize()-1];
    }

    for (int i = 0; i < dimXiV.getSize(); ++i) {
        VBarSouthBuffer[i] = VBAR[0][i];
        VBarNorthBuffer[i] = VBAR[dimEtaV.getSize()-1][i];
    }

    std::vector<size_t> start = {time, 0, 0};
    std::vector<size_t> count = {1, dimSRho.getSize(), dimEtaRho.getSize()};
    temp_west.putVar(start, count, TempWestBuffer.data());
    salt_west.putVar(start, count, SaltWestBuffer.data());
    temp_east.putVar(start, count, TempEastBuffer.data());
    salt_east.putVar(start, count, SaltEastBuffer.data());

    count = {1, dimSRho.getSize(), dimXiRho.getSize()};
    temp_south.putVar(start, count, TempSouthBuffer.data());
    salt_south.putVar(start, count, SaltSouthBuffer.data());
    temp_north.putVar(start, count, TempNorthBuffer.data());
    salt_north.putVar(start, count, SaltNorthBuffer.data());

    count = {1, dimSRho.getSize(), dimEtaU.getSize()};
    u_west.putVar(start, count, UWestBuffer.data());
    u_east.putVar(start, count, UEastBuffer.data());

    count = {1, dimSRho.getSize(), dimXiU.getSize()};
    u_south.putVar(start, count, USouthBuffer.data());
    u_north.putVar(start, count, UNorthBuffer.data());

    count = {1, dimSRho.getSize(), dimEtaV.getSize()};
    v_west.putVar(start, count, VWestBuffer.data());
    v_east.putVar(start, count, VEastBuffer.data());

    count = {1, dimSRho.getSize(), dimXiV.getSize()};
    v_south.putVar(start, count, VSouthBuffer.data());
    v_north.putVar(start, count, VNorthBuffer.data());

    count = {1, dimEtaRho.getSize()};
    zeta_west.putVar(start, count, ZetaWestBuffer.data());
    zeta_east.putVar(start, count, ZetaEastBuffer.data());

    count = {1, dimXiRho.getSize()};
    zeta_south.putVar(start, count, ZetaSouthBuffer.data());
    zeta_north.putVar(start, count, ZetaNorthBuffer.data());

    count = {1, dimEtaU.getSize()};
    ubar_west.putVar(start, count, UBarWestBuffer.data());
    ubar_east.putVar(start, count, UBarEastBuffer.data());

    count = {1, dimXiU.getSize()};
    ubar_south.putVar(start, count, UBarSouthBuffer.data());
    ubar_north.putVar(start, count, UBarNorthBuffer.data());

    count = {1, dimEtaV.getSize()};
    vbar_west.putVar(start, count, VBarWestBuffer.data());
    vbar_east.putVar(start, count, VBarEastBuffer.data());

    count = {1, dimXiV.getSize()};
    vbar_south.putVar(start, count, VBarSouthBuffer.data());
    vbar_north.putVar(start, count, VBarNorthBuffer.data());
}

std::string ROMSBoundary::getCurrentDateTime() const {
    std::time_t now = std::time(nullptr);
    char buf[sizeof "YYYY-MM-DD HH:MM:SS"];
    std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return buf;
}

void ROMSBoundary::setOceanTime(std::vector<double> ocean_time){
    OCEAN_TIME = ocean_time;
}

void ROMSBoundary::setZETA(std::vector<std::vector<double>> zeta) {
    ZETA = zeta;
}

void ROMSBoundary::setSALT(std::vector<std::vector<std::vector<double>>> salt) {
    SALT = salt;
}

void ROMSBoundary::setTEMP(std::vector<std::vector<std::vector<double>>> temp) {
    TEMP = temp;
}

void ROMSBoundary::setU(std::vector<std::vector<std::vector<double>>> u) {
    U = u;
}

void ROMSBoundary::setV(std::vector<std::vector<std::vector<double>>> v) {
    V = v;
}

void ROMSBoundary::setUBAR(std::vector<std::vector<double>> ubar) {
    UBAR = ubar;
}

void ROMSBoundary::setVBAR(std::vector<std::vector<double>> vbar) {
    VBAR = vbar;
}