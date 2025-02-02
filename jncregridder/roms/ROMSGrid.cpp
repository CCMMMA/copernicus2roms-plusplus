//
// Created by Ciro De Vita on 28/03/24.
//

#include "ROMSGrid.h"

ROMSGrid::ROMSGrid(const std::string& url) {
    try {
        ncDataset.open(url, netCDF::NcFile::read);

        // Find dimensions
        etaRho = ncDataset.getDim("eta_rho").getSize();
        xiRho = ncDataset.getDim("xi_rho").getSize();
        etaPsi = ncDataset.getDim("eta_psi").getSize();
        xiPsi = ncDataset.getDim("xi_psi").getSize();
        etaU = ncDataset.getDim("eta_u").getSize();
        xiU = ncDataset.getDim("xi_u").getSize();
        etaV = ncDataset.getDim("eta_v").getSize();
        xiV = ncDataset.getDim("xi_v").getSize();

        // Load other variables
        HC = loadVariable("hc");

        if (HC.empty()) {
            HC.resize(1);
            HC[0] = 5.0;

            TCLINE.resize(1);
            TCLINE[0] = 25.0;

            theta_s.resize(1);
            theta_s[0] = 3.0;

            theta_b.resize(1);
            theta_b[0] = 0.0;

            int N = 30;
            double ds = 1.0 / N;
            std::vector<double> lev(N);

            for (int i = 0; i < lev.size(); i++) {
                lev[i] = (1 + i) - 0.5;
            }

            s_rho.resize(lev.size());
            for (int i = 0; i < s_rho.size(); i++) {
                s_rho[i] = (lev[i] - N) * ds;
            }

            cs_r.resize(lev.size());
            s_w.resize(lev.size() + 1);
            s_w[0] = -1;

            cs_w.resize(lev.size() + 1);
            cs_w[0] = -1;

            double cff1 = (1 / sinh(theta_s[0]));
            double cff2 = (0.5 / tanh(0.5 * theta_s[0]));

            if (theta_s[0] > 0) {
                std::vector<double> pTheta(s_rho.size());
                std::vector<double> rTheta(s_rho.size());

                for (int i = 0; i < pTheta.size(); i++) {
                    pTheta[i] = sinh(theta_s[0] * s_rho[i]) * cff1;
                    rTheta[i] = tanh(theta_s[0] * (s_rho[i] + 0.5)) / (2.0 * tanh(0.5 * theta_s[0]) - 0.5);
                    cs_r[i] = (1.0 - theta_b[0]) * pTheta[i] + theta_b[0] * rTheta[i];
                    cs_w[i + 1] = (1.0 - theta_b[0]) * cff1 * sinh(theta_s[0] * s_w[i + 1]) + theta_b[0] * (cff2 * tanh(theta_s[0] * (s_w[i + 1] + 0.5)) - 0.5);
                }
            } else {
                for (int i = 0; i < s_rho.size(); i++) {
                    cs_r[i] = s_rho[i];
                    cs_w[i + 1] = s_w[i + 1];
                }
            }
        } else {
            theta_s = loadVariable("theta_s");
            theta_b = loadVariable("theta_b");
            s_rho = loadVariable("s_rho");
            cs_r = loadVariable("Cs_r");
            s_w = loadVariable("s_w");
            cs_w = loadVariable("Cs_w");
            TCLINE = loadVariable("Tcline");
        }

        ANGLE = loadVariable("angle");
        LATRHO = loadVariable("lat_rho");
        LONRHO = loadVariable("lon_rho");
        LATPSI = loadVariable("lat_psi");
        LONPSI = loadVariable("lon_psi");
        LATU = loadVariable("lat_u");
        LONU = loadVariable("lon_u");
        LATV = loadVariable("lat_v");
        LONV = loadVariable("lon_v");
        H = loadVariable("h");
        MASKRHO = loadVariable("mask_rho");
        MASKU = loadVariable("mask_u");
        MASKV = loadVariable("mask_v");
        Z = loadVariable("z");
    } catch (const netCDF::exceptions::NcException& e) {
        std::cerr << "NetCDF error: " << e.what() << std::endl;
    }
}

std::vector<double> ROMSGrid::loadVariable(const std::string& variable_name) {
    std::vector<double> variable_data;
    try {
        if (variable_name == "z") {
            std::vector<double> z(s_rho.size() * etaRho * xiRho, 0.0);

            double S = 0;
            size_t index = 0;
            for (size_t k = 0; k < s_rho.size(); ++k) {
                for (size_t j = 0; j < etaRho; ++j) {
                    for (size_t i = 0; i < xiRho; ++i) {
                        index = j * xiRho + i;
                        S = HC.front() * s_rho[k] + ((H[index] - HC.front()) * cs_r[k]);
                        z[k * etaRho * xiRho + j * xiRho + i] = S + zeta * (1 + (S / H[index]));
                    }
                }
            }
            variable_data = z;
        } else {
            auto variable = ncDataset.getVar(variable_name);
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
        }
    } catch (const netCDF::exceptions::NcException& e) {
        std::cerr << "NetCDF error: " << e.what() << ": " << variable_name << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error loading variable " << variable_name << ": " << e.what() << std::endl;
    }
    return variable_data;
}

std::string ROMSGrid::getURL() const {
    return url;
}

size_t ROMSGrid::getEtaRho() const {
    return etaRho;
}

size_t ROMSGrid::getXiRho() const {
    return xiRho;
}

size_t ROMSGrid::getEtaPsi() const {
    return etaPsi;
}

size_t ROMSGrid::getXiPsi() const {
    return xiPsi;
}

size_t ROMSGrid::getEtaU() const {
    return etaU;
}

size_t ROMSGrid::getXiU() const {
    return xiU;
}

size_t ROMSGrid::getEtaV() const {
    return etaV;
}

size_t ROMSGrid::getXiV() const {
    return xiV;
}

double ROMSGrid::getNoData() const {
    return no_data;
}

double ROMSGrid::getZeta() const {
    return zeta;
}

std::vector<double> ROMSGrid::getS_rho() const {
    return s_rho;
}

std::vector<double> ROMSGrid::getCs_r() const {
    return cs_r;
}

std::vector<double> ROMSGrid::getS_w() const {
    return s_w;
}

std::vector<double> ROMSGrid::getCs_w() const {
    return cs_w;
}

std::vector<double> ROMSGrid::getTheta_s() const {
    return theta_s;
}

std::vector<double> ROMSGrid::getTheta_b() const {
    return theta_b;
}

std::vector<double> ROMSGrid::getANGLE() const {
    return ANGLE;
}

std::vector<double> ROMSGrid::getLATRHO() const {
    return LATRHO;
}

std::vector<double> ROMSGrid::getLONRHO() const {
    return LONRHO;
}

std::vector<double> ROMSGrid::getLATPSI() const {
    return LATPSI;
}

std::vector<double> ROMSGrid::getLONPSI() const {
    return LONPSI;
}

std::vector<double> ROMSGrid::getLATU() const {
    return LATU;
}

std::vector<double> ROMSGrid::getLONU() const {
    return LONU;
}

std::vector<double> ROMSGrid::getLATV() const {
    return LATV;
}

std::vector<double> ROMSGrid::getLONV() const {
    return LONV;
}

std::vector<double> ROMSGrid::getH() const {
    return H;
}

std::vector<double> ROMSGrid::getHC() const {
    return HC;
}

std::vector<double> ROMSGrid::getMASKRHO() const {
    return MASKRHO;
}

std::vector<double> ROMSGrid::getMASKU() const {
    return MASKU;
}

std::vector<double> ROMSGrid::getMASKV() const {
    return MASKV;
}

std::vector<double> ROMSGrid::getZ() const {
    return Z;
}

std::vector<double> ROMSGrid::getTCLINE() const {
    return TCLINE;
}

void ROMSGrid::setNoData(double value) {
    no_data = value;
}

void ROMSGrid::setZeta(double value) {
    zeta = value;
}