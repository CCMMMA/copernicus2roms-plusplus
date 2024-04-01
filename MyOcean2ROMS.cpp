#include <iostream>
#include <string>
#include <mpi.h>
#include "jncregridder/roms/ROMSGrid.h"
#include "jncregridder/data/copernicus/CopernicusTem.h"
#include "jncregridder/data/copernicus/CopernicusSSH.h"
#include "jncregridder/data/copernicus/CopernicusCur.h"
#include "jncregridder/data/copernicus/CopernicusSal.h"
#include "jncregridder/util/BilinearInterpolator.h"
#include "jncregridder/util/BilinearInterpolator3D.h"
#include "jncregridder/roms/ROMSInit.h"
#include "jncregridder/roms/ROMSBoundary.h"

class MyOcean2ROMS {
public:
    MyOcean2ROMS(const std::string& gridPath, const std::string& dataPath, const std::string& ncepDate,
                 const std::string& initPath, const std::string& boundaryPath);

    double getModSimStartDate(const std::string &basicString);

    double calculate_julian_date(int year, int month, int day);
};

MyOcean2ROMS::MyOcean2ROMS(const std::string &gridPath, const std::string &dataPath, const std::string &ncepDate,
                           const std::string &initPath, const std::string &boundaryPath) {
    // Path to the ROMS grid
    const std::string& romsGridPath = gridPath;

    // Path to the MyOcean current copernicus-data
    std::string myOceanPathCur = dataPath + "/myoc_d00_" + ncepDate + "_cur.nc";

    // Path to the MyOcean temperature copernicus-data
    std::string myOceanPathTem = dataPath + "/myoc_d00_" + ncepDate + "_tem.nc";

    // Path to the MyOcean salinity copernicus-data
    std::string myOceanPathSal = dataPath + "/myoc_d00_" + ncepDate + "_sal.nc";

    // Path to the MyOcean sea surface height copernicus-data
    std::string myOceanPathSSH = dataPath + "/myoc_d00_" + ncepDate + "_ssh.nc";

    // Path to the output init file
    const std::string& romsInitPath = initPath;

    // Path to the output boundary file
    const std::string& romsBoundaryPath = boundaryPath;

    // Open ROMS grid copernicus-data
    ROMSGrid romsGrid(romsGridPath);

    // Get dimension size
    size_t etaRho = romsGrid.getEtaRho();
    size_t xiRho = romsGrid.getXiRho();
    size_t etaU = romsGrid.getEtaU();
    size_t xiU = romsGrid.getXiU();
    size_t etaV = romsGrid.getEtaV();
    size_t xiV = romsGrid.getXiV();

    std::cout << "Rho:\n\teta:" << etaRho << "\txi:" << xiRho << std::endl;
    std::cout << "U:\n\teta:" << etaU << "\txi:" << xiU << std::endl;
    std::cout << "V:\n\teta:" << etaV << "\txi:" << xiV << std::endl;

    // MASK at rho points
    auto maskRHO = romsGrid.getMASKRHO();
    std::vector<std::vector<int>> MASKRHO;
    MASKRHO.resize(etaRho);
    for (size_t i = 0; i < etaRho; ++i) {
        MASKRHO[i].resize(xiRho);
        for (size_t j = 0; j < xiRho; ++j) {
            size_t index = i * xiRho + j;
            MASKRHO[i][j] = int(maskRHO[index]);
        }
    }

    // MASK at u points
    auto maskU = romsGrid.getMASKU();
    std::vector<std::vector<int>> MASKU;
    MASKU.resize(etaU);
    for (size_t i = 0; i < etaU; ++i) {
        MASKU[i].resize(xiU);
        for (size_t j = 0; j < xiU; ++j) {
            size_t index = i * xiU + j;
            MASKU[i][j] = int(maskU[index]);
        }
    }

    // MASK at v points
    auto maskV = romsGrid.getMASKV();
    std::vector<std::vector<int>> MASKV;
    MASKV.resize(etaV);
    for (size_t i = 0; i < etaV; ++i) {
        MASKV[i].resize(xiV);
        for (size_t j = 0; j < xiV; ++j) {
            size_t index = i * xiV + j;
            MASKV[i][j] = int(maskV[index]);
        }
    }

    // LAT,LON at rho points
    auto latRHO = romsGrid.getLATRHO();
    std::vector<std::vector<double>> LATRHO;
    LATRHO.resize(etaRho);
    for (size_t i = 0; i < etaRho; ++i) {
        LATRHO[i].resize(xiRho);
        for (size_t j = 0; j < xiRho; ++j) {
            size_t index = i * xiRho + j;
            LATRHO[i][j] = latRHO[index];
        }
    }
    auto lonRHO = romsGrid.getLONRHO();
    std::vector<std::vector<double>> LONRHO;
    LONRHO.resize(etaRho);
    for (size_t i = 0; i < etaRho; ++i) {
        LONRHO[i].resize(xiRho);
        for (size_t j = 0; j < xiRho; ++j) {
            size_t index = i * xiRho + j;
            LONRHO[i][j] = lonRHO[index];
        }
    }

    // Z at rho/sigma points
    auto z = romsGrid.getZ();
    std::vector<std::vector<std::vector<double>>> romsZ;
    romsZ.resize(romsGrid.getS_rho().size());
    for (size_t t = 0; t < romsGrid.getS_rho().size(); ++t) {
        romsZ[t].resize(etaRho);
        for (size_t i = 0; i < etaRho; ++i) {
            romsZ[t][i].resize(xiRho);
            for (size_t j = 0; j < xiRho; ++j) {
                size_t index = t * etaRho * xiRho + i * xiRho + j;
                romsZ[t][i][j] = z[index];
            }
        }
    }

    // LAT,LON at u points
    auto latU = romsGrid.getLATU();
    std::vector<std::vector<double>> LATU;
    LATU.resize(etaU);
    for (size_t i = 0; i < etaU; ++i) {
        LATU[i].resize(xiU);
        for (size_t j = 0; j < xiU; ++j) {
            size_t index = i * xiU + j;
            LATU[i][j] = latU[index];
        }
    }
    auto lonU = romsGrid.getLONU();
    std::vector<std::vector<double>> LONU;
    LONU.resize(etaU);
    for (size_t i = 0; i < etaU; ++i) {
        LONU[i].resize(xiU);
        for (size_t j = 0; j < xiU; ++j) {
            size_t index = i * xiU + j;
            LONU[i][j] = lonU[index];
        }
    }

    // LAT,LON at v points
    auto latV = romsGrid.getLATV();
    std::vector<std::vector<double>> LATV;
    LATV.resize(etaV);
    for (size_t i = 0; i < etaV; ++i) {
        LATV[i].resize(xiV);
        for (size_t j = 0; j < xiV; ++j) {
            size_t index = i * xiV + j;
            LATV[i][j] = latV[index];
        }
    }
    auto lonV = romsGrid.getLONV();
    std::vector<std::vector<double>> LONV;
    LONV.resize(etaV);
    for (size_t i = 0; i < etaV; ++i) {
        LONV[i].resize(xiV);
        for (size_t j = 0; j < xiV; ++j) {
            size_t index = i * xiV + j;
            LONV[i][j] = lonV[index];
        }
    }

    // Get the angle between the xi axis and the real east
    auto angle = romsGrid.getANGLE();
    std::vector<std::vector<double>> ANGLE;
    ANGLE.resize(etaRho);
    for (size_t i = 0; i < etaRho; ++i) {
        ANGLE[i].resize(xiRho);
        for (size_t j = 0; j < xiRho; ++j) {
            size_t index = i * xiRho + j;
            ANGLE[i][j] = angle[index];
        }
    }

    std::cout << "MASKRHO:" << MASKRHO.size() << " " << MASKRHO[0].size() << std::endl;
    std::cout << "MASKU:" << MASKU.size() << " " << MASKU[0].size() << std::endl;
    std::cout << "MASKV:" << MASKV.size() << " " << MASKV[0].size() << std::endl;
    std::cout << "LATRHO:" << LATRHO.size() << " " << LATRHO[0].size() << std::endl;
    std::cout << "LONRHO:" << LONRHO.size() << " " << LONRHO[0].size() << std::endl;
    std::cout << "LATU:" << LATU.size() << " " << LATU[0].size() << std::endl;
    std::cout << "LONU:" << LONU.size() << " " << LONU[0].size() << std::endl;
    std::cout << "LATV:" << LATV.size() << " " << LATV[0].size() << std::endl;
    std::cout << "LONV:" << LONV.size() << " " << LONV[0].size() << std::endl;
    std::cout << "ANGLE:" << ANGLE.size() << " " << ANGLE[0].size() << std::endl;
    std::cout << "Z:" << romsZ.size() << " " << romsZ[0].size() << " " << romsZ[0][0].size() << std::endl;

    // Instantiate Copernicus data objects
    CopernicusTem dataTem(myOceanPathTem);
    CopernicusSal dataSal(myOceanPathSal);
    CopernicusSSH dataSSH(myOceanPathSSH);
    CopernicusCur dataCur(myOceanPathCur);

    auto LONXY = dataCur.getLON();
    auto LATXY = dataCur.getLAT();
    auto myOceanZ = dataCur.getZ();

    // Set the number of forcing time steps
    int forcingTimeSteps = dataCur.getTime();

    std::vector<double> oceanTime;
    std::vector<int> scrumTime;
    for (int t = 0; t < forcingTimeSteps; ++t) {
        // Set the value of each ocean time as delta form the simulation starting date
        oceanTime.push_back(getModSimStartDate(ncepDate) + t);
        // Set the scrum time
        scrumTime.push_back(t * 86400);
    }

    // Instantiate a ROMS init file
    ROMSInit romsInit(romsInitPath, romsGrid, forcingTimeSteps);
    romsInit.setOceanTime(oceanTime);
    romsInit.setScrumTime(scrumTime);
    romsInit.make();

    // Instantiate a ROMS boundary file
    ROMSBoundary romsBoundary(romsBoundaryPath, romsGrid, forcingTimeSteps);
    romsBoundary.setOceanTime(oceanTime);
    romsBoundary.make();

    for (int t = 0; t < forcingTimeSteps; ++t) {
        std::cout << "Time: " << t << " " << oceanTime[t] << std::endl;

        // 2D sea surface height
        auto valuesSSH = dataSSH.getZOS()[t];

        // 3D current U component
        auto valuesU = dataCur.getUO()[t];

        // 3D current V component
        auto valuesV = dataCur.getVO()[t];

        // 3D temperature
        auto valuesTem = dataTem.getTHETAO()[t];

        // 3D salinity
        auto valuesSal = dataSal.getSO()[t];

        // Create a 2D bilinear interpolator on Rho points
        BilinearInterpolator bilinearInterpolatorRho(LATXY, LONXY, LATRHO, LONRHO, MASKRHO);

        // Interpolate the SSH
        std::cout << "Interpolating SSH" << std::endl;
        auto SSH_ROMS = bilinearInterpolatorRho.interp(valuesSSH, 1e20, 1e37);

        // Create a 3D bilinear interpolator on Rho points
        BilinearInterpolator3D interpolator3DRho(LATXY, LONXY, myOceanZ, LATRHO, LONRHO, romsZ, MASKRHO, romsGrid);
        // Create a 3D bilinear interpolator on U points
        BilinearInterpolator3D interpolator3DU(LATXY, LONXY, myOceanZ, LATU, LONU, romsZ, MASKU, romsGrid);
        // Create a 3D bilinear interpolator on V points
        BilinearInterpolator3D interpolator3DV(LATXY, LONXY, myOceanZ, LATV, LONV, romsZ, MASKV, romsGrid);

        std::cout << "Interpolating SAL" << std::endl;
        auto SAL_ROMS = interpolator3DRho.interp(valuesSal, 1e20, 1e37);

        std::cout << "Interpolating TEMP" << std::endl;
        auto TEM_ROMS = interpolator3DRho.interp(valuesTem, 1e20, 1e37);

        std::cout << "Interpolating V" << std::endl;
        auto V_ROMS = interpolator3DV.interp(valuesV, 1e20, 1e37);

        std::cout << "Interpolating U" << std::endl;
        auto U_ROMS = interpolator3DU.interp(valuesU, 1e20, 1e37);

        std::cout << "Calculating UBAR" << std::endl;
        std::vector<std::vector<double>> UBAR(etaU, std::vector<double>(xiU, 1e37));
        for (int j = 0; j < etaU; j++) {
            for (int i = 0; i < xiU; i++) {
                if (MASKU[j][i] == 1) {
                    UBAR[j][i] = 0;
                    int count = 0;

                    for (int k=0; k<romsZ.size(); k++) {
                        if (U_ROMS[k][j][i] != 1e37) {
                            UBAR[j][i] += U_ROMS[k][j][i];
                            count++;
                        }
                    }

                    if (count > 0) {
                        UBAR[j][i] = UBAR[j][i]/count;
                    }
                } else {
                    UBAR[j][i] = 1e37;
                }
            }
        }

        std::cout << "Calculating VBAR" << std::endl;
        std::vector<std::vector<double>> VBAR(etaV, std::vector<double>(xiV, 1e37));
        for (int j = 0; j < etaV; j++) {
            for (int i = 0; i < xiV; i++) {
                if (MASKV[j][i] == 1) {
                    VBAR[j][i] = 0;
                    int count = 0;

                    for (int k=0; k<romsZ.size(); k++) {
                        if (V_ROMS[k][j][i] != 1e37) {
                            VBAR[j][i] += V_ROMS[k][j][i];
                            count++;
                        }
                    }

                    if (count > 0) {
                        VBAR[j][i] = VBAR[j][i]/count;
                    }
                } else {
                    VBAR[j][i] = 1e37;
                }
            }
        }

        std::cout << "Time: " << t << ": Saving init file..." << std::endl;
        romsInit.setZETA(SSH_ROMS);
        romsInit.setSALT(SAL_ROMS);
        romsInit.setTEMP(TEM_ROMS);
        romsInit.setV(V_ROMS);
        romsInit.setU(U_ROMS);
        romsInit.setUBAR(UBAR);
        romsInit.setVBAR(VBAR);
        romsInit.write(t);

        std::cout << "Time: " << t << ": Saving bry file..." << std::endl;
        romsBoundary.setZETA(SSH_ROMS);
        romsBoundary.setTEMP(TEM_ROMS);
        romsBoundary.setV(V_ROMS);
        romsBoundary.setU(U_ROMS);
        romsBoundary.setSALT(SAL_ROMS);
        romsBoundary.setUBAR(UBAR);
        romsBoundary.setVBAR(VBAR);
        romsBoundary.write(t);

        // break;
    }
}

double MyOcean2ROMS::getModSimStartDate(const std::string &ncepDate) {
    int year = std::stoi(ncepDate.substr(0, 4));
    int month = std::stoi(ncepDate.substr(4, 2));
    int day = std::stoi(ncepDate.substr(6, 2));

    double dSimStartDate = calculate_julian_date(year, month, day);
    double dModOffset = calculate_julian_date(1968, 5, 23);

    double dModSimStartDate = dSimStartDate - dModOffset;

    return dModSimStartDate;
}

double MyOcean2ROMS::calculate_julian_date(int year, int month, int day) {
    struct tm timeinfo = {};
    timeinfo.tm_year = year - 1900;
    timeinfo.tm_mon = month - 1;
    timeinfo.tm_mday = day;
    time_t rawtime = mktime(&timeinfo);
    return (rawtime / 86400.0) + 2440587.5;
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <gridPath> <dataPath> <ncepDate> <initPath> <boundaryPath>\n";
        return 1;
    }

    std::string gridPath = argv[1];
    std::string dataPath = argv[2];
    std::string ncepDate = argv[3];
    std::string initPath = argv[4];
    std::string boundaryPath = argv[5];

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << rank << " " << size << std::endl;

    /*
    auto start = std::chrono::high_resolution_clock::now();
    MyOcean2ROMS myOcean2Roms(gridPath, dataPath, ncepDate, initPath, boundaryPath);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    */

    return 0;
}