//
// Created by Ciro De Vita on 28/03/24.
//

#ifndef MYOCEAN2ROMS_COPERNICUSTEM_H
#define MYOCEAN2ROMS_COPERNICUSTEM_H

#include "CopernicusBase.h"

class CopernicusTem : public CopernicusBase {
public:
    CopernicusTem(const std::string& url);
    ~CopernicusTem();

    std::vector<std::vector<std::vector<std::vector<double>>>> getTHETAO() const;
    std::vector<std::vector<std::vector<double>>> getBOTTOMT() const;

private:
    std::vector<std::vector<std::vector<std::vector<double>>>> THETAO;
    std::vector<std::vector<std::vector<double>>> BOTTOMT;
};


#endif //MYOCEAN2ROMS_COPERNICUSTEM_H
