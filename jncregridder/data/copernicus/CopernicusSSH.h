//
// Created by Ciro De Vita on 28/03/24.
//

#ifndef MYOCEAN2ROMS_COPERNICUSSSH_H
#define MYOCEAN2ROMS_COPERNICUSSSH_H

#include "CopernicusBase.h"

class CopernicusSSH : public CopernicusBase {
public:
    CopernicusSSH(const std::string& url);
    virtual ~CopernicusSSH();

    std::vector<std::vector<std::vector<double>>> getZOS() const;

private:
    std::vector<std::vector<std::vector<double>>> ZOS;
};


#endif //MYOCEAN2ROMS_COPERNICUSSSH_H
