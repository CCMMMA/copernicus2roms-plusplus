//
// Created by Ciro De Vita on 28/03/24.
//

#ifndef MYOCEAN2ROMS_COPERNICUSCUR_H
#define MYOCEAN2ROMS_COPERNICUSCUR_H

#include "CopernicusBase.h"

class CopernicusCur : public CopernicusBase {
public:
    CopernicusCur(const std::string& url);
    virtual ~CopernicusCur();

    std::vector<std::vector<std::vector<std::vector<double>>>> getVO() const;
    std::vector<std::vector<std::vector<std::vector<double>>>> getUO() const;

private:
    std::vector<std::vector<std::vector<std::vector<double>>>> VO;
    std::vector<std::vector<std::vector<std::vector<double>>>> UO;
};


#endif //MYOCEAN2ROMS_COPERNICUSCUR_H
