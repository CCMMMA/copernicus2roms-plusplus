//
// Created by Ciro De Vita on 28/03/24.
//

#ifndef MYOCEAN2ROMS_COPERNICUSSAL_H
#define MYOCEAN2ROMS_COPERNICUSSAL_H

#include "CopernicusBase.h"

class CopernicusSal : public CopernicusBase {
public:
    CopernicusSal(const std::string& url);
    virtual ~CopernicusSal();

    std::vector<std::vector<std::vector<std::vector<double>>>> getSO() const;

private:
    std::vector<std::vector<std::vector<std::vector<double>>>> SO;
};


#endif //MYOCEAN2ROMS_COPERNICUSSAL_H
