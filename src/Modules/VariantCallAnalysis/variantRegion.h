#ifndef GEAR_VARIANTREGION_H
#define GEAR_VARIANTREGION_H
//
// Author(s) Pablo Galaviz (2020)
// e-mail  <pgalaviz@cmri.org.au>
//



//  This file is part of GEAR
//
//  GEAR is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  any later version.
//
//  GEAR is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GEAR.  If not, see <http://www.gnu.org/licenses/>.
//


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <genomeRegion.h>
#include "logger.h"
#include <map>
#include <string>
#include <vector>

namespace cmri{

    struct variantRegion : public genomeRegion {

        unsigned int total_bases=0;
        std::map<std::string,std::map<std::string,unsigned int>> mutations;


        std::string serialize() const override;
        void deserialize(const boost::property_tree::ptree &tree);

        inline void resetCount(){
            total_bases=0;
        }

        virtual bool operator==(const variantRegion &rhs) const {
            return genomeRegion::operator==(rhs) &&
                    total_bases == rhs.total_bases;
        }


        void operator+=(variantRegion &rhs)  {
            if(genomeRegion::operator==(rhs)){
                total_bases+= rhs.total_bases;
            }
            else{
                throw std::runtime_error("error regions are not match ");
            }
        }


    };

    typedef std::vector<variantRegion> vectorVariantRegion;
    typedef std::map<std::string,vectorVariantRegion> mapVectorVariantRegion;

    inline std::ostream& operator<<(std::ostream& result, const variantRegion& rhs)
    {
        result << rhs.serialize();
        return result;
    }


}



#endif //GEAR_VARIANTREGION_H
