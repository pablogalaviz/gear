#ifndef GEAR_MOTIFREGION_H
#define GEAR_MOTIFREGION_H
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
#include <genomeRegion.h>
#include "logger.h"
#include <map>
#include <string>
#include <vector>

namespace cmri{

    struct motifRegion : public genomeRegion {


        unsigned int reads_count=0;
        unsigned int total_bases=0;

        std::map<std::string,unsigned int> motifs;
        std::map<std::string,std::map<unsigned int,unsigned int>> motif_quality;
        std::map<std::string,unsigned int> regex;

        std::string serialize() const override;
         void deserialize(const boost::property_tree::ptree &tree);

        inline void resetCount(){
            for(auto &m : motifs){m.second=0;}
            for(auto &mq : motif_quality){
                for(auto &m : mq.second){m.second=0;}
            }
            for(auto &r : regex){r.second=0;}
            reads_count=0;
            total_bases=0;
        }

        virtual bool operator==(const motifRegion &rhs) const {
            return genomeRegion::operator==(rhs) &&
                   reads_count == rhs.reads_count &&
                   total_bases == rhs.total_bases
                   ;
        }


        void operator+=(motifRegion &rhs)  {
            if(genomeRegion::operator==(rhs)){
               for(auto &m : motifs){m.second+=rhs.motifs[m.first];}
                for(auto &mq : motif_quality){
                    for(auto &m : mq.second){
                        m.second+=rhs.motif_quality[mq.first][m.first];
                    }
                }
                for(auto &r : regex){r.second+=rhs.regex[r.first];}
                reads_count += rhs.reads_count;
                total_bases+= rhs.total_bases;
            }
            else{
                throw std::runtime_error("error regions are not match ");
            }
        }


    };

    inline std::ostream& operator<<(std::ostream& result, const motifRegion& rhs)
    {
        result << rhs.serialize();
        return result;
    }



    typedef std::vector<motifRegion> vectorMotifRegion;
    typedef std::map<std::string,vectorMotifRegion> mapVectorMotifRegion;


}

#endif //GEAR_MOTIFREGION_H
