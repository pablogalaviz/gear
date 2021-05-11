#ifndef GEAR_GENOMEREGION_H
#define GEAR_GENOMEREGION_H
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
#include <boost/lexical_cast.hpp>
#include <string>
#include <sstream>
#include "utils.h"

namespace cmri {

    struct genomeRegion {

        unsigned int start = 0;
        unsigned int end = 0;
        std::string name;

        inline bool intersect(const unsigned int location) const {
            return (start <= location && location <= end) || start==end;
        }


        inline bool intersect(const unsigned int other_start, const unsigned int other_end) const {
            if (start == end) { return true; }
            double l_other = static_cast<double>(other_start - start) / (end - start);
            double r_other = static_cast<double>(other_end - start) / (end - start);
            return l_other >= 0 && r_other <= 1; //fully inside region
            // any other case is outside.
        }



        virtual bool operator==(const genomeRegion &rhs) const {
            return start == rhs.start &&
                   end == rhs.end &&
                   name == rhs.name;
        }

        virtual bool operator!=(const genomeRegion &rhs) const {
            return !(rhs == *this);
        }

        virtual std::string serialize() const {
            std::stringstream result;
            result << "{";
            result << "\"start\":" << start << ",";
            result << "\"end\":" << end << ",";
            result << "\"name\":\"" << name << "\"";
            result << "}";
            return result.str();
        };

        virtual void deserialize(const boost::property_tree::ptree &tree) {
            try {
                start = tree.get<unsigned int>("start");
                end = tree.get<unsigned int>("end");
                name = tree.get<std::string>("name");
            }
            catch (std::exception &e) {
                LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
                LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
                exit(EIO);
            }
        }

    };

    inline std::ostream& operator<<(std::ostream& result, const genomeRegion& rhs)
    {
        result << rhs.serialize();
        return result;
    }




    template<class T>
    void deserialize(const boost::property_tree::ptree &tree,T &genome_region) {
        genome_region.deserialize(tree);
    }

    template<>
inline void deserialize(const boost::property_tree::ptree &tree,unsigned int &value) {
        value = tree.get_value<unsigned int>();
    }

    typedef std::vector<genomeRegion> vectorGenomeRegion;


    template<class T>
    void deserialize( const boost::property_tree::ptree &tree, std::vector<T> & result) {
            try {
                for(auto &t : tree) {
                    T new_item;
                    deserialize(t.second,new_item);
                    result.push_back(new_item);
                }
            }
            catch (std::exception &e) {
                LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
                LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
                exit(EIO);
            }

            }




    template<class T>
         void deserialize( const boost::property_tree::ptree &tree,std::map<std::string, T> &result) {

        try {
                for(auto &t : tree) {
                    T new_item;
                    deserialize(t.second,new_item);
                    result[t.first]=new_item;
                }
            }
            catch (std::exception &e) {
                LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
                LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
                exit(EIO);
            }


        }


    typedef std::vector<genomeRegion> vectorGenomeRegion;
    typedef std::map<std::string,vectorGenomeRegion> mapVectorGenomeRegion;


    template<class T>
    static std::map<std::string, T> tree2map(const boost::property_tree::ptree &tree,const std::string& name){

        std::map<std::string, T> result;
        try {
            for(auto &item : tree.get_child(name)){
                result[item.first]=boost::lexical_cast<T>(item.second.data());
            }
        }
        catch (std::exception &e) {
            LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
            LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
            exit(EIO);
        }

        return result;

    }



    }


#endif //GEAR_GENOMEREGION_H
