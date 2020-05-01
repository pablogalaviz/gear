#ifndef GEAR_MOTIF_H
#define GEAR_MOTIF_H
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

#include <string>
#include <map>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>
#include "logger.h"

namespace cmri{

    struct region_t {
        unsigned int start=0;
        unsigned int end=0;
        unsigned int count=0;
        std::string name;
        std::map<std::string,unsigned int> motifs;
        std::map<std::string,unsigned int> regex;

        inline bool intersect(const unsigned int other_start,const unsigned int other_end ){
            if(start==end){return true;}
            double l_other = static_cast<double>(other_start-start)/(end-start);
            double r_other = static_cast<double>(other_end-start)/(end-start);
            if(l_other >= 0 && r_other <= 1 ){return true;} //fully inside region
            return false; // any other case is outside.
        }

        inline void resetCount(){
            for(auto &m : motifs){m.second=0;}
            for(auto &r : regex){r.second=0;}
            count=0;
        }

        bool operator==(const region_t &rhs) const {
            return start == rhs.start &&
                   end == rhs.end &&
                   name == rhs.name &&
                   count == rhs.count;
        }

        bool operator!=(const region_t &rhs) const {
            return !(rhs == *this);
        }

        void operator+=(region_t &rhs)  {
            if(start == rhs.start &&
               end == rhs.end &&
               name == rhs.name){
               for(auto &m : motifs){m.second+=rhs.motifs[m.first];}
                for(auto &r : regex){r.second+=rhs.regex[r.first];}
                count += rhs.count;
            }
            else{
                throw std::runtime_error("error regions are not match ");
            }
        }


    };


typedef std::vector<region_t> region_list_t;
typedef std::map<std::string,region_list_t> chromosome_map_t;



    std::map<std::string, unsigned int> get_patterns(const boost::property_tree::ptree &tree,std::string name){
        std::map<std::string, unsigned int> result;
        try {
            for(auto &item : tree.get_child(name)){
                result[item.first]=boost::lexical_cast<unsigned int>(item.second.data());
            }
        }
        catch (std::exception &e) {
            LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
            LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
            exit(EIO);
        }

        return result;


    }


    region_list_t get_region_list(const boost::property_tree::ptree &tree){
        region_list_t result;
        try {
            for(auto &item : tree){
                region_t new_region;
                new_region.start = item.second.get<unsigned int>("start");
                new_region.end = item.second.get<unsigned int>("end");
                new_region.count = item.second.get<unsigned int>("count");
                new_region.name = item.second.get<std::string>("name");
                new_region.motifs = get_patterns(item.second,"motifs");
                new_region.regex = get_patterns(item.second,"regex");
                result.emplace_back(new_region);
            }
        }
        catch (std::exception &e) {
            LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
            LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
            exit(EIO);
        }

        return result;
    }


    chromosome_map_t deserialize(const std::string &file) {

        boost::property_tree::ptree motif_tree;
        try {
            boost::property_tree::read_json(file, motif_tree);
        }
        catch (std::exception &e) {
            LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
            LOGGER.error << "Exception parsing motif file: " << e.what() << std::endl;
            exit(EIO);
        }

        chromosome_map_t result;
        for(auto &item : motif_tree){
            result[item.first]= get_region_list(item.second);
        }


        return result;


    }


    void serialize(const std::string &file_name, const chromosome_map_t &data){

        std::ofstream file(file_name);
        file << "{";
        for(auto iter_data = data.begin(); iter_data != data.end(); ++iter_data){
            file << "\"" << iter_data->first << "\":[";
            for(int i =0; i < iter_data->second.size(); i++){
                auto v = iter_data->second[i];
                file << "{";
                file << "\"start\":" << v.start <<"," ;
                file << "\"end\":" << v.end <<",";
                file << "\"count\":" << v.count <<",";
                file << "\"name\":\"" << v.name <<"\"," ;
                file << "\"motifs\": {" ;
                for(auto motif_iter = v.motifs.begin(); motif_iter != v.motifs.end(); ++motif_iter){
                    file << "\""<< motif_iter->first <<"\":" << motif_iter->second << (std::next(motif_iter) != v.motifs.end()? ",":"") ;
                }
                file << "},\"regex\": {";
                for(auto regex_iter = v.regex.begin(); regex_iter != v.regex.end(); ++regex_iter){
                    file << "\"" << regex_iter->first << "\":" << regex_iter->second << (std::next(regex_iter) != v.regex.end() ? "," : "") ;
                }

                file << "}}" << (i < iter_data->second.size() - 1 ? "," : "");
            }
            file << "]" << (std::next(iter_data) != data.end() ? "," : "");
        }
        file << "}";

    }



}

#endif //GEAR_MOTIF_H
