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

#include "motifRegion.h"


std::string cmri::motifRegion::serialize() const {
    std::stringstream result;
    result << "{";
    result << "\"start\":" << start;
    result << ",\"end\":" << end;
    result << ",\"name\":\"" << name << "\"";
    result << ",\"count\":" << reads_count;
    result << ",\"total_bases\":" << total_bases;

    result << ",\"motifs\":{";
    for(auto iter_data = motifs.begin(); iter_data != motifs.end(); ++iter_data) {
        result << "\"" << iter_data->first <<"\" : " << cmri::serialize(iter_data->second);
        result << (std::next(iter_data) != motifs.end() ? "," : "");
    }
    result << "}";

    result << ",\"regex\":{";
    for(auto iter_data = regex.begin(); iter_data != regex.end(); ++iter_data) {
        result << "\"" << iter_data->first <<"\" : " << cmri::serialize(iter_data->second);
        result << (std::next(iter_data) != regex.end() ? "," : "");
    }
    result << "}";

    result << "}";
    return result.str();
}

void cmri::motifRegion::deserialize(const boost::property_tree::ptree &tree) {
    try {
        genomeRegion::deserialize(tree);
        reads_count = tree.get<unsigned int>("count");
        total_bases = tree.get<unsigned int>("total_bases");

        for(auto &item : tree.get_child("motifs")){
            std::map<unsigned int, unsigned int> values;
            for(auto &child_item : item.second){
                values[std::stoi(child_item.first)]=child_item.second.get_value<unsigned int>();
            }
            motifs[item.first]=values;
        }
        //TODO: try to make use of template
        //cmri::deserialize(tree.get_child("motif_quality"), motif_quality);
        for (auto &mq : motifs) {
            for (int i = 0; i < 100; i++) {
                mq.second[i] = 0;
            }
        }

        for(auto &item : tree.get_child("regex")){
            std::map<unsigned int, unsigned int> values;
            for(auto &child_item : item.second){
                values[std::stoi(child_item.first)]=child_item.second.get_value<unsigned int>();
            }
            regex[item.first]=values;
        }
        //TODO: try to make use of template
        //cmri::deserialize(tree.get_child("motif_quality"), motif_quality);
        for (auto &r : regex) {
            for (int i = 0; i < 100; i++) {
                r.second[i] = 0;
            }
        }


    }
    catch (std::exception &e) {
        LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
        LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
        exit(EIO);
    }
}

