//
// Created by pablo on 6/8/20.
//

#include "telomereRegion.h"

std::string cmri::telomere_signature_t::serialize() const {
    std::stringstream result;
    result << "{";
    result << "\"count\":" << count;
    result << ",\"histogram\":{";

    for(auto iter_data = histogram.begin(); iter_data != histogram.end(); ++iter_data) {
        result << "\"" << iter_data->first <<"\" : " << cmri::serialize(iter_data->second);
        result << (std::next(iter_data) != histogram.end() ? "," : "");
    }
    result << "}";

    result << ",\"context\":{";
    for(auto iter_data = context.begin(); iter_data != context.end(); ++iter_data) {
        result << "\"" << iter_data->first <<"\" : " << cmri::serialize(iter_data->second);
        result << (std::next(iter_data) != context.end() ? "," : "");
    }
    result << "}";

    result << "}";
    return result.str();
}


void cmri::telomere_signature_t::deserialize(const boost::property_tree::ptree &tree) {
    try {
        count = tree.get<unsigned int>("count");

        for(auto &item : tree.get_child("histogram")){
            for(auto &child_item : item.second){
                histogram[std::stoi(child_item.first)]=child_item.second.get_value<unsigned int>();
            }
        }

        for(auto &item : tree.get_child("context")){
            std::map<std::string, unsigned int> values;
            for(auto &child_item : item.second){
                context[child_item.first]=child_item.second.get_value<unsigned int>();
            }
        }

    }
    catch (std::exception &e) {
        LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
        LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
        exit(EIO);
    }
}



std::string cmri::telomereRegion::serialize() const {
    std::stringstream result;
    result << "{";
    result << "\"start\":" << start;
    result << ",\"end\":" << end;
    result << ",\"name\":\"" << name << "\"";
    result << ",\"total_bases\":" << total_bases;
    result << ",\"variants\":" << cmri::serialize(signature);
    result << "}";
    return result.str();
}

void cmri::telomereRegion::deserialize(const boost::property_tree::ptree &tree) {
    try {
        genomeRegion::deserialize(tree);
        total_bases = tree.get<unsigned int>("total_bases");
        cmri::deserialize(tree.get_child("signature"), signature);
/*
        for(auto &item : tree.get_child("signature")){
            std::map<std::string, telomere_signature_t> values;
            for(auto &child_item : item.second){
                values[std::stoi(child_item.first)]=child_item.second.get_value<unsigned int>();
            }
            histogram[item.first]=values;
        }
*/

    }
    catch (std::exception &e) {
        LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
        LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
        exit(EIO);
    }
}