//
// Created by pablo on 2/7/20.
//

#include "variantRegion.h"

std::string cmri::variantRegion::serialize() const {

    std::stringstream result;
    result << "{";
    result << "\"start\":" << start;
    result << ",\"end\":" << end;
    result << ",\"name\":\"" << name << "\"";
    result << ",\"total_bases\":" << total_bases;
    result << ",\"mutations\":" << cmri::serialize(mutations);

    result << "}";

    return result.str();

    return genomeRegion::serialize();
}

void cmri::variantRegion::deserialize(const boost::property_tree::ptree &tree) {
    try {
        genomeRegion::deserialize(tree);
        total_bases = tree.get<unsigned int>("total_bases");
        cmri::deserialize(tree.get_child("mutations"),mutations);
        //mutations = tree2map<unsigned int>(tree,"mutations");

    }
    catch (std::exception &e) {
        LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
        LOGGER.error << "Exception parsing json file: " << e.what() << std::endl;
        exit(EIO);
    }}
