#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "../../src/Modules/VariantCallAnalysis/variantRegion.h"


BOOST_AUTO_TEST_SUITE(variantRegionTest)


    BOOST_AUTO_TEST_CASE(testVariantRegion) {

        std::stringstream input_data;
        input_data << "{\"start\":9995,\"end\":11005,\"name\":\"p_telomere\",\"total_bases\":0,\"mutations\":{\"C>A\":{},\"C>G\":{},\"C>T\":{},\"CpG\":{},\"T>A\":{},\"T>C\":{},\"T>G\":{}}}";
        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        cmri::variantRegion variant_region;
        variant_region.deserialize(input_tree);

        std::string serialized_data = variant_region.serialize();
        BOOST_TEST_MESSAGE("serialized data " << serialized_data << "\ninput data      "<< input_data.str());

        std::stringstream output_data;
        output_data << serialized_data;
        boost::property_tree::ptree output_tree;
        boost::property_tree::read_json(output_data, output_tree);

        //test data
        BOOST_TEST(input_tree == output_tree);
        //test format
        BOOST_TEST(serialized_data == input_data.str());

    }


    BOOST_AUTO_TEST_CASE(testVariantRegionMap) {

        std::stringstream input_data;
        input_data << "{\"chr1\":[{\"start\":9995,\"end\":11005,\"name\":\"p_telomere\",\"total_bases\":0,\"mutations\":{\"C>A\":{},\"C>G\":{},\"C>T\":{},\"CpG\":{},\"T>A\":{},\"T>C\":{},\"T>G\":{}}},{\"start\":248945417,\"end\":248946427,\"name\":\"q_telomere\",\"total_bases\":0,\"mutations\":{\"C>A\":{},\"C>G\":{},\"C>T\":{},\"CpG\":{},\"T>A\":{},\"T>C\":{},\"T>G\":{}}},{\"start\":11006,\"end\":248945416,\"name\":\"inner_non_telomeric\",\"total_bases\":0,\"mutations\":{\"C>A\":{},\"C>G\":{},\"C>T\":{},\"CpG\":{},\"T>A\":{},\"T>C\":{},\"T>G\":{}}}]}";

        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        std::map<std::string,std::vector<cmri::variantRegion>> variant_region;
        cmri::deserialize(input_tree,variant_region);

        std::string serialized_data = serialize(variant_region);
        BOOST_TEST_MESSAGE("serialized data " << serialized_data << "\ninput data      "<< input_data.str());

        std::stringstream output_data;
        output_data << serialized_data;
        boost::property_tree::ptree output_tree;
        boost::property_tree::read_json(output_data, output_tree);

        //test data
        BOOST_TEST(input_tree == output_tree);
        //test format
        BOOST_TEST(serialized_data == input_data.str());

    }



BOOST_AUTO_TEST_SUITE_END()