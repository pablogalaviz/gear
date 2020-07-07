#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "genomeRegion.h"


BOOST_AUTO_TEST_SUITE(genomeRegionTest)


    BOOST_AUTO_TEST_CASE(testGenomeRegion) {

        std::stringstream input_data;
        input_data << "{\"start\":1,\"end\":10,\"name\":\"test\"}";
        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        cmri::genomeRegion genome_region;
        genome_region.deserialize(input_tree);

        std::string serialized_data = genome_region.serialize();
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


    BOOST_AUTO_TEST_CASE(testGenomeRegionList1) {

        std::stringstream input_data;
        input_data << "[{\"start\":1,\"end\":10,\"name\":\"element1\"}"
                ",{\"start\":2,\"end\":20,\"name\":\"element2\"}"
                ",{\"start\":3,\"end\":30,\"name\":\"element3\"}]";
        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        std::vector<cmri::genomeRegion> genome_region;
        cmri::deserialize(input_tree,genome_region);
        std::string serialized_data = cmri::serialize(genome_region);
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


    BOOST_AUTO_TEST_CASE(testGenomeRegionList2) {

        std::stringstream input_data;
        input_data << "[{\"item1\":{\"start\":1,\"end\":10,\"name\":\"element1\"},\"item2\":{\"start\":10,\"end\":20,\"name\":\"element2\"},\"item3\":{\"start\":20,\"end\":30,\"name\":\"element3\"}}"
                      ",{\"itemA\":{\"start\":10,\"end\":100,\"name\":\"elementA\"},\"itemB\":{\"start\":100,\"end\":200,\"name\":\"elementB\"}}]";
        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        std::vector<std::map<std::string,cmri::genomeRegion>> genome_region;
        cmri::deserialize(input_tree,genome_region);
        std::string serialized_data = cmri::serialize(genome_region);
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


    BOOST_AUTO_TEST_CASE(testGenomeRegionMap1) {

        std::stringstream input_data;
        input_data << "{\"item1\":{\"start\":1,\"end\":10,\"name\":\"element1\"}"
                ",\"item2\":{\"start\":2,\"end\":20,\"name\":\"element2\"}"
                ",\"item3\":{\"start\":2,\"end\":20,\"name\":\"element2\"}}";

        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        std::map<std::string,cmri::genomeRegion> genome_region;
        cmri::deserialize(input_tree,genome_region);

        std::string serialized_data = serialize(genome_region);
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

    BOOST_AUTO_TEST_CASE(testGenomeRegionMap2) {

        std::stringstream input_data;
        input_data << "{\"item1\":[{\"start\":1,\"end\":10,\"name\":\"element1\"},{\"start\":10,\"end\":20,\"name\":\"element2\"},{\"start\":20,\"end\":30,\"name\":\"element3\"}]"
                ",\"item2\":[{\"start\":1,\"end\":5,\"name\":\"element1\"},{\"start\":5,\"end\":10,\"name\":\"element2\"}]}";

        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        cmri::mapVectorGenomeRegion genome_region;
        cmri::deserialize(input_tree,genome_region);

        std::string serialized_data = serialize(genome_region);
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