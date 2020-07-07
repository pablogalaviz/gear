#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "../../src/Modules/MotifCount/motifRegion.h"


BOOST_AUTO_TEST_SUITE(motifRegionTest)


    BOOST_AUTO_TEST_CASE(testMotifRegion) {

        std::stringstream input_data;
        input_data << "{\"start\":9995,\"end\":11005,\"name\":\"p_telomere\",\"count\":0,\"total_bases\":0,\"motifs\":{\"CCCCAA\":0,\"CCCCTA\":0,\"CCCTCA\":0,\"TTATGG\":0},\"motif_quality\":{\"CCCCAA\":{},\"CCCCTA\":{},\"CCCTCA\":{},\"TTATGG\":{}},\"regex\":{\"(AGGTCA)(.{0})(AGGTCA)\":0,\"(AGGTCA)(.{0})(TGACCT)\":0,\"(GGGTCA)(.{0})(GGGTCA)\":0,\"(GGGTCA)(.{0})(TGACCC)\":0,\"(TGACCC)(.{0})(GGGTCA)\":0,\"(TGACCC)(.{0})(TGACCC)\":0,\"(TGACCT)(.{0})(AGGTCA)\":0,\"(TGACCT)(.{0})(TGACCT)\":0}}";
        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        cmri::motifRegion motif_region;
        motif_region.deserialize(input_tree);

        std::string serialized_data = motif_region.serialize();
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


    BOOST_AUTO_TEST_CASE(testMotifRegionMap) {

        std::stringstream input_data;
        input_data << "{\"chr1\":[{\"start\":9995,\"end\":11005,\"name\":\"p_telomere\",\"count\":0,\"total_bases\":0,\"motifs\":{\"CCCCAA\":0,\"CCCCTA\":0,\"CCCTCA\":0,\"TTATGG\":0},\"motif_quality\":{\"CCCCAA\":{},\"CCCCTA\":{},\"CCCTCA\":{},\"TTATGG\":{}},\"regex\":{\"(AGGTCA)(.{0})(AGGTCA)\":0,\"(AGGTCA)(.{0})(TGACCT)\":0,\"(GGGTCA)(.{0})(GGGTCA)\":0,\"(GGGTCA)(.{0})(TGACCC)\":0,\"(TGACCC)(.{0})(GGGTCA)\":0,\"(TGACCC)(.{0})(TGACCC)\":0,\"(TGACCT)(.{0})(AGGTCA)\":0,\"(TGACCT)(.{0})(TGACCT)\":0}},{\"start\":248945417,\"end\":248946427,\"name\":\"q_telomere\",\"count\":0,\"total_bases\":0,\"motifs\":{\"CCCCAA\":0,\"CCCCTA\":0,\"CCCTCA\":0,\"TTATGG\":0},\"motif_quality\":{\"CCCCAA\":{},\"CCCCTA\":{},\"CCCTCA\":{},\"TTATGG\":{}},\"regex\":{\"(AGGTCA)(.{0})(AGGTCA)\":0,\"(AGGTCA)(.{0})(TGACCT)\":0,\"(GGGTCA)(.{0})(GGGTCA)\":0,\"(GGGTCA)(.{0})(TGACCC)\":0,\"(TGACCC)(.{0})(GGGTCA)\":0,\"(TGACCC)(.{0})(TGACCC)\":0,\"(TGACCT)(.{0})(AGGTCA)\":0,\"(TGACCT)(.{0})(TGACCT)\":0}}],\"chr10\":[{\"start\":9995,\"end\":11005,\"name\":\"p_telomere\",\"count\":0,\"total_bases\":0,\"motifs\":{\"CCCCAA\":0,\"CCCCTA\":0,\"CCCTCA\":0,\"TTATGG\":0},\"motif_quality\":{\"CCCCAA\":{},\"CCCCTA\":{},\"CCCTCA\":{},\"TTATGG\":{}},\"regex\":{\"(AGGTCA)(.{0})(AGGTCA)\":0,\"(AGGTCA)(.{0})(TGACCT)\":0,\"(GGGTCA)(.{0})(GGGTCA)\":0,\"(GGGTCA)(.{0})(TGACCC)\":0,\"(TGACCC)(.{0})(GGGTCA)\":0,\"(TGACCC)(.{0})(TGACCC)\":0,\"(TGACCT)(.{0})(AGGTCA)\":0,\"(TGACCT)(.{0})(TGACCT)\":0}},{\"start\":133786417,\"end\":133787427,\"name\":\"q_telomere\",\"count\":0,\"total_bases\":0,\"motifs\":{\"CCCCAA\":0,\"CCCCTA\":0,\"CCCTCA\":0,\"TTATGG\":0},\"motif_quality\":{\"CCCCAA\":{},\"CCCCTA\":{},\"CCCTCA\":{},\"TTATGG\":{}},\"regex\":{\"(AGGTCA)(.{0})(AGGTCA)\":0,\"(AGGTCA)(.{0})(TGACCT)\":0,\"(GGGTCA)(.{0})(GGGTCA)\":0,\"(GGGTCA)(.{0})(TGACCC)\":0,\"(TGACCC)(.{0})(GGGTCA)\":0,\"(TGACCC)(.{0})(TGACCC)\":0,\"(TGACCT)(.{0})(AGGTCA)\":0,\"(TGACCT)(.{0})(TGACCT)\":0}}],\"chr11\":[{\"start\":59995,\"end\":61005,\"name\":\"p_telomere\",\"count\":0,\"total_bases\":0,\"motifs\":{\"CCCCAA\":0,\"CCCCTA\":0,\"CCCTCA\":0,\"TTATGG\":0},\"motif_quality\":{\"CCCCAA\":{},\"CCCCTA\":{},\"CCCTCA\":{},\"TTATGG\":{}},\"regex\":{\"(AGGTCA)(.{0})(AGGTCA)\":0,\"(AGGTCA)(.{0})(TGACCT)\":0,\"(GGGTCA)(.{0})(GGGTCA)\":0,\"(GGGTCA)(.{0})(TGACCC)\":0,\"(TGACCC)(.{0})(GGGTCA)\":0,\"(TGACCC)(.{0})(TGACCC)\":0,\"(TGACCT)(.{0})(AGGTCA)\":0,\"(TGACCT)(.{0})(TGACCT)\":0}},{\"start\":135075617,\"end\":135076627,\"name\":\"q_telomere\",\"count\":0,\"total_bases\":0,\"motifs\":{\"CCCCAA\":0,\"CCCCTA\":0,\"CCCTCA\":0,\"TTATGG\":0},\"motif_quality\":{\"CCCCAA\":{},\"CCCCTA\":{},\"CCCTCA\":{},\"TTATGG\":{}},\"regex\":{\"(AGGTCA)(.{0})(AGGTCA)\":0,\"(AGGTCA)(.{0})(TGACCT)\":0,\"(GGGTCA)(.{0})(GGGTCA)\":0,\"(GGGTCA)(.{0})(TGACCC)\":0,\"(TGACCC)(.{0})(GGGTCA)\":0,\"(TGACCC)(.{0})(TGACCC)\":0,\"(TGACCT)(.{0})(AGGTCA)\":0,\"(TGACCT)(.{0})(TGACCT)\":0}}]}";

        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(input_data, input_tree);

        std::map<std::string,std::vector<cmri::motifRegion>> motif_region;
        cmri::deserialize(input_tree, motif_region);

        std::string serialized_data = serialize(motif_region);
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