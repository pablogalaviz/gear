#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN  // in only one cpp file

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../../src/Modules/MotifCount/motifCount.h"
#include "options.h"


BOOST_AUTO_TEST_SUITE(motifCountTest)


    BOOST_AUTO_TEST_CASE(searchMotifTest) {

        std::string sequence = "TTAGGGCTTAGGGAAATTAGGGCCCTTAGGGactttAGGGTTaGGGTTaaCCC";
        std::string motif = "TTAGGG";
        int expected=4;
        BOOST_TEST(cmri::searchMotif(sequence, motif) == expected);

    }


    BOOST_AUTO_TEST_CASE(searchMotifWFilterTest) {

        std::string sequence = "TTAGGG"
                               "C"
                               "TTAGGG"
                               "AAA"
                               "TTAGGG"
                               "CCC"
                               "TTAGGG"
                               "ACT"
                               "TTAGGG"
                               "TTAGGG"
                               "TTAACCC";
        std::vector<int> quality={1,1,1,1,1,1
                                  ,7
                                  ,5,5,6,6,6,6
                                  ,3,3,2
                                  ,10,11,10,11,10,10
                                  ,5,3,5
                                  ,1,1,1,1,1,1
                                  ,9,8,9
                                  ,5,6,5,6,7,6
                                  ,13,13,14,14,12,13
                                  ,1,2,3,4,5,6,7};

        std::map<std::string,std::map<unsigned int,unsigned int>> motif_quality = {{"TTAGGG",{
            {1,0}
            ,{5,0}
            ,{6,0}
            ,{10,0}
            ,{13,0}
        }}};
        std::map<std::string,std::map<unsigned int,unsigned int>> expected = {{"TTAGGG",{
                         {1,2}
                        ,{5,2}
                        ,{6,0}
                        ,{10,1}
                        ,{13,1}
                }}};;
        cmri::searchMotifWFilter(sequence,quality,motif_quality);

        BOOST_TEST_MESSAGE("serialized data " << cmri::serialize(motif_quality));

        BOOST_TEST(motif_quality == expected);

    }


    std::string sample_regex_motif[] = {"(TTAGGG)(.{0})TTAGGG", "(TTAGGG)(.{1})TTAGGG", "(TTAGGG)(.{3})TTAGGG"};
    int sample_expected[] = {1, 1, 2};

    BOOST_DATA_TEST_CASE(searchRegexTest,
                         boost::unit_test::data::make(sample_regex_motif) ^ sample_expected
                         , regex_motif, expected) {
        std::string sequence = "TTAGGGCTTAGGGAAATTAGGGCCCTTAGGGACTTTAGGGTTAGGGTTAACCC";

        BOOST_TEST(cmri::searchRegex(sequence, regex_motif) == expected);

    }

    int sample_expected_consecutive[] = {1, 1, 3};
    BOOST_DATA_TEST_CASE(searchRegexConsecutiveTest,
                         boost::unit_test::data::make(sample_regex_motif) ^ sample_expected_consecutive
                         , regex_motif, expected) {
        std::string sequence = "TTAGGGCTTAGGGAAATTAGGGCCCTTAGGGACTTTAGGGTTAGGGTTAACCC";

        BOOST_TEST(cmri::searchRegexConsecutive(sequence, regex_motif) == expected);

    }

    std::string sample_input[] = {"data/input.fq",
                                        "data/input.bam"
    };

    std::string sample_result[] = {"data/result.json",
                                        "data/bam_result.json"
    };


    BOOST_DATA_TEST_CASE(processTest,boost::unit_test::data::make(sample_input)^sample_result,input_file,result_file){

        std::string motif_file = "data/input.json";
        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(motif_file, input_tree);
        std::map<std::string, std::vector<cmri::motifRegion>> motifs;
        cmri::deserialize(input_tree,motifs);

        boost::property_tree::ptree result_tree;
        boost::property_tree::read_json(result_file, result_tree);
        std::map<std::string, std::vector<cmri::motifRegion>> result;
        cmri::deserialize(result_tree,result);


        cmri::common_options_t common_options;
        common_options.input_file=input_file;
        cmri::motif_count_options_t motif_count_options;
        motif_count_options.quality_value=10;
        motif_count_options.quality_map=30;
        cmri::process(common_options,motif_count_options, motifs);

        BOOST_TEST(motifs == result);
    }


    BOOST_DATA_TEST_CASE(processMultiThreadingTest,boost::unit_test::data::make(sample_input)^sample_result,input_file,result_file){

        std::string motif_file = "data/input.json";

        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(motif_file, input_tree);
        std::map<std::string, std::vector<cmri::motifRegion>> motifs;
        cmri::deserialize(input_tree,motifs);

        boost::property_tree::ptree result_tree;
        boost::property_tree::read_json(result_file, result_tree);
        std::map<std::string, std::vector<cmri::motifRegion>> result;
        cmri::deserialize(result_tree,result);


        cmri::common_options_t common_options;
        common_options.input_file=input_file;
        common_options.threads = static_cast<int>(std::thread::hardware_concurrency());
        common_options.chunk_size = 1;

        cmri::motif_count_options_t motif_count_options;
        motif_count_options.quality_value=10;
        motif_count_options.quality_map=30;

        cmri::processMultiThreading(common_options,motif_count_options, motifs);

        BOOST_TEST(motifs == result);
    }


BOOST_AUTO_TEST_SUITE_END()