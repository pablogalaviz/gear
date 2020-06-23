#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN  // in only one cpp file

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "motifCount.h"


BOOST_AUTO_TEST_SUITE(motifCountTest)


    BOOST_DATA_TEST_CASE(searchMotifTest,
                         boost::unit_test::data::make({true, false}) ^ boost::unit_test::data::make({6, 4}), validate,
                         expected) {

        std::string sequence = "TTAGGGCTTAGGGAAATTAGGGCCCTTAGGGactttAGGGTTaGGGTTaaCCC";
        std::string motif = "TTAGGG";
        BOOST_TEST(cmri::searchMotif(sequence, motif, validate) == expected);

    }

    std::string sample_regex_motif[] = {"(TTAGGG)(.{0})TTAGGG", "(TTAGGG)(.{1})TTAGGG", "(TTAGGG)(.{3})TTAGGG",
                                        "(TTAGGG)(.{0})TTAGGG", "(TTAGGG)(.{1})TTAGGG", "(TTAGGG)(.{3})TTAGGG"};
    bool sample_validate[] = {true, true, true, false, false, false};
    int sample_expected[] = {1, 1, 2, 0, 1, 1};


    BOOST_DATA_TEST_CASE(searchRegexTest,
                         boost::unit_test::data::make(sample_regex_motif) ^ sample_expected ^
                         sample_validate, regex_motif, expected, validate) {
        std::string sequence = "TTAGGGCTTAGGGAAATTAGGGCCCTTAGGGactttAGGGTTaGGGTTaaCCC";

        BOOST_TEST(cmri::searchRegex(sequence, regex_motif, validate) == expected);

    }

    int sample_expected_consecutive[] = {1, 1, 3, 0, 1, 2};
    BOOST_DATA_TEST_CASE(searchRegexConsecutiveTest,
                         boost::unit_test::data::make(sample_regex_motif) ^ sample_expected_consecutive ^
                         sample_validate, regex_motif, expected, validate) {
        std::string sequence = "TTAGGGCTTAGGGAAATTAGGGCCCTTAGGGactttAGGGTTaGGGTTaaCCC";

        BOOST_TEST(cmri::searchRegexConsecutive(sequence, regex_motif, validate) == expected);

    }

    std::string sample_input[] = {"data/input.fq",
                                        "data/input.bam"
    };

    std::string sample_result[] = {"data/result.json",
                                        "data/bam_result.json"
    };


    BOOST_DATA_TEST_CASE(processTest,boost::unit_test::data::make(sample_input)^sample_result,input_file,result_file){

        std::string motif_file = "data/input.json";
        std::map<std::string, cmri::region_list_t> motifs = cmri::deserialize(motif_file);
        std::map<std::string, cmri::region_list_t> result = cmri::deserialize(result_file);
        cmri::options_t options;
        options.input_file=input_file;
        options.quality_value=10;
        options.quality_map=30;
        cmri::process(options, motifs);

        for(auto &m : motifs){
            BOOST_TEST(result[m.first]==m.second);
        }

        BOOST_TEST(motifs == result);
    }


    BOOST_DATA_TEST_CASE(processMultiThreadingTest,boost::unit_test::data::make(sample_input)^sample_result,input_file,result_file){

        std::string motif_file = "data/input.json";
        std::map<std::string, cmri::region_list_t> motifs = cmri::deserialize(motif_file);
        std::map<std::string, cmri::region_list_t> result = cmri::deserialize(result_file);
        cmri::options_t options;
        options.input_file=input_file;
        options.quality_value=10;
        options.quality_map=30;
        options.threads = static_cast<int>(std::thread::hardware_concurrency());
        options.chunk_size = 1;
        cmri::processMultiThreading(options, motifs);

        for(auto &m : motifs){
            BOOST_TEST(result[m.first]==m.second);
        }


        BOOST_TEST(motifs == result);
    }


BOOST_AUTO_TEST_SUITE_END()