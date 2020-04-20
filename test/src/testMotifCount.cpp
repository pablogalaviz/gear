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


BOOST_AUTO_TEST_SUITE_END()