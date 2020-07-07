#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../../src/Modules/TelomereAnalysis/sequenceNetwork.h"

BOOST_AUTO_TEST_SUITE(sequenceNetworkTest)

    std::string sample_sequence[] = {"TTAGGGTCAGGGTTATTCTTTAGGG"};

    BOOST_DATA_TEST_CASE(sequenceNetwrokTest,
                         boost::unit_test::data::make(sample_sequence), sequence) {

        cmri::sequenceNetwork sequence_network("data/transition_matrix.csv");

        sequence_network.nt_analysis(sequence,false);


    }

BOOST_AUTO_TEST_SUITE_END()
