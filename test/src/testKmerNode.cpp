#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include "../../src/Modules/TelomereAnalysis/kmerNode.h"


BOOST_AUTO_TEST_SUITE(kmerNodeTest)

    std::string sample_sequence[] = {"TTAGGG"};

    BOOST_DATA_TEST_CASE(newKmerNodeTest,
                         boost::unit_test::data::make(sample_sequence), sequence) {


        cmri::kmerNode root('R',0);
        for(int i =0 ; i < sequence.size(); i++){
            cmri::kmerNode new_node(sequence[i],i+1);
//            root.addChildren(new_node);
            BOOST_TEST_MESSAGE("kmer " << new_node.serialize()  );
        }
        BOOST_TEST_MESSAGE("kmer " << root.serialize()  );




    }

BOOST_AUTO_TEST_SUITE_END()