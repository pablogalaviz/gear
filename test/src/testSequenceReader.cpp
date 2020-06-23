#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include "sequenceReader.h"


BOOST_AUTO_TEST_SUITE(sequenceReaderTest)


    std::string sample_file_name[] = {"data/input.fq", "data/input.fq.gz","data/input.csv","data/input.bam"};


    BOOST_DATA_TEST_CASE(readerTest,
                         boost::unit_test::data::make(sample_file_name), file_name) {

        cmri::sequenceReader reader(file_name,10,30);

        int data_size = 0;
        int record = 0;
       cmri::read_item_t item;
       while(reader.get(item)){
           if(item.valid){
           data_size += item.sequence.size();
           //BOOST_TEST_MESSAGE("sequence " << item.sequence );
           record++;
           }
       }

       reader.close();

       BOOST_TEST( (record == reader.getCount()&&data_size > 0 ));

    }


BOOST_AUTO_TEST_SUITE_END()