#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "utils.h"

BOOST_AUTO_TEST_SUITE(utilsTest)

    std::string sample_file_name[] = {"input.bam", "input.csv", "input.csv.gz", "input.fasta", "input.fasta.gz",
                                      "input.fa", "input.fa.gz", "input.fastq", "input.fastq.gz", "input.fq",
                                      "input.fq.gz", "input.txt"
    };

    cmri::format_t sample_expected_format[] = {cmri::format_t::BAM, cmri::format_t::CSV, cmri::format_t::CSV,
                                               cmri::format_t::FASTA, cmri::format_t::FASTA, cmri::format_t::FASTA,
                                               cmri::format_t::FASTA, cmri::format_t::FASTQ, cmri::format_t::FASTQ,
                                               cmri::format_t::FASTQ, cmri::format_t::FASTQ, cmri::format_t::UNKNOWN
    };
    bool sample_compressed[] = {false, false, true, false, true, false, true, false, true, false,
                                true, false};


    BOOST_DATA_TEST_CASE(fileFormatTest,
                         boost::unit_test::data::make(sample_file_name) ^ sample_expected_format
                         , file_name, expected_format) {

        BOOST_TEST((cmri::file_format(file_name) & cmri::format_t::FILE_TYPE) == expected_format);

    }

    BOOST_DATA_TEST_CASE(compressFormatTest,
                         boost::unit_test::data::make(sample_file_name)  ^
                         sample_compressed, file_name, compressed) {

        BOOST_TEST( static_cast<bool>(cmri::file_format(file_name) & cmri::format_t::GZIP) == compressed);

    }


BOOST_AUTO_TEST_SUITE_END()