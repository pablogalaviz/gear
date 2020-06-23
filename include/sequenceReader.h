#ifndef GEAR_SEQUENCEREADER_H
#define GEAR_SEQUENCEREADER_H
//
// Author(s) Pablo Galaviz (2020)
// e-mail  <pgalaviz@cmri.org.au>
//



//  This file is part of GEAR
//
//  GEAR is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  any later version.
//
//  GEAR is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GEAR.  If not, see <http://www.gnu.org/licenses/>.
//


#include "csvParser.h"
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include "kseq.h"
#include <sam.h>
#include <string>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

namespace cmri {

    struct read_item_t {

        std::string sequence;
        unsigned int start = 0;
        unsigned int end = 0;
        std::string name;
        bool valid = true;

    };


    class sequenceReader {

        //Fastq/Fasta file handler.
        kseq_t *kseq;

        //Csv file handler.
        std::istream *instream;
        csvParser *csv_parser;
        std::ifstream file;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;

        //Sam/Bam file handler
        samFile *bam_file;
        bam_hdr_t *bam_header; //read header
        bam1_t *alignment; //initialize an alignment

        bool (sequenceReader::*getItem)(read_item_t &item);
        void (sequenceReader::*closeFile)();

        bool getFastxItem(read_item_t &item);
        void closeFastxFile();

        bool getCsvItem(read_item_t &item);
        void closeCsvFile();

        bool getBamItem(read_item_t &item);
        void closeBamFile();

        int count;
        int quality_value;
        int quality_map;
        int total_reads;

    public:
        sequenceReader(std::string input_file, int _quality_value, int _quality_map);

        ~sequenceReader() = default;

        inline bool get(read_item_t &item) {return (this->*getItem)(item);}
        inline void close() {(this->*closeFile)();}

        inline int getCount() const {
            return count;
        }
        inline int getTotalReads() const {
            return total_reads;
        }



    };

}
#endif //GEAR_SEQUENCEREADER_H
