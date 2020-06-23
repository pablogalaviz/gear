#include <utils.h>
#include <csvParser.h>
#include "sequenceReader.h"

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


cmri::sequenceReader::sequenceReader(std::string input_file, int _quality_value, int _quality_map) :
        quality_value(_quality_value), quality_map(_quality_map) {

    count = 0;
    total_reads= count_reads(input_file);

    auto format_flag = cmri::file_format(input_file);
    auto format = static_cast<cmri::format_t>(format_flag & cmri::format_t::FILE_TYPE);

    switch (format) {
        case cmri::format_t::FASTA:
        case cmri::format_t::FASTQ: {
            gzFile file = gzopen(input_file.c_str(), "r");
            kseq = kseq_init(file);
            getItem = &sequenceReader::getFastxItem;
            closeFile = &sequenceReader::closeFastxFile;
        }
            break;
        case cmri::format_t::CSV: {
            file.open(input_file, std::ios_base::in | std::ios_base::binary);
            if (format_flag & cmri::GZIP) { inbuf.push(boost::iostreams::gzip_decompressor()); }
            inbuf.push(file);
            //Convert streambuf to istream
            instream = new std::istream(&inbuf);
            inbuf.set_auto_close(false);
            std::string line;

            std::getline(*instream, line);
            csv_parser = new csvParser(line, std::vector<std::string>{"sequence"});

            getItem = &sequenceReader::getCsvItem;
            closeFile = &sequenceReader::closeCsvFile;
        }
            break;
        case cmri::format_t::BAM: {
            bam_file = hts_open(input_file.c_str(), "r");
            bam_header = sam_hdr_read(bam_file); //read header
            alignment = bam_init1(); //initialize an alignment
            getItem = &sequenceReader::getBamItem;
            closeFile = &sequenceReader::closeBamFile;
        }
            break;
        case cmri::format_t::UNKNOWN:
        default:
            cmri::LOGGER.error << "Unknown input file format " << std::endl;
            exit(EINVAL);
    }

}

bool cmri::sequenceReader::getFastxItem(read_item_t &item) {

    int l;
    if ((l = kseq_read(kseq)) >= 0) {
        item.sequence = kseq->seq.s;
        double mean_qv = 0;
        if (kseq->qual.l > 0) {
            for (int i = 0; i < kseq->qual.l; i++) {
                mean_qv += static_cast<int>(kseq->qual.s[i]) - 33;
            }
            mean_qv /= kseq->qual.l;
        }
        item.name = mean_qv < quality_value ? "qv_fail" : "unmapped";
        count++;
        item.valid= true;
        return true;
    }

    item.valid= false;
    return false;
}

void cmri::sequenceReader::closeFastxFile() {
    kseq_destroy(kseq);
}


bool cmri::sequenceReader::getCsvItem(read_item_t &item) {

    std::string line;
    if (std::getline(*instream, line)) {
        auto field = csv_parser->parseLine(line);
        item.sequence = field["sequence"];
        item.name = "unmapped";
        count++;
        item.valid= true;
        return true;
    }
    item.valid= false;
    return false;
}

void cmri::sequenceReader::closeCsvFile() {
    file.close();
}


bool cmri::sequenceReader::getBamItem(read_item_t &item) {

    if (sam_read1(bam_file, bam_header, alignment) > 0) {
        int chromosome_id = alignment->core.tid;

        if ((alignment->core.flag & BAM_FSECONDARY)
            || (alignment->core.flag & BAM_FDUP)
//            || (alignment->core.flag & BAM_FQCFAIL)
            || (alignment->core.flag & BAM_FSUPPLEMENTARY)
                ) {
            item.sequence = "";
            item.start = 0;
            item.name = "";
            item.end = 0;
            item.valid= false;
#ifdef DEBUG
            LOGGER.debug << "BAM_FSECONDARY,  BAM_FDUP, BAM_FQCFAIL, BAM_FSUPPLEMENTARY" << std::endl;
#endif
            return true;
        }

        //contig name (chromosome)
        std::string chromosome = chromosome_id >= 0 ? bam_header->target_name[chromosome_id] : "unmapped";
        auto query_name = bam_get_qname(alignment);
        auto quality = bam_get_seq(alignment); //quality string

        uint32_t len = alignment->core.l_qseq; //length of the read.
        int32_t start = alignment->core.pos + 1; //left most position of alignment in zero based coordinate (+1)
        int32_t end = start + len;

        std::string sequence = "";
        double mean_qv = 0;
        for (int i = 0; i < len; i++) {
            sequence += seq_nt16_str[bam_seqi(quality, i)]; //gets nucleotide id and converts them into IUPAC id.
            mean_qv += static_cast<int>(quality[i]);
        }
        mean_qv /= len;
        uint32_t mapping_quality = alignment->core.qual;

        if (chromosome_id >= 0) {
            if (mean_qv < quality_value) { chromosome = "qv_fail"; }
            else {
                if (mapping_quality < quality_map) { chromosome = "mapq_fail"; }
            }
        }
#ifdef DEBUG
        else {LOGGER.debug << "unmapped " << std::endl;}
#endif
        item.sequence = sequence;
        item.end = end;
        item.name = chromosome;
        item.start= start;
        item.valid= true;

        count++;
        return true;
    }

    item.valid= false;
    return false;
}

void cmri::sequenceReader::closeBamFile() {
    bam_destroy1(alignment);
    bam_hdr_destroy(bam_header);
    sam_close(bam_file);
}
