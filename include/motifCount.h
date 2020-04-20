#ifndef GEAR_MOTIFCOUNT_H
#define GEAR_MOTIFCOUNT_H

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

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include "csvParser.h"
#include "kseq.h"
#include <map>
#include "motif.h"
#include <regex>
#include <sam.h>
#include <string>
#include "utils.h"
#include <vector>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)


namespace cmri {


    int searchMotif(std::string sequence, const std::string &motif, bool validate) {

        if (validate) { for (auto &c: sequence) { c = toupper(c); }}
        int occurrences = 0;
        std::string::size_type start = 0;
        while ((start = sequence.find(motif, start)) != std::string::npos) {
            ++occurrences;
            // move start to find the rest of the string.
            start += motif.size();
        }
        return occurrences;
    }

    int searchRegex(std::string sequence, const std::string &regex, bool validate) {

        if (validate) { for (auto &c: sequence) { c = toupper(c); }}
        int occurrences = 0;

        std::regex basic_regex(regex);
        std::smatch match;

        while (std::regex_search(sequence, match, basic_regex)) {
            ++occurrences;
            // suffix to find the rest of the string.
            sequence = match.suffix().str();
        }

        return occurrences;
    }


    int searchRegexConsecutive(std::string sequence, const std::string &regex, bool validate) {

        if (validate) { for (auto &c: sequence) { c = toupper(c); }}
        int occurrences = 0;

        std::regex basic_regex(regex);
        std::smatch match;

        while (std::regex_search(sequence, match, basic_regex)) {
            ++occurrences;
            // suffix to find the rest of the string.
            auto m = match.position() + match.length() - 6;
            sequence = sequence.substr(m);
        }

        return occurrences;
    }


    void processFast(const std::string &input_file, std::map<std::string, region_list_t> &motif_map, bool validate,
                     int progress) {

        std::map<format_t, int> factor_map = {{FASTQ_GZ, 4},
                                              {FASTQ,    4,},
                                              {FASTA_GZ, 2},
                                              {FASTA,    2}};
        int factor = factor_map[file_format(input_file)];
        int lines_in_file = count_lines(input_file) / factor;
        LOGGER.info << "Processing: " << lines_in_file << " reads."<< std::endl;

        gzFile file = gzopen(input_file.c_str(), "r");
        kseq_t *record = kseq_init(file);
        int l;
        int count = 0;
        int next_log = progress;
        while ((l = kseq_read(record)) >= 0) {
            std::string sequence = record->seq.s;
            for (auto &item : motif_map["unmapped"]) {
                for (auto &m : item.motifs) {
                    m.second += searchMotif(sequence, m.first, validate);
                }
                for (auto &r : item.regex) {
                    r.second += searchRegex(sequence, r.first, validate);
                }
                count++;
                if (progress > 0 && count >= next_log) {
                    LOGGER.info << "Progress: " << count << " of " << lines_in_file << " "
                                << 100 * count / lines_in_file << "%" << std::endl;
                    next_log += progress;
                }
            }
        }
        LOGGER.info << "Total sequences analysed: " << count << std::endl;
    }


    void processCsv(const std::string &input_file, std::map<std::string, region_list_t> &motif_map, bool validate,
                    int progress,
                    bool decompress = false) {
        cmri::LOGGER.debug << "processCsv file " << std::endl;
        std::ifstream file(input_file, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        if (decompress) { inbuf.push(boost::iostreams::gzip_decompressor()); }
        inbuf.push(file);
        //Convert streambuf to istream
        std::istream instream(&inbuf);
        inbuf.set_auto_close(false);

        int lines_in_file = count_lines(input_file);
        LOGGER.info << "Processing: " << lines_in_file << " reads."<< std::endl;

        //Iterate lines
        std::string line;

        std::getline(instream, line);
        csvParser csv_parser(line, std::vector<std::string>{"sequence"});
        int count = 0;
        int next_log = progress;
        while (std::getline(instream, line)) {
            auto field = csv_parser.parseLine(line);
            for (auto &item : motif_map["unmapped"]) {
                for (auto &m : item.motifs) {
                    m.second += searchMotif(field["sequence"], m.first, validate);
                }
                for (auto &r : item.regex) {
                    r.second += searchRegex(field["sequence"], r.first, validate);
                }
            }
            count++;
            if (progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " of " << lines_in_file << " "
                            << 100 * count / lines_in_file << "%" << std::endl;
                next_log += progress;
            }
        }
        LOGGER.info << "Total sequences analysed: " << count << std::endl;
        //Cleanup
        file.close();
    }


    void
    processBam(const std::string &input_file, std::map<std::string, region_list_t> &motif_map,
               bool validate, int progress) {


        samFile *file = hts_open(input_file.c_str(), "r");
        bam_hdr_t *bam_header = sam_hdr_read(file); //read header
        bam1_t *alignment = bam_init1(); //initialize an alignment

        auto bam_index = sam_index_load(file, (input_file + ".bai").c_str());

        auto n_targets = bam_header->n_targets;
        int total_reads=0;
        for(int tid =0 ; tid < n_targets; tid++){
            uint64_t mapped;
            uint64_t unmapped;
            if (hts_idx_get_stat(bam_index, tid, &mapped, &unmapped)==0) {
                total_reads+=mapped+unmapped;
            }
        }

        LOGGER.info << "Processing: " << total_reads << " reads."<< std::endl;


        int count = 0;
        int next_log = progress;
        while (sam_read1(file, bam_header, alignment) > 0) {
            int chromosome_id = alignment->core.tid;

            //if (chromosome_id < 0) {continue;}


            if ((alignment->core.flag & BAM_FSECONDARY)
                || (alignment->core.flag & BAM_FDUP)
                || (alignment->core.flag & BAM_FQCFAIL)
                || (alignment->core.flag & BAM_FSUPPLEMENTARY)
                    ) { continue; }

            uint32_t mapping_quality = alignment->core.qual;
            /*
            //if(mapping_quality == 0){ continue;}
            LOGGER.debug << "mapQ " << mapping_quality << std::endl;

            if (alignment->core.flag != 0 && !(alignment->core.flag & BAM_FREVERSE)) {
                LOGGER.warning << "Flags " << bam_flag2str(alignment->core.flag) << std::endl;
            }
*/
            std::string chromosome =
                    chromosome_id >= 0 ? bam_header->target_name[chromosome_id] : "unknown"; //contig name (chromosome)

            auto query_name = bam_get_qname(alignment);

            auto quality = bam_get_seq(alignment); //quality string



            uint32_t len = alignment->core.l_qseq; //length of the read.
            int32_t start = alignment->core.pos + 1; //left most position of alignment in zero based coordianate (+1)
            int32_t end = start + len;

            std::string sequence = "";
            for (int i = 0; i < len; i++) {
                sequence += seq_nt16_str[bam_seqi(quality, i)]; //gets nucleotide id and converts them into IUPAC id.
                //LOGGER.debug << i << " - " << static_cast<int>(quality[i]) << std::endl;
            }


            if (motif_map.find(chromosome) != motif_map.end()) {
                for (auto &item : motif_map[chromosome]) {

                    if (!item.intersect(start, end)) { continue; }

                    for (auto &m : item.motifs) {
                        m.second += searchMotif(sequence, m.first, validate);
                    }
                    for (auto &r : item.regex) {
                        r.second += searchRegex(sequence, r.first, validate);
                    }
                }
            } else {
                if (motif_map.find("unmapped") != motif_map.end()) {
                    for (auto &item : motif_map["unmapped"]) {
                        //TODO:remove hack
                        int total_found=0;
                        for (auto &m : item.motifs) {
                            auto value = searchMotif(sequence, m.first, validate);
                            m.second += value;
                            total_found+=value;
                        }

                        if(total_found>0)
                        for (auto &r : item.regex) {
                            r.second += searchRegex(sequence, r.first, validate);
                        }
                    }
                }
            }
            count++;
            if (progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " of " << total_reads << " "
                            << 100.0 * count / total_reads << "%" << std::endl;
                next_log += progress;
            }

        }

        LOGGER.info << "Total sequences analysed: " << count << std::endl;
        bam_destroy1(alignment);
        sam_close(file);
    }


}

#endif //GEAR_MOTIFCOUNT_H
