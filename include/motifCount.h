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
#include <future>

KSEQ_INIT(gzFile, gzread)


namespace cmri {

    struct bam_seq_t {

        bam_seq_t(const std::string &sequence, unsigned int start, unsigned int anEnd, const std::string &name)
                : sequence(sequence), start(start), end(anEnd), name(name) {}

        std::string sequence;
        unsigned int start = 0;
        unsigned int end = 0;
        std::string name;
    };


    unsigned int searchMotif(std::string sequence, const std::string &motif, bool validate) {

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

    void asyncSearchMotif(const std::string &sequence, std::map<std::string, unsigned int> &motif, bool validate) {
        std::map<std::string, std::future<unsigned int> > future_motif_results;
        for (auto &m : motif) {
            future_motif_results[m.first] = std::async(std::launch::async, &searchMotif, sequence, m.first, validate);
        }
        for (auto &m : motif) {
            m.second += future_motif_results[m.first].get();
        }
    }

    unsigned int searchRegex(std::string sequence, const std::string &regex, bool validate) {

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

    void asyncSearchRegex(const std::string &sequence, std::map<std::string, unsigned int> &regex, bool validate) {
        std::map<std::string, std::future<unsigned int> > future_motif_results;
        for (auto &r : regex) {
            future_motif_results[r.first] = std::async(std::launch::async, &searchRegex, sequence, r.first, validate);
        }
        for (auto &r : regex) {
            r.second += future_motif_results[r.first].get();
        }
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


    void processFast(const options_t &options, std::map<std::string, region_list_t> &motif_map) {

        std::map<format_t, int> factor_map = {{FASTQ_GZ, 4},
                                              {FASTQ,    4,},
                                              {FASTA_GZ, 2},
                                              {FASTA,    2}};
        int factor = factor_map[file_format(options.input_file)];
        int lines_in_file = count_lines(options.input_file) / factor;
        LOGGER.info << "Processing: " << lines_in_file << " reads." << std::endl;

        gzFile file = gzopen(options.input_file.c_str(), "r");
        kseq_t *record = kseq_init(file);
        int l;
        int count = 0;
        int next_log = options.progress;
        while ((l = kseq_read(record)) >= 0) {
            std::string sequence = record->seq.s;
            for (auto &item : motif_map["unmapped"]) {
                for (auto &m : item.motifs) {
                    m.second += searchMotif(sequence, m.first, options.validate_sequence);
                }
                for (auto &r : item.regex) {
                    r.second += searchRegex(sequence, r.first, options.validate_sequence);
                }
            }

            count++;
            if (options.progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " of " << lines_in_file << " "
                            << 100 * count / lines_in_file << "%" << std::endl;
                next_log += options.progress;
            }

        }
        LOGGER.info << "Total sequences analysed: " << count << std::endl;

        kseq_destroy(record);
    }


    region_list_t
    processSequenceWorker(region_list_t region_list, std::vector<std::string> sequences, bool validate) {

        for (auto &item : region_list) {
            item.resetCount();
        }

        for (auto seq : sequences) {
            for (auto &item : region_list) {
                for (auto &m : item.motifs) {
                    m.second += searchMotif(seq, m.first, validate);
                }
                for (auto &r : item.regex) {
                    r.second += searchRegex(seq, r.first, validate);
                }
            }
        }

        return region_list;
    }

    void processFastMultiThreading(const options_t &options, std::map<std::string, region_list_t> &motif_map) {

        std::map<format_t, int> factor_map = {{FASTQ_GZ, 4},
                                              {FASTQ,    4,},
                                              {FASTA_GZ, 2},
                                              {FASTA,    2}};
        int factor = factor_map[file_format(options.input_file)];
        int lines_in_file = count_lines(options.input_file) / factor;
        LOGGER.info << "Processing: " << lines_in_file << " reads." << std::endl;

        std::vector<std::string> sequences[options.threads];


        gzFile file = gzopen(options.input_file.c_str(), "r");
        kseq_t *record = kseq_init(file);
        int l;
        int count = 0;
        int next_log = options.progress;
        while ((l = kseq_read(record)) >= 0) {
            std::string sequence = record->seq.s;

            sequences[count % options.threads].emplace_back(sequence);

            if (sequences[options.threads - 1].size() > options.chunk_size) {

                std::vector<std::future<region_list_t >> pool;
                for (auto &seq_data : sequences) {
                    pool.push_back(
                            std::async(std::launch::async, &processSequenceWorker, motif_map["unmapped"], seq_data,
                                       options.validate_sequence));
                    seq_data.clear();
                }
                for (auto &t : pool) {
                    t.wait();
                    int i = 0;
                    for (auto &result : t.get()) {
                        motif_map["unmapped"][i] += result;
                        i += 1;
                    }
                }
            }


            count++;
            if (options.progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " of " << lines_in_file << " "
                            << 100 * count / lines_in_file << "%" << std::endl;
                next_log += options.progress;
            }

        }
        LOGGER.info << "Total sequences analysed: " << count << std::endl;

        kseq_destroy(record);
    }


    void processCsvMultiThreading(const options_t &options, std::map<std::string, region_list_t> &motif_map, bool decompress = false) {

        cmri::LOGGER.debug << "processCsv file " << std::endl;
        std::ifstream file(options.input_file, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        if (decompress) { inbuf.push(boost::iostreams::gzip_decompressor()); }
        inbuf.push(file);
        //Convert streambuf to istream
        std::istream instream(&inbuf);
        inbuf.set_auto_close(false);

        int lines_in_file = count_lines(options.input_file);
        LOGGER.info << "Processing: " << lines_in_file << " reads." << std::endl;

        std::vector<std::string> sequences[options.threads];

        //Iterate lines
        std::string line;

        std::getline(instream, line);
        csvParser csv_parser(line, std::vector<std::string>{"sequence"});
        int count = 0;
        int next_log = options.progress;
        while (std::getline(instream, line)) {
            auto field = csv_parser.parseLine(line);
            std::string sequence = field["sequence"];

            sequences[count % options.threads].emplace_back(sequence);

            if (sequences[options.threads - 1].size() > options.chunk_size) {

                std::vector<std::future<region_list_t >> pool;
                for (auto &seq_data : sequences) {
                    pool.push_back(
                            std::async(std::launch::async, &processSequenceWorker, motif_map["unmapped"], seq_data,
                                       options.validate_sequence));
                    seq_data.clear();
                }
                for (auto &t : pool) {
                    t.wait();
                    int i = 0;
                    for (auto &result : t.get()) {
                        motif_map["unmapped"][i] += result;
                        i += 1;
                    }
                }
            }

            count++;
            if (options.progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " of " << lines_in_file << " "
                            << 100 * count / lines_in_file << "%" << std::endl;
                next_log += options.progress;
            }
        }
        LOGGER.info << "Total sequences analysed: " << count << std::endl;
        //Cleanup
        file.close();
    }



    void processCsv(const options_t &options, std::map<std::string, region_list_t> &motif_map,bool decompress = false) {
        cmri::LOGGER.debug << "processCsv file " << std::endl;
        std::ifstream file(options.input_file, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        if (decompress) { inbuf.push(boost::iostreams::gzip_decompressor()); }
        inbuf.push(file);
        //Convert streambuf to istream
        std::istream instream(&inbuf);
        inbuf.set_auto_close(false);

        int lines_in_file = count_lines(options.input_file);
        LOGGER.info << "Processing: " << lines_in_file << " reads." << std::endl;

        //Iterate lines
        std::string line;

        std::getline(instream, line);
        csvParser csv_parser(line, std::vector<std::string>{"sequence"});
        int count = 0;
        int next_log = options.progress;
        while (std::getline(instream, line)) {
            auto field = csv_parser.parseLine(line);
            for (auto &item : motif_map["unmapped"]) {
                for (auto &m : item.motifs) {
                    m.second += searchMotif(field["sequence"], m.first, options.validate_sequence);
                }
                for (auto &r : item.regex) {
                    r.second += searchRegex(field["sequence"], r.first, options.validate_sequence);
                }
            }


            count++;
            if (options.progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " of " << lines_in_file << " "
                            << 100 * count / lines_in_file << "%" << std::endl;
                next_log += options.progress;
            }
        }
        LOGGER.info << "Total sequences analysed: " << count << std::endl;
        //Cleanup
        file.close();
    }


    void processBam(const options_t &options, std::map<std::string, region_list_t> &motif_map) {


        samFile *file = hts_open(options.input_file.c_str(), "r");
        bam_hdr_t *bam_header = sam_hdr_read(file); //read header
        bam1_t *alignment = bam_init1(); //initialize an alignment

        std::string index_file_name = options.input_file + ".bai";
        open_file(index_file_name);
        auto bam_index = sam_index_load(file, index_file_name.c_str());

        auto n_targets = bam_header->n_targets;
        int total_reads = 0;
        for (int tid = 0; tid < n_targets; tid++) {
            uint64_t mapped;
            uint64_t unmapped;
            if (hts_idx_get_stat(bam_index, tid, &mapped, &unmapped) == 0) {
                total_reads += static_cast<int>(mapped + unmapped);
            }
        }

        LOGGER.info << "Processing: " << total_reads << " reads." << std::endl;

        int count = 0;
        int next_log = options.progress;
        while (sam_read1(file, bam_header, alignment) > 0) {
            int chromosome_id = alignment->core.tid;

            //if (chromosome_id < 0) {continue;}


            if ((alignment->core.flag & BAM_FSECONDARY)
                || (alignment->core.flag & BAM_FDUP)
                || (alignment->core.flag & BAM_FQCFAIL)
                || (alignment->core.flag & BAM_FSUPPLEMENTARY)
                    ) { continue; }

            /*
            if (alignment->core.flag != 0 && !(alignment->core.flag & BAM_FREVERSE)) {
                LOGGER.warning << "Flags " << bam_flag2str(alignment->core.flag) << std::endl;
            }
*/
            //contig name (chromosome)
            std::string chromosome = chromosome_id >= 0 ? bam_header->target_name[chromosome_id] : "unmapped";

            auto query_name = bam_get_qname(alignment);

            auto quality = bam_get_seq(alignment); //quality string


            uint32_t len = alignment->core.l_qseq; //length of the read.
            int32_t start = alignment->core.pos + 1; //left most position of alignment in zero based coordinate (+1)
            int32_t end = start + len;

            std::string sequence = "";
            double mean_qv=0;
            for (int i = 0; i < len; i++) {
                sequence += seq_nt16_str[bam_seqi(quality, i)]; //gets nucleotide id and converts them into IUPAC id.
                mean_qv+=static_cast<int>(quality[i]);
            }
            mean_qv/=len;
            uint32_t mapping_quality = alignment->core.qual;

            if(mapping_quality < options.quality_map || mean_qv < options.quality_value) {
                chromosome = "fail";
            }

            if (motif_map.find(chromosome) != motif_map.end()) {
                for (auto &item : motif_map[chromosome]) {

                    if (!item.intersect(start, end)) { continue; }

                    for (auto &m : item.motifs) {
                        m.second += searchMotif(sequence, m.first, options.validate_sequence);
                    }
                    for (auto &r : item.regex) {
                        r.second += searchRegex(sequence, r.first, options.validate_sequence);
                    }
                }
            } else {
                if (motif_map.find("other") != motif_map.end()) {
                    for (auto &item : motif_map["other"]) {
                        for (auto &m : item.motifs) {
                            m.second += searchMotif(sequence, m.first, options.validate_sequence);
                        }
                        for (auto &r : item.regex) {
                            r.second += searchRegex(sequence, r.first, options.validate_sequence);
                        }
                    }
                }
            }
            count++;
            if (options.progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " of " << total_reads << " "
                            << 100.0 * count / total_reads << "%" << std::endl;
                next_log += options.progress;
            }


        }

        LOGGER.info << "Total sequences analysed: " << count << std::endl;
        bam_destroy1(alignment);
        bam_hdr_destroy(bam_header);
        hts_idx_destroy(bam_index);
        sam_close(file);
    }


    std::map<std::string, region_list_t>
    processBamWorker(std::map<std::string, region_list_t> motif_map, std::vector<bam_seq_t> sequences, bool validate) {

        for (auto &item : motif_map) {
            for (auto &m : item.second) {
                m.resetCount();
            }
        }

        for (auto seq : sequences) {

            if (motif_map.find(seq.name) != motif_map.end()) {
                for (auto &item : motif_map[seq.name]) {

                    if (!item.intersect(seq.start, seq.end)) { continue; }

                    for (auto &m : item.motifs) {
                        m.second += searchMotif(seq.sequence, m.first, validate);
                    }
                    for (auto &r : item.regex) {
                        r.second += searchRegex(seq.sequence, r.first, validate);
                    }
                }
            } else {
                if (motif_map.find("other") != motif_map.end()) {
                    for (auto &item : motif_map["other"]) {
                        for (auto &m : item.motifs) {
                            m.second += searchMotif(seq.sequence, m.first, validate);
                        }
                        for (auto &r : item.regex) {
                            r.second += searchRegex(seq.sequence, r.first, validate);
                        }
                    }
                }
            }

        }

        return motif_map;
    }

    void processBamMultiThreading(const options_t &options, std::map<std::string, region_list_t> &motif_map) {


        samFile *file = hts_open(options.input_file.c_str(), "r");
        bam_hdr_t *bam_header = sam_hdr_read(file); //read header
        bam1_t *alignment = bam_init1(); //initialize an alignment

        std::string index_file_name = options.input_file + ".bai";
        open_file(index_file_name);
        auto bam_index = sam_index_load(file, index_file_name.c_str());

        auto n_targets = bam_header->n_targets;
        int total_reads = 0;
        for (int tid = 0; tid < n_targets; tid++) {
            uint64_t mapped;
            uint64_t unmapped;
            if (hts_idx_get_stat(bam_index, tid, &mapped, &unmapped) == 0) {
                total_reads += static_cast<int>(mapped + unmapped);
            }
        }

        LOGGER.info << "Processing: " << total_reads << " reads." << std::endl;

        std::vector<bam_seq_t> sequences[options.threads];

        int count = 0;
        int next_log = options.progress;
        while (sam_read1(file, bam_header, alignment) > 0) {
            int chromosome_id = alignment->core.tid;

            //if (chromosome_id < 0) {continue;}


            if ((alignment->core.flag & BAM_FSECONDARY)
                || (alignment->core.flag & BAM_FDUP)
                || (alignment->core.flag & BAM_FQCFAIL)
                || (alignment->core.flag & BAM_FSUPPLEMENTARY)
                    ) { continue; }

            //contig name (chromosome)
            std::string chromosome = chromosome_id >= 0 ? bam_header->target_name[chromosome_id] : "unknown";

            auto query_name = bam_get_qname(alignment);

            auto quality = bam_get_seq(alignment); //quality string

            uint32_t len = alignment->core.l_qseq; //length of the read.
            int32_t start = alignment->core.pos + 1; //left most position of alignment in zero based coordianate (+1)
            int32_t end = start + len;

            std::string sequence = "";
            double mean_qv=0;
            for (int i = 0; i < len; i++) {
                sequence += seq_nt16_str[bam_seqi(quality, i)]; //gets nucleotide id and converts them into IUPAC id.
                mean_qv+=static_cast<int>(quality[i]);
            }

            mean_qv/=len;
            uint32_t mapping_quality = alignment->core.qual;

            if(mapping_quality < options.quality_map || mean_qv < options.quality_value) {
                chromosome = "fail";
            }

            sequences[count % options.threads].emplace_back(bam_seq_t(sequence, start, end, chromosome));

            if (sequences[options.threads - 1].size() > options.chunk_size) {

                std::vector<std::future<std::map<std::string, region_list_t>  >> pool;
                for (auto &seq_data : sequences) {
                    pool.push_back(std::async(std::launch::async, &processBamWorker, motif_map, seq_data, options.validate_sequence));
                    seq_data.clear();
                }
                for (auto &t : pool) {
                    t.wait();
                    for (auto &result : t.get()) {
                        for (int i = 0; i < result.second.size(); i++) {
                            motif_map[result.first][i] += result.second[i];
                        }
                    }
                }

            }

            count++;
            if (options.progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " of " << total_reads << " "
                            << 100.0 * count / total_reads << "%" << std::endl;
                next_log += options.progress;
            }


        }


        std::vector<std::future<std::map<std::string, region_list_t>  >> pool;
        for (auto &seq_data : sequences) {
            pool.push_back(std::async(std::launch::async, &processBamWorker, motif_map, seq_data, options.validate_sequence));
        }
        for (auto &t : pool) {
            t.wait();
            for (auto &result : t.get()) {
                for (int i = 0; i < result.second.size(); i++) {
                    motif_map[result.first][i] += result.second[i];
                }
            }
        }


        LOGGER.info << "Total sequences analysed: " << count << std::endl;
        bam_destroy1(alignment);
        bam_hdr_destroy(bam_header);
        hts_idx_destroy(bam_index);
        sam_close(file);
    }


}

#endif //GEAR_MOTIFCOUNT_H
