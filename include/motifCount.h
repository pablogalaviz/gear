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


#include "motif.h"
#include <regex>
#include <future>
#include "utils.h"
#include "sequenceReader.h"



namespace cmri {


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


    std::map<std::string, region_list_t>
    processWorker(std::map<std::string, region_list_t> motif_map, std::vector<read_item_t> sequences, bool validate) {

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
                    item.count++;
                    item.total_bases+=seq.sequence.size();

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
                        item.count++;
                        item.total_bases+=seq.sequence.size();
                    }
                }
            }

        }

        return motif_map;
    }


    void process(const options_t &options, std::map<std::string, region_list_t> &motif_map) {

        cmri::sequenceReader reader(options.input_file,options.quality_value,options.quality_map);
        int total_reads = reader.getTotalReads();
        LOGGER.info << "Processing: " << total_reads << " reads." << std::endl;

        cmri::read_item_t read_item;
        int next_log = options.progress;
        while(reader.get(read_item)){
            if(!read_item.valid){continue;}

        if (motif_map.find(read_item.name) != motif_map.end()) {
                for (auto &item : motif_map[read_item.name]) {

                    if (!item.intersect(read_item.start, read_item.end)) { continue; }

                    for (auto &m : item.motifs) {
                        m.second += searchMotif(read_item.sequence, m.first, options.validate_sequence);
                    }
                    for (auto &r : item.regex) {
                        r.second += searchRegex(read_item.sequence, r.first, options.validate_sequence);
                    }
                    item.count++;
                    item.total_bases+=read_item.sequence.size();
                }
            } else {
                if (motif_map.find("other") != motif_map.end()) {
                    for (auto &item : motif_map["other"]) {
                        for (auto &m : item.motifs) {
                            m.second += searchMotif(read_item.sequence, m.first, options.validate_sequence);
                        }
                        for (auto &r : item.regex) {
                            r.second += searchRegex(read_item.sequence, r.first, options.validate_sequence);
                        }
                        item.count++;
                        item.total_bases+=read_item.sequence.size();
                    }
                }
            }

            if (options.progress > 0 && reader.getCount() >= next_log) {
                LOGGER.info << "Progress: " << reader.getCount() << " of " << total_reads << " "
                            << 100.0 * reader.getCount() / total_reads << "%" << std::endl;
                next_log += options.progress;
            }

        }
        LOGGER.info << "Total sequences analysed: " << reader.getCount() << std::endl;
        reader.close();
    }


    void processMultiThreading(const options_t &options, std::map<std::string, region_list_t> &motif_map) {

        cmri::sequenceReader reader(options.input_file,options.quality_value,options.quality_map);
        int total_reads = reader.getTotalReads();
        LOGGER.info << "Processing: " << total_reads << " reads." << std::endl;

        std::vector<read_item_t> sequences[options.threads];

        cmri::read_item_t read_item;
        int next_log = options.progress;
        while(reader.get(read_item)){
            if(!read_item.valid){continue;}

            sequences[reader.getCount() % options.threads].emplace_back(read_item);

            if (sequences[options.threads - 1].size() > options.chunk_size) {

                std::vector<std::future<std::map<std::string, region_list_t>  >> pool;
                for (auto &seq_data : sequences) {
                    pool.push_back(std::async(std::launch::async, &processWorker, motif_map, seq_data, options.validate_sequence));
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

            if (options.progress > 0 && reader.getCount() >= next_log) {
                LOGGER.info << "Progress: " << reader.getCount() << " of " << total_reads << " "
                            << 100.0 * reader.getCount() / total_reads << "%" << std::endl;
                next_log += options.progress;
            }


        }

        std::vector<std::future<std::map<std::string, region_list_t>  >> pool;
        for (auto &seq_data : sequences) {
            pool.push_back(std::async(std::launch::async, &processWorker, motif_map, seq_data, options.validate_sequence));
        }
        for (auto &t : pool) {
            t.wait();
            for (auto &result : t.get()) {
                for (int i = 0; i < result.second.size(); i++) {
                    motif_map[result.first][i] += result.second[i];
                }
            }
        }


        LOGGER.info << "Total sequences analysed: " << reader.getCount() << std::endl;
        reader.close();
    }


}

#endif //GEAR_MOTIFCOUNT_H
