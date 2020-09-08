
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

#include "motifCount.h"


unsigned int cmri::searchMotif(std::string sequence, const std::string &motif) {
    int occurrences = 0;
    std::string::size_type start = 0;
    while ((start = sequence.find(motif, start)) != std::string::npos) {
        ++occurrences;
        // move start to find the rest of the string.
        start += motif.size();
    }
    return occurrences;
}

void cmri::searchMotifWFilter(std::string sequence, std::vector<int> quality, std::map<std::string,std::map<unsigned int,unsigned int>> &motif_quality){

    for(auto & motif : motif_quality){
        int occurrences = 0;
        std::string::size_type start = 0;
        while ((start = sequence.find(motif.first, start)) != std::string::npos) {
            ++occurrences;
            int step = motif.first.size();
            int qv_mean = 0;
            for(int i=start;i < start+step; i++) {
                qv_mean += quality[i];
            }
            qv_mean= static_cast<int>(std::round(qv_mean/step));
            if(qv_mean <0 || qv_mean > 99){qv_mean=0;}
            motif.second[qv_mean]++;
            // move start to find the rest of the string.
            start += step;
        }

    }

}


unsigned int cmri::searchRegex(std::string sequence, const std::string &regex) {

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


int cmri::searchRegexConsecutive(std::string sequence, const std::string &regex) {

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


cmri::mapVectorMotifRegion
cmri::processWorker(mapVectorMotifRegion motif_map, const std::vector<cmri::read_item_t>& sequences, bool validate) {


    for (auto &item : motif_map) {
        for (auto &m : item.second) {
            m.resetCount();
        }
    }

    for (const auto& seq : sequences) {

        if (motif_map.find(seq.name) != motif_map.end()) {
            for (auto &item : motif_map[seq.name]) {

                if (!item.intersect(seq.start, seq.end)) { continue; }

                std::string sequence = seq.sequence;
                if (validate) { for (auto &c: sequence) { c = toupper(c); }}

                for (auto &m : item.motifs) {
                    m.second += searchMotif(sequence, m.first);
                }
                searchMotifWFilter(sequence,seq.qvalue,item.motif_quality);
                for (auto &r : item.regex) {
                    r.second += searchRegex(sequence, r.first);
                }
                item.reads_count++;
                item.total_bases+=sequence.size();

            }
        } else {
            if (motif_map.find("other") != motif_map.end()) {

                std::string sequence = seq.sequence;
                if (validate) { for (auto &c: sequence) { c = toupper(c); }}

                for (auto &item : motif_map["other"]) {
                    for (auto &m : item.motifs) {
                        m.second += searchMotif(sequence, m.first);
                    }
                    searchMotifWFilter(sequence,seq.qvalue,item.motif_quality);
                    for (auto &r : item.regex) {
                        r.second += searchRegex(sequence, r.first);
                    }
                    item.reads_count++;
                    item.total_bases+=sequence.size();
                }
            }
        }

    }

    return motif_map;
}


void cmri::process(const common_options_t &common_options,const motif_count_options_t &motif_count_options,mapVectorMotifRegion &motif_map) {

    cmri::sequenceReader reader(common_options.input_file,motif_count_options.quality_value,motif_count_options.quality_map);
    int total_reads = reader.getTotalReads();
    LOGGER.info << "Processing: " << total_reads << " reads." << std::endl;

    cmri::read_item_t read_item;
    int next_log = common_options.progress;
    while(reader.get(read_item)){
        if(!read_item.valid){continue;}

        if (motif_map.find(read_item.name) != motif_map.end()) {
            for (auto &item : motif_map[read_item.name]) {

                if (!item.intersect(read_item.start, read_item.end)) { continue; }

                std::string sequence = read_item.sequence;
                if (motif_count_options.validate_sequence) { for (auto &c: sequence) { c = toupper(c); }}


                for (auto &m : item.motifs) {
                    m.second += searchMotif(sequence, m.first);
                }
                searchMotifWFilter(sequence,read_item.qvalue,item.motif_quality);
                for (auto &r : item.regex) {
                    r.second += searchRegex(sequence, r.first);
                }
                item.reads_count++;
                item.total_bases+=read_item.sequence.size();
            }
        } else {
            if (motif_map.find("other") != motif_map.end()) {

                std::string sequence = read_item.sequence;
                if (motif_count_options.validate_sequence) { for (auto &c: sequence) { c = toupper(c); }}

                for (auto &item : motif_map["other"]) {
                    for (auto &m : item.motifs) {
                        m.second += searchMotif(sequence, m.first);
                    }
                    searchMotifWFilter(sequence,read_item.qvalue,item.motif_quality);
                    for (auto &r : item.regex) {
                        r.second += searchRegex(sequence, r.first);
                    }
                    item.reads_count++;
                    item.total_bases+=sequence.size();
                }
            }
        }

        if (common_options.progress > 0 && reader.getCount() >= next_log) {
            LOGGER.info << "Progress: " << reader.getCount() << " of " << total_reads << " "
                        << 100.0 * reader.getCount() / total_reads << "%" << std::endl;
            next_log += common_options.progress;
        }

    }
    LOGGER.info << "Total sequences analysed: " << reader.getCount() << std::endl;
    reader.close();
}


void cmri::processMultiThreading(const common_options_t &common_options,const motif_count_options_t &motif_count_options, mapVectorMotifRegion &motif_map) {

    cmri::sequenceReader reader(common_options.input_file,motif_count_options.quality_value,motif_count_options.quality_map);
    int total_reads = reader.getTotalReads();
    LOGGER.info << "Processing: " << total_reads << " reads." << std::endl;

    std::vector<read_item_t> sequences[common_options.threads];

    cmri::read_item_t read_item;
    int next_log = common_options.progress;
    while(reader.get(read_item)){
        if(!read_item.valid){continue;}

        sequences[reader.getCount() % common_options.threads].emplace_back(read_item);

        if (sequences[common_options.threads - 1].size() > common_options.chunk_size) {

            std::vector<std::future<mapVectorMotifRegion  >> pool;
            for (auto &seq_data : sequences) {
                pool.push_back(std::async(std::launch::async, &processWorker, motif_map, seq_data, motif_count_options.validate_sequence));
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

        if (common_options.progress > 0 && reader.getCount() >= next_log) {
            LOGGER.info << "Progress: " << reader.getCount() << " of " << total_reads << " "
                        << 100.0 * reader.getCount() / total_reads << "%" << std::endl;
            next_log += common_options.progress;
        }


    }

    std::vector<std::future<mapVectorMotifRegion>> pool;
    for (auto &seq_data : sequences) {
        pool.push_back(std::async(std::launch::async, &processWorker, motif_map, seq_data, motif_count_options.validate_sequence));
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




void cmri::mainMotifCount(const common_options_t &common_options, const motif_count_options_t &motif_count_options){


    boost::property_tree::ptree input_tree;
    boost::property_tree::read_json(motif_count_options.motif_file, input_tree);
    mapVectorMotifRegion motifs;
    cmri::deserialize(input_tree,motifs);

    if (motifs.empty()) {
        cmri::LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
        cmri::LOGGER.error << "Invalid argument. No valid motif list found!" << std::endl;
        exit (EINVAL);
    }


    if (common_options.threads > 1) {
        cmri::processMultiThreading(common_options, motif_count_options, motifs);
    } else {
        cmri::process(common_options,motif_count_options, motifs);
    }

    std::ofstream file(common_options.output_path + "/output.json");
    file << cmri::serialize(motifs);
    file.close();

}
