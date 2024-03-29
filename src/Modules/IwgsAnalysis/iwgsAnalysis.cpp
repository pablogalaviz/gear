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

#include <boost/filesystem.hpp>
#include "iwgsAnalysis.h"
#include "sequenceClassification.h"
#include <zlib.h>
#include <boost/property_tree/json_parser.hpp>
#include <future>
#include <random>
#include <iterator>
#include <algorithm>


std::pair<std::string,std::string> cmri::motifFilterWorker(const std::vector<std::pair<my_kseq_t,my_kseq_t>> &sequences,
                                                           iwgsAnalysis iwgs ){

    std::stringstream result1;
    std::stringstream result2;

    for(auto &item : sequences) {
        if (iwgs.consecutive_count(item.first.name, item.first.seq, item.first.qual)  ||
            iwgs.consecutive_count(item.second.name, item.second.seq, item.second.qual) ) {
            result1 << "@" << item.first.name << " " << item.first.comment << std::endl;
            result1 << item.first.seq << std::endl;
            result1 << "+" << std::endl;
            result1 << item.first.qual << std::endl;

            result2 << "@" << item.second.name << " " << item.second.comment << std::endl;
            result2 << item.second.seq << std::endl;
            result2 << "+" << std::endl;
            result2 << item.second.qual << std::endl;

        }
    }

    return std::make_pair(result1.str(),result2.str());
}


std::map<char,double> cmri::qv_histogram(const std::string &sequence, const std::string &quality) {

    std::map<char,double> result{{'A',0},{'C',0},{'G',0},{'T',0},{'N',0}};

    int qv_th=20;

    for(int i=0; i < sequence.size(); i++){
        int qv = (int)quality[sequence[i]]-33;
        if( qv < qv_th) { result[sequence[i]]++; }
    }

    return result;
}

std::pair<std::string,std::string> cmri::qvFilterWorker(const std::vector<std::pair<my_kseq_t,my_kseq_t>> &sequences){

    std::stringstream result1;
    std::stringstream result2;

    int count_th=30;

    for(auto &item : sequences) {

        std::map<char,double> seq1_hist = qv_histogram(item.first.seq,item.first.qual);
        std::map<char,double> seq2_hist = qv_histogram(item.second.seq,item.second.qual);

        bool low_qv=false;
        std::stringstream qv_comment1;
        for(auto & sh : seq1_hist){
            if(sh.second > count_th){
                low_qv=true;
                qv_comment1 << "|" << sh.first << '-' << sh.second;
            }
        }

        std::stringstream qv_comment2;
        for(auto & sh : seq2_hist){
            if(sh.second > count_th){
                low_qv=true;
                qv_comment2 << "|" << sh.first << '-' << sh.second;
            }
        }

        if(low_qv){
            result1 << "@" << item.first.name << " " << item.first.comment <<":" << qv_comment1.str() << std::endl;
            result1 << item.first.seq << std::endl;
            result1 << "+" << std::endl;
            result1 << item.first.qual << std::endl;

            result2 << "@" << item.second.name << " " << item.second.comment <<":" << qv_comment2.str()  << std::endl;
            result2 << item.second.seq << std::endl;
            result2 << "+" << std::endl;
            result2 << item.second.qual << std::endl;
        }


    }

    return std::make_pair(result1.str(),result2.str());
}

void cmri::mainIwgsAnalysis(const common_options_t &common_options,
                              const iwgs_analysis_options_t &iwgs_options) {

        gzFile ifile1 = gzopen(common_options.input_file.c_str(), "r");
        kseq_t *kseq1 = kseq_init(ifile1);

        gzFile ifile2 = gzopen(iwgs_options.input_file.c_str(), "r");
        kseq_t *kseq2 = kseq_init(ifile2);

        int l1,l2;
        int count =0;

        boost::filesystem::path output1(common_options.output_path);
        output1.append(boost::filesystem::basename(common_options.input_file));
        std::ofstream ofile1(output1.string(), std::ios_base::out);

        boost::filesystem::path output2(common_options.output_path);
        output2.append(boost::filesystem::basename(iwgs_options.input_file));
        std::ofstream ofile2(output2.string(), std::ios_base::out);

        boost::property_tree::ptree input_tree;
        boost::property_tree::read_json(iwgs_options.variants_file, input_tree);

        std::vector<std::string> variants;
        for(auto &item : input_tree.get_child("variants")){
            variants.push_back(item.second.get_value<std::string>());
        }

        iwgsAnalysis iwgs("TTAGGG", variants, iwgs_options.count_filter_threshold);

        std::vector< std::pair<my_kseq_t,my_kseq_t> > sequences[common_options.threads];

        int next_log = common_options.progress;
        while ((l1 = kseq_read(kseq1)) >= 0 && (l2 = kseq_read(kseq2)) >= 0 ) {

            my_kseq_t seq1, seq2;
            seq1.comment = kseq1->comment.s;
            seq1.name = kseq1->name.s;
            seq1.seq = kseq1->seq.s;
            seq1.qual = kseq1->qual.s;

            seq2.comment = kseq2->comment.s;
            seq2.name = kseq2->name.s;
            seq2.seq = kseq2->seq.s;
            seq2.qual = kseq2->qual.s;

            sequences[count % common_options.threads].emplace_back(std::make_pair(seq1,seq2));

            if (sequences[common_options.threads - 1].size() > common_options.chunk_size) {
                std::vector<std::future< std::pair<std::string,std::string>  >> pool;

                for (auto &seq_data : sequences) {
                    pool.push_back(std::async(std::launch::async, &motifFilterWorker, seq_data, iwgs));
                    seq_data.clear();
                }
                for (auto &t : pool) {
                    t.wait();
                    auto result = t.get();
                    ofile1 << result.first;
                    ofile1.flush();
                    ofile2 << result.second;
                    ofile2.flush();
                }


            }

            count++;

            if (common_options.progress > 0 && count >= next_log) {
                LOGGER.info << "Progress: " << count << " reads " << std::endl;
                next_log += common_options.progress;
            }


            /*
           if(count > 10){
                break;
            }
            */
        }

        LOGGER.info << "Total number of sequences: " << count << std::endl;

        ofile1.close();
        ofile2.close();

}


void cmri::mainQVSelector(const common_options_t &common_options, const iwgs_analysis_options_t &iwgs_options) {

    gzFile ifile1 = gzopen(common_options.input_file.c_str(), "r");
    kseq_t *kseq1 = kseq_init(ifile1);

    gzFile ifile2 = gzopen(iwgs_options.input_file.c_str(), "r");
    kseq_t *kseq2 = kseq_init(ifile2);

    int l1,l2;
    int count =0;

    boost::filesystem::path output1(common_options.output_path);
    output1.append(boost::filesystem::basename(common_options.input_file));
    std::ofstream ofile1(output1.string(), std::ios_base::out);

    boost::filesystem::path output2(common_options.output_path);
    output2.append(boost::filesystem::basename(iwgs_options.input_file));
    std::ofstream ofile2(output2.string(), std::ios_base::out);

    boost::property_tree::ptree input_tree;
    boost::property_tree::read_json(iwgs_options.variants_file, input_tree);

    std::vector< std::pair<my_kseq_t,my_kseq_t> > sequences[common_options.threads];

    int next_log = common_options.progress;
    while ((l1 = kseq_read(kseq1)) >= 0 && (l2 = kseq_read(kseq2)) >= 0 ) {

        my_kseq_t seq1, seq2;
        seq1.comment = kseq1->comment.s;
        seq1.name = kseq1->name.s;
        seq1.seq = kseq1->seq.s;
        seq1.qual = kseq1->qual.s;

        seq2.comment = kseq2->comment.s;
        seq2.name = kseq2->name.s;
        seq2.seq = kseq2->seq.s;
        seq2.qual = kseq2->qual.s;

        sequences[count % common_options.threads].emplace_back(std::make_pair(seq1,seq2));

        if (sequences[common_options.threads - 1].size() > common_options.chunk_size) {
            std::vector<std::future< std::pair<std::string,std::string>  >> pool;

            for (auto &seq_data : sequences) {
                pool.push_back(std::async(std::launch::async, &qvFilterWorker, seq_data));
                seq_data.clear();
            }
            for (auto &t : pool) {
                t.wait();
                auto result = t.get();
                ofile1 << result.first;
                ofile1.flush();
                ofile2 << result.second;
                ofile2.flush();
            }


        }

        count++;

        if (common_options.progress > 0 && count >= next_log) {
            LOGGER.info << "Progress: " << count << " reads " << std::endl;
            next_log += common_options.progress;
        }

    }

    LOGGER.info << "Total number of sequences: " << count << std::endl;

    ofile1.close();
    ofile2.close();



}

void cmri::mainRandomSelector(const common_options_t &common_options) {

    gzFile ifile1 = gzopen(common_options.input_file.c_str(), "r");
    kseq_t *kseq1 = kseq_init(ifile1);

    int l1;
    int count =0;

    boost::filesystem::path output1(common_options.output_path);
    output1.append(boost::filesystem::basename(common_options.input_file));
    std::ofstream ofile1(output1.string(), std::ios_base::out);

    double read_size = 36444356; //count_reads(common_options.input_file);
    double sample_size = 15000;
    double probability = RAND_MAX*sample_size/read_size;
    /* initialize random seed: */
    srand (time(NULL));

    int sample_count=0;

    int next_log = common_options.progress;
    while ((l1 = kseq_read(kseq1)) >= 0  ) {

        if(rand() < probability) {

            ofile1 << "@" << kseq1->name.s << " " << kseq1->comment.s << std::endl;
            ofile1 << kseq1->seq.s << std::endl;
            ofile1 << "+" << std::endl;
            ofile1 << kseq1->qual.s << std::endl;
            ofile1.flush();
            sample_count++;
            if(sample_count>sample_size){
                break;
            }
        }

        count++;
        if (common_options.progress > 0 && count >= next_log) {
            LOGGER.info << "Progress: " << count << " reads " << std::endl;
            next_log += common_options.progress;
        }

    }

    LOGGER.info << "Total number of sequences: " << count << std::endl;

    ofile1.close();

}


bool cmri::iwgsAnalysis::consecutive_count(const std::string &name, const std::string &sequence, const std::string &quality){


    int forward_consecutive=0;
    int current =0;

    std::string::size_type start = 0;
    while (start < sequence.size()) {

        bool found=false;
        for(auto &motif : variants_forward){
            if(sequence.substr(start,motif.size())==motif){
                current++;
                start += motif.size();
                found=true;
                break;
            }
        }
        if(!found) {
            forward_consecutive = std::max(forward_consecutive, current);
            current = 0;
            start++;
        }

    }


    int reverse_consecutive=0;
    current =0;

    start = 0;
    while (start < sequence.size()) {

        bool found=false;
        for(auto &motif : variants_reverse){
            if(sequence.substr(start,motif.size())==motif){
                current++;
                start += motif.size();
                found=true;
                break;
            }
        }
        if(!found) {
            reverse_consecutive = std::max(reverse_consecutive, current);
            current = 0;
            start++;
        }

    }

    /*
    if(forward_consecutive > 5 or reverse_consecutive > 5) {
        LOGGER.debug << "sequence id " << name << std::endl;
        LOGGER.debug << "seq " << sequence << std::endl;
        LOGGER.debug << "quality " << quality << std::endl;
        LOGGER.debug << "current: " << current << " forward_consecutive " << forward_consecutive
                     << " reverse_consecutive " << reverse_consecutive << std::endl;
    }
    */

    return std::max(forward_consecutive,reverse_consecutive) >= count_filter_threshold;

}



cmri::iwgsAnalysis::iwgsAnalysis(const std::string &motif, const std::vector<std::string> &variants, int countFilterThreshold) :
        motif_forward(motif)
        ,variants_forward(variants)
        ,count_filter_threshold(countFilterThreshold) {

    motif_reverse = reverse_complement(motif);

    for(auto &variants : variants_forward){
        variants_reverse.push_back(reverse_complement(variants));
    }


}
