//
// Created by pablo on 30/10/20.
//

#ifndef GEAR_TELOMEREMUTATIONS_H
#define GEAR_TELOMEREMUTATIONS_H

#include <regex>
#include "options.h"
#include "utils.h"

namespace cmri {

    struct sbs_t {
        int pos;
        char value;
        int qv;
        std::string serialize() const {
            std::stringstream result;
            result << "{";
            result << "\"pos\":" << pos << ",";
            result << "\"qv\":" << qv << ",";
            result << "\"value\":\"" << value << "\"";
            result << "}";
            return result.str();
        };
    };

    inline std::ostream& operator<<(std::ostream& result, const sbs_t& rhs)
    {
        result << rhs.serialize();
        return result;
    }

    struct indel_t{
        std::string kind;
        int pos;
        std::string seq;
        double mean_qv;
        std::string serialize() const {
            std::stringstream result;
            result << "{";
            result << "\"pos\":" << pos << ",";
            result << "\"mean_qv\":" << mean_qv << ",";
            result << "\"seq\":\"" << seq << "\",";
            result << "\"kind\":\"" << kind << "\"";
            result << "}";
            return result.str();
        };

    };

    inline std::ostream& operator<<(std::ostream& result, const indel_t& rhs)
    {
        result << rhs.serialize();
        return result;
    }


    struct mutations_t {
        std::string name;
        std::string comment;
        int rs;
        int re;
        int qs;
        int qe;
        int seq_len;
        int mapq;
        int mlen;
        int blen;
        double score=0;
        std::string seq;
        std::string qv;
        std::map<int, std::vector<sbs_t>> sbs;
        std::map<int, std::vector<indel_t>> indels;

        std::map<std::string,int> count;
        std::map<std::string,int> variants;
        std::vector<int> sequence_position;

        double mean_ins_size=0;
        double mean_ins_qv=0;
        int ins_count=0;

        double mean_del_size=0;
        double mean_del_qv=0;
        int del_count=0;
        int reverse;

        std::string cs_str;


        std::string serialize() const {
            std::stringstream result;
            result << "{";
            result << "\"rs\":" << rs << ",";
            result << "\"re\":" << re << ",";
            result << "\"qs\":" << qs << ",";
            result << "\"qe\":" << qe << ",";
            result << "\"seq_len\":" << seq_len << ",";
            result << "\"mapq\":" << mapq << ",";
            result << "\"mlen\":" << mlen << ",";
            result << "\"blen\":" << blen << ",";
            result << "\"score\":" << score << ",";
            result << "\"reverse\":" << reverse << ",";

            /*
            for(auto &item : count){
                result << "\"" <<item.first <<"\":" << item.second << ",";
            }
            */
            for(auto &item : variants){
                result << "\"" <<item.first <<"\":" << item.second << ",";
            }

            result << "\"mean_ins_size\":" << mean_ins_size << ",";
            result << "\"mean_ins_qv\":" << mean_ins_qv << ",";
            result << "\"ins_count\":" << ins_count << ",";

            result << "\"mean_del_size\":" << mean_del_size << ",";
            result << "\"mean_del_qv\":" << mean_del_qv << ",";
            result << "\"del_count\":" << del_count << ",";

            result << "\"sbs\":" << ::cmri::serialize(sbs) << ",";
            result << "\"indels\":" << ::cmri::serialize(indels) << ",";

            result << "\"seq\":\"" << seq << "\",";
            result << "\"qv\":\"" << qv << "\",";
            result << "\"cs_str\":\"" << cs_str << "\",";
            result << "\"comment\":\"" << comment << "\",";
            result << "\"name\":\"" << name << "\"";
            result << "}";
            return result.str();
        };

        void find_mutations(){

            for(auto &item : sbs){

                if(item.second.size()>1){
                    continue;
                }

                char motif[] = {'T','T','A','G','G','G'};
                for(auto &s : item.second){
                    motif[s.pos]=s.value;
                }

                std::string mutation = "";
                for(auto &c : motif){
                    mutation+=c;
                }

                if(count.find(mutation) == count.end()){
                    count[mutation]=1;
                }
                else{
                    count[mutation]+=1;
                }
                sequence_position.push_back(item.first);

            }

            for(auto &item : indels){

                for(auto &i : item.second){
                    if(i.kind == "ins"){
                        ins_count++;
                        mean_ins_size+=i.seq.size();
                        mean_ins_qv+=i.mean_qv;
//                        if(i.seq.size()==1){

//                        }

                    }
                    else{
                        del_count++;
                        mean_del_size+=i.seq.size();
                        mean_del_qv+=i.mean_qv;
                    }
                    LOGGER.debug << i.serialize() << std::endl;
                }

                if (ins_count > 0){
                    mean_ins_qv/= ins_count;
                    mean_ins_size/=ins_count;
                }
                if (del_count > 0){
                    mean_del_qv/= del_count;
                    mean_del_size/=del_count;
                }

            }

        }


    };

    inline std::ostream& operator<<(std::ostream& result, const mutations_t& rhs)
    {
        result << rhs.serialize();
        return result;
    }


inline   std::map<std::string,int> find_variants(std::string sequence){
        std::map<std::string,int> result;

        std::regex basic_regex("TTAGGG|ATAGGG|GTAGGG|CTAGGG|TAAGGG|TGAGGG|TCAGGG|TTCGGG|TTTGGG|TTGGGG|TTACGG|TTATGG|TTAAGG|TTAGCG|TTAGTG|TTAGAG|TTAGGC|TTAGGT|TTAGGA");
        std::smatch match;

        while (std::regex_search(sequence, match, basic_regex)) {
            // suffix to find the rest of the string.
            std::string prefix = match.prefix().str();
            if (!prefix.empty()) {
                if (result.find(prefix) == result.end()) {
                    result[prefix] = 1;
                } else {
                    result[prefix] += 1;
                }
            }
            sequence = match.suffix().str();
        }
/*
        for(auto & item : result){
            LOGGER.debug << item.first << " : " << item.second << std::endl;
        }
*/
        return result;
    }


    int mainTelomereMutations(common_options_t common_options, telomere_mutations_options_t telomere_mutation_options );


}


#endif //GEAR_TELOMEREMUTATIONS_H
