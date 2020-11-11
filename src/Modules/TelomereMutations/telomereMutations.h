//
// Created by pablo on 30/10/20.
//

#ifndef GEAR_TELOMEREMUTATIONS_H
#define GEAR_TELOMEREMUTATIONS_H

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
        int rs;
        int re;
        int mapq;
        int mlen;
        int blen;
        double score=0;
        std::map<int, std::vector<sbs_t>> sbs;
        std::map<int, std::vector<indel_t>> indels;

        std::map<std::string,int> count;
        std::vector<int> sequence_position;

        double mean_ins_size=0;
        double mean_ins_qv=0;
        int ins_count=0;

        double mean_del_size=0;
        double mean_del_qv=0;
        int del_count=0;

        std::string serialize() const {
            std::stringstream result;
            result << "{";
            result << "\"rs\":" << rs << ",";
            result << "\"re\":" << re << ",";
            result << "\"mapq\":" << mapq << ",";
            result << "\"mlen\":" << mlen << ",";
            result << "\"blen\":" << blen << ",";
            result << "\"score\":" << score << ",";
            result << "\"count\":" << cmri::serialize(count) << ",";

            result << "\"mean_ins_size\":" << mean_ins_size << ",";
            result << "\"mean_ins_qv\":" << mean_ins_qv << ",";
            result << "\"ins_count\":" << ins_count << ",";

            result << "\"mean_del_size\":" << mean_del_size << ",";
            result << "\"mean_del_qv\":" << mean_del_qv << ",";
            result << "\"del_count\":" << del_count << ",";

//            result << "\"sbs\":" << ::cmri::serialize(sbs) << ",";
//            result << "\"indels\":" << ::cmri::serialize(indels) << ",";

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


    int mainTelomereMutations(common_options_t common_options, telomere_mutations_options_t telomere_mutation_options );

   std::vector<std::pair<char,std::string>> parseTag(std::string tag);

}


#endif //GEAR_TELOMEREMUTATIONS_H
