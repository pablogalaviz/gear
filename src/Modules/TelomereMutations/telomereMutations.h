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
        std::map<int, std::vector<sbs_t>> sbs;
        std::map<int, std::vector<indel_t>> indels;

        std::string serialize() const {
            std::stringstream result;
            result << "{";
            result << "\"rs\":" << rs << ",";
            result << "\"re\":" << re << ",";
            result << "\"mapq\":" << mapq << ",";

            result << "\"sbs\":" << ::cmri::serialize(sbs) << ",";

            result << "\"name\":\"" << name << "\"";
            result << "}";
            return result.str();
        };

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
