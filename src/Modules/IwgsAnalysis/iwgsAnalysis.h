
#ifndef GEAR_IWGSANALYSIS_H
#define GEAR_IWGSANALYSIS_H

#include "options.h"
#include "sequenceReader.h"

namespace cmri {

    struct my_kseq_t{
        std::string name;
        std::string seq;
        std::string qual;
        std::string comment;
    };


    class iwgsAnalysis {

        std::string motif_forward;
        std::string motif_reverse;

        std::vector<std::string> variants_forward;
        std::vector<std::string> variants_reverse;

        int count_filter_threshold;
    public:
        iwgsAnalysis(const std::string &motif, const std::vector<std::string> &variants, int countFilterThreshold);
        bool consecutive_count(const std::string &name, const std::string &sequence, const std::string &quality);

        bool count_filter(std::string name, std::string sequence, std::string quality);

    };

    std::pair<std::string,std::string> filterWorker(const std::vector<std::pair<my_kseq_t,my_kseq_t>> &sequences,
                                                    iwgsAnalysis iwgs );

    void mainIwgsAnalysis(const common_options_t &common_options,
                                 const iwgs_analysis_options_t &iwgs_options);
}


#endif //GEAR_IWGSANALYSIS_H
