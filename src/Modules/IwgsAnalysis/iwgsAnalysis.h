
#ifndef GEAR_IWGSANALYSIS_H
#define GEAR_IWGSANALYSIS_H

#include "options.h"

namespace cmri {

    class iwgsAnalysis {

        std::string motif_forward;
        std::string motif_reverse;

        std::vector<std::string> variants_forward;
        std::vector<std::string> variants_reverse;

        int count_filter_threshold;
    public:
        iwgsAnalysis(const std::string &motif, const std::vector<std::string> &variants, int countFilterThreshold);


        bool count_filter(std::string name, std::string sequence, std::string quality);

    };



    void mainIwgsAnalysis(const common_options_t &common_options,
                                 const iwgs_analysis_options_t &iwgs_options);
}

#endif //GEAR_IWGSANALYSIS_H
