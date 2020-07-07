//
// Created by pablo on 26/6/20.
//

#ifndef GEAR_VARIANTCALLANALYSIS_H
#define GEAR_VARIANTCALLANALYSIS_H

#include <string>
#include "options.h"

namespace cmri {

    struct region_field {
        std::string name;
        std::string chromosome;
        int start;
        int end;

    };

    void mainVariantCallAnalysis(common_options_t common_options, variant_call_analysis_options_t variant_call_analysis_options);



}
#endif //GEAR_VARIANTCALLANALYSIS_H