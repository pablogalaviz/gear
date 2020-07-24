//
// Created by pablo on 26/6/20.
//

#ifndef GEAR_VARIANTCALLANALYSIS_H
#define GEAR_VARIANTCALLANALYSIS_H

#include "faidx.h"
#include <string>
#include "options.h"
#include "vcf.h"

namespace cmri {

    struct vc_record_t {
        std::string reference;
        std::vector<std::string> variant;
    };



    std::string get_context(bcf1_t *vcf_record, faidx_t *ref_file_index,std::string chromosome,bool reverse,int left_padding, int right_padding);

    bool is_reverse(const std::string &allele, const int &variant_type,int indx, bool reverse);

    void mainVariantCallAnalysis(common_options_t common_options, variant_call_analysis_options_t variant_call_analysis_options);

}
#endif //GEAR_VARIANTCALLANALYSIS_H
