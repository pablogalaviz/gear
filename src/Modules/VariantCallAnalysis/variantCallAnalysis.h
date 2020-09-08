//
// Created by pablo on 26/6/20.
//

#ifndef GEAR_VARIANTCALLANALYSIS_H
#define GEAR_VARIANTCALLANALYSIS_H

#include "faidx.h"
#include <string>
#include "options.h"
#include "vcf.h"
#include "variantCallRecord.h"

namespace cmri {


    std::string get_context(faidx_t *ref_file_index,std::string chromosome,int position,int size, int left_padding, int right_padding);

    bool is_reverse(const std::string &allele, const int &variant_type,int indx, bool reverse);

    bool is_reverse(const std::string &allele, const int &variant_type);

    std::string get_key_from_snp(faidx_t *ref_file_index, std::string reference, std::string alternate, variantCallRecord variant_call_record);
    std::string get_key_from_mnp(std::string reference, std::string alternate, variantCallRecord variant_call_record);
    std::string get_key_from_insertion(faidx_t *ref_file_index, std::string reference, std::string alternate, variantCallRecord variant_call_record);
    std::string get_key_from_deletion(faidx_t *ref_file_index, std::string reference, std::string alternate, variantCallRecord variant_call_record);

    void mainVariantCallAnalysis(common_options_t common_options, variant_call_analysis_options_t variant_call_analysis_options);


}
#endif //GEAR_VARIANTCALLANALYSIS_H
