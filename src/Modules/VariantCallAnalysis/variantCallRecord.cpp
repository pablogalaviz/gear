//
// Created by pablo on 24/7/20.
//

#include "logger.h"
#include "variantCallRecord.h"

cmri::variantCallRecord::variantCallRecord(bcf_hdr_t *vcf_header, bcf1_t *vcf_record) {

    chromosome = bcf_hdr_id2name(vcf_header, vcf_record->rid);
    position = vcf_record->pos;
    auto filter_value = bcf_has_filter(vcf_header, vcf_record, "PASS");
    if (filter_value == 1) {
        filter = "PASS";
    } else {
        if (filter_value == 0) {
            filter = "FAIL";
        } else {
            filter = "NP";
        }
    }

    bcf_unpack(vcf_record, BCF_UN_ALL);

    int32_t *gt_arr = NULL, ngt_arr = 0;
    int number_of_genotypes = bcf_get_genotypes(vcf_header, vcf_record, &gt_arr, &ngt_arr);
    if (number_of_genotypes <= 0) {
        valid = false;// GT not present
        free(gt_arr);
        return;
    }

    int number_of_samples = bcf_hdr_nsamples(vcf_header);
    int max_ploidy = number_of_genotypes / number_of_samples;


    for (int i = 0; i < number_of_samples; i++) {
        int32_t *ptr = gt_arr + i * max_ploidy;
        variant_t new_variant;
        new_variant.is_variant= false;
        for (int j = 0; j < max_ploidy; j++) {
            // if true, the sample has smaller ploidy
            if (ptr[j] == bcf_int32_vector_end) { break; }

            // missing allele
            if (bcf_gt_is_missing(ptr[j])) { continue; }

            // the VCF 0-based allele index
            int allele_index = bcf_gt_allele(ptr[j]);
            std::string allele = vcf_record->d.allele[allele_index];

            if (j == 0) {
                new_variant.reference = allele;
            } else {
                new_variant.variants.push_back(allele);
                if(new_variant.reference != allele){new_variant.is_variant=true;}
            }

        }
        samples[vcf_header->samples[i]]=new_variant;

    }

    free(gt_arr);

    valid=true;
}









