//
// Created by pablo on 26/6/20.
//

#include "../MotifCount/motifRegion.h"
#include "variantCallAnalysis.h"
#include "variantRegion.h"
#include <sequenceReader.h>
#include "kseq.h"


std::string
cmri::get_context(bcf1_t *vcf_record, faidx_t *ref_file_index, std::string chromosome, bool reverse, int left_padding,
                  int right_padding) {

    int len;
    std::stringstream region;
    int pos = vcf_record->pos + 1;
    region << chromosome << ":" << pos - left_padding << "-" << pos + right_padding;
    std::string result = fai_fetch(ref_file_index, region.str().c_str(), &len);
    if (reverse) {
        result = reverse_complement(result);
    } else {
        for (auto &c: result) { c = toupper(c); }
    }
    return result;
}

bool cmri::is_reverse(const std::string &allele, const int &variant_type, int index, bool reverse) {
    if (variant_type == VCF_SNP && index == 0) { return allele == "G" || allele == "A"; }
    if (variant_type == VCF_MNP && index == 0) {
        return allele == "GT"
               || allele == "GG"
               || allele == "AG"
               || allele == "GA"
               || allele == "CA"
               || allele == "AA";
    }
    if (variant_type == VCF_INDEL) {
        if (allele.size() > 1)//deletion
        {
            return allele[1] == 'G' || allele[1] == 'A';
        }
        if (allele.size() > 1 && index == 1)//insertion
        {
            return allele[0] == 'G' || allele[0] == 'A';
        }
    }

    return reverse;

}


void cmri::mainVariantCallAnalysis(common_options_t common_options,
                                   variant_call_analysis_options_t variant_call_analysis_options) {


    boost::property_tree::ptree region_data;
    boost::property_tree::read_json(variant_call_analysis_options.regions, region_data);

    std::map<std::string, std::vector<variantRegion>> regions;
    cmri::deserialize(region_data, regions);

    cmri::LOGGER.debug << "Input: " << serialize(regions) << std::endl;

    if (regions.empty()) {
        cmri::LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
        cmri::LOGGER.error << "Invalid argument. No valid motif list found!" << std::endl;
        exit(EINVAL);
    }

    faidx_t *ref_file_index = fai_load(variant_call_analysis_options.reference.c_str());


    htsFile *vcf_file = vcf_open(common_options.input_file.c_str(), "r");
    bcf_hdr_t *vcf_header = vcf_hdr_read(vcf_file);

    bcf1_t *vcf_record = bcf_init1();

    int number_of_samples = bcf_hdr_nsamples(vcf_header);


    for (auto &item_vector: regions) {
        for (auto &item : item_vector.second) {
            for (auto &item_mutation : item.mutations) {
                for (int i = 0; i < number_of_samples; i++) {
                    item_mutation.second[vcf_header->samples[i]] = 0;
                }
            }
        }
    }


    std::map<int, std::string> type_map{
            {0,  "VCF_REF"},
            {1,  "VCF_SNP"},
            {2,  "VCF_MNP"},
            {4,  "VCF_INDEL"},
            {8,  "VCF_OTHER"},
            {16, "VCF_BND"},
            {32, "VCF_OVERLAP"}
    };


    int count = 0;
    while (bcf_read(vcf_file, vcf_header, vcf_record) == 0) {

        int variant_type = bcf_get_variant_types(vcf_record);
        if (variant_type != VCF_SNP && variant_type != VCF_MNP && variant_type != VCF_INDEL) { continue; }

        std::string chromosome = bcf_hdr_id2name(vcf_header, vcf_record->rid);
        if (regions.find(chromosome) == regions.end()) {
            LOGGER.warning << "regions.find(chromosome) == regions.end()" << std::endl;
            continue; }

        for (auto &item : regions[chromosome]) {

            if (!item.intersect(vcf_record->pos)) { continue; }

            std::string filter_str;
            auto filter = bcf_has_filter(vcf_header, vcf_record, "PASS");
            if (filter == 1) { filter_str = "PASS"; }
            else {
                if (filter == 0) { filter_str = "FAIL"; }
                else { filter_str = "NP"; }
            }


            bcf_unpack(vcf_record, BCF_UN_ALL);

            int32_t *gt_arr = NULL, ngt_arr = 0;
            int number_of_genotypes = bcf_get_genotypes(vcf_header, vcf_record, &gt_arr, &ngt_arr);
            if (number_of_genotypes <= 0) {
                free(gt_arr);
                LOGGER.warning << "number_of_genotypes <= 0"<<std::endl;
                continue;
            } // GT not present

            int max_ploidy = number_of_genotypes / number_of_samples;
            for (int i = 0; i < number_of_samples; i++) {
                int32_t *ptr = gt_arr + i * max_ploidy;
                bool reverse = false;
                std::vector<std::string> alleles;
                for (int j = 0; j < max_ploidy; j++) {
                    // if true, the sample has smaller ploidy
                    if (ptr[j] == bcf_int32_vector_end) { break; }

                    // missing allele
                    if (bcf_gt_is_missing(ptr[j])) { LOGGER.warning << "bcf_gt_is_missing(ptr[j])"<<std::endl; continue; }

                    // the VCF 0-based allele index
                    int allele_index = bcf_gt_allele(ptr[j]);

                    std::string allele = vcf_record->d.allele[allele_index];

                    reverse = is_reverse(allele, variant_type, j, reverse);

                    alleles.push_back(allele);
                }

                if (alleles.size() < 2) { LOGGER.warning << "alleles.size() < 2"<<std::endl; continue; }

                std::stringstream key_stream;
                for (int j = 0; j < alleles.size(); j++) {
                    if (reverse) {
                        alleles[j] = reverse_complement(alleles[j]);
                    }
                    key_stream << alleles[j] << (j < max_ploidy - 1 ? ">" : "");
                }


                if (variant_type == VCF_SNP || variant_type == VCF_MNP) {
                    std::string context = variant_type == VCF_SNP ? "_" +
                                                                    get_context(vcf_record, ref_file_index, chromosome,
                                                                                reverse, 1, 1) : "";
                    std::string mutation_key = key_stream.str() + context + ":" + filter_str;

                    if (item.mutations.find(mutation_key) == item.mutations.end()) { continue; }

                    item.mutations[mutation_key][vcf_header->samples[i]] += 1;

                } else {

                    bool insertion = alleles[0].size() < alleles[1].size();
                    std::string ID_type = insertion ? "INS" : "DEL";

                    for (int repeats = 1; repeats < 6; repeats++) {
                        if ((insertion && alleles[1].size() == repeats + 1) ||
                            (!insertion && alleles[0].size() == repeats + 1)) {
                            std::string bp = insertion ? alleles[1].substr(1, repeats + 1) : alleles[0].substr(1,
                                                                                                               repeats +
                                                                                                               1);
                            std::vector<std::string> mutation_keys;
                            std::string context = get_context(vcf_record, ref_file_index, chromosome, reverse,
                                                              -alleles[0].size(), repeats * 20);

                            std::string repeat_key;
                            if (repeats == 1) { repeat_key = "_" + bp + "_" + std::to_string(repeats) + "_"; }
                            else {
                                if (repeats < 5) { repeat_key = "_repeats_" + std::to_string(repeats) + "_"; }
                                else { repeat_key = "_repeats_5+_"; }
                            }
                            int k = 0;
                            for (int j = 0; j < context.size(); j += repeats) {
                                std::string count_key = k < 5 ? std::to_string(k++) : "5+";
                                std::string new_key =
                                        ID_type + repeat_key + count_key +
                                        ":" + filter_str;
                                std::string context_bp = context.substr(j, repeats);
                                if (context_bp != bp) {
                                    mutation_keys.push_back(new_key);
                                    break;
                                }
                            }

                            for (auto mutation_key : mutation_keys) {
                                if (item.mutations.find(mutation_key) == item.mutations.end()) { continue; }
                                item.mutations[mutation_key][vcf_header->samples[i]] += 1;
                            }

                            int pos = vcf_record->pos + 1;
                            LOGGER.debug << chromosome << ":" << pos << (reverse ? " reverse" : "") << std::endl;
                            LOGGER.debug << key_stream.str() << std::endl;
                            LOGGER.debug << context << std::endl;

                            //micro-homologies
                            if (repeats > 1 && !insertion) {
                                for (int size = 1; size < repeats; size++) {

                                    std::string left_context = get_context(vcf_record, ref_file_index, chromosome,
                                                                           reverse,
                                                                           repeats - size-1, 0);
                                    std::string right_context = get_context(vcf_record, ref_file_index, chromosome,
                                                                            reverse,
                                                                            -repeats - 1, 2 * repeats - size);

                                    std::string mutation_str;

                                    if (left_context == bp.substr(size, repeats - size) ||
                                        right_context == bp.substr(0, repeats - size)) {

                                        std::string mutation_key = "DEL_MH_"+(repeats <5?std::to_string(repeats):"5+")+"_"+(repeats-size <5?std::to_string(repeats-size):"5+")+":"+filter_str;
                                        LOGGER.debug << mutation_key << std::endl;
                                        LOGGER.debug << left_context << "." << bp << "." << right_context << std::endl;
                                        LOGGER.debug << bp.substr(size, repeats - size) <<" - " << bp.substr(0, repeats - size) << std::endl;

                                        if (item.mutations.find(mutation_key) != item.mutations.end()) {
                                            item.mutations[mutation_key][vcf_header->samples[i]] += 1;
                                            break;
                                        }
                                    }
                                }

                            }

                            continue;
                        }


                    }
                }


            }

            item.total_bases++;
            cmri::LOGGER.debug
                    << chromosome << " "
                    << item.start << " "
                    << vcf_record->pos << " "
                    << item.end << " "
                    << item.name << " "
                    << item.total_bases << " "
                    << vcf_record->n_sample << " "
                    << vcf_record->n_allele << " "
                    << std::endl;

        }

        //if (count > 10) { break; }
        count++;
    }

    bcf_hdr_destroy(vcf_header);
    vcf_close(vcf_file);

    std::ofstream file(common_options.output_path + "/output.json");
    file << cmri::serialize(regions);
    file.close();


}
