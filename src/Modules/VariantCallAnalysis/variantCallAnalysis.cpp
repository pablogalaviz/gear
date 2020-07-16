//
// Created by pablo on 26/6/20.
//

#include "../MotifCount/motifRegion.h"
#include "faidx.h"
#include "variantCallAnalysis.h"
#include "vcf.h"
#include "variantRegion.h"
#include <zlib.h>
#include <sequenceReader.h>
#include "kseq.h"

std::string get_snp_key() {
    std::string result;


    return result;
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

        if (bcf_get_variant_types(vcf_record) != VCF_SNP && bcf_get_variant_types(vcf_record) != VCF_MNP) { continue; }

        std::string chromosome = bcf_hdr_id2name(vcf_header, vcf_record->rid);
        if (regions.find(chromosome) == regions.end()) { continue; }

        for (auto &item : regions[chromosome]) {

            if (!item.intersect(vcf_record->pos)) { continue; }

            bcf_unpack(vcf_record, BCF_UN_ALL);

            if (bcf_get_variant_types(vcf_record) == VCF_MNP) {
                cmri::LOGGER.debug << "variant type: " << type_map[bcf_get_variant_types(vcf_record)] << std::endl;
            }
            /*
            LOGGER.debug << "cvf_record allele ";
            for (int i = 0; i < vcf_record->d.m_allele; i++) {
                LOGGER.debug << vcf_record->d.allele[i] << " ";
            }
            LOGGER.debug << std::endl;
*/
            int32_t *gt_arr = NULL, ngt_arr = 0;
            int number_of_genotypes = bcf_get_genotypes(vcf_header, vcf_record, &gt_arr, &ngt_arr);
            if (number_of_genotypes <= 0) {
                free(gt_arr);
                continue;
            } // GT not present


            int max_ploidy = number_of_genotypes / number_of_samples;
            for (int i = 0; i < number_of_samples; i++) {
                int32_t *ptr = gt_arr + i * max_ploidy;
                std::stringstream key_stream;
                bool reverse = false;
                for (int j = 0; j < max_ploidy; j++) {
                    // if true, the sample has smaller ploidy
                    if (ptr[j] == bcf_int32_vector_end) break;

                    // missing allele
                    if (bcf_gt_is_missing(ptr[j])) continue;

                    // the VCF 0-based allele index
                    int allele_index = bcf_gt_allele(ptr[j]);

                    // is phased?
                    int is_phased = bcf_gt_is_phased(ptr[j]);

                    std::string allele = vcf_record->d.allele[allele_index];
                    if (j == 0 && (
                            allele == "G"
                            || allele == "A"
                            || allele == "GT"
                            || allele == "GG"
                            || allele == "AG"
                            || allele == "GA"
                            || allele == "CA"
                            || allele == "AA"
                    )) {
                        reverse = true;
                    }
                    if (reverse) {
                        allele = reverse_complement(allele);
                    }

                    key_stream << allele << (j < max_ploidy - 1 ? ">" : "");
                }

                int len;
                std::stringstream region;
                int pos = vcf_record->pos + 1;
                region << chromosome << ":" << pos - 1 << "-" << pos + 1;
                std::string seq = fai_fetch(ref_file_index, region.str().c_str(), &len);
                if (reverse) {
                    LOGGER.debug << "to reverse seq " << seq << std::endl;
                    seq = reverse_complement(seq);
                }

                LOGGER.debug << "region " << region.str() << std::endl;
                LOGGER.debug << "len " << len << std::endl;
                LOGGER.debug << "mutation key " << key_stream.str() << std::endl;
                LOGGER.debug << "subtype " << seq << std::endl;

                std::string mutation_key =
                        bcf_get_variant_types(vcf_record) == VCF_SNP ? key_stream.str()+"_"+ seq
                                                                     : key_stream.str();

                if (item.mutations.find(mutation_key) == item.mutations.end()) { continue; }
                item.mutations[mutation_key][vcf_header->samples[i]] += 1;
                if(bcf_get_variant_types(vcf_record) == VCF_MNP ) {
                    cmri::LOGGER.debug << vcf_header->samples[i] << " " << key_stream.str() << std::endl;
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
