//
// Created by pablo on 26/6/20.
//

#include "../MotifCount/motifRegion.h"
#include "variantCallAnalysis.h"
#include "variantRegion.h"
#include <sequenceReader.h>
#include "kseq.h"
#include "variantCallRecord.h"


std::string
cmri::get_context(faidx_t *ref_file_index, std::string chromosome, int position, int size, int left_padding,
                  int right_padding) {
    int len;
    std::string region = chromosome + ":"
                         + std::to_string(position - left_padding + 1) + "-"
                         + std::to_string(position + size + right_padding);
    std::string result = fai_fetch(ref_file_index, region.c_str(), &len);
    for (auto &c: result) { c = toupper(c); }
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

bool cmri::is_reverse(const std::string &allele, const int &variant_type) {
    if (variant_type == VCF_SNP) { return allele == "G" || allele == "A"; }
    if (variant_type == VCF_MNP) {
        return allele == "GT"
               || allele == "GG"
               || allele == "AG"
               || allele == "GA"
               || allele == "CA"
               || allele == "AA";
    }
    return false;
}

std::string cmri::get_key_from_snp(faidx_t *ref_file_index, std::string reference, std::string alternate,
                                   variantCallRecord variant_call_record) {

    std::string result;
    std::string context = get_context(ref_file_index, variant_call_record.getChromosome(),
                                      variant_call_record.getPosition(), 1, 1, 1);

    bool reverse = reference == "G" || reference == "A";
    if (reverse) {
        result = reverse_complement(reference) + ">" + reverse_complement(alternate);
        result += "_" + reverse_complement(context);
    } else {
        result = reference + ">" + alternate;
        result += "_" + context;
    }
    result += ":" + variant_call_record.getFilter();

    return result;

}

std::string cmri::get_key_from_mnp(std::string reference, std::string alternate,
                                   variantCallRecord variant_call_record) {

    std::string result;

    bool reverse = reference == "G" || reference == "A";
    if (reverse) {
        result = reverse_complement(reference) + ">" + reverse_complement(alternate);
    } else {
        result = reference + ">" + alternate;
    }
    result += ":" + variant_call_record.getFilter();
    return result;
}


std::string cmri::get_key_from_insertion(faidx_t *ref_file_index, std::string reference, std::string alternate,
                                         variantCallRecord variant_call_record) {

    std::string result;
    int ins_size = alternate.size() - reference.size();

    for (int repeats = 1; repeats < 6; repeats++) {
        if (ins_size == repeats) {
            std::string bp = alternate.substr(reference.size(), ins_size);
            bool reverse = (bp[0] == 'G' || bp[0] == 'A') && repeats == 1;

            std::string left_context = get_context(ref_file_index, variant_call_record.getChromosome(),
                                                   variant_call_record.getPosition(), reference.size(),
                                                   repeats * 6, 0);
            std::string right_context = get_context(ref_file_index, variant_call_record.getChromosome(),
                                                    variant_call_record.getPosition() + reference.size(), 0,
                                                    0, repeats * 6);
            if (reverse) {
                bp = reverse_complement(bp);
                left_context = reverse_sequence(left_context);
                right_context = reverse_sequence(right_context);
            }

            int matches = 0;
            for (int i = left_context.size() - repeats; i > 0; i -= repeats) {
                if (bp != left_context.substr(i, repeats)) { break; }
                matches++;
            }

            for (int i = 0; i < right_context.size(); i += repeats) {
                if (bp != right_context.substr(i, repeats)) { break; }
                matches++;
            }


            std::string repeats_str = repeats < 5 ? std::to_string(repeats) : "5+";
            std::string matches_str = matches < 5 ? std::to_string(matches) : "5+";

            result = "INS_" + (repeats == 1 ? bp : "repeats") + "_" + repeats_str + "_" + matches_str + ":" +
                     variant_call_record.getFilter();

            LOGGER.debug << variant_call_record.getChromosome() << ":"
                         << variant_call_record.getPosition() << std::endl;
            LOGGER.debug << left_context << "." << bp << "." << right_context << std::endl;
            LOGGER.debug << "result: " << result << std::endl;

            break;
        }
    }

    return result;
}

std::string cmri::get_key_from_deletion(faidx_t *ref_file_index, std::string reference, std::string alternate,
                                        variantCallRecord variant_call_record) {

    std::string result;
    int deletion_size = reference.size() - alternate.size();


    std::string bp = reference.substr(alternate.size(), deletion_size);
    bool reverse = (bp[0] == 'G' || bp[0] == 'A') && deletion_size == 1;

    std::string left_context = get_context(ref_file_index, variant_call_record.getChromosome(),
                                           variant_call_record.getPosition(), alternate.size(),
                                           deletion_size * 6, 0);
    std::string right_context = get_context(ref_file_index, variant_call_record.getChromosome(),
                                            variant_call_record.getPosition() + reference.size(), 0,
                                            0, deletion_size * 6);


    if (reverse) {
        bp = reverse_complement(bp);
        left_context = reverse_sequence(left_context);
        right_context = reverse_sequence(right_context);
    }

    int matches = 0;
    for (int repeat_size = 0; repeat_size < right_context.size(); repeat_size += deletion_size) {
        if (bp != right_context.substr(repeat_size, deletion_size)) { break; }
        matches++;
    }

    for (int repeat_size = left_context.size() - deletion_size; repeat_size > 0; repeat_size -= deletion_size) {
        if (bp != left_context.substr(repeat_size, deletion_size)) { break; }
        matches++;
    }

    int homology = 0;
    if (deletion_size > 1 && matches == 0) {
        int right_homology = 0;
        for (int homology_size = 1; homology_size < deletion_size ; homology_size++) {
            if (bp.substr(0, homology_size) != right_context.substr(0, homology_size)) { break; }
            right_homology = homology_size;
        }
        int left_homology = 0;
        for (int homology_size = 1; homology_size < deletion_size ; homology_size++) {
            if (bp.substr(bp.size() - homology_size, homology_size) !=
                left_context.substr(left_context.size() - homology_size, homology_size)) { break; }
            left_homology = homology_size;
        }
        homology = std::max(left_homology, right_homology);
    }


    std::string repeats_str = deletion_size < 5 ? std::to_string(deletion_size) : "5+";

    if (homology > 0) {
        std::string homology_str = homology < 5 ? std::to_string(homology) : "5+";
        result = "DEL_MH_" + repeats_str + "_" + homology_str + ":" +
                 variant_call_record.getFilter();
    } else {
        std::string matches_str = matches < 5 ? std::to_string(matches) : "5+";
        result = "DEL_" + (deletion_size == 1 ? bp : "repeats") + "_" + repeats_str + "_" + matches_str + ":" +
                 variant_call_record.getFilter();
    }

    LOGGER.debug << variant_call_record.getChromosome() << ":"
                 << variant_call_record.getPosition() << std::endl;
    LOGGER.debug << left_context << "." << bp << "." << right_context << std::endl;
    LOGGER.debug << "result: " << result << std::endl;
    LOGGER.debug << std::endl;

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


    int count = 0;
    while (bcf_read(vcf_file, vcf_header, vcf_record) == 0) {

        int variant_type = bcf_get_variant_types(vcf_record);
        if (variant_type != VCF_SNP && variant_type != VCF_MNP && variant_type != VCF_INDEL) { continue; }

        variantCallRecord variant_call_record(vcf_header, vcf_record);
        if (!variant_call_record.isValid()) { continue; }

        std::string chromosome = variant_call_record.getChromosome();

        if (regions.find(chromosome) == regions.end()) {
            chromosome="other";
        }

        for (auto &item : regions[chromosome]) {
            if (!item.intersect(variant_call_record.getPosition())) { continue; }

            for (auto &sample : variant_call_record.getSamples()) {
                for (auto &alternate : sample.second.variants) {
                    if (!sample.second.is_variant) { continue; }

                    std::string mutation_key;

                    if (variant_type == VCF_SNP) {
                        mutation_key = get_key_from_snp(ref_file_index, sample.second.reference, alternate,
                                                        variant_call_record);
                    } else {
                        if (variant_type == VCF_MNP) {
                            mutation_key = get_key_from_mnp(sample.second.reference, alternate,
                                                            variant_call_record);
                        } else {
                            auto reference = sample.second.reference;
                            bool insertion = reference.size() < alternate.size();

                            if (insertion) {
                                mutation_key = get_key_from_insertion(ref_file_index, sample.second.reference,
                                                                      alternate,
                                                                      variant_call_record);
                            } else {
                                mutation_key = get_key_from_deletion(ref_file_index, sample.second.reference, alternate,
                                                                     variant_call_record);
                            }

                        }
                    }

                    if (item.mutations.find(mutation_key) == item.mutations.end()) { continue; }
                    item.mutations[mutation_key][sample.first] += 1;
                    item.total_bases++;
                    count++;
                }
            }

        }
    }

    LOGGER.info << "Total number of mutations " << count << std::endl;

    bcf_hdr_destroy(vcf_header);
    vcf_close(vcf_file);

    std::ofstream file(common_options.output_path + "/output.json");
    file << cmri::serialize(regions);
    file.close();


}
