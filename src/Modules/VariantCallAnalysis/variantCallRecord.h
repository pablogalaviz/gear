//
// Created by pablo on 24/7/20.
//

#ifndef GEAR_VARIANTCALLRECORD_H
#define GEAR_VARIANTCALLRECORD_H

#include <string>
#include <vcf.h>
#include <vector>
namespace cmri {

    struct variant_t {
        std::string reference;
        std::vector<std::string> variants;
        bool is_variant;
    };


    class variantCallRecord {

        std::string chromosome;
        int position;
        std::string filter;
        bool valid;
        std::map<std::string ,variant_t> samples;

    public:

        variantCallRecord(bcf_hdr_t *vcf_header, bcf1_t *vcf_record);

        const std::string &getChromosome() const {
            return chromosome;
        }

        int getPosition() const {
            return position;
        }

        const std::string &getFilter() const {
            return filter;
        }

        bool isValid() const {
            return valid;
        }

        const std::map<std::string, variant_t> &getSamples() const {
            return samples;
        }

    };

}

#endif //GEAR_VARIANTCALLRECORD_H
