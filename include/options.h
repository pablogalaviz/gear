//
// Created by pablo on 26/6/20.
//

#ifndef GEAR_OPTIONS_H
#define GEAR_OPTIONS_H

#include <string>
#include <thread>
#include "utils.h"

namespace cmri {

    enum class task_t : int {
        MotifCount, VariantCallAnalysis, GenomeAnalysis, IwgsAnalysis, TelomereMutations, RandomSelector, QVSelector
    };

    static std::map<std::string,task_t> str2task {
            {"MotifCount",task_t::MotifCount}
            ,{"VariantCallAnalysis",task_t::VariantCallAnalysis}
            , {"GenomeAnalysis",task_t::GenomeAnalysis}
            , {"IwgsAnalysis",task_t::IwgsAnalysis}
            , {"TelomereMutations",task_t::TelomereMutations}
            , {"RandomSelector",task_t::RandomSelector}
            , {"QVSelector",task_t::QVSelector}
    };

    struct common_options_t {
        bool backup;
        std::string output_path;
        std::string input_file;
        int progress = 0;
        int chunk_size = 1;
        int threads = 1;

        void validate() {
            int max_threads = static_cast<int>(std::thread::hardware_concurrency());
            threads = cmri::clip(threads, 1, max_threads);
            progress = std::max(progress, 0);
            //test if files exists.
            cmri::open_file(input_file, "expecting input file").close();
        }
    };


    struct motif_count_options_t {
        std::string motif_file;
        int quality_value = 0;
        int quality_map = 0;
        bool validate_sequence;

        void validate() {
            cmri::open_file(motif_file, "expecting morif file");
            quality_value = cmri::clip(quality_value, 0, 92);
            quality_map = cmri::clip(quality_map, 0, 254);
        }

    };


    struct variant_call_analysis_options_t {
        std::string regions;
        std::string reference;

        void validate() const {
            cmri::open_file(regions, "expecting region description file.").close();
            cmri::open_file(reference, "expecting reference file.").close();
            cmri::open_file(reference + ".fai", "expecting reference index file.").close();

        }
    };

    struct genome_analysis_options_t {
        std::string regions;
        std::string target_file;
        bool validate_sequence;

        void validate() const {
            cmri::open_file(regions, "expecting region description file.").close();
            cmri::open_file(target_file, "expecting fasta telomere reference file.").close();
        }
    };


    struct iwgs_analysis_options_t {
        std::string input_file;
        std::string variants_file;
        int count_filter_threshold;

        void validate() {
            //test if files exists.
            cmri::open_file(input_file, "expecting input file").close();
            cmri::open_file(variants_file, "expecting variants file").close();
            count_filter_threshold = std::max(1, count_filter_threshold);
        }


    };


    struct telomere_mutations_options_t {

        std::string target_file;
        std::string query_file;
        size_t trimming_window_mean;
        int trimming_threshold;


        void validate() {
            //test if files exists.
            cmri::open_file(target_file, "expecting target file").close();
            cmri::open_file(query_file, "expecting query file").close();
            trimming_window_mean=clip(trimming_window_mean,2ul,100ul);
            trimming_threshold=clip(trimming_threshold,0,100);
        }


    };


}


#endif //GEAR_OPTIONS_H
