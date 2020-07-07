//
// Created by pablo on 26/6/20.
//

#ifndef GEAR_OPTIONS_H
#define GEAR_OPTIONS_H

#include <string>
#include <thread>
#include "utils.h"

namespace cmri {

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
            cmri::open_file(input_file,"expecting input file").close();
        }
    };


    struct motif_count_options_t {
        std::string motif_file;
        int quality_value = 0;
        int quality_map = 0;
        bool validate_sequence;

        void validate() {
            cmri::open_file(motif_file,"expecting morif file");
            quality_value = cmri::clip(quality_value, 0, 92);
            quality_map = cmri::clip(quality_map, 0, 254);
        }

    };


    struct variant_call_analysis_options_t {
        std::string regions;

        void validate() const {
            cmri::open_file(regions,"expecting region description file.").close();
        }
    };

}


#endif //GEAR_OPTIONS_H