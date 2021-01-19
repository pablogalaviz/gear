#ifndef GEAR_GENOMEANALYSIS_H
#define GEAR_GENOMEANALYSIS_H

//
// Author(s) Pablo Galaviz (2020)
// e-mail  <pgalaviz@cmri.org.au>
//



//  This file is part of GEAR
//
//  GEAR is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  any later version.
//
//  GEAR is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GEAR.  If not, see <http://www.gnu.org/licenses/>.
//

#include <options.h>
#include "minimap.h"

namespace cmri {


    enum class region_state_t : int { COMPLETE=0, EXTEND=1, REDUCE=2 };

    struct bed_region_t {
        std::string chrom;
        unsigned long chromStart=0;
        unsigned long chromEnd=0;
        std::string name="";
        unsigned int score=0;
        char strand='+';
        unsigned long thickStart=0;
        unsigned long thickEnd=0;

        region_state_t state;
        bool is_telomere;

    };

    class genomeAnalysis {

        mm_idxopt_t iopt;
        mm_mapopt_t mopt;
        int n_threads = 3;
        mm_idx_reader_t *index_reader;
        mm_idx_t *mi;
        std::string target_file;

        cmri::bed_region_t find_telomere(const std::string &sequence, const cmri::bed_region_t &previous_region, const std::string::size_type &start, const std::string::size_type &end);

    public:
        genomeAnalysis(const genome_analysis_options_t &options);



        std::vector<cmri::bed_region_t> process(const std::string sequence,const std::string chrom);

    };

    void mainGenomeAnalysis(const common_options_t &common_options, const genome_analysis_options_t &telomere_options);

}


#endif //GEAR_GENOMEANALYSIS_H
