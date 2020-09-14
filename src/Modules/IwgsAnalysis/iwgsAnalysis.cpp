
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

#include "iwgsAnalysis.h"
#include <faidx.h>

void cmri::mainIwgsAnalysis(const common_options_t &common_options,
                              const iwgs_analysis_options_t &iwgs_options) {


    cmri::open_file(common_options.input_file + ".fai", "expecting index file.").close();

    faidx_t *ref_file_index = fai_load(common_options.input_file.c_str());

    faidx_t *ref_file_index_pair;

    if(iwgs_options.pair_ended) {
        cmri::open_file(iwgs_options.input_file + ".fai", "expecting index file.").close();
        ref_file_index_pair = fai_load(iwgs_options.input_file.c_str());
    }




}