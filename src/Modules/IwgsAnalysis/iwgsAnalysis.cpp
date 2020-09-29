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
#include "../MotifCount/motifCount.h"
#include "sequenceClassification.h"
#include <zlib.h>


void cmri::mainIwgsAnalysis(const common_options_t &common_options,
                              const iwgs_analysis_options_t &iwgs_options) {

        cmri::open_file(common_options.input_file, "expecting index file.").close();
        gzFile file = gzopen(common_options.input_file.c_str(), "r");
        kseq_t *kseq = kseq_init(file);
        int l;
        int count =0;


        std::vector<std::string> variants;
        variants.push_back("TCAGGG");
        variants.push_back("TGAGGG");

        iwgsAnalysis iwgs("TTAGGG", variants, 5);

        while ((l = kseq_read(kseq)) >= 0) {

            if(iwgs.count_filter(kseq->name.s, kseq->seq.s, kseq->qual.s)){
                break;

            }
            count++;
            if(count > 10000){
            }

        }

        LOGGER.debug << "count: " << count << std::endl;


}


bool cmri::iwgsAnalysis::count_filter(std::string name, std::string sequence, std::string quality) {


    if(searchMotif(sequence, motif_forward) > count_filter_threshold){

        LOGGER.debug << "forward process " << name << std::endl;
        LOGGER.debug << "seq " << sequence << std::endl;

        return true;
    }



    if(searchMotif(sequence, motif_reverse)  > count_filter_threshold){

        sequenceClassification sequence_classification(variants_reverse);

        LOGGER.debug << "reverse process " << name << std::endl;
        LOGGER.debug << "seq " << sequence << std::endl;
        LOGGER.debug << "quality " << quality << std::endl;

        std::string::size_type previous = 0;
        std::string::size_type start = 0;
        while ((start = sequence.find(motif_reverse, start)) != std::string::npos) {
            // move start to find the rest of the string.

            std::string::size_type segment_size = start-previous;
            if(segment_size > 0){
                std::string segment = sequence.substr(previous,segment_size);
                std::string segment_qv = quality.substr(previous,segment_size);
                LOGGER.debug << "segment " <<  segment << " " <<previous << ":" << start << std::endl;
                sequence_classification.append_segment(segment,segment_qv);
            }

            sequence_classification.append_motif(quality.substr(start,motif_reverse.size()));

            start += motif_reverse.size();
            previous = start;
        }

        return true;
    }

    return false;

}

cmri::iwgsAnalysis::iwgsAnalysis(const std::string &motif, const std::vector<std::string> &variants, int countFilterThreshold) :
        motif_forward(motif)
        ,variants_forward(variants)
        ,count_filter_threshold(countFilterThreshold) {

    motif_reverse = reverse_complement(motif);

    for(auto &variants : variants_forward){
        variants_reverse.push_back(reverse_complement(variants));
    }


}
