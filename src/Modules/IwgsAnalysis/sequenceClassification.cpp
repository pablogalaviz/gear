//
// Created by Pablo Galaviz on 22/9/20.
//

#include "sequenceClassification.h"

#include <utility>
#include <logger.h>

sequenceClassification::sequenceClassification(std::vector<std::string> variants) : variants(std::move(variants)) {

    confidence=0;
}

void sequenceClassification::append_motif(std::string quality) {

    double mean_qv = calculate_mean_qv(quality);

    sequenceMotif new_main_motif;
    new_main_motif.kind = 0;
    new_main_motif.qv=mean_qv;
    sequence_motifs.push_back(new_main_motif);

}

void sequenceClassification::append_segment(std::string segment, std::string quality) {


}
