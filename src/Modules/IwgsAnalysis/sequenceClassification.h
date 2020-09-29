//
// Created by Pablo Galaviz on 22/9/20.
//

#ifndef GEAR_SEQUENCECLASSIFICATION_H
#define GEAR_SEQUENCECLASSIFICATION_H


#include <vector>
#include <string>
#include <cmath>

struct sequenceMotif{
    int kind;
    double qv;
};



class sequenceClassification {
public:
    explicit sequenceClassification(std::vector<std::string> variants);

    void append_motif(std::string quality);

    void append_segment(std::string segment,std::string quality);


private:

    std::vector<std::string> variants;

    std::vector<sequenceMotif> sequence_motifs;

    int confidence;

};


inline double calculate_mean_qv(const std::string &quality){
    double mean_qv=0;
    for(auto &c : quality) {
        mean_qv +=  1- std::pow(10.0,-(static_cast<int>(c) - 33)/10);
    }
    mean_qv/= quality.size();
    return mean_qv;
}


#endif //GEAR_SEQUENCECLASSIFICATION_H
