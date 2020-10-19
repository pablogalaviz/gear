#ifndef GEAR_MOTIFCOUNT_H
#define GEAR_MOTIFCOUNT_H

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


#include "motifRegion.h"
#include <regex>
#include <future>
#include "utils.h"
#include "sequenceReader.h"
#include "options.h"


namespace cmri {

    //given a DNA sequence string and a motif count the number of occurrences of the given motif in the string
    unsigned int searchMotif(std::string sequence, const std::string &motif);

    //given a DNA seqience string and a vector of phred quality values (per bp) creates a motif occurrence histogram of quality values.
    void searchMotifWFilter(std::string sequence, std::vector<int> quality, std::map<std::string,std::map<unsigned int,unsigned int>> &motif_quality);

    //given a DNA sequence string and a regular expression count the number of occurrences of the given regular expression in the string
    unsigned int searchRegex(std::string sequence, const std::string &regex);
    void searchRegex(std::string sequence, std::vector<int> quality,
                     std::map<std::string, std::map<unsigned int, unsigned int>> &regex_quality);

    //given a DNA sequence string and a regular expression count the number of CONSECUTIVE occurrences of the given regular expression in the string
    int searchRegexConsecutive(std::string sequence, const std::string &regex);

    mapVectorMotifRegion
    processWorker(mapVectorMotifRegion motif_map, const std::vector<read_item_t>& sequences, bool validate);

    void process(const common_options_t &common_options,const motif_count_options_t &motif_count_options, mapVectorMotifRegion &motif_map);

    void processMultiThreading(const common_options_t &common_options,const motif_count_options_t &motif_count_options, mapVectorMotifRegion &motif_map);

    void mainMotifCount(const common_options_t &common_options, const motif_count_options_t &motif_count_options);


}

#endif //GEAR_MOTIFCOUNT_H
