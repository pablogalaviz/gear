#ifndef GEAR_SEQUENCENETWORK_H
#define GEAR_SEQUENCENETWORK_H

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

#include "kmerNode.h"

namespace cmri {

    typedef std::tuple<char, int, double> nodeValue;
    typedef std::vector<nodeValue> transitionVector;
    typedef std::map<nodeIndex, transitionVector> transitionMatrix;


    class sequenceNetwork {

        transitionMatrix tMatrix;
        transitionMatrix reverse_tMatrix;


        void makeNetwork(std::string segment, std::shared_ptr<kmerNode> &node, transitionMatrix &tMatrix,
                         std::shared_ptr<kmerNode> &bestNode);

        std::vector<nodeIndex> parseHeader(std::string line);

        std::pair<nodeIndex, std::vector<double> > parseRow(std::string line);

        transitionMatrix readTransitionMatrix(std::string path);

        std::string serializeTransitionMatrix(transitionMatrix m);

        std::string serializeTransitionVector(transitionVector v);


    public :

        sequenceNetwork(std::string filename);

        ~sequenceNetwork() {};

        void nt_analysis(std::string segment, bool reverse);

    };

    static std::map<std::string, telomere_t> map2telomere= {
             {"AA",TTAGGG}
            ,{"AC",NT}
            ,{"AG",TTGGGG}
            ,{"AT",NT}
            ,{"AN",NT}
            ,{"CA",TCAGGG}
            ,{"CC",TCAGGG}
            ,{"CG",TCAGGG}
            ,{"CT",TCAGGG}
            ,{"CN",TCAGGG}
            ,{"GA",TGAGGG}
            ,{"GC",TGAGGG}
            ,{"GG",TGAGGG}
            ,{"GT",TGAGGG}
            ,{"GN",TGAGGG}
            ,{"TA",TTAGGG}
            ,{"TC",TTAGGG}
            ,{"TG",TTGGGG}
            ,{"TT",TTAGGG}
            ,{"TN",TTAGGG}
            ,{"NA",TTAGGG}
            ,{"NC",NT}
            ,{"NG",TTGGGG}
            ,{"NT",NT}
            ,{"NN",NT}
    };


    static telomere_t bestMatch(std::string target)
    {
        return map2telomere[target.substr(1,2)];
    }



}

#endif //GEAR_SEQUENCENETWORK_H
