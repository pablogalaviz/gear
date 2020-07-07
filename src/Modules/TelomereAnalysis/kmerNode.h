
#ifndef GEAR_KMERNODE_H
#define GEAR_KMERNODE_H

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


#include <memory>
#include <vector>
#include <map>
#include <string>

namespace cmri {

    enum base_state_t : int {
        MATCH = 0, INSERTION = 1, MISMATCH = 2, DELETION = 3, UNKNOWN = 4
    };

    enum telomere_t : int {
        TTAGGG = 1, TCAGGG = 2, TGAGGG = 3, TTGGGG = 4, NT = 0
    };

    typedef std::pair<char, int> nodeIndex;

    static std::map <base_state_t, std::string> state_map = {
            {MATCH,     "MATCH"},
            {INSERTION, "INSERTION"},
            {MISMATCH,  "MISMATCH"},
            {DELETION,  "DELETION"},
            {UNKNOWN,   "UNKNOWN"}
    };


    struct BasePairSequence {

        std::string seq = "";
        std::string matching = "NNNNNN";
        double score;
        double id;

        std::vector<int> base_state();

        telomere_t likely_match;


    };

    class kmerNode : public std::enable_shared_from_this<kmerNode> {

        nodeIndex index;
        char observed;
        double cost = 0;

        std::vector <std::shared_ptr<kmerNode>> children;
        std::shared_ptr <kmerNode> parent;
        base_state_t state = UNKNOWN;

    public:

        kmerNode(char n, int l) : index(n, l), observed(n) {};

        kmerNode(std::shared_ptr <kmerNode> *other);

        ~kmerNode() {};

        void addChildren(std::shared_ptr <kmerNode> &c) {

            children.push_back(c);
            c->setParent(shared_from_this());

        }

        void setParent(const std::shared_ptr <kmerNode> &p) {
            parent = p;
        }

        inline void setCost(double c) { cost = c; }

        inline void setState(base_state_t s) { state = s; }

        inline void setObserved(char t) { observed = t; }

        inline int getLevel() { return index.second; }

        inline double getCost() { return cost; }

        inline char getName() { return index.first; }

        inline std::string getState(){return state_map[state];}

        inline nodeIndex getIndex() { return index; }

        std::string serialize();

        std::string inverseSerialize();

        void getBPS(BasePairSequence &bp);

    };
}

#endif //GEAR_KMERNODE_H
