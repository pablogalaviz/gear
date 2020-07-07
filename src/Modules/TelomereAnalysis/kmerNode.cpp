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
#include <sstream>

cmri::kmerNode::kmerNode(std::shared_ptr<kmerNode>  *other){
    index = (*other)->index;
    observed = (*other)->observed;
    cost = (*other)->cost;
    children = (*other)->children;
    parent = (*other)->parent;
    state = (*other)->state;
}


std::string cmri::kmerNode::serialize(){

    std::stringstream result;

    result << "{ \"name\" : \"" << index.first << "\","
           << "\"level\" : " << index.second << ","
           << "\"observed\" : \"" << observed << "\","
           << "\"state\" : \"" << getState() << "\","
           << "\"cost\" : " << cost << ","
           << "\"children\" : [";
    for(int i =0; i < children.size(); i++){
        result << children[i]->serialize() << ",";
    }
    if( children.size() > 0)
        result.seekp(-1, std::ios_base::end);;

    result << "]}";

    return result.str();


}

std::string cmri::kmerNode::inverseSerialize(){

    std::stringstream result;

    result << "{ \"name\" : \"" << index.first << "\","
           << "\"level\" : " << index.second << ","
           << "\"observed\" : \"" << observed << "\","
           << "\"state\" : \"" << getState() << "\","
           << "\"cost\" : " << cost ;
    if(parent != nullptr ){
        if(parent->parent != nullptr ){
            result << ",\"parent\" : " << parent->inverseSerialize();
        }
    }
    result << "}";

    return result.str();


}


void cmri::kmerNode::getBPS(BasePairSequence &bp){

    bp.seq.insert(0,1,observed);
    bp.matching.replace(index.second-1,1,1,index.first);

    if (parent != nullptr)
    {
        if (parent->parent != nullptr)
        {
            parent->getBPS(bp);
        }
    }
}
