
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

#include <algorithm>
#include <logger.h>
#include "sequenceNetwork.h"


cmri::sequenceNetwork::sequenceNetwork(std::string filename) {

    tMatrix = readTransitionMatrix(filename);
    reverse_tMatrix.clear();

    std::map<char,char> base_pair_map = {{'A','T'},{'C','G'},{'G','C'},{'T','A'},{'R','R'}};

    for(auto &item : tMatrix){
        nodeIndex node_index(base_pair_map[item.first.first],item.first.second);
        transitionVector transition_vector;

        for(int i=0; i < item.second.size(); i++){
            nodeValue nv = item.second[i];
            transition_vector.push_back(nodeValue(base_pair_map[std::get<0>(nv)],std::get<1>(nv),std::get<2>(nv)));
        }
        reverse_tMatrix[node_index] = transition_vector;

    }

    //LOGGER.debug <<  serializeTransitionMatrix(reverse_tMatrix) << std::endl;

}



void cmri::sequenceNetwork::makeNetwork(std::string segment, std::shared_ptr<kmerNode> &node, transitionMatrix &tMatrix, std::shared_ptr<kmerNode> &bestNode)
{

    char current = segment[0];
    std::shared_ptr<kmerNode> nn = std::make_shared<kmerNode>(current, node->getLevel() + 1);

    transitionVector node2nn = tMatrix[node->getIndex()];

    if (nn->getLevel() <= 6 && node->getCost() >= -2)
    {
        bool valid_transition = false;
        for (int i = 0; i < node2nn.size(); i++)
        {
            nodeValue v = node2nn[i];
            if (std::get<2>(v) > 0 && std::get<0>(v) == nn->getName() && std::get<1>(v) == nn->getLevel())
            {
                valid_transition = true;
                double cost = node->getCost() + std::get<2>(v);
                nn->setCost(cost);
                nn->setState(MATCH);
                node->addChildren(nn);

                if (segment.size() > 1)
                    makeNetwork(segment.substr(1, segment.size()), nn, tMatrix, bestNode);

                if (nn->getCost() >= bestNode->getCost() && nn->getLevel() >= bestNode->getLevel() && nn->getLevel()<= 6 )
                {
                    bestNode = nn;
                }
            }
        }
        if (!valid_transition)
        {
            transitionVector shift2node = tMatrix[nn->getIndex()];
            for (int i = 0; i < shift2node.size(); i++)
            {
                nodeValue v = shift2node[i];
                std::shared_ptr<kmerNode> mm = std::make_shared<kmerNode>(std::get<0>(v), std::get<1>(v));
                double cost = node->getCost() + std::get<2>(v);
                mm->setCost(cost);
                mm->setObserved(nn->getName());
                if(mm->getLevel() < nn->getLevel()){
                    mm->setState(INSERTION);
                }
                else{
                    if(mm->getLevel() == nn->getLevel()){
                        mm->setState(MISMATCH);
                    }
                    else{
                        mm->setState(DELETION);
                    }
                }

                if (segment.size() > 1)
                    makeNetwork(segment.substr(1, segment.size()), mm, tMatrix, bestNode);

                node->addChildren(mm);
                if (mm->getCost() >= bestNode->getCost() && mm->getLevel() > bestNode->getLevel())
                {
                    bestNode = mm;
                }
            }
        }
    }
    else{
        LOGGER.warning << "level: " << nn->getLevel() << " cost: " <<  node->getCost() << std::endl;
    }
}

void cmri::sequenceNetwork::nt_analysis(std::string segment, bool reverse)
{
    if (reverse)
    {
        std::reverse(segment.begin(), segment.end());
    }

    std::string reconstructed="";

    while(segment.size()>3){

        std::shared_ptr<kmerNode> root = std::make_shared<kmerNode>('R', 0);
        std::shared_ptr<kmerNode> bestNode = std::make_shared<kmerNode>('R', 0);
        bestNode->setCost(-1);
        if (reverse)
        {
            makeNetwork(segment, root, reverse_tMatrix, bestNode);
        }
        else
        {
            makeNetwork(segment, root, tMatrix, bestNode);
        }

        LOGGER.warning << root->serialize() << std::endl;

        LOGGER.warning << bestNode->inverseSerialize() << std::endl;

        BasePairSequence bp;
        bp.score = bestNode->getCost();


        bestNode->getBPS(bp);

        LOGGER.info << "Original:" << segment << std::endl;
        LOGGER.info << "sequence:" << bp.seq << std::endl;
        LOGGER.info << "match:" << bp.matching << std::endl;
        LOGGER.info << "type:" << bestMatch(bp.matching) << std::endl;
        LOGGER.info << "score:" << bp.score << std::endl;
        LOGGER.info << "id:" << bp.score+bestMatch(bp.matching) << std::endl;

        reconstructed += bp.matching;

        if(bp.seq.size() <= segment.size()){
            segment = segment.substr(bp.seq.size(), segment.size());
        }

    };

    reconstructed+=segment;

    LOGGER.info << "result: " << reconstructed << std::endl;

}

std::vector< cmri::nodeIndex > cmri::sequenceNetwork::parseHeader(std::string header){

    std::vector< nodeIndex > result;
    std::stringstream ss(header);
    std::string line;
    while (std::getline(ss, line,','))
    {
        if(line.size() > 1){
            if(isdigit(line[1]) && (line[0] == 'A' || line[0] == 'C' || line[0] == 'G' || line[0] == 'T' || line[0] == 'R'))
            {
                int level = line[1]-'0';
                result.push_back( nodeIndex(line[0],level) );
            }
        }
    }

    return result;

}

std::pair<cmri::nodeIndex, std::vector<double>> cmri::sequenceNetwork::parseRow(std::string line)
{

    std::pair<nodeIndex, std::vector<double>> result;

    std::stringstream ss(line);
    std::string field;

    std::getline(ss, field, ',');

    if (isdigit(field[1]) && (field[0] == 'A' || field[0] == 'C' || field[0] == 'G' || field[0] == 'T' || field[0]=='R'))
    {
        int level = line[1] - '0';
        result.first = nodeIndex(field[0], level);
    }
    else
    {
        LOGGER.error << "Transition matrix invalid format. Expecting first column Node-Level pair. Got: " << field << std::endl;
        exit(-1);
    }

    std::vector<double> values;
    while (std::getline(ss, field, ','))
    {
        try
        {
            values.push_back(std::stod(field));
        }
        catch (...)
        {
            LOGGER.warning << "invalid value on transition matrix, expecting floating-point number. Got: " << field << std::endl;
            values.push_back(0);
        }
    }

    result.second = values;

    return result;
}

cmri::transitionMatrix cmri::sequenceNetwork::readTransitionMatrix(std::string file_name)
{

    transitionMatrix result;

    std::ifstream ifs(file_name.c_str());

    if (ifs.is_open())
    {
        std::string line;

        std::getline(ifs, line);

        std::vector<nodeIndex> header = parseHeader(line);

        while (std::getline(ifs, line))
        {
            std::pair< nodeIndex, std::vector<double> > row = parseRow(line);
            if(row.second.size() == header.size())
            {
                transitionVector value;
                for(int i=0; i < header.size(); i++){
                    if(row.second[i] != 0){
                        std::tuple<char,int,double>  v = std::make_tuple(header[i].first,header[i].second,row.second[i]);
                        value.push_back(v);
                    }
                }
                result[row.first] = value;
            }
            else{
                LOGGER.warning << "Invalid transition matrix data, expecting "<< header.size() << " columns. Got: " << row.second.size() << std::endl;
            }

        }

        ifs.close();
    }

    else
    {
        LOGGER.error << "Can not open file: " << file_name << std::endl;
        exit(-1);
    }

    return result;
}

std::string cmri::sequenceNetwork::serializeTransitionMatrix(transitionMatrix matrix)
{

    std::stringstream result;

    result << "{";

    for (auto &v : matrix)
    {
        result << "\"" << v.first.first << v.first.second << "\" : " << serializeTransitionVector(v.second) << ",";
    }

    if(matrix.size()>0)
        result.seekp(-1, std::ios_base::end);

    result << "}";

    return result.str();
}

std::string cmri::sequenceNetwork::serializeTransitionVector(transitionVector v)
{

    std::stringstream result;

    result << "[";

    for (int i = 0; i < v.size(); i++)
    {
        std::tuple<char, int, double> k = v[i];
        result << "{\"" << std::get<0>(k) << std::get<1>(k) << "\":" << std::get<2>(k) << "},";
    }
    if (v.size() > 0)
        result.seekp(-1, std::ios_base::end);

    result << "]";


    return result.str();
}
