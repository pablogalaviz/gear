
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <faidx.h>
#include <bedWriter.h>
#include "genomeAnalysis.h"
#include "telomereRegion.h"

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


void cmri::mainGenomeAnalysis(const common_options_t &common_options,
                              const genome_analysis_options_t &telomere_options) {

    boost::property_tree::ptree region_data;
    boost::property_tree::read_json(telomere_options.regions, region_data);

    std::map<std::string, std::vector<telomereRegion>> regions;
    cmri::deserialize(region_data, regions);

    if (regions.empty()) {
        cmri::LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
        cmri::LOGGER.error << "Invalid argument. No valid region list found!" << std::endl;
        exit(EINVAL);
    }

    std::string index_file = common_options.input_file + ".fai";
    cmri::open_file(index_file, "expecting index file.").close();

    faidx_t *ref_file_index = fai_load(common_options.input_file.c_str());

    bedWriter bed_writer(common_options.output_path + "/output.bed");


    for (auto &chromosome_item : regions) {
        std::string chromosome = chromosome_item.first;

        cmri::LOGGER.info << "Procesing contig: " << chromosome << std::endl;

        for (auto &region_item : chromosome_item.second) {
            if (region_item.start == region_item.end) { continue; }

            int len;
            std::string region = chromosome + ":"
                                 + std::to_string(region_item.start) + "-"
                                 + std::to_string(region_item.end);
            std::string sequence = fai_fetch(ref_file_index, region.c_str(), &len);
            region_item.total_bases = sequence.size();

            if(telomere_options.validate) {
                for (auto &c: sequence) { c = toupper(c); }
            }

            int variant_index=0;
            for (auto &variant : region_item.variants) {
                variant_index++;
                std::string::size_type start = 0;
                start = sequence.find(variant.first, start);
                if(start==std::string::npos){continue;}

                auto variant_size =  variant.first.size();
                std::string::size_type current_position = start+variant_size;
                std::string::size_type previous_position = start;
                int count =1;
                int acc=0;
                while ((current_position = sequence.find(variant.first, current_position)) != std::string::npos) {
                    count++;
                    if(current_position != previous_position + variant_size){

                        std::stringstream context;
                        std::string left_context= start>=1 ? sequence.substr(start-1,1) : "";
                        std::string right_context= previous_position + variant_size+1 < sequence.size() ? sequence.substr(previous_position + variant_size+1,1) : "";
                        context << left_context << "_" << right_context;
                        int consecutive = sequence.substr(start,previous_position + variant_size-start).size()/variant_size;
                        acc+=consecutive;
                        variant.second.histogram[consecutive]+=1;
                        variant.second.context[context.str()]+=1;

                        bed_writer.append(chromosome,region_item.start+start-1,region_item.start+previous_position + variant_size-1, variant.first, consecutive, true, 0, 0,255*(variant_index % 2), 255 * ((variant_index + 1) % 2), 51 * (variant_index % 5));

                        //LOGGER.debug <<  chromosome <<":" << region_item.start+start  <<"-"  << region_item.start+previous_position + variant_size << std::endl;
                        //LOGGER.debug <<  "consecutive:" << consecutive  <<" context: "  << context.str() << std::endl;
                        //LOGGER.debug <<  "seq:" << sequence.substr(start,previous_position + variant_size-start) << std::endl;
                        start=current_position;

                    }

                    previous_position = current_position;
                    current_position+=variant_size;
                }

                if(start!=previous_position+variant_size)
                {
                    std::stringstream context;
                    std::string left_context= start>=1 ? sequence.substr(start-1,1) : "";
                    std::string right_context= previous_position + variant_size+1 < sequence.size() ? sequence.substr(previous_position + variant_size+1,1) : "";
                    context << left_context << "_" << right_context;
                    int consecutive = sequence.substr(start,previous_position + variant_size-start).size()/variant_size;
                    acc+=consecutive;
                    variant.second.histogram[consecutive]+=1;
                    variant.second.context[context.str()]+=1;
                }

                assert(count == acc);


                variant.second.count=count;

            }



        }
    }


    std::ofstream file(common_options.output_path + "/output.json");
    file << cmri::serialize(regions);
    file.close();


}
