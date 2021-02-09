
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

#include <math.h>

void cmri::mainGenomeAnalysis(const common_options_t &common_options,
                              const genome_analysis_options_t &telomere_options) {

    genomeAnalysis ga(telomere_options);

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

            if (telomere_options.validate_sequence) {
                for (auto &c: sequence) { c = toupper(c); }
            }

            int variant_index = 0;
            for (auto &variant : region_item.variants) {
                variant_index++;
                std::string::size_type start = 0;
                start = sequence.find(variant.first, start);
                if (start == std::string::npos) { continue; }

                auto variant_size = variant.first.size();
                std::string::size_type current_position = start + variant_size;
                std::string::size_type previous_position = start;
                int count = 1;
                int acc = 0;
                while ((current_position = sequence.find(variant.first, current_position)) != std::string::npos) {
                    count++;
                    if (current_position != previous_position + variant_size) {

                        std::stringstream context;
                        std::string left_context = start >= 1 ? sequence.substr(start - 1, 1) : "";
                        std::string right_context =
                                previous_position + variant_size + 1 < sequence.size() ? sequence.substr(
                                        previous_position + variant_size + 1, 1) : "";
                        context << left_context << "_" << right_context;
                        int consecutive =
                                sequence.substr(start, previous_position + variant_size - start).size() / variant_size;
                        acc += consecutive;
                        variant.second.histogram[consecutive] += 1;
                        variant.second.context[context.str()] += 1;

                        bed_writer.append(chromosome, region_item.start + start - 1,
                                          region_item.start + previous_position + variant_size - 1, variant.first,
                                          consecutive, true, 0, 0, 255 * (variant_index % 2),
                                          255 * ((variant_index + 1) % 2), 51 * (variant_index % 5));

                        //LOGGER.debug <<  chromosome <<":" << region_item.start+start  <<"-"  << region_item.start+previous_position + variant_size << std::endl;
                        //LOGGER.debug <<  "consecutive:" << consecutive  <<" context: "  << context.str() << std::endl;
                        //LOGGER.debug <<  "seq:" << sequence.substr(start,previous_position + variant_size-start) << std::endl;
                        start = current_position;

                    }

                    previous_position = current_position;
                    current_position += variant_size;
                }

                if (start != previous_position + variant_size) {
                    std::stringstream context;
                    std::string left_context = start >= 1 ? sequence.substr(start - 1, 1) : "";
                    std::string right_context =
                            previous_position + variant_size + 1 < sequence.size() ? sequence.substr(
                                    previous_position + variant_size + 1, 1) : "";
                    context << left_context << "_" << right_context;
                    int consecutive =
                            sequence.substr(start, previous_position + variant_size - start).size() / variant_size;
                    acc += consecutive;
                    variant.second.histogram[consecutive] += 1;
                    variant.second.context[context.str()] += 1;
                }

                assert(count == acc);


                variant.second.count = count;

            }


            auto bed_regions = ga.process(sequence, chromosome);

            std::ofstream file(common_options.output_path + "/output_t.bed", std::ofstream::app);
            for (auto &item : bed_regions) {
                file << item.chrom << "\t"
                     << item.chromStart << "\t"
                     << item.chromEnd << "\t"
                     << item.name << "\t"
                     << item.score << "\t"
                     << item.strand << "\t"
                     << item.thickStart << "\t"
                     << item.thickEnd // << "\t"
                     << std::endl;
            }
            file.close();


        }
    }


    std::ofstream file(common_options.output_path + "/output.json");
    file << cmri::serialize(regions);
    file.close();


}

cmri::genomeAnalysis::genomeAnalysis(const genome_analysis_options_t &options) {


    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR; // perform alignment
    mopt.flag |= MM_F_OUT_CG; // perform alignment

    // open index reader
    target_file = options.target_file;


}


cmri::bed_region_t
cmri::genomeAnalysis::find_telomere(const std::string &sequence, const cmri::bed_region_t &previous_region,
                                    const std::string::size_type &start, const std::string::size_type &end) {
    cmri::bed_region_t result;
    result.chrom = previous_region.chrom;
    result.score = 0;

    // open index reader
    index_reader = mm_idx_reader_open(target_file.c_str(), &iopt, 0);

    while ((mi = mm_idx_reader_read(index_reader, n_threads)) != 0) { // traverse each part of the index
        mm_mapopt_update(&mopt,
                         mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread

        mm_reg1_t *reg;
        int j, i, n_reg;

        reg = mm_map(mi, sequence.size(), sequence.c_str(), &n_reg, tbuf, &mopt, 0); // get all hits for the query

        for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
            mm_reg1_t *r = &reg[j];
            assert(r->p); // with MM_F_CIGAR, this should not be NULL

            LOGGER.debug << "len:"
                         << sequence.size() << "\tqs:"
                         << start << "\tqe:"
                         << start + sequence.size() << "\t"
                         << "+-"[r->rev] << std::endl;

            LOGGER.debug << mi->seq[r->rid].name << "\tlen:"
                         << mi->seq[r->rid].len << "\trs:"
                         << r->qs << "\tre:"
                         << r->qe << "\tmlen:"
                         << r->mlen << "\tblen:"
                         << r->blen << "\tmapq:"
                         << r->mapq << "\t cg ";


            void *km = nullptr;
            char *cs_str = NULL;
            int max_len = 0;
            int n_cs = mm_gen_cs(km, &cs_str, &max_len, mi, r, sequence.c_str(), true);

            LOGGER.debug << cs_str << std::endl;
            LOGGER.debug << sequence << std::endl;

            unsigned int score = static_cast<unsigned int>(std::round(
                    ((static_cast<double>(r->blen + r->mlen) / sequence.size() + r->mapq / 60.0) / 3.0) * 1000));

            if (result.score < score) {
                result.name = "Telomere";
                result.score = score;
                result.chromStart = start;
                result.chromEnd = start + sequence.size();
                result.thickStart = start + r->qs;
                result.thickEnd = start + r->qe;
                result.strand = "+-"[r->rev];
            }

        }

        if (result.score > 600) {
            if (previous_region.state == region_state_t::EXTEND) {
                if (previous_region.score < result.score) { result.state = region_state_t::EXTEND; }
                else {
                    if (previous_region.score > result.score) { result.state = region_state_t::REDUCE; }
                    else {
                        result.state = region_state_t::COMPLETE;
                    }
                }
            } else {
                if (previous_region.state == region_state_t::REDUCE) {
                    if (previous_region.score < result.score) { result.state = region_state_t::REDUCE; }
                    else {
                        if (previous_region.score > result.score) { result.state = region_state_t::EXTEND; }
                        else {
                            result.state = region_state_t::COMPLETE;
                        }
                    }
                }

            }

        } else {
            result.state = region_state_t::COMPLETE;
            result.is_telomere = false;
        }


    }


    return result;
}

std::vector<cmri::bed_region_t> cmri::genomeAnalysis::process(const std::string sequence, const std::string chrom) {

    std::vector<bed_region_t> result;

    std::string::size_type start = 0;
    std::string::size_type end = std::min(100UL, sequence.size() - 1);
    std::string::size_type size = end - start;

    bed_region_t previous_region;
    previous_region.state = region_state_t::EXTEND;
    previous_region.is_telomere = true;
    previous_region.chrom = chrom;

    while (end < sequence.size()) {
        std::string subsequence = sequence.substr(start, size);

        bed_region_t region = find_telomere(subsequence, previous_region, start, end);

        auto delta_size = region.chromEnd - region.chromStart;

        if (region.state == region_state_t::EXTEND) {

            size = std::min<std::string::size_type>(sequence.size() - start, 2 * delta_size);
            if (size == 0) {
                region.state = region_state_t::COMPLETE;
            } else {
                end = region.chromStart + size;
                previous_region = region;
                continue;
            }
        }

        if (region.state == region_state_t::REDUCE) {
            auto s = 2.0*delta_size / 3.0;
            size = std::min(sequence.size() - start, static_cast<std::string::size_type>(std::round(s)));
            if (size == 0) {
                region.state = region_state_t::COMPLETE;
            } else {
                end = region.chromStart + size;
                previous_region = region;
                continue;
            }
        }

        if (region.state == region_state_t::COMPLETE) {
            if (region.is_telomere) {
                result.push_back(region);
            }
            start = region.chromEnd;
            size = std::min(100UL, sequence.size() - start);
            end = start + size;
            previous_region = region;
        } else {
            LOGGER.error << "Unknown state found in:" << __FILE__ << ":" << __LINE__ << std::endl;
            exit(-1);
        }

    }

    return result;

}




/*
 *
std::vector<cmri::bed_region_t> cmri::genomeAnalysis::process(const std::string sequence,const std::string chrom) {

    std::vector<bed_region_t> result;

    // open index reader
    index_reader=mm_idx_reader_open(target_file.c_str(), &iopt, 0);

    while ((mi = mm_idx_reader_read(index_reader, n_threads)) != 0) { // traverse each part of the index
        mm_mapopt_update(&mopt,
                         mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread

        unsigned long start =0;
        unsigned long size = std::min(100UL,sequence.size()-start);

        while(start+size <= sequence.size())
        { // each kseq_read() call reads one query sequence
            mm_reg1_t *reg;
            int j, i, n_reg;
            std::string subsequence = sequence.substr(start, size);
            reg = mm_map(mi, size, subsequence.c_str(), &n_reg, tbuf, &mopt, 0); // get all hits for the query

            bed_region_t region;
            region.chrom=chrom;
            region.score=0;
            for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                mm_reg1_t *r = &reg[j];
                assert(r->p); // with MM_F_CIGAR, this should not be NULL

                LOGGER.debug << "len:"
                             << size << "\tqs:"
                             << start << "\tqe:"
                             << start+size << "\t"
                             << "+-"[r->rev] << std::endl;

                LOGGER.debug << mi->seq[r->rid].name << "\tlen:"
                             << mi->seq[r->rid].len << "\trs:"
                             << r->qs << "\tre:"
                             << r->qe << "\tmlen:"
                             << r->mlen << "\tblen:"
                             << r->blen << "\tmapq:"
                             << r->mapq << "\t cg ";


                void *km = nullptr;
                char *cs_str = NULL;
                int max_len = 0;
                int n_cs = mm_gen_cs(km, &cs_str, &max_len, mi, r, subsequence.c_str(), true);

                LOGGER.debug << cs_str << std::endl;
                LOGGER.debug << subsequence << std::endl;

                unsigned int score = static_cast<unsigned int>(std::round((( static_cast<double>(r->blen+ r->mlen)/size + r->mapq/60.0)/3.0)*1000));

                if(region.score < score){
                    region.name = "Telomere";
                    region.score = score;
                    region.chromStart = start;
                    region.chromEnd = start + size;
                    region.thickStart = start + r->qs;
                    region.thickEnd = start + r->qe;
                    region.strand = "+-"[r->rev];
                }

            }

            if(region.score > 0) {
                result.push_back(region);
            }

            size = std::min(150UL,sequence.size()-start);
            start += std::min(50UL,sequence.size()-start);


        }

    }

    return result;

}
 *
 */