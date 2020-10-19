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
#include <boost/program_options.hpp>
#include <chrono>

#include "csvParser.h"
#include "Modules/MotifCount/motifCount.h"
#include "Modules/VariantCallAnalysis/variantCallAnalysis.h"
#include "Modules/GenomeAnalysis/genomeAnalysis.h"
#include "utils.h"
#include "options.h"
#include "Modules/IwgsAnalysis/iwgsAnalysis.h"


int main(const int ac, char *av[]) {


    try {

        auto start_time = std::chrono::system_clock::now();

        std::string task{};
        std::string parameters;
        boost::program_options::options_description genericOptions(
                "GEAR, Genomic sEquence AnalyzeR.  \nAllowed common:");
        genericOptions.add_options()
                ("debug,d", "Shows debug messages in log")
                ("help,h", "Shows a help message")
                ("task", boost::program_options::value<std::string>(&task), "Perform one of the following tasks: [MotifCount, GenomeAnalysis, VariantCallAnalysis, IwgsAnalysis]")
                ("parameters,p", boost::program_options::value<std::string>(&parameters), "Parameters file")
                ("silent,s", "Shows only errors");


        cmri::common_options_t common;
        boost::program_options::options_description commonOptions("Common Options:");
        commonOptions.add_options()
                ("common.backup", boost::program_options::value<bool>(&common.backup)->default_value(true), "Create a backup of previous output")
                ("common.chunk_size", boost::program_options::value<int>(&common.chunk_size)->default_value(10000), "Size of the reading chuck")
                ("common.input_file,i", boost::program_options::value<std::string>(&common.input_file), "Input file")
                ("common.output_path,o", boost::program_options::value<std::string>(&common.output_path)->default_value("output"), "Output directory name")
                ("common.progress", boost::program_options::value<int>(&common.progress)->default_value(0), "Show progress message every X records (0 - off)")
                ("common.threads", boost::program_options::value<int>(&common.threads)->default_value(1), "Number of threads")
                ;

        cmri::motif_count_options_t motif_count;
        boost::program_options::options_description motifCountOptions("Motif Count Options:");
        motifCountOptions.add_options()
                ("motif_count.motifs", boost::program_options::value<std::string>(&motif_count.motif_file), "Motif per region definition in json format")
                ("motif_count.quality_value", boost::program_options::value<int>(&motif_count.quality_value)->default_value(0), "Mean base quality threshold")
                ("motif_count.quality_map", boost::program_options::value<int>(&motif_count.quality_map)->default_value(0), "Quality Mapping threshold")
                ("motif_count.validate", boost::program_options::value<bool>(&motif_count.validate_sequence)->default_value(false), "Validate sequences (slow)")
        ;

        cmri::variant_call_analysis_options_t variant_call_analysis;
        boost::program_options::options_description variantCallAnalysisOptions("Variant Call Analysis Options:");
        variantCallAnalysisOptions.add_options()
                ("variant_call_analysis.regions", boost::program_options::value<std::string>(&variant_call_analysis.regions), "Variants per region definition in json format")
                ("variant_call_analysis.reference", boost::program_options::value<std::string>(&variant_call_analysis.reference), "Fasta reference (required index).")
                ;

        cmri::genome_analysis_options_t genome_analysis;
        boost::program_options::options_description genomeAnalysisOptions("Genome Analysis Options:");
        genomeAnalysisOptions.add_options()
                ("genome_analysis.regions", boost::program_options::value<std::string>(&genome_analysis.regions), "Chromosome region definition in json format")
                ;

        cmri::iwgs_analysis_options_t iwgs_analysis;
        boost::program_options::options_description iwgsAnalysisOptions("Illumina WGS Analysis Options:");
        genomeAnalysisOptions.add_options()
                ("iwgs_analysis.input_file", boost::program_options::value<std::string>(&iwgs_analysis.input_file), "Second fastq input file")
                ("iwgs_analysis.variants_file", boost::program_options::value<std::string>(&iwgs_analysis.variants_file), "List of variants")
                ;

        boost::program_options::positional_options_description positional;
        positional.add("task", 1);

        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions)
        .add(commonOptions)
        .add(motifCountOptions)
        .add(variantCallAnalysisOptions)
        .add(genomeAnalysisOptions)
        .add(iwgsAnalysisOptions)
                ;

        boost::program_options::options_description configFileOptions;
        configFileOptions.add(commonOptions)
        .add(motifCountOptions)
        .add(variantCallAnalysisOptions)
        .add(genomeAnalysisOptions)
        .add(iwgsAnalysisOptions)
                ;

        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::command_line_parser(ac, av).options(cmdlineOptions).positional(positional).run(),vm);
        boost::program_options::notify(vm);

        if (vm.count("help") || vm.count("task") == 0 ) {
            if (vm.count("task") == 0 && vm.count("help") == 0)
                std::cout << "MISSING TASK OPTION!" << std::endl;
            std::cerr << cmdlineOptions << std::endl;
            return 0;
        }

        if (vm.count("parameters"))
        {
            auto par_file = cmri::open_file(parameters,"expecting parameters file");
            store(parse_config_file(par_file, configFileOptions), vm);
            notify(vm);
        }


        cmri::create_output_directory(common.output_path, common.backup);
        cmri::log_command(common.output_path, ac, av);
        std::string log_file = common.output_path + "/output.log";

        bool debug = vm.count("debug");
        bool silent = vm.count("silent");
        if (silent)
            cmri::LOGGER.init(cmri::log_t::ERROR, log_file);
        else {
            if (debug)
                cmri::LOGGER.init(cmri::log_t::DBG, log_file);
            else
                cmri::LOGGER.init(cmri::log_t::INFO, log_file);
        }

        cmri::welcome("GEAR, Genomic sEquence AnalyzeR.");

        cmri::show_options(vm);
        common.validate();

        if(task == "MotifCount") {
            motif_count.validate();
            cmri::mainMotifCount(common,motif_count);
        }
        else {
            if (task == "VariantCallAnalysis") {
                variant_call_analysis.validate();
                cmri::mainVariantCallAnalysis(common, variant_call_analysis);
            }else{
                if (task == "GenomeAnalysis") {
                    genome_analysis.validate();
                    cmri::mainGenomeAnalysis(common, genome_analysis);
                }else{
                    if (task == "IwgsAnalysis") {
                        iwgs_analysis.validate();
                        cmri::mainIwgsAnalysis(common, iwgs_analysis);
                    }else{
                        cmri::LOGGER.error << "Unknown task: " << task;
                        std::cerr << cmdlineOptions << std::endl;
                    }
                }
            }
        }

        cmri::goodbye(start_time);

    }
    catch (std::exception &e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return -1;
    }
    catch (...) {
        std::cerr << "Exception of unknown type!" << std::endl;
        return -1;
    }


    return 0;
}
