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
#include "motifCount.h"
#include "utils.h"


int main(const int ac, char *av[]) {


    try {

        auto start_time = std::chrono::system_clock::now();

        boost::program_options::options_description genericOptions(
                "GEAR, Genomic sEquence AnalyzeR.  \nAllowed options:");

        std::string output_path;
        std::string motif_file;
        cmri::options_t options;
        bool backup;
        genericOptions.add_options()
                ("backup,b", boost::program_options::value<bool>(&backup)->default_value(true),"Create a backup of previous output")
                ("chunk_size,c", boost::program_options::value<int>(&options.chunk_size)->default_value(10000), "Size of the reading chuck")
                ("debug,d", "Shows debug messages in log")
                ("help,h", "Shows a help message")
                ("input,i", boost::program_options::value<std::string>(&options.input_file), "Input file")
                ("motifs,m", boost::program_options::value<std::string>(&motif_file),"Motif per region definition in json format")
                ("output,o", boost::program_options::value<std::string>(&output_path)->default_value("output"),"Output directory name")
                ("progress,p", boost::program_options::value<int>(&options.progress)->default_value(0), "Show progress message every X records (0 - off)")
                ("quality_value", boost::program_options::value<int>(&options.quality_value)->default_value(0), "Mean base quality threshold")
                ("quality_map", boost::program_options::value<int>(&options.quality_map)->default_value(0), "Quality Mapping threshold")
                ("silent,s", "Shows only errors")
                ("threads,t", boost::program_options::value<int>(&options.threads)->default_value(1),"Number of threads")
                ("validate,v", boost::program_options::value<bool>(&options.validate_sequence)->default_value(false),"Validate sequences (slow)")
                ;

        boost::program_options::options_description cmdlineOptions;
        cmdlineOptions.add(genericOptions);

        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::command_line_parser(ac, av).options(cmdlineOptions).run(),
                                      vm);
        boost::program_options::notify(vm);

        cmri::create_output_directory(output_path, backup);
        cmri::log_command(output_path, ac, av);
        std::string log_file = output_path + "/output.log";

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

        if (vm.count("help") || vm.count("input") == 0 || vm.count("motifs") == 0) {
            if (vm.count("input") == 0 && vm.count("help") == 0)
                std::cout << "MISSING INPUT FILE!" << std::endl;
            if (vm.count("motifs") == 0 && vm.count("help") == 0)
                std::cout << "MISSING MOTIF FILE!" << std::endl;

            std::cerr << cmdlineOptions << std::endl;
            return 0;
        }


        cmri::show_options(vm);

        auto m_file = cmri::open_file(motif_file);
        std::map<std::string, cmri::region_list_t> motifs = cmri::deserialize(motif_file);

        options.validate();


        if (motifs.empty()) {
            cmri::LOGGER.error << "From: " << __FILE__ << ":" << __LINE__ << std::endl;
            cmri::LOGGER.error << "Invalid argument. No valid motif list found!" << std::endl;
            return (EINVAL);
        }


        std::string output_file = output_path + "/output.csv";
        cmri::format_t format = cmri::file_format(options.input_file);
        switch (format) {

            case cmri::UNKNOWN:
                cmri::LOGGER.error << "Unknown input file format " << std::endl;
                return EINVAL;
                break;
            case cmri::FASTA:
            case cmri::FASTQ:
            case cmri::FASTA_GZ:
            case cmri::FASTQ_GZ:
                if(options.threads >1) {
                    cmri::processFastMultiThreading(options, motifs);
                }
                else{
                    cmri::processFast(options, motifs);
                }
                break;
            case cmri::CSV:
                if(options.threads >1) {
                    cmri::processCsvMultiThreading(options, motifs);
                }
                else{
                    cmri::processCsv(options, motifs);
                }
                break;
            case cmri::CSV_GZ:
                if(options.threads >1) {
                    cmri::processCsvMultiThreading(options, motifs,true);
                }
                else{
                    cmri::processCsv(options, motifs, true);
                }
                break;
            case cmri::BAM:
                if(options.threads >1) {
                    cmri::processBamMultiThreading(options, motifs);
                }
                else{
                    cmri::processBam(options, motifs);
                }
                break;
        }


        cmri::serialize(output_path + "/output.json", motifs);

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
