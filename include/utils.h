#ifndef GEAR_UTILS_H
#define GEAR_UTILS_H

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



#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <chrono>
#include "logger.h"
#include <stdlib.h>


namespace cmri {


    enum format_t : int {
        UNKNOWN = 0, FASTA = 1, FASTA_GZ = 2, FASTQ = 3, FASTQ_GZ = 4, CSV = 5, CSV_GZ = 6, BAM = 7
    };

    inline format_t file_format(std::string name) {

        std::vector<std::string> name_vec;
        boost::algorithm::split(name_vec, name, boost::is_any_of("."));
        size_t name_size = name_vec.size();
        if (name_size >= 1) {
            std::string ext1 = name_vec[name_size - 1];
            bool compressed = false;
            if (ext1 == "gzip" || ext1 == "gz") {
                if (name_size == 2) { return format_t::UNKNOWN; }
                ext1 = name_vec[name_size - 2];
                compressed = true;
            }

            if (ext1 == "bam") { return format_t::BAM; }
            if (ext1 == "fastq" || ext1 == "fq") {
                if (compressed) { return format_t::FASTQ_GZ; }
                return format_t::FASTQ;
            }
            if (ext1 == "fasta" || ext1 == "fa") {
                if (compressed) { return format_t::FASTA_GZ; }
                return format_t::FASTA;
            }
            if (ext1 == "csv") {
                if (compressed) { return format_t::CSV_GZ; }
                return format_t::CSV;
            }

        }
        return format_t::UNKNOWN;
    }


    template<class T>
    bool is_type(const boost::any &operand) {
        return operand.type() == typeid(T);
    }

    static std::ifstream open_file(const std::string &file_name) {
        std::ifstream ifs(file_name.c_str());
        if (!ifs) {
            LOGGER.error << "No such file or directory: " << file_name << std::endl;
            exit(ENOENT);
        }
        return ifs;
    }


    void show_options(boost::program_options::variables_map vm);

    void welcome(const std::string &code_name);

    std::string get_time_str(long value, const std::string& unit);

    void goodbye(std::chrono::system_clock::time_point start);

    inline bool is_a_valid_sequence(const std::string &sequence) {
        if (sequence.size() == 0) {
            return false;
        }

        for (int i = 0; i < sequence.size(); i++) {
            char base = sequence[i];
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
                return false;
            }
        }
        return true;
    }

    template<typename T>
    inline T clip(const T &n, const T &lower, const T &upper) {
        return std::max(lower, std::min(n, upper));
    }

    inline void create_output_directory(std::string path, bool backup) {

        int sys_out;
        if (backup) {
            std::string rmdir = "rm -rf " + path + "_prev";
            sys_out = system(rmdir.c_str());
            if (sys_out) {
                LOGGER.debug << "rm -rf return value " << sys_out << std::endl;
            }
            std::string mvdir = "mv -f " + path + " " + path + "_prev 2>/dev/null";
            sys_out = system(mvdir.c_str());
            if (sys_out) {
                LOGGER.debug << "mv -f return value " << sys_out << std::endl;
            }
        } else {
            std::string rmdir = "rm -rf " + path;
            sys_out = system(rmdir.c_str());
        }

        std::string mkdir = "mkdir " + path;
        sys_out = system(mkdir.c_str());

    }

    inline void log_command(std::string path, const int ac, char *av[]) {

        std::string cmd = "echo '#!/bin/bash' > " + path + "/command.sh";
        int sys_out = system(cmd.c_str());

        cmd = "echo cd `pwd` >> " + path + "/command.sh";
        sys_out = system(cmd.c_str());

        std::stringstream param;

        for (int i = 0; i < ac; i++)

            param << av[i] << " ";

        cmd = "echo " + param.str() + " >> " + path + "/command.sh;";
        sys_out = system(cmd.c_str());


        cmd = "chmod +x " + path + "/command.sh";
        sys_out = system(cmd.c_str());


    }


    inline int count_lines(const std::string &file_name) {

        int result = 0;
        format_t format = file_format(file_name);
        switch (format) {
            case FASTA :
            case FASTQ :
            case CSV :
                {
                std::ifstream file(file_name);
                result = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
                }
                break;
            case FASTA_GZ :
            case FASTQ_GZ :
            case CSV_GZ :
                {
                std::ifstream file(file_name, std::ios_base::in | std::ios_base::binary);
                boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
                inbuf.push(boost::iostreams::gzip_decompressor());
                inbuf.push(file);
                //Convert streambuf to istream
                std::istream instream(&inbuf);
                result=std::count(std::istreambuf_iterator<char>(instream), std::istreambuf_iterator<char>(), '\n');

                }
                break;
        }
        return result;
    }


}


#endif //GEAR_UTILS_H
