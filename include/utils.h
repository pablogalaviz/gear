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
#include <sam.h>
#include <vector>

namespace cmri {


    enum format_t : int {
        UNKNOWN = 0b0000'0000, FASTA = 0b0000'0001, CSV = 0b0000'0010, FASTQ = 0b0000'0100, BAM = 0b0000'1000, GZIP = 0b0001'0000 , FILE_TYPE = 0b0000'1111, FASTA_GZ = 0b0001'0001, CSV_GZ = 0b0001'0010, FASTQ_GZ = 0b0001'0100
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

    static std::ifstream open_file(const std::string &file_name, const std::string &error_message) {
        std::ifstream ifs(file_name.c_str());
        if (!ifs) {
            LOGGER.error << "No such file or directory: " << file_name << std::endl;
            LOGGER.error << error_message << std::endl;
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

    inline int get_total_reads(samFile *bam_file, bam_hdr_t *bam_header ){

        std::string index_file_name = bam_file->fn;
        index_file_name += ".bai";
        open_file(index_file_name,"expecting bam index").close(); // check if file exists
        auto bam_index = sam_index_load(bam_file, index_file_name.c_str());
        auto n_targets = bam_header->n_targets;
        int result = 0;
        for (int tid = 0; tid < n_targets; tid++) {
            uint64_t mapped;
            uint64_t unmapped;
            if (hts_idx_get_stat(bam_index, tid, &mapped, &unmapped) == 0) {
                result += static_cast<int>(mapped + unmapped);
            }
        }

        hts_idx_destroy(bam_index);
        return result;
    }

    inline int count_reads(const std::string &file_name) {

        int result = 0;
        int format = file_format(file_name);

        char token = (format&format_t::FASTA)  ? '>' :'\n';
        int factor = (format&format_t::FASTQ ) ? 4 : 1;

        if(format == BAM) {
            samFile *bam_file = hts_open(file_name.c_str(), "r");
            bam_hdr_t *bam_header = sam_hdr_read(bam_file); //read header
            return get_total_reads(bam_file, bam_header);
        }
        if( format&format_t::GZIP&(format_t::FASTA|format_t::FASTQ|format_t::CSV) ){

            std::ifstream file(file_name, std::ios_base::in | std::ios_base::binary);
            boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
            inbuf.push(boost::iostreams::gzip_decompressor());
            inbuf.push(file);
            //Convert streambuf to istream
            std::istream instream(&inbuf);
            return std::count(std::istreambuf_iterator<char>(instream), std::istreambuf_iterator<char>(), token)/factor;
        }
        if (format & (format_t::FASTA | format_t::FASTQ | format_t::CSV)) {
            std::ifstream file(file_name);
            return std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), token) /
                   factor;
        }
        return result;
    }


}


#endif //GEAR_UTILS_H
