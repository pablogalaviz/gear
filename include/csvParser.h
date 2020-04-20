#ifndef GEAR_CSVPARSER_H
#define GEAR_CSVPARSER_H

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

#include <string>
#include <utility>
#include <vector>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include "logger.h"

namespace cmri {

    class csvParser {
        std::vector<std::string> header;
        bool raise_warning;
        std::string separator;
    public:

        explicit csvParser(const std::string &header_str,
                           const std::vector<std::string> &required = std::vector<std::string>(),
                           const std::string separator = ",", const bool &raise_warning = false
        ) : separator(separator), raise_warning(raise_warning) {
            boost::algorithm::split(header, header_str, boost::is_any_of(separator));
            if (header.size() == 0 || header_str == "") {
                cmri::LOGGER.error << "From " << __FILE__ << ":" << __LINE__ << std::endl;
                cmri::LOGGER.error << "Invalid argument: Empty csv header" << std::endl;
                exit(EINVAL);
            }
            if (required.size() > 0) {
                for (auto &r : required) {
                    bool found = false;
                    for (auto &h : header) {
                        if (h == r) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        cmri::LOGGER.error << "From " << __FILE__ << ":" << __LINE__ << std::endl;
                        cmri::LOGGER.error << "Missing required field <" << r << "> in csv file." << std::endl;
                        exit(EINVAL);
                    }
                }


            }


        }

        virtual ~csvParser() = default;

        std::map<std::string, std::string> parseLine(const std::string &line) const {

            std::map<std::string, std::string> result;

            std::vector<std::string> fields;
            boost::algorithm::split(fields, line, boost::is_any_of(separator));
            if (fields.size() != header.size() && raise_warning) {
                LOGGER.warning << "Expecting line of size " << header.size() << " but got " << fields.size()
                               << " fields. ";
                LOGGER.warning << "Possible corrupted parsing from " << __FILE__ << ":" << __LINE__ << " fields. ";
            }
            int parse_lines = std::min(fields.size(), header.size());

            //#pragma omp parallel for
            for (int i = 0; i < parse_lines; i++) {
                result[header[i]] = fields[i];
            }

            return result;
        }
    };

}


#endif //GEAR_CSVPARSER_H
