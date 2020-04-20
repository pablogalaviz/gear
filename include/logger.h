#ifndef GEAR_LOGGER_H
#define GEAR_LOGGER_H

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
#include <map>
#include <fstream>
#include <iomanip>
#include <iostream>

#define LOGGER Logger::getInstance()

namespace cmri
{

    enum class log_t : int { DBG=0, INFO=1, WARNING=2, ERROR=3 };

    class Logger
    {
    public :

        static Logger& getInstance()
        {
            static Logger instance;
            return instance;
        }

    private :

        Logger(){};

        class LoggerStream
        {

            bool m_new_line;
            log_t m_level;
            log_t m_current;
            std::string m_label;

            std::map<log_t,std::string>color_map = {
                    {log_t::DBG  , "\033[1;32m"},
                    {log_t::INFO   , "\033[1;94m"},
                    {log_t::WARNING, "\033[1;33m"},
                    {log_t::ERROR  , "\033[1;31m"}};

        public:

            std::ofstream file;

            LoggerStream(){};

            ~LoggerStream(){file.close ();};

            void init(const log_t level,
                      const std::string log_file,
                      const log_t current,
                      const std::string label)
            {
                m_new_line = true;
                m_current=current;
                m_label = label;
                m_level=level;
                file.open (log_file, std::ios::out);
            };

            template<class T>
            LoggerStream&  operator<< (const T out)
            {
                if(m_level <= m_current)
                {
                    std::stringstream ss;
                    if(m_new_line)
                    {
                        std::time_t t = std::time(nullptr);
                        std::tm tm = *std::localtime(&t);
                        ss << color_map[m_current]
                           << m_label
                           << std::put_time(&tm, "[%F %T] | ")
                           << "\033[0m";
                        m_new_line = false;
                    }
                    ss << out;
                    file << ss.str();
                    file.flush();
                    {
                        if(m_current == log_t::ERROR)
                            std::cerr << ss.str();
                        else
                            std::cout << ss.str();
                    }

                }
                return *this;
            }

            LoggerStream& operator<< (std::ostream& (*pfun)(std::ostream&))
            {
                if(m_level <= m_current)
                {
                    pfun(file);
                    if(m_current == log_t::ERROR)
                        pfun(std::cerr);
                    else
                        pfun(std::cout);
                    m_new_line=true;
                }
                return *this;
            }

        };

        int m_id;

    public:

        Logger(Logger const&) = delete;
        void operator=(Logger const&)  = delete;

        inline void init(const log_t level,const std::string& log_file)
        {
            debug.init(level, log_file, log_t::DBG, "Debug   ");
            info.init(level, log_file, log_t::INFO, "Info    ");
            warning.init(level, log_file, log_t::WARNING, "Warning ");
            error.init(level, log_file, log_t::ERROR, "Error   ");
        }

        LoggerStream debug;
        LoggerStream info;
        LoggerStream warning;
        LoggerStream error;

    };


}

#endif //GEAR_LOGGER_H
