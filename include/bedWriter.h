//
// Created by pablo on 12/8/20.
//

#ifndef GEAR_BEDWRITER_H
#define GEAR_BEDWRITER_H

#include <fstream>

namespace cmri {

    class bedWriter {

        std::ofstream file;

    public:
        bedWriter(std::string file_name);

        void append(std::string chromosome, int start, int end, std::string _name, int score =0, bool strand=true, int thickStart=0, int thickEnd=0, int r=255, int g=0,int b=0);

        virtual ~bedWriter();

    };



}

#endif //GEAR_BEDWRITER_H
