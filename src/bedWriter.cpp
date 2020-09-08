//
// Created by pablo on 12/8/20.
//

#include "bedWriter.h"

cmri::bedWriter::bedWriter(std::string file_name) {

    file.open(file_name);

    file << "#chrom"
         << "\tchromStart"
         << "\tchromEnd"
         << "\tname"
         << "\tscore"
         << "\tstrand"
         << "\tthickStart"
         << "\tthickEnd"
         << "\titemRgb"
         << std::endl;

}

cmri::bedWriter::~bedWriter() {
    file.close();
}

void cmri::bedWriter::append(std::string chromosome, int start, int end, std::string name, int score, bool strand, int thickStart,
                             int thickEnd, int r, int g,int b) {

    file << chromosome
         << "\t" << start
         << "\t" << end
         << "\t" << name
         << "\t" << score
         << "\t" << (strand ? "+" : "-")
         << "\t" << (thickStart > 0 ? thickStart : start)
         << "\t" << (thickEnd > 0 ? thickEnd : end)
         << "\t" << r << "," << g << "," << b
         << std::endl;

}
