#ifndef _BED_WRITER_HPP
#define _BED_WRITER_HPP

#include <htslib/sam.h>
#include <fstream>


class bed_writer {

    public:

    bed_writer(std::string& outfile_str): 
        outfile_str(outfile_str),
        outfile(outfile_str) {
    }

    void write_record_str(const std::string& record) {
        outfile << record << "\n";
    }

    void write_record_items (

        const std::string& chrom,
        const unsigned long chromStart,
        const unsigned long chromEnd,
        const std::string& name,
        const int score,
        const std::string& strand) {

        std::string out_str = chrom + "\t" + 
            std::to_string(chromStart) + "\t" +
            std::to_string(chromEnd) + "\t" +
            name + "\t" +
            std::to_string(score) + "\t" +
            strand;

        outfile << out_str << "\n";
    }
        
    private:

    std::string outfile_str;
    std::ofstream outfile;    

};

#endif
