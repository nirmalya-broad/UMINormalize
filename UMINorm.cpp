#include <iostream>
#include <string>
#include "bam_reader.hpp"
#include "bam_writer.hpp"
#include "bam_record.hpp"
#include <vector>

void dump_sorted_records (std::vector<bam_record> sorted_vec, std::string temp_outfile, bam_hdr_t* lhdr) {
    bam_writer writer(temp_outfile, lhdr);
    for(const bam_record& lrecord: sorted_vec) {
        writer.write_record(lrecord.full_rec);
    }
}


std::string get_temp_file(std::string lstr, int count) {

    std::regex reg("(\\.[s|b]am$)");
    std::string lcount = std::to_string(count);
    std::string res = std::regex_replace(lstr, reg, "_temp_" + lcount + "$1");
    return res;
}


int main(int argc, char** argv) {
    std::string infile_str(argv[1]);
    std::string outfile_str(argv[2]);
    
    unsigned long size_lim = 100000000;

    bam_reader obj(infile_str);
    bam_hdr_t* lhdr = obj.get_sam_header();


    for(unsigned int temp_count = 1; ; temp_count++) {
        std::string ret_str;
        std::vector<bam_record> brvec;
        unsigned long used_size = 0;
        bam_record next_rec;

        while(!(ret_str = obj.read_record(next_rec)).empty()) {
            const int lsize = next_rec.get_size();
            used_size += lsize;
            brvec.push_back(next_rec);
            if (used_size > size_lim) {
                // Create a new file
                std::sort(brvec.begin(), brvec.end());
                // Dump the sorted records to the temp file
                std::string temp_str = get_temp_file(outfile_str, temp_count);
                std::cout << "Dumping data to file: " << temp_str << "\n";
                std::cout << "Vector size: " << brvec.size() << "\n";
                dump_sorted_records(brvec, temp_str, lhdr); 
                break;
            }
        }
        if (ret_str.empty()) {
            break;
        }
    }

    return 0;
}
