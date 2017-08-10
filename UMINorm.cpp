#include <iostream>
#include <string>
#include <vector>
#include <queue>

#include "bam_reader.hpp"
#include "bam_writer.hpp"
#include "bam_record.hpp"

std::string get_temp_file(std::string lstr, unsigned int count) {

    std::regex reg("(\\.[s|b]am$)");
    std::string lcount = std::to_string(count);
    std::string res = std::regex_replace(lstr, reg, "_temp_" + lcount + "$1");
    return res;
}

void dump_sorted_records (std::vector<bam_record> brvec, std::string outfile_str, 
        unsigned int temp_count, bam_hdr_t* lhdr) {
    std::sort(brvec.begin(), brvec.end(), compare_bam_less());
    std::string temp_str = get_temp_file(outfile_str, temp_count);
    bam_writer writer(temp_str, lhdr);
    std::cout << "Dumping data to file: " << temp_str << "\n";
    std::cout << "Vector size: " << brvec.size() << "\n";
    for(const bam_record& lrecord: brvec) {
        writer.write_record(lrecord.full_rec);
    }
}


int main(int argc, char** argv) {
    std::string infile_str(argv[1]);
    std::string outfile_str(argv[2]);
    
    unsigned long size_lim = 100000000;

    bam_reader obj(infile_str);
    bam_hdr_t* lhdr = obj.get_sam_header();

    unsigned int temp_count = 0;
    while(true) {
        std::string ret_str;
        std::vector<bam_record> brvec;
        unsigned long used_size = 0;
        bam_record next_rec;

        while(!(ret_str = obj.read_record(next_rec)).empty()) {
            const int lsize = next_rec.get_size();
            used_size += lsize;
            brvec.push_back(next_rec);
            if (used_size > size_lim) {
                temp_count++;
                dump_sorted_records(brvec, outfile_str, temp_count, lhdr); 
                break;
            }
        }
        
        if (ret_str.empty()) {
            if (used_size > 0) {
                temp_count++;
                dump_sorted_records(brvec, outfile_str, temp_count, lhdr);
            }
            break;
        }
    }


    // k-way merging and writing to the outfile

    std::map<unsigned int, bam_reader> reader_map;
    std::priority_queue<bam_record, std::vector<bam_record>, compare_bam_greater> bam_pq;

    for (unsigned int j = 1; j <= temp_count; j++) {
        std::string temp_str = get_temp_file(outfile_str, j);
        std::cout << "Opening tempfile for reading: " << outfile_str << "\n";
        //reader_map.emplace(j, temp_str);
        reader_map.emplace(j, temp_str);
        // Get the first read; it is expected that the first read would 
        // be useful.
        bam_record lrec;
        reader_map[j].read_record(lrec);
        lrec.reader_index = j;
        bam_pq.push(lrec);
    } 

    std::cout << "Size of priority queue: " << bam_pq.size() << "\n";

    // Create the final writer
    // Get the header from the first object and use it for creating the 
    // header of output
    //bam_hdr_t* lhdr = reader_map[1].get_sam_header();
    bam_writer writer(outfile_str, lhdr);
    //std::set<unsigned int> empty_readers;

    while(!bam_pq.empty()) {
        bam_record lrec = bam_pq.top();
        unsigned int reader_index = lrec.reader_index;
        bam_pq.pop();
        // Write it to the outfile
        writer.write_record(lrec.full_rec);
        // get the index of the lrec and get one from that reader
        // 
        bam_record lrec_new;
        std::string ret = reader_map[reader_index].read_record(lrec_new);
        if (!ret.empty()) {
            // transfer the reader_index
            lrec_new.reader_index = reader_index;
            bam_pq.push(lrec_new);
        }
    }

    return 0;
}
