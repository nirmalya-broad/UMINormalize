#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <experimental/filesystem>
#include <boost/program_options.hpp>

namespace fs = std::experimental::filesystem;
namespace po = boost::program_options;

#include "bam_reader.hpp"
#include "bam_writer.hpp"
#include "bam_record.hpp"

class uminorm {
    private:
        std::string infile_str;
        std::string outfile_str;
        std::string tempdir_str;
        std::string prefix_str;
        std::string coll_str;
        po::options_description desc;
        bam_hdr_t* lhdr = NULL;
    public:
        std::string get_temp_file(std::string lstr, unsigned int count); 
        void dump_sorted_records (std::vector<bam_record> brvec, std::string outfile_str,
            unsigned int temp_count, bam_hdr_t* lhdr);
        bool will_break_feature(bam_record first_rec, bam_record last_rec, bam_record this_rec);
        bool will_break_coordinate(bam_record first_rec, bam_record last_rec, bam_record this_rec);
        bool will_break(bam_record first_record, bam_record last_record, bam_record this_record, std::string coll_type);
        void write_collapse(std::vector<bam_record>& local_vec, std::ofstream& coll_writer);
        int split_files();
        void sort_and_merge_files(int split_count);
        void main_func();
        bool parse_args(int argc, char* argv[]);
        void print_help();
        std::string get_outfile_suffix_path(std::string suf);
        void initialize();
        
};

std::string uminorm::get_outfile_suffix_path(std::string suf) {
    fs::path lpath(outfile_str);
    std::string lfilename = lpath.filename();
    std::regex reg("(\\.[s|b]am$)");
    std::string new_name = std::regex_replace(lfilename, reg, suf);
    std::string res = tempdir_str + "/" + new_name;
    return res;

}

std::string uminorm::get_temp_file(std::string lstr, unsigned int count) {

    fs::path lpath(lstr);
    std::string lfilename = lpath.filename();
    std::regex reg("(\\.[s|b]am$)");
    std::string lcount = std::to_string(count);
    std::string new_name = std::regex_replace(lfilename, reg, "_temp_" + lcount + "$1");
    std::string res = tempdir_str + "/" + new_name;
    return res;
}

void uminorm::dump_sorted_records (std::vector<bam_record> brvec, 
        std::string outfile_str, unsigned int temp_count, bam_hdr_t* lhdr) {
    std::sort(brvec.begin(), brvec.end(), compare_bam_less());
    std::string temp_str = get_temp_file(outfile_str, temp_count);
    bam_writer writer(temp_str, lhdr);
    std::cout << "Dumping data to file: " << temp_str << "\n";
    std::cout << "Vector size: " << brvec.size() << "\n";
    for(const bam_record& lrecord: brvec) {
        writer.write_record(lrecord.full_rec);
    }
}

bool uminorm::will_break_feature(bam_record first_rec, bam_record last_rec, bam_record this_rec) {
    // Compare between last_rec and this_rec; in some cases first rec and
    // last_rec would be identical.

    // We should change this gap
    int brake_gap = 500;
    if (last_rec.ref_name_id != this_rec.ref_name_id) {
        return true;
    } else if (0 != strcmp(last_rec.umi, this_rec.umi)) {
        return true;
    } else if (last_rec.strand != this_rec.strand) {
        return true;
    } else if ((this_rec.start_pos - last_rec.start_pos) > brake_gap) {
        return true;
    } else {
        return false;
    }

}


bool uminorm::will_break_coordinate(bam_record first_rec, bam_record last_rec, bam_record this_rec) {
    if (last_rec.ref_name_id != this_rec.ref_name_id) {
        return true;
    } else if (0 != strcmp(last_rec.umi, this_rec.umi)) {
        return true;
    } else if (last_rec.strand != this_rec.strand) {
        return true;
    } else {
        return false;
    }
} 


bool uminorm::will_break(bam_record first_record, bam_record last_record, bam_record this_record, std::string coll_type) {
    if (coll_type.compare("feature")) {
        return will_break_feature(first_record, last_record, this_record);
    } else if (coll_type.compare("coordinate")) {
        return will_break_coordinate(first_record, last_record, this_record);
    } else {
        std::string throw_msg = "Illegal umi brake option: " + coll_type;
        throw std::runtime_error(throw_msg);
    }
}


void uminorm::write_collapse(std::vector<bam_record>& local_vec, std::ofstream& coll_writer) {
    // Get the last record, specifically the name of the query

    bam_record& last_rec = local_vec.back();
    std::string qname_str(last_rec.qname);
    unsigned long startPos = local_vec.front().start_pos;
    unsigned long endPos = local_vec.back().end_pos;
    int totalReads = local_vec.size();
    int totalGap = endPos - startPos + 1;
    std::string strand_str(1, local_vec.front().strand);
    coll_writer << "representative read: " << qname_str << " total_reads: " << totalReads << " gap: " << totalGap << " strand: " << strand_str << " start_pos: " << startPos << " end_pos: " << endPos << "\n";
     coll_writer << "------------------------------------\n";
    for (auto& lrec : local_vec) {
        std::string full_rec_str(lrec.full_rec);
        coll_writer << full_rec_str << "\n";
    }
    coll_writer << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

}

void initialize() {

}
int uminorm::split_files() {

    unsigned long size_lim = 100000000;

    bam_reader obj(infile_str);
    lhdr = obj.get_sam_header_dup();

    // if the tempdir does not exists, create it
    fs::path tempdir_path(tempdir_str);
    if (!fs::exists(tempdir_path)) {
        fs::create_directories(tempdir_path);
    }

    unsigned int temp_count = 0;
    while(true) {
        std::string ret_str;
        std::vector<bam_record> brvec;
        unsigned long used_size = 0;
        bam_record next_rec;

        while(!(ret_str = obj.read_record(next_rec)).empty()) {
            const int lsize = next_rec.get_size();
            //std::cout << "lsize: " << lsize << " used_size: " << used_size << "\n";
            used_size += lsize;
            brvec.push_back(std::move(next_rec));
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

    return temp_count;
}

void uminorm::sort_and_merge_files(int split_count) {

    std::map<unsigned int, bam_reader> reader_map;
    std::priority_queue<bam_record, std::vector<bam_record>, compare_bam_greater> bam_pq;

    for (unsigned int j = 1; j <= split_count; j++) {
        std::string temp_str = get_temp_file(outfile_str, j);
        std::cout << "Opening tempfile for reading: " << temp_str << "\n";
        //reader_map.emplace(j, temp_str);
        reader_map.emplace(j, temp_str);
        // Get the first read; it is expected that the first read would 
        // be useful.
        bam_record lrec;
        reader_map[j].read_record(lrec);
        lrec.reader_index = j;
        bam_pq.push(std::move(lrec));
    } 

    std::cout << "Size of priority queue: " << bam_pq.size() << "\n";

    // Create the final writer
    // Get the header from the first object and use it for creating the 
    // header of output
    //bam_hdr_t* lhdr = reader_map[1].get_sam_header();
    //std::set<unsigned int> empty_readers;
    std::string outfile_log_str = get_outfile_suffix_path("_log.txt");
    std::string sorted_sam_str = get_outfile_suffix_path("_sorted.sam");
    bam_writer writer(outfile_str, lhdr);
    std::cout << "sorted_sam_str: " << sorted_sam_str << "\n";
    bam_writer writer_sorted(sorted_sam_str, lhdr);

    bam_record first_record;
    bam_record last_record;   
    std::vector<bam_record> local_vec;
    std::ofstream outfile_log(outfile_log_str);
 
    bool fresh_start = true;
    while(!bam_pq.empty()) {
        bam_record lrec = bam_pq.top();
        //std::cout << "full_rec: " << std::string(lrec.full_rec) << ", is_mapped: " << std::boolalpha << lrec.is_mapped << "\n";
        if (lrec.is_mapped) {
            if (fresh_start) {
                // If this is the first record, please record all the information.
                first_record = lrec;
                last_record = lrec;
                fresh_start = false;
                local_vec.push_back(lrec);
            } else {
                bool break_status = will_break(first_record, last_record, lrec, coll_str);
                if (!break_status) {
                    last_record = lrec;
                    local_vec.push_back(lrec);
                } else {
                    writer.write_record(last_record.full_rec);
                    write_collapse(local_vec, outfile_log);
                    local_vec.clear();
                    first_record = lrec;
                    last_record = lrec;
                    local_vec.push_back(lrec);
                }
            }
        }
        writer_sorted.write_record(lrec.full_rec);
        unsigned int reader_index = lrec.reader_index;
        bam_pq.pop();
        // Write it to the outfile
        // get the index of the lrec and get one from that reader
        // 
        bam_record lrec_new;
        std::string ret = reader_map[reader_index].read_record(lrec_new);
        if (!ret.empty()) {
            // transfer the reader_index
            lrec_new.reader_index = reader_index;
            bam_pq.push(std::move(lrec_new));
        }
    }
}

void uminorm::main_func() {
    int split_count = split_files();
    sort_and_merge_files(split_count);
}

void uminorm::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: umi_norm -i <infile> -o <outfile> -t <tempdir> -c <collapse_type>"
    "\n\n";
}

bool uminorm::parse_args(int argc, char* argv[]) {

    bool all_set = true;
    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "input sam/bam file.")
        ("outfile,o", po::value<std::string>(&outfile_str), "output file.")
        ("tempdir,t", po::value<std::string>(&tempdir_str), "temporary directory.")
        ("collapse_type,c", po::value<std::string>(&coll_str), "type of collapse.")
        ;

        po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        return 0;
    } else {
    }

    if (vm.count("infile")) {
        std::cout << "Infile is set to: " << infile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: infile is not set.\n";
    }

    if (vm.count("outfile")) {
        std::cout << "Outfile is set to " << outfile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: outfile is not set.\n";
    }

    if (vm.count("tempdir")) {
        std::cout << "Tempdir is set to " << tempdir_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: tempdir is not set.\n";
    }

    if (vm.count("collapse_type")) {
        std::cout << "Collapse_type is set to " << coll_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: Collapse_type is not set.\n";
    }

    return all_set;

}

int main(int argc, char** argv) {

    uminorm uno;
    bool all_set = true;
    try {
        all_set = uno.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
    } catch(...) {
        return 0;
    }
    
    if (!all_set) {
        uno.print_help();
        return 0;
    }

    try {
        uno.main_func();
    } catch(const std::runtime_error& e) {
        std::cerr << "error: " << e.what() << "\n";
    }

    return 0;
}
