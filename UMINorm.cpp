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



class args_c {
    public:
        po::options_description desc;
        std::string infile_str;
        std::string outdir_str;
        std::string prefix_str;
        std::string coll_str;
        bool parse_args(int argc, char* argv[]); 
        void print_help();
};

class uminorm {
    private:
        std::string infile_str;
        std::string outfile_str;
        std::string outfile_all_str;
        std::string outdir_str;
        std::string logdir_str;
        std::string prefix_str;
        std::string coll_str;
        int total_split_count = 0;
        bam_reader obj;
        unsigned long size_lim = 200000000;
        int brake_gap = 500;
        bam_hdr_t* lhdr = NULL;
    public:
        uminorm(args_c args_o);
        std::string get_temp_file(unsigned int count); 
        void dump_sorted_records (std::vector<bam_record> brvec, 
            unsigned int temp_count, bam_hdr_t* lhdr);
        bool will_break_feature(bam_record first_rec, bam_record last_rec, bam_record this_rec);
        bool will_break_coordinate(bam_record first_rec, bam_record last_rec, bam_record this_rec);
        bool will_break(bam_record first_record, bam_record last_record, bam_record this_record, std::string coll_type);
        void write_collapse(std::vector<bam_record>& local_vec, std::ofstream& coll_writer);
        void split_n_sort_files();
        void merge_files();
        void main_func();
        bool parse_args(int argc, char* argv[]);
        void print_help();
        std::string get_outfile_suffix_path(std::string suf);
        void initialize();
        void clean();
        
};

uminorm::uminorm(args_c args_o)
    : infile_str(args_o.infile_str),
    outdir_str(args_o.outdir_str),
    prefix_str(args_o.prefix_str),
    coll_str(args_o.coll_str),
    obj(infile_str) {}

std::string uminorm::get_outfile_suffix_path(std::string suf) {
    std::string res = logdir_str + "/" + prefix_str + suf;
    return res;

}

std::string uminorm::get_temp_file(unsigned int count) {

    std::string res = logdir_str + "/" + prefix_str + "_" + std::to_string(count) + ".bam";
    return res;
}

void uminorm::dump_sorted_records (std::vector<bam_record> brvec, 
        unsigned int temp_count, bam_hdr_t* lhdr) {
    std::sort(brvec.begin(), brvec.end(), compare_bam_less());
    std::string temp_str = get_temp_file(temp_count);
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

void uminorm::initialize() {
    fs::path outdir_path(outdir_str);
    if (!fs::exists(outdir_path)) {
        fs::create_directories(outdir_path);
    }

    lhdr = obj.get_sam_header();
    // Outdir would be those place dedicated specifically for UMI
    outfile_str = outdir_str + "/" + prefix_str + "_u.bam";
    logdir_str = outdir_str + "/logdir";
    
    fs::path logdir_path (logdir_str);
    if (!fs::exists(logdir_path)) {
        fs::create_directories(logdir_path);
    }
}

void uminorm::split_n_sort_files() {


    // We duplicate the sam header to keep it alive after this function exits.
    // ideally we have to have a copy constructor for bam_reader.

    // if the outdir does not exists, create it

    unsigned int split_count = 0;
    while(true) {
        std::string ret_str;
        std::vector<bam_record> brvec;
        unsigned long used_size = 0;
        bam_record next_rec;

        while(!(ret_str = obj.read_record(next_rec)).empty()) {
            if (next_rec.is_mapped) {
                const int lsize = next_rec.get_size();
                //std::cout << "lsize: " << lsize << " used_size: " << used_size << "\n";
                used_size += lsize;
                brvec.push_back(std::move(next_rec));
                if (used_size > size_lim) {
                    split_count++;
                    dump_sorted_records(brvec, split_count, lhdr); 
                    break;
                }
            }
        }
        
        if (ret_str.empty()) {
            if (used_size > 0) {
                split_count++;
                dump_sorted_records(brvec, split_count, lhdr);
            }
            break;
        }
    }
    total_split_count = split_count;

}

void uminorm::merge_files() {

    std::map<unsigned int, bam_reader> reader_map;
    std::priority_queue<bam_record, std::vector<bam_record>, compare_bam_greater> bam_pq;

    for (unsigned int j = 1; j <= total_split_count; j++) {
        std::string temp_str = get_temp_file(j);
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
    std::string sorted_bam_str = get_outfile_suffix_path("_sorted.bam");
    bam_writer writer(outfile_str, lhdr);
    std::cout << "sorted_sam_str: " << sorted_bam_str << "\n";
    bam_writer writer_sorted(sorted_bam_str, lhdr);

    bam_record first_record;
    bam_record last_record;   
    std::vector<bam_record> local_vec;
    std::ofstream outfile_log(outfile_log_str);
 
    bool fresh_start = true;
    unsigned long lcount = 0;
    while(!bam_pq.empty()) {
        bam_record lrec = bam_pq.top();
        //std::cout << "full_rec: " << std::string(lrec.full_rec) << ", is_mapped: " << std::boolalpha << lrec.is_mapped << "\n";
        if (lrec.is_mapped) {
            lcount++;
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
            writer_sorted.write_record(lrec.full_rec);
        }
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
    std::cout << "mapped_count: " << lcount << "\n";
}

void uminorm::clean() {
    for (unsigned int j = 1; j <= total_split_count; j++) {
        std::string temp_str = get_temp_file(j);
        
        fs::path temp_path(temp_str);
        bool n = fs::remove(temp_path);
        std::cout << "Deleted: " << temp_str << "\n";
    }
}

void uminorm::main_func() {
    split_n_sort_files();
    merge_files();
}

void args_c::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: umi_norm -i <infile> -o <outdir> -p <prefix> -c <collapse_type>"
    "\n\n";
}

bool args_c::parse_args(int argc, char* argv[]) {

    bool all_set = true;
    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "Input sam/bam file.")
        ("prefix,p", po::value<std::string>(&prefix_str), "Prefix.")
        ("outdir,o", po::value<std::string>(&outdir_str), "Output directory.")
        ("collapse_type,c", po::value<std::string>(&coll_str), "Type of collapse.")
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

    if (vm.count("outdir")) {
        std::cout << "Outdir is set to " << outdir_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: outdir is not set.\n";
    }

    if (vm.count("prefix")) {
        std::cout << "Prefix is set to " << prefix_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: prefix is not set.\n";
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
    
    args_c args_o;
    bool all_set = true;
    try {
        all_set = args_o.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
    } catch(...) {
        return 0;
    }
    
    if (!all_set) {
        args_o.print_help();
        return 0;
    }

    uminorm uno(args_o);
    try {
        uno.initialize();
        uno.main_func();
        uno.clean();
    } catch(const std::runtime_error& e) {
        std::cerr << "error: " << e.what() << "\n";
    }

    return 0;
}
