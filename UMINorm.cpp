#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <chrono>
#include <random>
#include <experimental/filesystem>
//#include <filesystem>
#include <boost/program_options.hpp>

namespace fs = std::experimental::filesystem;
namespace po = boost::program_options;

#include "bam_reader.hpp"
#include "bam_writer.hpp"
#include "bam_record.hpp"
#include "bed_writer.hpp"

class args_c {
    public:
        po::options_description desc;
        std::string infile_str;
        std::string outdir_str;
        std::string prefix_str;
        std::string coll_str;
        unsigned int size_lim_M;
        bool parse_args(int argc, char* argv[]); 
        void print_help();
};

class uminorm {
    private:
        std::string infile_str;
        std::string outfile_str;
        std::string bedfile_str;
        std::string gapfile_str;
        std::string outfile_all_str;
        std::string outdir_str;
        std::string logdir_str;
        std::string prefix_str;
        std::string coll_str;
        int total_split_count = 0;
        bam_reader obj;
        unsigned int size_lim_M;
        unsigned long size_lim;
        int brake_gap = 500;
        bam_hdr_t* lhdr = NULL;
        unsigned seed = 100;
        std::default_random_engine generator;
    public:
        uminorm(args_c args_o);
        int get_rand_pos(int vec_size);
        std::string get_temp_file(unsigned int count); 
        void dump_sorted_records (std::vector<bam_record> brvec, 
            unsigned int temp_count, bam_hdr_t* lhdr);
        bool will_break_feature(bam_record first_rec, bam_record last_rec, bam_record this_rec);
        bool will_break_coordinate(bam_record first_rec, bam_record last_rec, bam_record this_rec);
        bool will_break(bam_record first_record, bam_record last_record, bam_record this_record, std::string coll_type);
        void write_collapse(std::vector<bam_record>& local_vec, std::ofstream& coll_writer, std::ofstream& coll_len,  int final_pos);
        void split_n_sort_files();
        void merge_files();
        void main_func();
        bool parse_args(int argc, char* argv[]);
        void print_help();
        std::string get_outfile_suffix_path(std::string suf);
        void initialize();
        void clean();
        std::string get_bed_str(std::vector<bam_record> local_vec);
        void throw_ineq_exception(std::string first_str, std::string sec_str);
        void throw_neg_execption(long lvar);

        std::string get_gap_str(const bam_record& first_rec, 
            const bam_record& last_rec);

        bool will_write_gap_coordinate(bam_record last_rec, 
            bam_record this_rec);

};

uminorm::uminorm(args_c args_o)
    : infile_str(args_o.infile_str),
    outdir_str(args_o.outdir_str),
    prefix_str(args_o.prefix_str),
    coll_str(args_o.coll_str),
    size_lim_M(args_o.size_lim_M),
    obj(infile_str),
    generator(seed) {
        size_lim = size_lim_M * 1000000;    
    }

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

bool uminorm::will_break_coordinate(bam_record first_rec, bam_record last_rec, bam_record this_rec) {
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

bool uminorm::will_write_gap_coordinate(bam_record last_rec, bam_record this_rec) {
    // Compare between last_rec and this_rec; in some cases first rec and
    // last_rec would be identical.

    // We should change this gap
    if (last_rec.ref_name_id != this_rec.ref_name_id) {
        return false;
    } else if (0 != strcmp(last_rec.umi, this_rec.umi)) {
        return false;
    } else if (last_rec.strand != this_rec.strand) {
        return false;
    } else {
        return true;
    }
}

bool uminorm::will_break_feature(bam_record first_rec, bam_record last_rec, bam_record this_rec) {
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
    if (0 == coll_type.compare("feature")) {
        return will_break_feature(first_record, last_record, this_record);
    } else if (0 == coll_type.compare("coordinate")) {
        return will_break_coordinate(first_record, last_record, this_record);
    } else {
        std::string throw_msg = "Illegal umi brake option: " + coll_type;
        throw std::runtime_error(throw_msg);
    }
}

void uminorm::write_collapse(std::vector<bam_record>& local_vec, std::ofstream& coll_writer, std::ofstream& coll_len, int final_pos) {
    // Get the last record, specifically the name of the query

    bam_record& last_rec = local_vec.back();
    std::string qname_str(last_rec.qname);
    unsigned long startPos = local_vec.front().start_pos;
    unsigned long endPos = local_vec.back().end_pos;
    int totalReads = local_vec.size();
    int totalGap = endPos - startPos + 1;
    std::string strand_str(1, local_vec.front().strand);
    coll_writer << "representative read: " << qname_str << " total_reads: " << totalReads << " gap: " << totalGap << " final_pos: " << final_pos << " strand: " << strand_str << " start_pos: " << startPos << " end_pos: " << endPos << "\n";
    coll_len << totalReads << "\n";
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
    bedfile_str = outdir_str + "/" + prefix_str + ".bed";
    gapfile_str = outdir_str + "/" + prefix_str + "_gap.txt";
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

    unsigned long read_counter = 0;
    unsigned int split_count = 0;
    while(true) {
        std::string ret_str;
        std::vector<bam_record> brvec;
        unsigned long used_size = 0;
        bam_record next_rec;

        while(!(ret_str = obj.read_record(next_rec)).empty()) {
            read_counter++;
            if (read_counter %100000 == 0) {
                std::cout << "The value of read_counter: " << std::to_string(read_counter) << "\n";
                std::cout << "ret_str: ---" << ret_str << "---\n";
            }
            if (next_rec.is_mapped) {
                const int lsize = next_rec.get_size();
                //std::cout << "lsize: " << lsize << " used_size: " << used_size << "\n";
                used_size += lsize;
                brvec.push_back(std::move(next_rec));
                if (used_size > size_lim) {
                    split_count++;
                    dump_sorted_records(brvec, split_count, lhdr); 
                    std::cout << "split_count: " << std::to_string(split_count) << "\n";
                    break;
                }
            }
        }
        std::cout << "Reached out of the while loop" << "\n"; 
        if (ret_str.empty()) {
            if (used_size > 0) {
                split_count++;
                dump_sorted_records(brvec, split_count, lhdr);
                    std::cout << "split_count: " << std::to_string(split_count) << "\n";
            }
            break;
        }
    }
    total_split_count = split_count;
    std::cout << "Reached end of split and sort" << "\n"; 

}

int uminorm::get_rand_pos(int vec_size) {
    if (vec_size == 1) {
        return 0;
    } else {
        int vec_size_t = vec_size - 1;
        std::uniform_int_distribution<int> distribution(0, vec_size_t);
        return distribution(generator);
        
    }
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
    std::string coll_len_str = get_outfile_suffix_path("_coll_len.txt");
    bam_writer writer(outfile_str, lhdr);

    bed_writer bwriter(bedfile_str);

    std::ofstream gwriter(gapfile_str);

    std::cout << "sorted_sam_str: " << sorted_bam_str << "\n";
    bam_writer writer_sorted(sorted_bam_str, lhdr);

    bam_record first_record;
    bam_record last_record;   
    std::vector<bam_record> local_vec;
    std::ofstream outfile_log(outfile_log_str);
    std::ofstream coll_len(coll_len_str);
 
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

                // Get the gap between lrec and last_record
                // Get the starting position between lrec and last_record
                // Here last_record is actually the first record which is 
                // on the left side on the coordinate and lrec is the 
                // second record which is on the right side of the coordinate.

                if (will_write_gap_coordinate(last_record, lrec)) {
                    std::string gap_str = get_gap_str(last_record, lrec);
                    gwriter << gap_str << "\n";
                }
                
                bool break_status = will_break(first_record, last_record, \
                    lrec, coll_str);
                if (!break_status) {
                    last_record = lrec;
                    local_vec.push_back(lrec);
                } else {
                    // Write bed information for the umi chain
                    // Get the corresponding bed information
                    std::string bed_str = get_bed_str(local_vec);
                    bwriter.write_record_str(bed_str);

                    int rand_pos = get_rand_pos(local_vec.size());
                    bam_record final_rec = local_vec.at(rand_pos); 
                    writer.write_record(final_rec.full_rec);
                    write_collapse(local_vec, outfile_log, coll_len, rand_pos);
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

    // Check if local_vec has something; if yes push the last read
    if (!local_vec.empty()) {
        
        // Get the corresponding bed information
        std::string bed_str = get_bed_str(local_vec);
        bwriter.write_record_str(bed_str);

        // This works as the last break point
        int rand_pos = get_rand_pos(local_vec.size());
        bam_record final_rec = local_vec.at(rand_pos);
        writer.write_record(final_rec.full_rec);
        write_collapse(local_vec, outfile_log, coll_len, rand_pos);
        local_vec.clear();

    }
    std::cout << "mapped_count: " << lcount << "\n";
}

// This would take the first and the last bam record from local_vec. We
// assume that they represent the two ends of a isolated UMI chain.

void uminorm::throw_ineq_exception(std::string first_str, std::string sec_str) {

    if (first_str != sec_str) {
        std::string throw_msg = "First and second strings are not equal. \
            First string: " + first_str + ", sec_str: " + sec_str;
        throw std::runtime_error(throw_msg);
    }
}

void uminorm::throw_neg_execption(long lvar) {

    if (lvar < 0) {
        std::string throw_msg = "Unexpected negative value. Value: " + \
            std::to_string(lvar);
        throw std::runtime_error(throw_msg);
    }
}

std::string uminorm::get_bed_str(std::vector<bam_record> local_vec) {
    bam_record& first_rec = local_vec.front();
    bam_record& last_rec = local_vec.back();

    char first_strand = first_rec.strand;
    char last_strand = last_rec.strand;
    std::string first_strand_s(1, first_strand);
    std::string last_strand_s(1, last_strand);
    throw_ineq_exception(first_strand_s, last_strand_s);

    std::string first_umi_s(first_rec.umi);
    std::string last_umi_s(last_rec.umi);
    throw_ineq_exception(first_umi_s, last_umi_s);

    unsigned long startPos = first_rec.start_pos;
    unsigned long endPos = last_rec.end_pos;
    long totalGap = endPos - startPos + 1;
    throw_neg_execption(totalGap);
    std::string startPos_s = std::to_string(startPos);
    std::string endPos_s = std::to_string(endPos);

    // We need strand, stratPos, endPos and umi string to create a bed record.
    int first_ref_name_id = first_rec.ref_name_id;
    int last_ref_name_id = last_rec.ref_name_id;

    std::string first_ref_name_s = std::to_string(first_ref_name_id);
    std::string last_ref_name_s = std::to_string(last_ref_name_id);

    throw_ineq_exception(first_ref_name_s, last_ref_name_s);   

    // We shall have to get the tid of one of the two alignments
    const char* refarr = sam_hdr_tid2name(lhdr, first_ref_name_id);
    std::string lref_s(refarr); 
    
    std::string lname = first_umi_s + "_" + startPos_s + "_" + endPos_s +\
        "_" + std::to_string(totalGap);
    int lscore = 0;
    std::string lscore_s = std::to_string(lscore);
    std::string bed_str = lref_s + "\t" + 
        startPos_s + "\t" + 
        endPos_s + "\t" + 
        lname + "\t" + 
        lscore_s + "\t" +
        first_strand;
        
    return (bed_str);
}

std::string uminorm::get_gap_str(const bam_record& first_rec, 
    const bam_record& last_rec) {

    // Check that refname, strand and UMI for next_rec and first_rec are 
    char first_strand = first_rec.strand;
    char last_strand = last_rec.strand;
    std::string first_strand_s(1, first_strand);
    std::string last_strand_s(1, last_strand);
    throw_ineq_exception(first_strand_s, last_strand_s);

    std::string first_umi_s(first_rec.umi);
    std::string last_umi_s(last_rec.umi);
    throw_ineq_exception(first_umi_s, last_umi_s);

    int first_ref_name_id = first_rec.ref_name_id;
    int last_ref_name_id = last_rec.ref_name_id;

    std::string first_ref_name_s = std::to_string(first_ref_name_id);
    std::string last_ref_name_s = std::to_string(last_ref_name_id);

    throw_ineq_exception(first_ref_name_s, last_ref_name_s);   

    // We shall have to get the tid of one of the two alignments
    const char* refarr = sam_hdr_tid2name(lhdr, first_ref_name_id);
    std::string lref_s(refarr); 
    
    unsigned long first_start_pos = first_rec.start_pos;
    unsigned long last_start_pos = last_rec.start_pos;
    std::string first_start_pos_s = std::to_string(first_start_pos);
    std::string last_start_pos_s = std::to_string(last_start_pos);

    long gap_two_recs = last_start_pos - first_start_pos;
    throw_neg_execption(gap_two_recs);
    std::string gap_two_recs_str = std::to_string(gap_two_recs);

    std::string gap_str = first_umi_s + "\t" +
        first_strand_s + "\t" +
        gap_two_recs_str + "\t" +
        first_ref_name_s + "\t" +
        first_start_pos_s + "\t" +
        last_start_pos_s;

    return (gap_str);
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
        ("size_lim_M,s", po::value(&size_lim_M)->default_value(200),
            "Size of memory in megabyte")
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

    std::cout << "size_lim_M is set to " << std::to_string(size_lim_M) << "\n";
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
