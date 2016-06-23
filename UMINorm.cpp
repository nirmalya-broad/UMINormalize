#include <stxxl/sorter>

#include <limits>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <boost/program_options.hpp>

#include <regex>

namespace po = boost::program_options;


/* Parses from the simple input that Jonathan provides */

struct jp_record {

	long long start_c;
	long long end_c;
	char strand;
	char is_paired;
	std::string umi_str;
	std::string rname;


	jp_record(long long _start_c, long long _end_c, char _strand, 
		char _is_paired, std::string _umi_str, std::string _rname):
		start_c(_start_c), 
		end_c(_end_c), 
		strand(_strand), 
		is_paired(_is_paired), 
		umi_str(_umi_str),
		rname(_rname)
	{}

	jp_record(const jp_record& jpr):
		start_c(jpr.start_c),
        end_c(jpr.end_c),
        strand(jpr.strand),
        is_paired(jpr.is_paired),
        umi_str(jpr.umi_str),
        rname(jpr.rname)
	{}

	jp_record() {}

};

struct jc_Comparator {

	bool operator () (const jp_record& a, const jp_record& b) const {
		if (a.start_c < b.start_c) {
        	return true;
       }

		return false;
	}

	jp_record min_value() const {
        jp_record dummy;
        dummy.start_c = std::numeric_limits<long long>::min();

        return dummy;
    }

	jp_record max_value() const {

        jp_record dummy;
        dummy.start_c = std::numeric_limits<long long>::max();

        return dummy;

    }


	
};


struct jpr_Comparator {
	bool operator () (const jp_record& a, const jp_record& b) const {
		int umi_c = a.umi_str.compare(b.umi_str);
		if (umi_c < 0) {
			return true;
		} else if (umi_c == 0) {
			if (a.strand < b.strand) {
				return true;
			} else if (a.strand == b.strand) {
				if (a.start_c < b.start_c) {
					return true;
				} else if (a.start_c == b.start_c) {
					if (a.rname.compare(b.rname) < 0) {
						return true;
					}
				}
			}
		}
		return false;
	}

	jp_record min_value() const {
		jp_record dummy;
		dummy.start_c = std::numeric_limits<long long>::min();
		dummy.strand = '+';
		dummy.umi_str = "AAAAAA";

		return dummy;
	}

	jp_record max_value() const {

		jp_record dummy;
		dummy.start_c = std::numeric_limits<long long>::max();
		dummy.strand = '-';
		dummy.umi_str = "TTTTTT";

		return dummy;

	}
};


class cl_args {
	
	public:
	void print_help(); 
	bool parse_args(int argc, char* argv[]);
	po::options_description desc;

	int lmem_size;
    int lgap_cluster;
    std::string infile_str;
    std::string outfile_str;
    std::string cluster_file_str;
	std::string sorted_file_str;
	
};

class UmiNorm {
	public:
	int parse_jpr_file();
	void sort_jp_records();
	void print_jp_records();
	UmiNorm(cl_args l_args);
	
	stxxl::sorter<jp_record, jpr_Comparator, 1 * 1024 * 1024> jpr_sorter_vec;
	std::map<std::string, int> umiCountMap;
	std::vector<int> gap_vec;


	private:
	int lmem_size;
	int lgap_cluster;
	std::string infile_str;
	std::string outfile_str;
	std::string cluster_file_str;
	std::string sorted_file_str;

};

UmiNorm::UmiNorm(cl_args l_args):
	lmem_size(l_args.lmem_size),
	lgap_cluster(l_args.lgap_cluster),
	infile_str(l_args.infile_str),
	outfile_str(l_args.outfile_str),
	cluster_file_str(l_args.cluster_file_str),
	sorted_file_str(l_args.sorted_file_str),
	jpr_sorter_vec(jpr_Comparator(), l_args.lmem_size * 1024 * 1024)
{}		



void cl_args::print_help() {
    std::cout << desc << "\n";
	std::cout << "Usage: uminorm -i <infile>"
			" -o <outfile> -a <all_out>\n\n";
}

bool cl_args::parse_args(int argc, char* argv[]) {
	 bool all_set = true;
    desc.add_options()
        ("help,h", "produce help mesage")
        ("infile,i", po::value<std::string>(&infile_str), "Input file")
        ("outfile,o", po::value<std::string>(&cluster_file_str), "Contains normalized read clusters")
		("all_out,a", po::value<std::string>(&outfile_str), "Contains detailed sorted reads")
		("sorted_out,s", po::value<std::string>(&sorted_file_str), "Contains sorted clusters")
		("memory_size,m", po::value<int>(&lmem_size)->default_value(256), "Optional/Memory (RAM) in use in MB")
		("gap_cluster,g", po::value<int>(&lgap_cluster)->default_value(500), "Optional/minium gap between two reads in an UMI cluster")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
		//print_help();
        return 0;
    } else {
		//all_set =false;
	}

    if (!vm.count("infile")) {
        std::cout << "Infile not set.\n";
        all_set = false;
		
    } else {
		
        std::cout << "Infile is set to: " << infile_str << ".\n";
    }

    if (!vm.count("outfile")) {
        std::cout << "outfile not set.\n";
        all_set = false;
    } else {
        std::cout << "outfile is set to: " << cluster_file_str << ".\n";
    }

	if (!vm.count("all_out")) {
        std::cout << "all_out not set.\n";
        all_set = false;
	} else {
		std::cout << "all_out is set to: " << outfile_str << ".\n";
	}

	if (!vm.count("sorted_out")) {
        std::cout << "sorted_out not set.\n";
        all_set = false;
    } else {
        std::cout << "sorted_out is set to: " << sorted_file_str << ".\n";
    }


	return all_set;
}


int UmiNorm::parse_jpr_file() {
	std::regex expr("(\\d+)\\.\\s+(\\d+)\\.\\s+(\\S)\\s+(\\w)\\s+(\\S+)\\:(\\w{6})");
	std::smatch what;
	
	std::string lstr;
	std::ifstream infile(infile_str.c_str());

	while (std::getline(infile, lstr)) {
	
		bool res = std::regex_search(lstr, what, expr);
		//std::cout << lstr << "\n";

		if (res) {
			std::string start_str = what[1];
			std::string end_str = what[2];
			std::string strand_str = what[3];
			std::string is_paired_str = what[4];
			std::string rname = what[5];
			std::string umi_str = what[6];
		
			long long lstart = std::stoll(start_str);
			long long lend = std::stoll(end_str);
			char lstrand = '.';
			if (0 == strand_str.compare("+")) {
				lstrand = '+';
			} else if (0 == strand_str.compare("-")) {
				lstrand = '-';
			}
			
			// Junk
			char is_paired = 'J';

			if (0 == is_paired_str.compare("Y") || 0 == is_paired_str.compare("y")) {
				is_paired = 'Y';
			} else if (0 == is_paired_str.compare("N") || 0 == is_paired_str.compare("n")){
				is_paired = 'N';
			}

			//std::cout << lstart << " " << lend << " " << lstrand << " " << is_paired 
			//	<< " " << umi_str << "\n";

			jp_record ljp_record = {
				lstart,
				lend,
				lstrand,
				is_paired,
				umi_str,
				rname
			};
			
			jpr_sorter_vec.push(ljp_record);
			umiCountMap[umi_str]++;
		}
	}

	return 0;
}

void UmiNorm::sort_jp_records() {
	jpr_sorter_vec.sort();
}


void UmiNorm::print_jp_records() {

	std::ofstream outfile(outfile_str);
	std::ofstream cluster_file(cluster_file_str);
	std::ofstream sorted_file(sorted_file_str);

	std::string last_umi = "";
	char last_strand = '.';
	long long last_start = 0;

	jp_record first_record;
	jp_record last_record;

	
	/* willdump is set to true when there is an end of the existing 
	 * and start of a new cluster. This would happen when the
     * stream of reads change UMI, or strand or enter to a new
     * segment.
	*/
	bool willDump = false;
	std::string last_delim = "";
	int mem_count = 0;


	stxxl::sorter<jp_record, jc_Comparator, 1 * 1024 * 1024> cluster_vec(jc_Comparator(), (lmem_size/2) * 1024 * 1024);

	while (!jpr_sorter_vec.empty()) {
		if ((*jpr_sorter_vec).umi_str.compare(last_umi) != 0) {
			if (last_umi.compare("") != 0) {
				outfile << "----------------------------------------\n";
				last_delim =  "----------------------------------------\n";
				willDump = true;
			} else {
				first_record = *jpr_sorter_vec;
			}
		} else if ((*jpr_sorter_vec).strand != last_strand) {
			outfile << "............\n";
			last_delim = "............\n";
			willDump = true;
		}
			
		
		if ((*jpr_sorter_vec).strand == last_strand 
			&& (*jpr_sorter_vec).umi_str.compare(last_umi) == 0) {
			int lgap = (*jpr_sorter_vec).start_c - last_start;
			gap_vec.push_back(lgap);
			if (lgap > 500) {
				//outfile  << ">>>>> Lgap: " << lgap << "\n";
				outfile  << ">>>>>\n";
				last_delim = ">>>>>\n";
				willDump = true;
			}
		}
		
		//std::cout << (*jpr_sorter_vec).rname << " " << willDump << "\n";

		if (willDump) {
			std::string combined_name;
			if (first_record.rname.compare(last_record.rname) != 0) {
				combined_name = first_record.rname +  
				 "--" + last_record.rname + "--" + std::to_string(mem_count);
			} else {
				combined_name = first_record.rname + "--" + std::to_string(mem_count);
			}
			cluster_file << first_record.umi_str << " "
				<< first_record.strand << " "
				<< first_record.start_c << " "
				<< last_record.end_c << " "
				<< combined_name << "\n";


			jp_record lc_record = {
                first_record.start_c,
                last_record.end_c,
                first_record.strand,
                first_record.is_paired,
                first_record.umi_str,
                combined_name
            };
			cluster_vec.push(lc_record);

			cluster_file << last_delim;
			first_record = *jpr_sorter_vec;
			mem_count =0;
			willDump = false;
			
		}

		last_umi = (*jpr_sorter_vec).umi_str;
		last_strand = (*jpr_sorter_vec).strand;
		last_start = (*jpr_sorter_vec).start_c;

        outfile << (*jpr_sorter_vec).umi_str << " " 
			<< (*jpr_sorter_vec).strand << " "
			<< (*jpr_sorter_vec).start_c 
			<< " " << (*jpr_sorter_vec).end_c 
			<< " " << (*jpr_sorter_vec).rname << "\n";

		last_record = *jpr_sorter_vec;

        ++jpr_sorter_vec;
		mem_count++;

    }
	
	const std::string combined_name = first_record.rname
                + "--" + last_record.rname;
        cluster_file << first_record.umi_str << " "
        << first_record.strand << " "
        << first_record.start_c << " "
        << last_record.end_c << " "
        << combined_name << "\n";

	jp_record lc_record = {
        first_record.start_c,
       last_record.end_c,
       first_record.strand,
       first_record.is_paired,
       first_record.umi_str,
       combined_name
    };
	cluster_vec.push(lc_record);
	cluster_vec.sort();
	
	while (!cluster_vec.empty()) {
		sorted_file << (*cluster_vec).umi_str << " "
            << (*cluster_vec).strand << " "
            << (*cluster_vec).start_c
            << " " << (*cluster_vec).end_c
            << " " << (*cluster_vec).rname << "\n";
		++cluster_vec;
	}


		
	outfile.close();
	cluster_file.close();
	sorted_file.close();

}


int main(int argc, char* argv[]) {



	bool all_set = true;
	cl_args l_args;

	try {
		all_set = l_args.parse_args(argc, argv);	
	} catch(std::exception& e) {
		std::cerr << "error: " << e.what() << "\n";
        //lbs.print_help();
        return 1;

	} catch (...) {
		return 0;
	}

 	if (!all_set) {
		l_args.print_help();
		return 0;
	}

	UmiNorm umi_obj(l_args);

	umi_obj.parse_jpr_file();
	umi_obj.sort_jp_records();
	umi_obj.print_jp_records();

	return 0;

}

