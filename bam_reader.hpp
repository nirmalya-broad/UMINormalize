#ifndef _BAM_READER_HPP
#define _BAM_READER_HPP

#include <cstring>
#include <cmath>
#include <regex>
#include <htslib/sam.h>
#include "bam_record.hpp"

class bam_reader {
    public:

    bam_reader() = default;

    bam_reader(std::string& infile_str) {
        const char* format = NULL;
        if (has_suffix(infile_str, "sam")) {
            format = "r";
        } else if (has_suffix(infile_str, "bam")) {
            format = "rb";
        } else {
            std::string lstr = "File with illegal suffix: " + infile_str + "\n";
            throw std::runtime_error(lstr);
        }

        const char* infile_cstr = infile_str.c_str();
        if (!(fp = sam_open(infile_cstr, format))) {
            is_good = false;
            throw std::runtime_error("Error in sam_open");
        }

        lhdr = sam_hdr_read(fp);
        lread = bam_init1();

    }

    bam_hdr_t* get_sam_header() {
        return lhdr;
    }

    char* get_refname(bam1_t* lread) {
        int32_t core_tid = (lread -> core).tid;
        char* refname = NULL;
        if (core_tid == -1) {
            //std::string lqname = std::string(bam_get_qname(lread));
            //std::cout << "Negative core_tid with qname: " << lqname << "\n";
            std::string retval = "*";
            refname = (char*)retval.c_str();
            return refname;
           
        } else {
            refname = (lhdr -> target_name)[core_tid];
        }
        if (NULL == refname) {
            throw std::runtime_error("NULL refname.");
        }
        return refname;
    }


    std::string read_record(bam_record& bam_rec) {
        int ret_val = -1;
        if ((ret_val = sam_read1(fp, lhdr, lread)) >= 0) {
            std::string next_record = bam_format1(lhdr, lread);
            // Get the other things
            uint16_t flag = (lread -> core).flag;
            bool is_mapped = flag & BAM_FMUNMAP;
            // Process the read and return
            // Get the query_name
            char* qname = bam_get_qname(lread);
            std::string qname_str(qname);
            // Get the ref_name
            char* ref_name = get_refname(lread);
            std::string ref_name_str(ref_name);

            // Get umi_str
            // Get start_pos
            // Get umi string; basically extract the umi_XXXXXX part from the read
            std::regex expr("^\\S+?umi_(\\w+).*$");
            std::smatch what;
            std::string umi_str;
            bool res = std::regex_search(qname_str, what, expr);

            if (res) {
                umi_str = what[1];
            } else {
                std::string err_str = "umi str not found, qname: " + qname_str;
                throw std::runtime_error(err_str);
            }


            // Get start_pos
            unsigned long start_pos = (lread -> core).pos;

            // Get end_pos
            unsigned int q_len = (lread -> core).l_qseq;
            unsigned long end_pos = start_pos + q_len -1;

            // Get strand 
            char lstrand = '.';
            if (bam_is_rev(lread)) {
                lstrand = '-';
            } else {
                lstrand = '+';
            }
            
            // Now fill all the fields of bam_record
            bam_rec.ref_name = ref_name_str;
            bam_rec.qname = qname_str;
            bam_rec.start_pos = start_pos;
            bam_rec.end_pos = end_pos;
            bam_rec.strand = lstrand;
            bam_rec.umi = umi_str;
            bam_rec.full_rec = next_record;    
            bam_rec.is_mapped = is_mapped;
             
            return next_record;
        }
        return std::string();
    }

    bool has_suffix(const std::string &str, const std::string &suf)
    {
        return str.size() >= suf.size() &&
           str.compare(str.size() - suf.size(), suf.size(), suf) == 0;
    }


    // Copied from bam.c of samtools
    std::string bam_format1(const bam_hdr_t *header, const bam1_t *b)
    {
        kstring_t str;
        str.l = str.m = 0; str.s = NULL;
        if (sam_format1(header, b, &str) < 0) {
            free(str.s);
            str.s = NULL;
            return std::string();
        }
        std::string res_str(str.s);
        free(str.s);
        return res_str;
    }

    ~bam_reader() {
        bam_destroy1(lread);
        bam_hdr_destroy(lhdr);   
        sam_close(fp); // clean up 
    }

    bool good() {
        return is_good; 
    }

    private:
    std::string infile_str;    
    htsFile *fp = NULL;
    bam_hdr_t *lhdr = NULL;
    bam1_t* lread = NULL;
    bool is_good = true;


};
#endif


