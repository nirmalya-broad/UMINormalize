#ifndef _BAM_READER_HPP
#define _BAM_READER_HPP

#include <cstring>
#include <cmath>
#include <regex>
#include <htslib/sam.h>
#include <regex.h>
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

   char* regex_match_cstr(const char* regex_str, char* full_rec,
        int match_len) {
        regex_t regex;
        int nmatch = 2;
        regmatch_t pmatch[2];

        int reti = regcomp(&regex, regex_str, REG_EXTENDED);
        if (reti) {
            fprintf(stderr, "Could not compile regex\n");
            return NULL;
        } 


        reti = regexec(&regex, full_rec, 2, pmatch, 0);
        if (!reti) {
            char* result = new char[match_len + 1];
            int len = pmatch[1].rm_eo - pmatch[1].rm_so;
            memcpy(result, full_rec + pmatch[1].rm_so, len);
            result[len] = 0;
            regfree(&regex);
            return result;
        }
        return NULL;
    }
 

    std::string read_record(bam_record& bam_rec) {
        int ret_val = -1;
        // Return value of sam_read1:
        // 0 if successful; otherwise negative
        if ((ret_val = sam_read1(fp, lhdr, lread)) >= 0) {
            // Get the other things
            char* next_record = bam_format1(lhdr, lread);
            uint16_t flag = (lread -> core).flag;
            bool is_mapped = flag & BAM_FMUNMAP;
            // Process the read and return
            // Get the query_name
            char* qname_src = bam_get_qname(lread);
            int qname_len = strlen(qname_src) + 1;
            char* qname = new char[qname_len];
            memcpy(qname, qname_src, qname_len);
            // Get the ref_name

            // Get umi_str
            // Get start_pos
            // Get umi string; extract the umi_XXXXXX part from the read
            // Note the length of UMI is 6 base pair, which is indicated here
            const char* regex_str = "^\\S+?umi_(\\w{6}).*$";
            int match_len = 6;
            char* umi_str = regex_match_cstr(regex_str, qname, match_len);
            if (!umi_str) {
                std::string err_str = "umi str not found, qname: " + std::string(qname);
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
            
            int reader_index = -1;
            int ref_name_id = (lread -> core).tid;
            std::string next_rec_str(next_record);
            bam_record bam_rec_temp(is_mapped, ref_name_id, umi_str, lstrand,
                start_pos, end_pos, qname, next_record, reader_index); 
            bam_rec = std::move(bam_rec_temp);

            delete[] umi_str;
            free(next_record);
            delete[] qname;
             
            return next_rec_str;
        }
        return std::string();
    }

    bool has_suffix(const std::string &str, const std::string &suf)
    {
        return str.size() >= suf.size() &&
           str.compare(str.size() - suf.size(), suf.size(), suf) == 0;
    }


    // Copied from bam.c of samtools
    char* bam_format1(const bam_hdr_t *header, const bam1_t *b) {
        kstring_t str;
        str.l = str.m = 0; str.s = NULL;
        if (sam_format1(header, b, &str) < 0) {
            free(str.s);
            str.s = NULL;
        }
        return str.s;
    }

    ~bam_reader() {
        bam_destroy1(lread);
        bam_hdr_destroy(lhdr);   
        sam_close(fp); // clean up 
    }


    private:
    std::string infile_str;    
    htsFile *fp = NULL;
    bam_hdr_t *lhdr = NULL;
    bam1_t* lread = NULL;


};
#endif


