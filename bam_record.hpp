#ifndef _BAM_RECORD_HPP
#define _BAM_RECORD_HPP

#include <iostream>

class bam_record {

    public:
    bool is_mapped;
    std::string ref_name;
    std::string umi;
    char strand;
    unsigned long start_pos;
    unsigned long end_pos;
    std::string qname;
    std::string full_rec;
    unsigned int reader_index;


    unsigned int get_size() {
        unsigned int lsize = sizeof(bool) + sizeof(char) + 
            2 * sizeof(unsigned long) + ref_name.size() + umi.size() + 
            qname.size() + full_rec.size();
        return lsize;
    }
};

struct compare_bam_less {

    bool operator()(const bam_record& a, const bam_record& b)   {
        int ref_c = a.ref_name.compare(b.ref_name);
        if (ref_c < 0) {
            return true;
        } else if (0 == ref_c) {
            int umi_c = a.umi.compare(b.umi);
            if (umi_c < 0) {
                return true;
            } else if (umi_c == 0) {
                if (a.strand < b.strand) {
                    return true;
                } else if (a.strand == b.strand) {
                    if (a.start_pos < b.start_pos) {
                        return true;
                    } else if (a.start_pos == b.start_pos) {
                        int qname_c = a.qname.compare(b.qname);
                        if (qname_c < 0) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }
};

struct compare_bam_greater {

    bool operator()(const bam_record& a, const bam_record& b)   {
        int ref_c = a.ref_name.compare(b.ref_name);
        if (ref_c > 0) {
            return true;
        } else if (0 == ref_c) {
            int umi_c = a.umi.compare(b.umi);
            if (umi_c > 0) {
                return true;
            } else if (umi_c == 0) {
                if (a.strand > b.strand) {
                    return true;
                } else if (a.strand == b.strand) {
                    if (a.start_pos > b.start_pos) {
                        return true;
                    } else if (a.start_pos == b.start_pos) {
                        int qname_c = a.qname.compare(b.qname);
                        if (qname_c > 0) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }
};

#endif
