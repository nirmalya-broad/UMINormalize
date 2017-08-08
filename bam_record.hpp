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

    bool operator < (const bam_record &b) const {
        int ref_c = ref_name.compare(b.ref_name);
        if (ref_c < 0) {
            return true;
        } else if (0 == ref_c) {
            int umi_c = umi.compare(b.umi);
            if (umi_c < 0) {
                return true;
            } else if (umi_c == 0) {
                if (strand < b.strand) {
                    return true;
                } else if (strand == b.strand) {
                    if (start_pos < b.start_pos) {
                        return true;
                    } else if (start_pos == b.start_pos) {
                        int qname_c = qname.compare(b.qname);
                        if (qname_c < 0) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    unsigned int get_size() {
        unsigned int lsize = sizeof(bool) + sizeof(char) + 
            2 * sizeof(unsigned long) + ref_name.size() + umi.size() + 
            qname.size() + full_rec.size();
        return lsize;
    }
};

#endif
