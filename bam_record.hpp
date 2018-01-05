#ifndef _BAM_RECORD_HPP
#define _BAM_RECORD_HPP

#include <iostream>

class bam_record {

    public:
    bool is_mapped;
    int ref_name_id = -1;
    char* umi;
    char strand;
    unsigned long start_pos;
    unsigned long end_pos;
    char* qname;
    char* full_rec;
    int reader_index = -1;

    unsigned int get_size() {
        unsigned int lsize = sizeof(bool) + sizeof(char) + 
            2 * sizeof(unsigned long) + 2 * sizeof (int) + 
            strlen(umi) + strlen(qname) + strlen(full_rec) + 3;
        return lsize;
    }

    ~bam_record() {
        if (umi != nullptr) {
            delete[] umi;
        }
        if (qname != nullptr) {
            delete[] qname;
        }
        if (full_rec != nullptr) {
            delete[] full_rec;
        }
    }

    bam_record() {
        umi = nullptr;
        qname = nullptr;
        full_rec = nullptr;
    }

    // Constructor
    bam_record (bool _is_mapped, int _ref_name_id, char* _umi, char _strand,
            unsigned long _start_pos, unsigned long _end_pos, char* _qname, 
            char* _full_rec, int _reader_index) {
        is_mapped = _is_mapped;
        ref_name_id = _ref_name_id;
        strand = _strand;
        start_pos = _start_pos;
        end_pos = _end_pos;
        reader_index = _reader_index;
 
        size_t umi_len  = strlen(_umi) + 1;
        umi = new char[umi_len];
        memcpy(umi, _umi, umi_len);
 
        size_t qname_len = strlen(_qname) + 1;
        qname = new char[qname_len];
        memcpy(qname, _qname, qname_len);

        size_t full_rec_len = strlen(_full_rec) + 1;
        full_rec = new char[full_rec_len];
        memcpy(full_rec, _full_rec, full_rec_len);

        //std::cout << std::to_string(ref_name_id) << " " << std::string(umi) << " " << 
        //    std::string(1, strand) << " " << std::to_string(start_pos) << " " << 
        //   std::string(qname) << " " << std::string(full_rec) << "\n";

    }

    // Copy constructor
    bam_record(const bam_record& that) {
        is_mapped = that.is_mapped;
        ref_name_id = that.ref_name_id;
        strand = that.strand;
        start_pos = that.start_pos;
        end_pos = that.end_pos;
        reader_index = that.reader_index;
 
        size_t umi_len  = strlen(that.umi) + 1;
        umi = new char[umi_len];
        memcpy(umi, that.umi, umi_len);
 
        size_t qname_len = strlen(that.qname) + 1;
        qname = new char[qname_len];
        memcpy(qname, that.qname, qname_len);

        size_t full_rec_len = strlen(that.full_rec) + 1;
        full_rec = new char[full_rec_len];
        memcpy(full_rec, that.full_rec, full_rec_len);        
    }

    // Copy assignment operator
    bam_record& operator=(const bam_record& that) {
        bam_record temp(that);
        swap(*this, temp);
        return *this;
    }

    friend void swap(bam_record& a, bam_record& b) {
        std::swap(a.is_mapped, b.is_mapped);
        std::swap(a.ref_name_id, b.ref_name_id);
        std::swap(a.strand, b.strand);
        std::swap(a.start_pos, b.start_pos);
        std::swap(a.end_pos, b.end_pos);
        std::swap(a.reader_index, b.reader_index);
        std::swap(a.qname, b.qname);
        std::swap(a.umi, b.umi);
        std::swap(a.full_rec, b.full_rec);
    }

    // Move constructor
    bam_record(bam_record&& that)   
    {
        // First copy the primitive types
        is_mapped = that.is_mapped;
        ref_name_id = that.ref_name_id;
        strand = that.strand;
        start_pos = that.start_pos;
        end_pos = that.end_pos;
        reader_index = that.reader_index;
        umi = that.umi;
        that.umi = nullptr;
        qname = that.qname;
        that.qname = nullptr;
        full_rec = that.full_rec;
        that.full_rec = nullptr;
    }

    // Move assignment operator
    bam_record& operator=(bam_record&& that)   
    {
        // First copy the primitive types
        is_mapped = that.is_mapped;
        ref_name_id = that.ref_name_id;
        strand = that.strand;
        start_pos = that.start_pos;
        end_pos = that.end_pos;
        reader_index = that.reader_index;
        umi = that.umi;
        that.umi = nullptr;
        qname = that.qname;
        that.qname = nullptr;
        full_rec = that.full_rec;
        that.full_rec = nullptr;
        return *this;
    }

};

struct compare_bam_less {

    bool operator()(const bam_record& a, const bam_record& b)   {
        if (a.ref_name_id < b.ref_name_id) {
            return true;
        } else if (a.ref_name_id == b.ref_name_id) {
            int umi_c = strcmp(a.umi, b.umi);
            if (umi_c < 0) {
                return true;
            } else if (umi_c == 0) {
                if (a.strand < b.strand) {
                    return true;
                } else if (a.strand == b.strand) {
                    if (a.start_pos < b.start_pos) {
                        return true;
                    } else if (a.start_pos == b.start_pos) {
                        int qname_c = strcmp(a.qname, b.qname);
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
        if (a.ref_name_id > b.ref_name_id) {
            return true;
        } else if (a.ref_name_id == b.ref_name_id) {
            // In reality one might have to do only this
            // comparison, so the string comparison between
            // a.qname and b.qname may not be required
            // except extreme cases.
            int umi_c = strcmp(a.umi, b.umi);
            if (umi_c > 0) {
                return true;
            } else if (umi_c == 0) {
                if (a.strand > b.strand) {
                    return true;
                } else if (a.strand == b.strand) {
                    if (a.start_pos > b.start_pos) {
                        return true;
                    } else if (a.start_pos == b.start_pos) {
                        int qname_c = strcmp(a.qname, b.qname);
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
