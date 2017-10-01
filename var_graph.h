#include <cstdio>
#include <unordered_map>
#include <stdint.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "kmer_graph.h"

#define qK 8
#define K (4*qK) //multiple of 4

uint32_t nhash[4] = { 371158683,398345360,160141404,125238179 };

struct acgt4;
typedef struct acgt4 acgt4_t;
struct kmer_key;
typedef struct kmer_key kmer_key_t;
struct kmer_value;
typedef struct kmer_value kmer_value_t;
struct var_value;
typedef struct var_value var_value_t;

const char i2c[4] = { 'A','C','G','T' };
const uint8_t c2i[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15, 0,15, 1, 15,15,15, 2, 15,15,15,15, 15,15,15,15,
    15,15,15,15,  3,15,15,15, 15,15,15,15, 15,15,15,15,
    15, 0,15, 1, 15,15,15, 2, 15,15,15,15, 15,15,15,15,
    15,15,15,15,  3,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};  // c2i['A'] == 0, 'C' = 1, 'G' = 2, 'T' = 3

struct acgt4 {
    uint8_t n0 : 2;
    uint8_t n1 : 2;
    uint8_t n2 : 2;
    uint8_t n3 : 2;
    inline uint8_t at(uint8_t i) const {
        return ((i > 1) ? (i == 3 ? n3 : n2) : (i == 0 ? n0 : n1));
    }
    inline void set(uint8_t i, uint8_t val) {
        if (i > 1) {
            if (i == 3) {
                n3 = val;
            }
            else {
                n2 = val;
            }
        }
        else {
            if (i == 0) {
                n0 = val;
            }
            else {
                n1 = val;
            }
        }
    }
};

struct kmer_key {
    acgt4_t acgt4s[qK];

    inline uint8_t acgt_at(uint32_t i) const { return acgt4s[i / 4].at(i % 4); }
    inline void acgt_set(uint32_t i, uint8_t val) { acgt4s[i / 4].set(i % 4, val); }

    std::string to_string() const {
        std::string s;
        for (int i = 0; i < K; ++i) {
            switch (acgt_at(i)) {
            case 0:
                s += 'A';
                break;
            case 1:
                s += 'C';
                break;
            case 2:
                s += 'G';
                break;
            case 3:
                s += 'T';
            }
        }
        return s;
    }

    bool operator==(const struct kmer_key &other) const {
        for (int32_t i = 0; i < K; ++i)
            if (acgt_at(i) != other.acgt_at(i)) return false;
        return true;
    }
};

struct kmer_value {
    uint32_t acgt_prev : 4;
    uint8_t flag : 2;
    uint32_t count : 26;
    uint32_t acgt_next : 4;
    uint32_t path_length : 28;
    kmer_value() {
        acgt_prev = 0000; //binary
        acgt_next = 0000; //binary
        count = 0;
        path_length = 0;
        flag = 0;
    }
};

struct var_value {

    uint32_t a_next_count : 16;
    uint32_t c_next_count : 16;
    uint32_t g_next_count : 16;
    uint32_t t_next_count : 16;

    var_value() {
        a_next_count = 0;
        c_next_count = 0;
        g_next_count = 0;
        t_next_count = 0;
    }
    //need to decide what to include for this
};

struct kmer_key_hasher {
    std::size_t operator()(const kmer_key& k) const {
        std::size_t seed = 0;
        for (int32_t i = 0; i < K; ++i) {
            seed ^= (nhash[k.acgt_at(i)] + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        }
        return seed;
    }
};



class varGraph {
private:
    // Map to store graph and pointer to reference
    std::unordered_map<kmer_key_t, var_value_t, kmer_key_hasher> vmap;
    kmerGraph* refGraph;

public:
    std::unordered_map<kmer_key_t, var_value_t, kmer_key_hasher>::iterator vmap_it_t;
    bool cflag = false;

    // Effects: Different Constructors for varGraph. First one still needs to be finished
    varGraph(const char* ref_graph_name, const char* var_graph_name);
    varGraph(const char* ref_graph_name);
    varGraph(kmerGraph* r);

    // Effects: Deletes refGraph if it was created by new
    ~varGraph();

    // Effects: Reads from fastq files. 
    // If kmer meets quality standard and isn't found in reference, inserts into hashtable
    void fq_read(std::string filename); 

    // Effects: Inserts key if its quality is sufficient
    void insert(std::string seq, std::string qual);

    //Effects: Removes kmers from graph which only appear once
    void error_prune();

    // Effects: Calculates quality of input string. Returns true if it meets quality threshold
    bool quality(std::string q);

    // Effects: Checks whether key is found in reference
    bool is_variation(std::string kmer);

    // Effects: Converts string to key. Can take pointer or string input
    kmer_key_t string_to_key(const char* String);
    kmer_key_t string_to_key(const std::string & String);

    // Requires: Key in hashtable
    // Effects: Outputs key
    void print_key(kmer_key_t Key);
    
    // Effects: Prints map in key value pairs. Value printed is count of kmer
    void print_map();

    // Effects: Returns Total number of nodes as well as how many are unique vs non-unique
    void nodes();

    // Effects: Stores vmap to output file
    bool store_to_file(std::string filename);
};