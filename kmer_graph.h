#include <cstdio>
#include <unordered_map>
#include <stdint.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#define qK 8
#define K (4*qK) //multiple of 4

// Hash function 
uint32_t nhash[4] = { 371158683,398345360,160141404,125238179 };

// Map entry types
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

struct kmer_key_hasher {
    std::size_t operator()(const kmer_key& k) const {
        std::size_t seed = 0;
        for (int32_t i = 0; i < K; ++i) {
            seed ^= (nhash[k.acgt_at(i)] + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        }
        return seed;
    }
};

class kmerGraph {
private:
    // Map to store graph
    std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher> kmap;

public:
    typedef std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher>::iterator kmap_it_t;

    //Modifies: kmap
    //Effects: clears kmap
    void clear_map();

    // Requires: String of pointer to start of string
    // Effects: Inserts key value pair into hash table
    int32_t insert(const char* s);
    int32_t insert(const std::string& s);

    // Effects: Returns true if key is found in hash table
    bool has_key(const std::string s);
    bool has_key(kmer_key_t k);


    // Effects: Returns iterator to key if found or end of map
    kmap_it_t find(const std::string s);
    kmap_it_t find(kmer_key_t k);

    // Requires: Key present in hashtable
    // Effects: Returns value associated with a given key
    kmer_value_t value_at(const std::string s);
    kmer_value_t value_at(kmer_key_t k);

    // Effects: Prints key at iterator position
    void print_key(kmap_it_t Position);

    // Effects: Prints key to standard output
    const void print_key(kmer_key_t Key);

    // Effects: Prints out last base of key at iterator position
    void print_last_k(kmap_it_t Position);

    // Requires: Key present in hashtable
    // Effects: Outputs hashtable in key value pairs. 
    // Can output count, prev base, or next base for value. 
    void print_map();

    // Effects: Returns complemantary base to input base
    char comp_base(char c);

    // Effects: Returns complemantary sequence to input sequence
    std::string string_complement(std::string s);

    // Modifies: cin and cout
    // Effects: Returns string read from user inputed file. Can cout sequence if prompted.
    std::string string_reader();

    // Modifies: Output file
    // Effects: Creates file with user given name or writes to existing file. 
    // Writes in hashtable converted to binary. 
    bool store_to_file(std::string filename);

    // Modifies: Cin, kmap
    // Effects: Loads graph from user inputed file into hashtable
    bool load_from_file();

    // Effects: Returns kmer with the largest count
    void largest_count();

    // Effects: Returns Total number of nodes as well as how many are unique vs non-unique
    void nodes();

    //Longest Path

    // Effects: Checks if current node has exactly one child
    bool one_child_next(kmap_it_t CurrentNode);

    // Effects: Returns length of path ending in current node
    int32_t path_length(kmap_it_t &CurrentNode);

    // Effects: Returns next kmer in path
    kmer_key_t next_kmer(kmap_it_t CurrentNode);
    
    // Effects: Changes path length
    void change_path_length(kmap_it_t CurrentNode, uint32_t NewLength);

    // Effects: Finds the path from current node to the end
    uint32_t path_length_calculation(kmap_it_t CurrentNode);

    uint8_t loop_solver(kmap_it_t CurrentNode);

    // Effects: Finds the longest path in map
    int32_t longest_unique_path();

    // Effects: Prints out longest path
    void printing_out_paths_kmer();
    void printing_out_paths_string();

    // Effects: Converts input string into a key type
    const kmer_key_t string_to_key(const char* String);
    const kmer_key_t string_to_key(const std::string & String);

    // Effects: Prints out length of path starting from input kmer
    void printing_NodePath_string(std::string InputKmer);

};