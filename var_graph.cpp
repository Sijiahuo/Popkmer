#include <cstdio>
#include <unordered_map>
#include <stdint.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "kmer_graph.h"
#include "var_graph.h"

#define qK 8
#define K (4*qK) //multiple of 4

    varGraph::varGraph(const char* ref_graph_name, const char* var_graph_name) {
        refGraph = new kmerGraph();
        cflag = true;
        //refGraph->loadFromFile;
        //loadFromFile(var_graph_name); //these two lines give me errors so i just commented them out.
    }


    varGraph::varGraph(const char* ref_graph_name) {
        refGraph = new kmerGraph();
        cflag = true;
    }

    varGraph::varGraph(kmerGraph* r) {
        refGraph = r;
        cflag = false;
        // kmerGraph g;
        // g.loadFromFile
    }

    varGraph::~varGraph() {
        // if it was created with new (probably need a flag)
        if (cflag) { 
            delete refGraph; 
        }
    }

    void varGraph::fq_read(std::string filename) {
        std::ifstream ins;
        ins.open(filename);

        // if invalid file is entered
        while (!ins.is_open()) {
            std::cin.clear();
            std::string errorName;
            getline(std::cin, errorName);
            std::cout << "Invalid filename" << std::endl << "Enter filename: ";
            std::cin >> filename;
            ins.open(filename);
        }

        std::string identifier;
        std::string sequence;
        std::string plus;
        std::string quality;
        std::string kmer = "";
        std::string qual = "";
        //const int SIZE = 150;

        while (ins >> identifier >> sequence >> plus >> quality) {
            insert(sequence, quality);
        }
    }

    void varGraph::insert(std::string seq, std::string qual) {
        uint64_t SIZE = seq.length();
        int i;
        int j;
        std::string kmer = "";
        std::string kmer_qual = "";
        for (i = 0; i < SIZE - K + 1; i++) {
            kmer = "";
            kmer_qual = "";
            for (j = i; j < i + K; j++) {
                kmer += seq[j];
                kmer_qual += qual[j];
            }
            if (is_variation(kmer) && quality(kmer_qual)) {  //&& quality(kmer_qual)
                var_value& var_v = vmap[string_to_key(kmer)];
                if (c2i[seq[j]] == 0) {
                    ++var_v.a_next_count;
                }
                else if (c2i[seq[j]] == 1) {
                    ++var_v.c_next_count;
                }
                else if (c2i[seq[j]] == 2) {
                    ++var_v.g_next_count;
                }
                else if (c2i[seq[j]] == 3) {
                    ++var_v.t_next_count;
                }
            }
        }
    }

    bool varGraph::quality(std::string q) {
        double error_score = 0;
        for (int i = 0; i < q.length(); i++) {
            error_score += pow(10, -(((double)q[i]) - 33) / 10);
            //cout << error_score;
        }
        //cout << error_score;
        if (error_score >= 2) {//function that finds quality of line
            return false;
        }
        return true;
    }

    bool varGraph::is_variation(std::string kmer) {
        kmer_key_t kmer_key;
        kmer_key = string_to_key(kmer);
        if (refGraph->has_key(kmer_key)) {
            return false;
        }
        return true;
    }

    void varGraph::error_prune() {
        vmap_it_t = vmap.begin();
        while (vmap_it_t != vmap.end()) {
            const var_value_t& ivalue = vmap_it_t->second;
            uint32_t vcount = (uint32_t)ivalue.a_next_count + (uint32_t)ivalue.c_next_count + (uint32_t)ivalue.g_next_count + (uint32_t)ivalue.t_next_count;
            //kmer_key_t ikey = vmap_it_t->first;
            //if (vcount == 1) { vmap.erase(ikey); }
            if (vcount == 1) {
                std::cerr << "****" << vmap_it_t->first.to_string() << std::endl;
                vmap_it_t = vmap.erase(vmap_it_t);
            }
            else {
                vmap_it_t++;
            }
            //if ( vmap_it_t != vmap.end() ) vmap_it_t++;
        }
    }

    kmer_key_t varGraph::string_to_key(const char* String) {
        kmer_key_t NewKey;

        for (int j = 0; j < K; ++j) {
            NewKey.acgt_set(j, c2i[String[j]]);
        }
        return NewKey;
    }

    kmer_key_t varGraph::string_to_key(const std::string & String) {
        return(string_to_key(String.c_str()));
    }

    void varGraph::print_key(kmer_key_t Key) {
        for (int32_t i = 0; i<K; ++i) {
            std::cout << i2c[Key.acgt_at(i)];
        }
        std::cout << " ";
    }

    void varGraph::print_map() {
        std::unordered_map<kmer_key_t, var_value_t, kmer_key_hasher> ::iterator its = vmap.begin();
        while (its != vmap.end()) {
            const kmer_key_t& ikey = its->first;
            for (int i = 0; i < K; i++) {
                std::cout << i2c[ikey.acgt_at(i)];
            }
            const var_value_t& ivalue = its->second;
            uint32_t vcount = (uint32_t)ivalue.a_next_count + (uint32_t)ivalue.c_next_count + (uint32_t)ivalue.g_next_count + (uint32_t)ivalue.t_next_count;
            std::cout << " :: " << vcount << " ";
            its++;
        }
    }

    void varGraph::nodes() {
        // declare iterator
        std::unordered_map<kmer_key_t, var_value_t, kmer_key_hasher> ::iterator its = vmap.begin();
        int non_unique_nodes = 0;

        while (its != vmap.end()) {
            const var_value_t& ivalue = its->second;
            uint32_t vcount = (uint32_t)ivalue.a_next_count + (uint32_t)ivalue.c_next_count + (uint32_t)ivalue.g_next_count + (uint32_t)ivalue.t_next_count;
            if (vcount != 1) {
                non_unique_nodes++;
            }
            its++;
        }
        std::cout << "Total nodes: " << vmap.size() << std::endl << "Number of unique nodes: " 
            << vmap.size() - non_unique_nodes << std::endl
            << "Number of non-unique nodes: " << non_unique_nodes << std::endl;
    }

    bool varGraph::store_to_file(std::string filename) {
        std::ofstream outs;

        outs.open(filename, std::ios::binary);
        while (!outs.is_open()) {
            std::cin.clear();
            std::string errorName;
            getline(std::cin, errorName);
            std::cout << "Invalid filename" << std::endl << "Enter filename: ";
            std::cin >> filename;
            outs.open(filename, std::ios::binary);
        }

        //Writes out kmer size and size of map
        uint32_t kmer_k = K;
        uint64_t sz_map = vmap.size();
        outs.write((char*)&kmer_k, sizeof(uint32_t));
        outs.write((char*)&sz_map, sizeof(uint64_t));

        //write out dictionary in format keyvalue
        for (vmap_it_t = vmap.begin(); vmap_it_t != vmap.end(); ++vmap_it_t) {
            outs.write((char*)&(vmap_it_t->first), sizeof(kmer_key_t));
            outs.write((char*)&(vmap_it_t->second), sizeof(var_value_t));
        }
        outs.close();
        return true;
    }
