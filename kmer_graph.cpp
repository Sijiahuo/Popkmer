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


void kmerGraph::clear_map() {
    kmap.clear();
}
    
    int32_t kmerGraph::insert(const char* s) {

        kmer_key_t kmer_map_key;
        uint32_t prev; //number between 0-3
        uint32_t next; //number between 0-3
        int32_t uniq_words = 0;

        // to insert first k-mer
        int i = 0;
        kmer_map_key.acgt_set(i, c2i[s[i]]);

        for (i = 0; i < K; i++) {
            if (s[i] == '\0') {
                return 0;
            }
            else if (c2i[s[i]] > 3) {//(s[i] == 'N') {
                while (s[i] != '\0' && (c2i[s[i]] > 3)) {//(s[i] == 'N') {) {
                    i++;
                }
                return insert(&s[i]);
            }
            kmer_map_key.acgt_set(i, c2i[s[i]]);
        }

        kmap_it_t it = kmap.find(kmer_map_key);
        if (it == kmap.end()) { //check to see if it is unique
            ++uniq_words;
        }

        kmer_value& v = kmap[kmer_map_key];
        ++(v.count);
        if ((s[i] != '\0')) {
            v.acgt_next = (v.acgt_next | (0x01 << c2i[s[i]]));
        }

        prev = c2i[s[0]];

        while (s[i] != '\0') {
            bool flag2 = true;

            if (c2i[s[i]] > 3) {//(s[i] == 'N') {
                return uniq_words + insert(&s[i]);
            }

            for (int j = 0; j < K - 1; ++j) {
                kmer_map_key.acgt_set(j, kmer_map_key.acgt_at(j + 1));
            }

            kmer_map_key.acgt_set(K - 1, c2i[s[i]]);

            if (flag2) {
                it = kmap.find(kmer_map_key);
                if (it == kmap.end()) { //check to see if it is unique
                    ++uniq_words;
                }

                if (s[i + 1] != '\0') { //initialize/update "next" value
                    next = c2i[s[i + 1]];
                }
                else {
                    next = 4; //won't update acgt_next
                }

                kmer_value& v2 = kmap[kmer_map_key];
                //update count,acgt_prev,acgt_next
                ++(v2.count);
                if (prev < 4) v2.acgt_prev = (v2.acgt_prev | (0x01 << prev));
                if (next < 4) v2.acgt_next = (v2.acgt_next | (0x01 << next));


                prev = kmer_map_key.acgt_at(0); //update "prev" value
            }
            i++;
        }


        return uniq_words;
    }

    int32_t kmerGraph::insert(const std::string& s) { //calls insert function if passed a string
        return insert(s.c_str());
    }

    // basic query functions
    bool kmerGraph::has_key(kmer_key_t k) { return (kmap.find(k) != kmap.end()); }
    bool kmerGraph::has_key(const std::string s) {
        kmer_key_t key;
        for (int i = 0; i < K; i++) {
            key.acgt_set(i, c2i[s[i]]);
        }
        return (has_key(key));
    }

    kmerGraph::kmap_it_t kmerGraph::find(kmer_key_t k) { return kmap.find(k); }
    kmerGraph::kmap_it_t kmerGraph::find(const std::string s) {
        kmer_key_t key;
        for (int i = 0; i < K; i++) {
            key.acgt_set(i, c2i[s[i]]);
        }
        return (find(key));
    }

    kmer_value_t kmerGraph::value_at(kmer_key_t k) { return kmap[k]; }
    kmer_value_t kmerGraph::value_at(const std::string s) {
        kmer_key_t key;
        for (int i = 0; i < K; i++) {
            key.acgt_set(i, c2i[s[i]]);
        }
        return (value_at(key));
    }


    void kmerGraph::print_map() {
        int choice;
        std::cout << "Which value component would you like to print?" << std::endl << "1) Count"
            << std::endl << "2) Previous base" << std::endl << "3) Next base" << std::endl;
        std::cin >> choice;

        std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher> ::iterator its = kmap.begin();
        while (its != kmap.end()) {
            const kmer_key_t& ikey = its->first;
            for (int i = 0; i < K; i++) {
                std::cout << i2c[ikey.acgt_at(i)];
            }

            const kmer_value_t& ivalue = its->second;
            if (choice == 1) {
                std::cout << " :: " << (uint32_t)ivalue.count << " ";
            }
            else if (choice == 2) {
                char base;
                if ((uint32_t)ivalue.acgt_prev == 0) {
                    base = 'Z';
                }
                else if ((uint32_t)ivalue.acgt_prev == 1) {
                    base = 'A';
                }
                else if ((uint32_t)ivalue.acgt_prev == 2) {
                    base = 'C';
                }
                else if ((uint32_t)ivalue.acgt_prev == 4) {
                    base = 'G';
                }
                else if ((uint32_t)ivalue.acgt_prev == 8) {
                    base = 'T';
                }
                std::cout << " :: " << base << " ";
            }
            else if (choice == 3) {
                char base;
                if ((uint32_t)ivalue.acgt_next == 0) {
                    base = 'Z';
                }
                else if ((uint32_t)ivalue.acgt_next == 1) {
                    base = 'A';
                }
                else if ((uint32_t)ivalue.acgt_next == 2) {
                    base = 'C';
                }
                else if ((uint32_t)ivalue.acgt_next == 4) {
                    base = 'G';
                }
                else if ((uint32_t)ivalue.acgt_next == 8) {
                    base = 'T';
                }
                std::cout << " :: " << base << " ";
            }
            its++;
        }
    }

    char kmerGraph::comp_base(char c) {
        if (c == 'A' || c == 'a') return 'T';
        else if (c == 'T' || c == 't') return 'A';
        else if (c == 'C' || c == 'c') return 'G';
        else if (c == 'G' || c == 'g') return 'C';
        else return c;
    }

    std::string kmerGraph::string_complement(std::string s) {
        //complement string
        std::string c = s;

        for (int i = 0; i < s.length(); i++) c[i] = comp_base(s[s.length() - 1 - i]);
        return c;
    }

    ////////////////////////////Longest Path///////////////////////////////////////

    //Declare variables.

    int32_t LongestPathLength = 1;      //If there're no unique_path exist.
    std::vector<kmerGraph::kmap_it_t> Heads; //Store the heads of the longest paths. Save as iterators;
    bool loop = false;
    kmerGraph::kmap_it_t InfiniteLoopNode;

    //Check whether the current node has exactly one child
    bool kmerGraph::one_child_next(kmerGraph::kmap_it_t CurrentNode) {
        const kmer_value_t& itvalue = CurrentNode->second;
        uint32_t nextn = itvalue.acgt_next;
        return ((nextn == 1) || (nextn == 2) || (nextn == 4) || (nextn == 8));
    }

    int32_t kmerGraph::path_length(kmerGraph::kmap_it_t &CurrentNode) {
        const kmer_value_t& itvalue = CurrentNode->second;
        uint32_t pathlength = itvalue.path_length;
        return(pathlength);
    }

    // Return the Next Kmer
    kmer_key_t kmerGraph::next_kmer(kmerGraph::kmap_it_t CurrentNode) {

        kmer_key_t NextNodes;

        //const kmer_key_t& currentkey=CurrentNode->first;
        const kmer_value_t& NextKmer = CurrentNode->second;

        for (int j = 1; j < K; ++j) {
            NextNodes.acgt_set(j - 1, CurrentNode->first.acgt_at(j));
        }

        if (NextKmer.acgt_next == 8) {
            NextNodes.acgt_set(K - 1, 3);
        }
        else if (NextKmer.acgt_next == 4) {
            NextNodes.acgt_set(K - 1, 2);
        }
        else if (NextKmer.acgt_next == 2) {
            NextNodes.acgt_set(K - 1, 1);
        }
        else if (NextKmer.acgt_next == 1) {
            NextNodes.acgt_set(K - 1, 0);
        }
        return NextNodes;
    }

    //Change the length of the path
    void kmerGraph::change_path_length(kmerGraph::kmap_it_t CurrentNode, uint32_t NewLength) {
        const kmerGraph::kmap_it_t& pathlengthiter = CurrentNode;
        pathlengthiter->second.path_length = NewLength;
    }

    //A recursive function to find the path from current node to the end

    uint32_t kmerGraph::path_length_calculation(kmerGraph::kmap_it_t CurrentNode) {
        if (!one_child_next(CurrentNode)) {
            change_path_length(CurrentNode, 1);
            return 1;
        }
        else if (path_length(CurrentNode) >= 1) {
            return (path_length(CurrentNode));
        } //use original path length
        else {
            change_path_length(CurrentNode, -2);
            kmer_key_t Next_Kmer = next_kmer(CurrentNode);
            kmap_it_t NextKey = find(Next_Kmer);   //Searching for the next iterator according to the key
            if (path_length(NextKey) == -2) {
                loop = true;
                InfiniteLoopNode = NextKey;
                change_path_length(CurrentNode, 1);
                return 1;
            }

            if (NextKey == kmap.end()) {
                fprintf(stderr, "Cannot find next kmer %s from %s\n", Next_Kmer.to_string().c_str(), CurrentNode->first.to_string().c_str());
                abort();
            }

            change_path_length(CurrentNode, path_length_calculation(NextKey) + 1);
            if (path_length(CurrentNode) == LongestPathLength) Heads.push_back(CurrentNode);

            if (path_length(CurrentNode)>LongestPathLength) {
                LongestPathLength = path_length(CurrentNode);
                Heads.clear(); //Delete the original heads
                Heads.push_back(CurrentNode);
            }

            return (path_length(CurrentNode));
        }
    }

    uint8_t count37 = 0; //Everytime use LoopSolver, remember to set the count to 0.
    uint8_t kmerGraph::loop_solver(kmerGraph::kmap_it_t CurrentNode) {
        kmer_key_t Next_Kmer = next_kmer(CurrentNode);
        kmap_it_t  NextKey = find(Next_Kmer);   //Searching for the next iterator according to the key
        if (NextKey != InfiniteLoopNode) {
            count37++;
            change_path_length(CurrentNode, loop_solver(NextKey));
            return path_length(CurrentNode);
        }

        else {
            count37++;
            change_path_length(CurrentNode, count37);
            return count37;
        }
    }

    //The real function to find the longest path
    int32_t kmerGraph::longest_unique_path() {
        if (Heads.size() != 0) {
            Heads.clear();
        }

        kmerGraph::kmap_it_t Current = kmerGraph::kmap.begin(); //Start to loop over the dictionary.

        while (Current != kmap.end()) {
            if (path_length(Current) == -1) { //Not Visited
                Current->second.path_length = path_length_calculation(Current);
                if (loop == true) {
                    count37 = 0;
                    loop_solver(InfiniteLoopNode);
                    loop = false;
                }
            }
            Current++;
        }
        return LongestPathLength;
    }

    //Print out the key of a iterator
    void kmerGraph::print_key(kmap_it_t Position) {
        const kmer_key_t & ikey = Position->first;
        for (int32_t i = 0; i<K; ++i) {
            std::cout << i2c[ikey.acgt_at(i)];
        }
    }

    //Print out the last node of a iterator
    void kmerGraph::print_last_k(kmerGraph::kmap_it_t Position) {
        const kmer_key_t & ikey = Position->first;
        std::cout << i2c[ikey.acgt_at(K - 1)];

    }

    //Output the entire paths of the longest unique paths
    void kmerGraph::printing_out_paths_kmer() {
        if (LongestPathLength == 1) {
            printf("\nThere is no unique path with K equals to %d\n", (int)K);

        }
        else {
            printf("\nThe length of the unique path with K equals to %d is %d\n", (int)K, (int)LongestPathLength);
            std::cout << "The number of such path is " << Heads.size() << std::endl; //The number of paths
            std::cout << "The path(s) is(are):" << std::endl;

            for (std::vector<kmap_it_t>::iterator it = Heads.begin(); it != Heads.end(); ++it) {
                kmap_it_t TempIt = *it;
                int LoopCount = -1; //Count of the loop not if any

                while (path_length(TempIt) != 1 && (LoopCount>2 | LoopCount == -1)) { //The current node has decendents
                    print_key(TempIt); //Print out current kmer
                    uint8_t CurrentPathLength = path_length(TempIt); //Calculate the pathlength of currentnode
                    std::cout << " ";
                    TempIt = find(next_kmer(TempIt));
                    uint8_t NextPathLength = path_length(TempIt); //Take a look at the pathlength of next node
                    if (NextPathLength == CurrentPathLength) {
                        if (LoopCount == -1) {
                            LoopCount = CurrentPathLength;
                        }
                        else if (LoopCount>0) {
                            LoopCount -= 1; //Reduce by one.
                        }
                    } //Mark the length of the loop
                }
                print_key(TempIt); //Printout the ending node of the path
                std::cout << std::endl; //One entire path has been printed out.
            }
        }
    }

    void kmerGraph::printing_out_paths_string() {
        if (LongestPathLength == 1) {
            printf("\nThere is no unique path with K equals to %d\n", (int)K);

        }
        else {
            int stringlength = LongestPathLength + K - 1; //Count the length of the string
            printf("\nThe string length of the unique path with K equals to %d is %d\n", (int)K, (int)stringlength);

            std::cout << "The number of such path is " << Heads.size() << std::endl; //The number of paths
            std::cout << "The path(s) is(are): " << std::endl;
            for (std::vector<kmap_it_t>::iterator it = Heads.begin(); it != Heads.end(); ++it) {
                print_key(*it); //The whole string of the head node
                int LoopCount = -1; //Count the number of nodes in an infinite loop

                kmap_it_t TempIt = *it;

                while (path_length(TempIt) != 1 && (LoopCount>1 | LoopCount == -1)) { //The current node has decendents
                    uint8_t CurrentPathLength = path_length(TempIt); //Calculate the pathlength of currentnode
                    TempIt = find(next_kmer(TempIt));
                    print_last_k(TempIt); //Print out next kmer
                    uint8_t NextPathLength = path_length(TempIt); //Take a look at the pathlength of next node

                    if (NextPathLength == CurrentPathLength) {
                        if (LoopCount == -1) {
                            LoopCount = CurrentPathLength - 1;
                        }
                        else if (LoopCount>0) {
                            LoopCount -= 1; //Reduce by one.
                        }
                    }//Mark the length of the loop
                }
                std::cout << std::endl; //One entire path has been printed out.
            }
        }
    }

    const kmer_key_t kmerGraph::string_to_key(const char* String) {
        kmer_key_t NewKey;

        for (int j = 0; j < K; ++j) {
            NewKey.acgt_set(j, c2i[String[j]]);
        }
        return NewKey;
    }

    const kmer_key_t kmerGraph::string_to_key(const std::string & String) {
        return(string_to_key(String.c_str()));
    }

    const void kmerGraph::print_key(kmer_key_t Key) {
        for (int32_t i = 0; i<K; ++i) {
            std::cout << i2c[Key.acgt_at(i)];
        }
    }

    void kmerGraph::printing_NodePath_string(std::string InputKmer) {

        kmer_key_t Input = string_to_key(InputKmer); //Change string to key
        kmap_it_t StartingIterator = find(Input); //Search for the interator of the key


        if (path_length(StartingIterator) == 1) {
            std::cout << "\nThere is no path starting from current kmer." << std::endl;
        }
        else {
            int stringlength = path_length(StartingIterator) + K - 1; //Count the length of the string

            printf("The string length of the path starting from %s is %d\n", InputKmer.c_str(), stringlength);

            std::cout << "The path is: " << std::endl;

            print_key(StartingIterator); //The whole string of the head node
            int LoopCount = -1; //Count the number of nodes in an infinite loop

            kmap_it_t TempIt = StartingIterator;

            while (path_length(TempIt) != 1 && (LoopCount>1 | LoopCount == -1)) { //The current node has decendents
                uint8_t CurrentPathLength = path_length(TempIt); //Calculate the pathlength of currentnode
                TempIt = find(next_kmer(TempIt));
                print_last_k(TempIt); //Print out next kmer
                uint8_t NextPathLength = path_length(TempIt); //Take a look at the pathlength of next node

                if (NextPathLength == CurrentPathLength) {
                    if (LoopCount == -1) {
                        LoopCount = CurrentPathLength - 1;
                    }
                    else if (LoopCount>0) {
                        LoopCount -= 1; //Reduce by one.
                    }

                }//Mark the length of the loop


            }
            std::cout << std::endl; //One entire path has been printed out.

        }

    }

    std::string kmerGraph::string_reader() {
        std::ifstream ins;
        std::string filename;
        std::string fileSequence;
        std::string baseSequence = "";

        //Open file (however many tries that takes)
        std::cout << "Enter filename: ";
        std::cin >> filename;
        ins.open(filename);
        while (!ins.is_open()) {
            std::cin.clear();
            std::string errorName;
            getline(std::cin, errorName);
            std::cout << "Invalid filename" << std::endl << "Enter filename: ";
            std::cin >> filename;
            ins.open(filename);
        }


        //Reads first line (starts with non-sequence characters
        // that need to be thrown out)
        ins >> fileSequence;
        for (int i = 0; i < fileSequence.length(); i++) {
            if (fileSequence[i] == 'N' || fileSequence[i] == 'A' ||
                fileSequence[i] == 'G' || fileSequence[i] == 'C' ||
                fileSequence[i] == 'T') {
                baseSequence += fileSequence[i];
            }
        }

        // Reads in remainder of sequence
        while (ins >> fileSequence) {
            baseSequence += fileSequence;
        }

        // Prints full sequence if prompted
        std::string userInput;
        std::cout << "Print Sequence (Yes or No)? ";
        std::cin >> userInput;
        if (userInput == "Yes" || userInput == "yes") {
            std::cout << baseSequence << std::endl;
        }

        return baseSequence;
    }


    bool kmerGraph::store_to_file(std::string filename) {
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
        uint64_t sz_map = kmerGraph::kmap.size();
        outs.write((char*)&kmer_k, sizeof(uint32_t));
        outs.write((char*)&sz_map, sizeof(uint64_t));

        //write out dictionary in format keyvalue
        for (kmerGraph::kmap_it_t it = kmerGraph::kmap.begin(); it != kmap.end(); ++it) {
            outs.write((char*)&(it->first), sizeof(kmer_key_t));
            outs.write((char*)&(it->second), sizeof(kmer_value_t));
        }
        outs.close();
        return true;
    }

    bool kmerGraph::load_from_file() {
        std::ifstream ins;
        std::string filename;

        //Open file (however many tries that takes)
        std::cout << "Enter filename: ";
        std::cin >> filename;
        ins.open(filename, std::ios::binary);
        while (!ins.is_open()) {
            std::cin.clear();
            std::string errorName;
            getline(std::cin, errorName);
            std::cout << "Invalid filename" << std::endl << "Enter filename: ";
            std::cin >> filename;
            ins.open(filename, std::ios::binary);
        }

        // Read in kmer length and size of map
        uint32_t kmer_k;
        uint64_t sz_map;
        ins.read((char*)&kmer_k, sizeof(uint32_t));
        ins.read((char*)&sz_map, sizeof(uint64_t));

        if (kmer_k != K) {
            std::cout << "Wrong K value";
            //not sure what this print statement does
            //fprintf(stderr, "The K value does not match %u vs %u", kmer_k, K);
            return false;
        }

        kmer_key_t mykey;
        kmer_value_t myvalue;

        for (uint64_t i = 0; i < sz_map; ++i) {
            ins.read((char*)&mykey, sizeof(kmer_key_t));
            ins.read((char*)&myvalue, sizeof(kmer_value_t));
            kmap.emplace(mykey, myvalue);
        }
        return true;
    }

    void kmerGraph::largest_count() {
        // declare iterator
        std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher> ::iterator its = kmap.begin();

        // set variables to hold key & value info
        const kmer_key_t& ikey = its->first;
        const kmer_value_t& ivalue = its->second;

        //variables to hold largest valued key/its value
        std::string skey = "";
        uint32_t value = (uint32_t)ivalue.count;

        for (int i = 0; i < K; i++) {
            skey += i2c[ikey.acgt_at(i)];
        }

        while (its != kmap.end()) {
            const kmer_value_t& ivalue = its->second;
            if ((uint32_t)ivalue.count > value) {
                // replaces old key with new one
                const kmer_key_t& ikey = its->first;
                skey = "";
                for (int i = 0; i < K; i++) {
                    skey += i2c[ikey.acgt_at(i)];
                }

                // replaces old value with new one
                value = (uint32_t)ivalue.count;
            }
            else if ((uint32_t)ivalue.count == value) {
                // Adds new key to old one
                const kmer_key_t& ikey = its->first;
                skey += ", ";
                for (int i = 0; i < K; i++) {
                    skey += i2c[ikey.acgt_at(i)];
                }
            }

            its++;
        }
        std::cout << "The largest kmer value pair is " << skey << " :: " << value << std::endl;
    }


    void kmerGraph::nodes() {
        // declare iterator
        std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher> ::iterator its = kmap.begin();
        int non_unique_nodes = 0;

        while (its != kmap.end()) {
            const kmer_value_t& ivalue = its->second;
            if ((uint32_t)ivalue.count != 1) {
                non_unique_nodes++;
            }
            its++;
        }
        std::cout << "Total nodes: " << kmap.size() << std::endl << "Number of unique nodes: " << kmap.size() - non_unique_nodes << std::endl
            << "Number of non-unique nodes: " << non_unique_nodes << std::endl;
    }
