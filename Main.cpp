#include <cstdio>
#include <unordered_map>
#include <stdint.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;
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

    string to_string() const {
        string s;
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



class kmerGraph {
private:
    std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher> kmap;

public:
    typedef std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher>::iterator kmap_it_t;

    // to construct the hashtable
    int32_t insert(const char* s) {//, uint32_t prev, uint32_t next){

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

    int32_t insert(const std::string& s) { //calls insert function if passed a string
        return insert(s.c_str());
    }

   // basic query functions
    bool has_key(kmer_key_t k) { return (kmap.find(k) != kmap.end()); }
    bool has_key(const string s) {
        kmer_key_t key;
        for (int i = 0; i < K; i++) {
            key.acgt_set(i, c2i[s[i]]);
        }
        return (has_key(key));
    }

    kmap_it_t find(kmer_key_t k) { return kmap.find(k); }
    kmap_it_t find(const string s) {
        kmer_key_t key;
        for (int i = 0; i < K; i++) {
            key.acgt_set(i, c2i[s[i]]);
        }
        return (find(key));
    }

    kmer_value_t value_at(kmer_key_t k) { return kmap[k]; }
    kmer_value_t value_at(const string s) {
        kmer_key_t key;
        for (int i = 0; i < K; i++) {
            key.acgt_set(i, c2i[s[i]]);
        }
        return (value_at(key));
    }


    void printMap() {
        std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher> ::iterator its = kmap.begin();
        while (its != kmap.end()) {
            const kmer_key_t& ikey = its->first;
            for (int i = 0; i < K; i++) {
                cout << i2c[ikey.acgt_at(i)];
            }
            const kmer_value_t& ivalue = its->second;
            cout << " :: " << (uint32_t)ivalue.count << " ";
            //cout << " :: " << (uint32_t)ivalue.acgt_prev << " ";
            //cout << " :: " << (uint32_t)ivalue.acgt_next << " ";
            its++;
        }
    }

    char compBase(char c) {
        if (c == 'A' || c == 'a') return 'T';
        else if (c == 'T' || c == 't') return 'A';
        else if (c == 'C' || c == 'c') return 'G';
        else if (c == 'G' || c == 'g') return 'C';
        else return c;
    }

    std::string stringComplement(std::string s) {
        //complement string
        std::string c = s;

        for (int i = 0; i < s.length(); i++) c[i] = compBase(s[s.length() - 1 - i]);
        return c;
    }

    ////////////////////////////Longest Path///////////////////////////////////////

    //Declare variables.

    int32_t LongestPathLength = 1;      //If there're no unique_path exist.
    vector<kmap_it_t> Heads;            //Store the heads of the longest paths. Save as iterators;
    bool loop = false;
    kmap_it_t InfiniteLoopNode;

    //Check whether the current node has exactly one child
    bool OneChildNext(kmap_it_t CurrentNode) {
        const kmer_value_t& itvalue = CurrentNode->second;
        uint32_t nextn = itvalue.acgt_next;
        return ((nextn == 1) || (nextn == 2) || (nextn == 4) || (nextn == 8));
    }

    int32_t PathLength(kmap_it_t &CurrentNode) {
        const kmer_value_t& itvalue = CurrentNode->second;
        uint32_t pathlength = itvalue.path_length;
        return(pathlength);
    }

    // Return the Next Kmer
    kmer_key_t NextKmer(kmap_it_t CurrentNode) {

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
    void ChangePathLength(kmap_it_t CurrentNode, uint32_t NewLength) {
        const kmap_it_t& pathlengthiter = CurrentNode;
        pathlengthiter->second.path_length = NewLength;
    }

    //A recursive function to find the path from current node to the end

    uint32_t PathLengthCalculation(kmap_it_t CurrentNode) {
        if (!OneChildNext(CurrentNode)) {
            ChangePathLength(CurrentNode, 1);
            return 1;
        }
        else if (PathLength(CurrentNode) >= 1) {
            return (PathLength(CurrentNode));
        } //use original path length
        else {
            ChangePathLength(CurrentNode, -2);
            kmer_key_t Next_Kmer = NextKmer(CurrentNode);
            kmap_it_t NextKey = find(Next_Kmer);   //Searching for the next iterator according to the key
            if (PathLength(NextKey) == -2) {
                loop = true;
                InfiniteLoopNode = NextKey;
                ChangePathLength(CurrentNode, 1);
                return 1;
            }

            if (NextKey == kmap.end()) {
                fprintf(stderr, "Cannot find next kmer %s from %s\n", Next_Kmer.to_string().c_str(), CurrentNode->first.to_string().c_str());
                abort();
            }

            ChangePathLength(CurrentNode, PathLengthCalculation(NextKey) + 1);
            if (PathLength(CurrentNode) == LongestPathLength) Heads.push_back(CurrentNode);

            if (PathLength(CurrentNode)>LongestPathLength) {
                LongestPathLength = PathLength(CurrentNode);
                Heads.clear(); //Delete the original heads
                Heads.push_back(CurrentNode);
            }

            return (PathLength(CurrentNode));
        }
    }

    uint8_t count = 0; //Everytime use LoopSolver, remember to set the count to 0.
    uint8_t LoopSolver(kmap_it_t CurrentNode) {
        kmer_key_t Next_Kmer = NextKmer(CurrentNode);
        kmap_it_t  NextKey = find(Next_Kmer);   //Searching for the next iterator according to the key
        if (NextKey != InfiniteLoopNode) {
            count++;
            ChangePathLength(CurrentNode, LoopSolver(NextKey));
            return PathLength(CurrentNode);
        }

        else {
            count++;
            ChangePathLength(CurrentNode, count);
            return count;
        }
    }

    //The real function to find the longest path
    int32_t longest_unique_path() {
        if (Heads.size() != 0) {
            Heads.clear();
        }

        kmap_it_t Current = kmap.begin(); //Start to loop over the dictionary.

        while (Current != kmap.end()) {
            if (PathLength(Current) == -1) { //Not Visited
                Current->second.path_length = PathLengthCalculation(Current);
                if (loop == true) {
                    count = 0;
                    LoopSolver(InfiniteLoopNode);
                    loop = false;
                }
            }
            Current++;
        }
        return LongestPathLength;
    }

    //Print out the key of a iterator
    void PrintKey(kmap_it_t Position) {
        const kmer_key_t & ikey = Position->first;
        for (int32_t i = 0; i<K; ++i) {
            cout << i2c[ikey.acgt_at(i)];
        }
    }

    //Print out the last node of a iterator
    void PrintLastK(kmap_it_t Position) {
        const kmer_key_t & ikey = Position->first;
        cout << i2c[ikey.acgt_at(K - 1)];

    }

    //Output the entire paths of the longest unique paths
    void printing_out_paths_kmer() {
        if (LongestPathLength == 1) {
            printf("\nThere is no unique path with K equals to %d\n", (int)K);

        }
        else {
            printf("\nThe length of the unique path with K equals to %d is %d\n", (int)K, (int)LongestPathLength);
            cout << "The number of such path is " << Heads.size() << endl; //The number of paths
            cout << "The path(s) is(are):" << endl;

            for (vector<kmap_it_t>::iterator it = Heads.begin(); it != Heads.end(); ++it) {
                kmap_it_t TempIt = *it;
                int LoopCount = -1; //Count of the loop not if any

                while (PathLength(TempIt) != 1 && (LoopCount>2 | LoopCount == -1)) { //The current node has decendents
                    PrintKey(TempIt); //Print out current kmer
                    uint8_t CurrentPathLength = PathLength(TempIt); //Calculate the pathlength of currentnode
                    cout << " ";
                    TempIt = find(NextKmer(TempIt));
                    uint8_t NextPathLength = PathLength(TempIt); //Take a look at the pathlength of next node
                    if (NextPathLength == CurrentPathLength) {
                        if (LoopCount == -1) {
                            LoopCount = CurrentPathLength;
                        }
                        else if (LoopCount>0) {
                            LoopCount -= 1; //Reduce by one.
                        }
                    } //Mark the length of the loop
                }
                PrintKey(TempIt); //Printout the ending node of the path
                cout << endl; //One entire path has been printed out.
            }
        }
    }

    void printing_out_paths_string() {
        if (LongestPathLength == 1) {
            printf("\nThere is no unique path with K equals to %d\n", (int)K);

        }
        else {
            int stringlength = LongestPathLength + K - 1; //Count the length of the string
            printf("\nThe string length of the unique path with K equals to %d is %d\n", (int)K, (int)stringlength);

            cout << "The number of such path is " << Heads.size() << endl; //The number of paths
            cout << "The path(s) is(are): " << endl;
            for (vector<kmap_it_t>::iterator it = Heads.begin(); it != Heads.end(); ++it) {
                PrintKey(*it); //The whole string of the head node
                int LoopCount = -1; //Count the number of nodes in an infinite loop

                kmap_it_t TempIt = *it;

                while (PathLength(TempIt) != 1 && (LoopCount>1 | LoopCount == -1)) { //The current node has decendents
                    uint8_t CurrentPathLength = PathLength(TempIt); //Calculate the pathlength of currentnode
                    TempIt = find(NextKmer(TempIt));
                    PrintLastK(TempIt); //Print out next kmer
                    uint8_t NextPathLength = PathLength(TempIt); //Take a look at the pathlength of next node

                    if (NextPathLength == CurrentPathLength) {
                        if (LoopCount == -1) {
                            LoopCount = CurrentPathLength - 1;
                        }
                        else if (LoopCount>0) {
                            LoopCount -= 1; //Reduce by one.
                        }
                    }//Mark the length of the loop
                }
                cout << endl; //One entire path has been printed out.
            }
        }
    }

    const kmer_key_t StringToKey(const char* String) {
        kmer_key_t NewKey;

        for (int j = 0; j < K; ++j) {
            NewKey.acgt_set(j, c2i[String[j]]);
        }
        return NewKey;
    }

    const kmer_key_t StringToKey(const string & String) {
        return(StringToKey(String.c_str()));
    }

    const void PrintKey(kmer_key_t Key) {
        for (int32_t i = 0; i<K; ++i) {
            cout << i2c[Key.acgt_at(i)];
        }
    }

    //Print out the determinate path of an arbitrary node

    void printing_NodePath_string(string InputKmer) {

        kmer_key_t Input = StringToKey(InputKmer); //Change string to key
        kmap_it_t StartingIterator = find(Input); //Search for the interator of the key


        if (PathLength(StartingIterator) == 1) {
            cout << "\nThere is no path starting from current kmer." << endl;
        }
        else {
            int stringlength = PathLength(StartingIterator) + K - 1; //Count the length of the string

            printf("The string length of the path starting from %s is %d\n", InputKmer.c_str(), stringlength);

            //cout << "The string length of the path is " << (int)stringlength << endl;
            cout << "The path is: " << endl;

            PrintKey(StartingIterator); //The whole string of the head node
            int LoopCount = -1; //Count the number of nodes in an infinite loop

            kmap_it_t TempIt = StartingIterator;

            while (PathLength(TempIt) != 1 && (LoopCount>1 | LoopCount == -1)) { //The current node has decendents
                uint8_t CurrentPathLength = PathLength(TempIt); //Calculate the pathlength of currentnode
                TempIt = find(NextKmer(TempIt));
                PrintLastK(TempIt); //Print out next kmer
                uint8_t NextPathLength = PathLength(TempIt); //Take a look at the pathlength of next node

                if (NextPathLength == CurrentPathLength) {
                    if (LoopCount == -1) {
                        LoopCount = CurrentPathLength - 1;
                    }
                    else if (LoopCount>0) {
                        LoopCount -= 1; //Reduce by one.
                    }

                }//Mark the length of the loop


            }
            cout << endl; //One entire path has been printed out.

        }

    }



    string stringReader() {
        ifstream ins;
        string filename;
        string fileSequence;
        string baseSequence = "";

        //Open file (however many tries that takes)
        cout << "Enter filename: ";
        cin >> filename;
        ins.open(filename);
        while (!ins.is_open()) {
            cin.clear();
            string errorName;
            getline(cin, errorName);
            cout << "Invalid filename" << endl << "Enter filename: ";
            cin >> filename;
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
        string userInput;
        cout << "Print Sequence (Yes or No)? ";
        cin >> userInput;
        if (userInput == "Yes" || userInput == "yes") {
            cout << baseSequence << endl;
        }

        return baseSequence;
    }


    bool storeToFile(string filename) {
        ofstream outs;

        outs.open(filename, ios::binary);
        while (!outs.is_open()) {
            cin.clear();
            string errorName;
            getline(cin, errorName);
            cout << "Invalid filename" << endl << "Enter filename: ";
            cin >> filename;
            outs.open(filename, ios::binary);
        }

        //Writes out kmer size and size of map
        uint32_t kmer_k = K;
        uint64_t sz_map = kmap.size();
        outs.write((char*)&kmer_k, sizeof(uint32_t));
        outs.write((char*)&sz_map, sizeof(uint64_t));

        //write out dictionary in format keyvalue
        for (kmap_it_t it = kmap.begin(); it != kmap.end(); ++it) {
            outs.write((char*)&(it->first), sizeof(kmer_key_t));
            outs.write((char*)&(it->second), sizeof(kmer_value_t));
        }
        outs.close();
        return true;
    }

    bool loadFromFile() {
        ifstream ins;
        string filename;

        //Open file (however many tries that takes)
        cout << "Enter filename: ";
        cin >> filename;
        ins.open(filename, ios::binary);
        while (!ins.is_open()) {
            cin.clear();
            string errorName;
            getline(cin, errorName);
            cout << "Invalid filename" << endl << "Enter filename: ";
            cin >> filename;
            ins.open(filename, ios::binary);
        }

        // Read in kmer length and size of map
        uint32_t kmer_k;
        uint64_t sz_map;
        ins.read((char*)&kmer_k, sizeof(uint32_t));
        ins.read((char*)&sz_map, sizeof(uint64_t));

        if (kmer_k != K) {
            cout << "Wrong K value";
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

    void largestCount() {
        // declare iterator
        std::unordered_map<kmer_key_t, kmer_value_t, kmer_key_hasher> ::iterator its = kmap.begin();

        // set variables to hold key & value info
        const kmer_key_t& ikey = its->first;
        const kmer_value_t& ivalue = its->second;

        //variables to hold largest valued key/its value
        string skey = "";
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

            its++;
        }
        cout << "The largest kmer value pair is " << skey << " :: " << value << endl;
    }


    void nodes() {
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
        cout << "Total nodes: " << kmap.size() << endl << "Number of unique nodes: " << kmap.size() - non_unique_nodes << endl
            << "Number of non-unique nodes: " << non_unique_nodes << endl;
    }
};







class varGraph {
private:
    std::unordered_map<kmer_key_t, var_value_t, kmer_key_hasher> vmap;
    kmerGraph* refGraph;

public:
    std::unordered_map<kmer_key_t, var_value_t, kmer_key_hasher>::iterator vmap_it_t;
    bool cflag = false;


    //Different Constructors for varGraph. We should choose which one we want to use
    varGraph(const char* ref_graph_name, const char* var_graph_name) {
        refGraph = new kmerGraph();
        cflag = true;
        //refGraph->loadFromFile;
        //loadFromFile(var_graph_name); //these two lines give me errors so i just commented them out.
    }


    varGraph(const char* ref_graph_name) {
        refGraph = new kmerGraph();
        cflag = true;
    }

    varGraph(kmerGraph* r) {
        refGraph = r;
        cflag = false;
        // kmerGraph g;
        // g.loadFromFile
    }

    ~varGraph() {
        // if it was created with new (probably need a flag)
        if (cflag) { delete refGraph; }
    }

    varGraph() {

    }


    void fq_read(string filename) {
        ifstream ins;
        ins.open(filename);

        // if invalid file is entered
        while (!ins.is_open()) {
            cin.clear();
            string errorName;
            getline(cin, errorName);
            cout << "Invalid filename" << endl << "Enter filename: ";
            cin >> filename;
            ins.open(filename);
        }

        string identifier;
        string sequence;
        string plus;
        string quality;
        string kmer = "";
        string qual = "";
        //const int SIZE = 150;

        while (ins >> identifier >> sequence >> plus >> quality) {
            insert(sequence, quality);
        }
    }

    void insert(string seq, string qual) {
        uint64_t SIZE = seq.length();
        int i;
        int j;
        string kmer = "";
        string kmer_qual = "";
        for (i = 0; i < SIZE - K + 1; i++) {
            kmer = "";
            kmer_qual = "";
            for (j = i; j < i + K; j++) {
                kmer += seq[j];
                kmer_qual += qual[j];
            }
            if (isVariation(kmer) && quality(kmer_qual)) {
                var_value& var_v = vmap[StringToKey(kmer)];
                if (c2i[seq[j]] == 0) {
                    ++var_v.a_next_count;
                }
                else if (c2i[seq[j]] == 1) {
                    ++var_v.c_next_count;
                }
                else if (c2i[seq[i]] == 2) {
                    ++var_v.g_next_count;
                }
                else if (c2i[seq[i]] == 3) {
                    ++var_v.t_next_count;
                }
            }
        }
    }



    bool quality(string q) {
        double error_score = 0;
        for (int i = 0; i < q.length(); i++) {
            error_score += pow(10, -(((double)q[i]) - 33) / 10);
            //cout << error_score;
        }
        //cout << error_score;
        if (error_score >= K) {//function that finds quality of line
            return false;
        }
        return true;
    }

    bool isVariation(string kmer) {
        kmer_key_t kmer_key;
        kmer_key = StringToKey(kmer);
        if (refGraph->has_key(kmer_key)) {
            return false;
        }
        return true;
    }




    kmer_key_t StringToKey(const char* String) {
        kmer_key_t NewKey;

        for (int j = 0; j < K; ++j) {
            NewKey.acgt_set(j, c2i[String[j]]);
        }
        return NewKey;
    }

    kmer_key_t StringToKey(const string & String) {
        return(StringToKey(String.c_str()));
    }

    void PrintKey(kmer_key_t Key) {
        for (int32_t i = 0; i<K; ++i) {
            cout << i2c[Key.acgt_at(i)];
        }
        cout << " ";
    }

    void printMap() {
        std::unordered_map<kmer_key_t, var_value_t, kmer_key_hasher> ::iterator its = vmap.begin();
        while (its != vmap.end()) {
            const kmer_key_t& ikey = its->first;
            for (int i = 0; i < K; i++) {
                cout << i2c[ikey.acgt_at(i)];
            }
            const var_value_t& ivalue = its->second;
            uint32_t vcount = (uint32_t)ivalue.a_next_count + (uint32_t)ivalue.c_next_count + (uint32_t)ivalue.g_next_count + (uint32_t)ivalue.t_next_count;
            cout << " :: " << vcount << " ";
            its++;
        }
    }

    void nodes() {
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
        cout << "Total nodes: " << vmap.size() << endl << "Number of unique nodes: " << vmap.size() - non_unique_nodes << endl
            << "Number of non-unique nodes: " << non_unique_nodes << endl;
    }




};










int main() {
    kmerGraph kG;
    varGraph vG(&kG);
    string s = kG.stringReader();
    kG.insert(s);
    vG.fq_read("Text.txt");
    vG.nodes();

    return 0;
}

