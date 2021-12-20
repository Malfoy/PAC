#ifndef BESTH
#define BESTH

#include "Bloom.h"
#include "ExponentialBloom.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 



using namespace std;



//TODO SERIALIZATION
template <class T>
class Best{
public:
    vector<Bloom*> leaf_filters;
    ExponentialBloom<T>* trunk;
    uint current_level;
    uint K;
    uint64_t offsetUpdatekmer;
    uint64_t leaf_filters_size;
    uint64_t trunk_size;
    uint number_hash_function;
    uint64_t nb_insertions;



    Best(const uint64_t Itrunk_size,const uint64_t Ileaf_filters_size,const uint Inumber_hash_function,const uint Ik){
        K=Ik;
        trunk_size=Itrunk_size;
        leaf_filters_size=Ileaf_filters_size;
        number_hash_function=Inumber_hash_function;
        trunk=new ExponentialBloom<T>(trunk_size,number_hash_function);
        current_level=0;
        offsetUpdatekmer=1;
        offsetUpdatekmer<<=2*K;
        nb_insertions=0;
        leaf_filters.push_back(new Bloom(leaf_filters_size,number_hash_function));
    }



    ~Best(){
        delete trunk;
        for(uint i=0;i<leaf_filters.size();++i){
            delete leaf_filters[i];
        }

    }


    //LOW LEVEL FUNCTIONS
    void insert_key(const uint64_t key);
    void change_level();
    void check_key_leaf(const uint64_t key,const uint level)const;
    uint64_t rcb(uint64_t min)const;
    void update_kmer(uint64_t& min, char nuc)const;
    void update_kmer_RC(uint64_t& min, char nuc)const;
    void insert_last_leaf_trunk();
    void load(zstr::ifstream* out);

    //HIGH LEVEL FUNCTIONS
    void insert_sequence(const string& reference) ;
    void insert_file(const string filename);
    void insert_file_of_file(const string filename);
    vector<T> query_key(const uint64_t key);
    vector<uint> query_sequence(const string& reference);
    void get_stats()const;
    void serialize(zstr::ostream* out)const;
    void optimize();
};

template class Best<uint8_t>;
template class Best<uint16_t>;
template class Best<uint32_t>;
#endif