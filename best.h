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
    string prefix;
    bool filter;
    Bloom* protection;
    ExponentialBloom<T>* trunk;
    ExponentialBloom<T>* reverse_trunk;
    uint K;
    uint64_t offsetUpdatekmer;
    uint64_t leaf_filters_size;
    uint64_t trunk_size;
    uint number_hash_function;
    uint64_t nb_insertions;



    Best(const uint64_t Itrunk_size,const uint64_t Ileaf_filters_size,const uint Inumber_hash_function,const uint Ik,const string Iprefix){
        K=Ik;
        prefix=Iprefix;
        trunk_size=Itrunk_size;
        leaf_filters_size=Ileaf_filters_size;
        number_hash_function=Inumber_hash_function;
        trunk=NULL;
        reverse_trunk=NULL;
        offsetUpdatekmer=1;
        offsetUpdatekmer<<=2*K;
        nb_insertions=0;
    }



    ~Best(){
        if(trunk!=NULL){delete trunk;}
        if(reverse_trunk!=NULL){delete reverse_trunk;}
        if(filter){delete protection;}
        for(uint i=0;i<leaf_filters.size();++i){
            delete leaf_filters[i];
        }
    }


    //LOW LEVEL FUNCTIONS
    void insert_key(const uint64_t key,uint level);
    void change_level();
    void check_key_leaf(const uint64_t key,const uint level)const;
    uint64_t rcb(uint64_t min)const;
    void update_kmer(uint64_t& min, char nuc)const;
    void update_kmer_RC(uint64_t& min, char nuc)const;
    void insert_last_leaf_trunk(uint level,ExponentialBloom<T>* EB);
    void load(zstr::ifstream* out);

    //HIGH LEVEL FUNCTIONS
    void insert_sequence(const string& reference) ;
    void insert_file(const string filename);
    void insert_file_of_file(const string filename);
    vector<T> query_key(const uint64_t key);
    vector<uint> query_sequence(const string& reference);
    void construct_reverse_trunk();
    void construct_trunk();
    void get_stats()const;
    void serialize(zstr::ostream* out)const;
    void optimize();
    void optimize(uint i);
    void dump(uint i);
    void add_leaf();
    void free_ram();
};

template class Best<uint8_t>;
template class Best<uint16_t>;
template class Best<uint32_t>;
#endif