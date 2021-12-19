#ifndef BCARDIH
#define BCARDIH

#include "Bloom.h"
#include "ExponentialBloom.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 



using namespace std;

template <class T>
class Bcardi{
public:
    vector<Bloom*> big_filters;
    ExponentialBloom<T>* little_filters;
    uint current_level;
    uint K;
    uint64_t offsetUpdatekmer;

    Bcardi(const uint number_big_filter,const uint64_t bigest_filter_size,const uint64_t little_filter_size,const uint number_hash_function,const uint Ik){
        K=Ik;
        little_filters=new ExponentialBloom<T>(little_filter_size,number_hash_function);
        for(uint i=0;i<number_big_filter;++i){
            big_filters.push_back(new Bloom(bigest_filter_size,number_hash_function));
        }
        current_level=0;
        offsetUpdatekmer=1;
        offsetUpdatekmer<<=2*K;
    }

    ~Bcardi(){
        delete little_filters;
        for(uint i=0;i<big_filters.size();++i){
            delete big_filters[i];
        }

    }
    //LOW LEVEL FUNCTIONS
    void insert_key(const uint64_t key);
    void change_level();
    void check_key_level(const uint64_t key,const uint level)const;
    uint64_t rcb(uint64_t min)const;
    void update_kmer(uint64_t& min, char nuc)const;
    void update_kmer_RC(uint64_t& min, char nuc)const;

    vector<uint64_t> get_cardinalities()const;
    //HIGH LEVEL FUNCTIONS
    void insert_sequence(const string& reference) ;
    void insert_file(const string filename);
    void insert_file_of_file(const string filename);
};

#endif