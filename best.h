#ifndef BESTH
#define BESTH



#include "Bloom.h"
#include "ExponentialBloom.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include "BitMagic/src/bm.h"
#include "BitMagic/src/bmserial.h"
#include "BitMagic/src/bmundef.h"



using namespace std;



template <class T>
class Bloom;



template <class T>
class Best{
public:
    vector<Bloom<T>*> leaf_filters;
    string prefix;
    bool filter;
    ExponentialBloom<T>* trunk;
    ExponentialBloom<T>* reverse_trunk;
    uint K;
    uint64_t size;
    uint number_hash_function;
    uint64_t number_bit_set;
    uint64_t number_bit_set_abt;
    uint64_t disk_space_used;
    zstr::ofstream* out;
    bool write;



    Best(const uint64_t Isize,const uint Inumber_hash_function,const uint Ik,const string Iprefix, bool Iwrite=true){
        K=Ik;
        prefix=Iprefix;
        write=Iwrite;
        if(write){
	        out=new zstr::ofstream(Iprefix,ios::trunc);
        }
        size=Isize-1;
        number_hash_function=Inumber_hash_function;
        trunk=NULL;
        reverse_trunk=NULL;
        number_bit_set_abt=number_bit_set=disk_space_used=0;
    }



    ~Best(){
        if(trunk!=NULL){delete trunk;}
        if(reverse_trunk!=NULL){delete reverse_trunk;}
        for(uint i=0;i<leaf_filters.size();++i){
            delete leaf_filters[i];
        }
        if(write){
	        delete out;
        }else{
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
    void insert_leaf_trunk(uint indice,uint level,ExponentialBloom<T>* EB);
    void load(uint64_t leaf_number,bool double_index);

    //HIGH LEVEL FUNCTIONS
    void insert_sequence(const string& reference) ;
    void insert_file(const string filename);
    void insert_file_of_file(const string filename);
    vector<T> query_key(const uint64_t key);
    vector<uint> query_sequence(const string& reference);
    void construct_reverse_trunk();
    uint64_t construct_trunk();
    void get_stats()const;
    void serialize();
    void optimize();
    void optimize(uint i);
    void dump(uint i,bm::serializer<bm::bvector<> >& bvs);
    void add_leaf();
    void free_ram();
    void load_bf(uint64_t leaf_number);
    T query_key_max(const uint64_t key);
    T query_key_min(const uint64_t key);

};



template class Best<uint8_t>;
template class Best<uint16_t>;
template class Best<uint32_t>;



#endif
