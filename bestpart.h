#ifndef BESTPARTH
#define BESTPARTH

#include "Bloom.h"
#include "ExponentialBloom.h"
#include "best.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include <omp.h> 



using namespace std;





template <class T>
class BestPart{
public:
    vector< Best<T>* > buckets;
    omp_lock_t* mutex_array;
    uint K;
    uint64_t offsetUpdatekmer;
    uint64_t offsetUpdateminimizer;
    uint64_t leaf_filters_size;
    uint64_t trunk_size;
    uint64_t leaf_number;
    uint number_hash_function;
    uint small_minimizer_size;
    uint large_minimizer_size;
    uint bucket_number;
    uint large_minimizer_number;
    bool filter;
    bool hot;
    bool use_double_index;
    string w_dir;



    BestPart(const uint64_t Itrunk_size,const uint64_t Ileaf_filters_size,const uint Inumber_hash_function,const uint Ik, bool Ifilter,bool Ihot,const string Iwdir, bool Idouble,uint bucketing){
        K=Ik;
        use_double_index=Idouble;
        w_dir=Iwdir;
        hot=Ihot;
        filter=Ifilter;
        leaf_number=0;
        small_minimizer_size=bucketing;
        bucket_number=1<<(2*small_minimizer_size);
        large_minimizer_size=small_minimizer_size+3;
        large_minimizer_number=1<<(2*large_minimizer_size);
        trunk_size=Itrunk_size;
        leaf_filters_size=Ileaf_filters_size;
        number_hash_function=Inumber_hash_function;
        buckets.resize(bucket_number,NULL);
        mutex_array=new omp_lock_t[bucket_number];
        for(uint32_t i=0;i<bucket_number;++i){
            buckets[i]=new Best<T>(trunk_size/bucket_number,leaf_filters_size/bucket_number,number_hash_function,K,"Leon"+to_string(i));
            omp_init_lock(&mutex_array[i]);
        }
        offsetUpdatekmer=1;
        offsetUpdatekmer<<=2*K;
        offsetUpdateminimizer=1;
        offsetUpdateminimizer<<=(2*large_minimizer_size);
    }



    ~BestPart(){
        for(uint32_t i=0;i<bucket_number;++i){
            delete buckets[i];
        }
        delete [] mutex_array;
    }


    //LOW LEVEL FUNCTIONS
    void change_level();
    void check_key_leaf(const uint64_t key,const uint level)const;
    uint64_t rcb(uint64_t min,uint64_t n)const;
    void updateK(uint64_t& min, char nuc)const;
    void updateM(uint64_t& min, char nuc)const;
    void updateRCK(uint64_t& min, char nuc)const;
    void updateRCM(uint64_t& min, char nuc)const;
    void insert_last_leaf_trunk();
    uint64_t regular_minimizer_pos(uint64_t seq, uint64_t& position);
    uint64_t canonize(uint64_t x, uint64_t n);
    void insert_keys(const vector<uint64_t>& key,uint minimizer,uint level,Bloom* unique_filter);
    vector<T> query_keys(const vector<uint64_t>& key,uint minimizer);
    vector<pair<vector<uint64_t>,uint64_t> > get_super_kmers(const string& ref);
    //HIGH LEVEL FUNCTIONS
    void insert_sequence(const string& reference,uint level,Bloom* unique_filter) ;
    void insert_file(const string& filename,uint level);
    void insert_file_of_file(const string& filename);
    vector<uint> query_key(const uint64_t key);
    vector<uint32_t> query_sequence(const string& reference);
    void get_stats()const;
    void optimize();
    void query_file(const string& reference,const string& output);
    void serialize()const;
    void load(const string& existing_index);
    void index();
    void double_index();
    void add_leaf();
    uint64_t nb_bit_set()const;
    uint64_t nb_insertions()const ;
};

template class BestPart<uint8_t>;
template class BestPart<uint16_t>;
template class BestPart<uint32_t>;

#endif
