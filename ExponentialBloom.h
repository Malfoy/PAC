#ifndef EXPONENTIALBLOOMH
#define EXPONENTIALBLOOMH

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include "omp.h"



using namespace std;


class ExponentialBloom{
public:
    uint64_t size;
    uint number_hash_functions;
    uint64_t number_nonempty_cell;
    vector<uint8_t> filter;
    vector<uint64_t> level_bit_set;
    omp_lock_t lock[1024];

    void insert_key(uint64_t key,uint level);
    //return min level
    uint8_t check_key(uint64_t key)const;
    //return check key of a giver bloom
    bool check_key(uint64_t key,uint level)const;
    uint64_t get_cardinality(uint level)const;

    ExponentialBloom(const uint64_t Isize, const uint Inumber_hash_functions){
        size=Isize;
        number_hash_functions=Inumber_hash_functions;
        level_bit_set.resize(256,0);
        filter.resize(size,0);
        for(uint32_t i=0; i<1024;++i) {
            omp_init_lock(&lock[i]);
        }
    }
};
#endif