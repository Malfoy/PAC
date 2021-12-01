#ifndef BLOOMH
#define BLOOMH
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include "omp.h"



using namespace std;



class Bloom{
public:
    uint64_t size;
    uint number_hash_functions;
    uint64_t number_bit_set;
    vector<bool> filter;
    omp_lock_t lock[1024];

    void insert_key(uint64_t key);
    bool check_key(uint64_t key)const;
    uint64_t get_cardinality()const;
    
    Bloom(const uint64_t Isize, uint Inumber_hash_functions){
        size=Isize;
        number_hash_functions=Inumber_hash_functions;
        number_bit_set=0;
        filter.resize(size,false);
        for(uint32_t i=0; i<1024;++i) {
            omp_init_lock(&lock[i]);
        }
    }
};

#endif