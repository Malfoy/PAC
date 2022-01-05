#ifndef BLOOMH
#define BLOOMH
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include "omp.h"
#include "BitMagic/src/bm.h"



using namespace std;





class Bloom{
public:
//TODO POWER OF TWO SIZE
    uint64_t size;
    uint number_hash_functions;
    uint64_t number_bit_set;
    bm::bvector<> * BV;

    void insert_key(uint64_t key);
    void optimize(){
        BV->optimize(NULL,bm::bvector<>::opt_compress);
    }
    bool check_key(uint64_t key)const;
    uint64_t get_cardinality()const;
    
    Bloom(const uint64_t Isize, uint Inumber_hash_functions){
        size=Isize;
        number_hash_functions=Inumber_hash_functions;
        number_bit_set=0;
        BV=new bm::bvector<>(size,bm::BM_GAP);
    }

    ~Bloom(){
        delete BV;
    }
};

#endif