#ifndef BLOOMH
#define BLOOMH
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include "omp.h"
#include "utils.h"
#include "BitMagic/src/bm.h"
#include "BitMagic/src/bmserial.h"
#include "BitMagic/src/bmundef.h"



using namespace std;





class Bloom{
public:
//TODO POWER OF TWO SIZE
    string filename;
    bool available;
    uint64_t size;
    uint number_hash_functions;
    bm::bvector<> * BV;

    void insert_key(uint64_t key);
    void optimize(){
        BV->optimize(NULL,bm::bvector<>::opt_compress);
    }
    bool check_key(uint64_t key);
    uint64_t get_cardinality()const;
    uint64_t dump_disk(bm::serializer<bm::bvector<> >& bvs);
    void load_disk();
    void free_ram();
    
    Bloom(const uint64_t Isize, uint Inumber_hash_functions,const string& Ifilename,bool Iavailable=true){
        filename=Ifilename;
        size=Isize-1;
        number_hash_functions=Inumber_hash_functions;
        BV=new bm::bvector<>(Isize,bm::BM_GAP);
        available=Iavailable;
        
    }

    ~Bloom(){
        delete BV;
    }
};

#endif
