#ifndef BLOOMH
#define BLOOMH
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include "omp.h"
#include "utils.h"
#include "best.h"
#include "BitMagic/src/bm.h"
#include "BitMagic/src/bmserial.h"
#include "BitMagic/src/bmundef.h"



using namespace std;



template <class T>
class Best;



template <class T>
class Bloom{
public:
    bm::bvector<> * BV;
    Best<T>* father;
    
    void insert_key(uint64_t key);
    void optimize(){
        BV->optimize(NULL,bm::bvector<>::opt_compress);
    }
    bool check_key(uint64_t key);
    uint64_t get_cardinality()const;
    uint64_t dump_disk(bm::serializer<bm::bvector<> >& bvs,zstr::ofstream* out,uint32_t i);
    void load_disk(zstr::ifstream* in);
    void free_ram();
    Bloom(Best<T>* Ifather);
    ~Bloom(){ delete BV;}
};



template class Bloom<uint8_t>;
template class Bloom<uint16_t>;
template class Bloom<uint32_t>;



#endif
