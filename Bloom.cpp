#include "Bloom.h"
#include "utils.h"
#include "BitMagic/src/bm.h"



void Bloom::insert_key(uint64_t key){
    for(uint64_t i=0; i<number_hash_functions;++i){
        uint64_t h=hash_family(key,i)%size;//TODO SIZE POWER OF TWO
        if(filter[h]==false){
            number_bit_set++;
            filter[h]=true;
        }
    }
}



bool Bloom::check_key(uint64_t key)const{
    for(uint64_t i=0; i<number_hash_functions;++i){
        uint64_t h=hash_family(key,i)%size;//TODO SIZE POWER OF TWO
        if(filter[h]==false){
            return false;
        }
    }
    return true;
}


uint64_t Bloom::get_cardinality()const{
    return number_bit_set/number_hash_functions;
}