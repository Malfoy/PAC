#include "ExponentialBloom.h"
#include "utils.h"


    
template <class T>
void ExponentialBloom<T>::insert_key(uint64_t key,uint level){
    for(uint64_t i=0; i<number_hash_functions;++i){
        uint64_t h=hash_family(key,i)&size;
        if((*filter)[h]==level){
            (*filter)[h]++;
        }
    }
}



template <class T>
T ExponentialBloom<T>::check_key(uint64_t key)const{
    T result=-1;
    for(uint64_t i=0; i<number_hash_functions;++i){
        uint64_t h=hash_family(key,i)&size;
        if((*filter)[h]<result){
            result=(*filter)[h];
        }
    }
    return result;
}



//return check key of a giver bloom
template <class T>
bool ExponentialBloom<T>::check_key(uint64_t key,uint level)const{
    for(uint64_t i=0; i<number_hash_functions;++i){
        uint64_t h=hash_family(key,i)&size;
        if((*filter)[h]<level){
            return false;
        }
    }
    return true;
}

