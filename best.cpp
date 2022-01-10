#include "Bloom.h"
#include "ExponentialBloom.h"
#include "utils.h"
#include "best.h"
#include "BitMagic/src/bmserial.h"
#include "BitMagic/src/bmundef.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include <filesystem>



using namespace std;
using namespace filesystem;



template <class T>
void Best<T>::construct_reverse_trunk(){
    reverse_trunk=new ExponentialBloom<T>(trunk_size,number_hash_function);
    uint64_t lfs(leaf_filters.size()-1);
    for(int i = leaf_filters.size()-1; i>=0; i--){
        insert_leaf_trunk(i,lfs-i,trunk);//TODO TO CHECK
        //~ insert_last_leaf_trunk(i,reverse_trunk);
        
    }
}



template <class T>
void Best<T>::construct_trunk(){
    trunk=new ExponentialBloom<T>(trunk_size,number_hash_function);
    uint64_t lfs(leaf_filters.size()-1);
    for(uint i = 0; i<leaf_filters.size(); i++){
        //~ insert_last_leaf_trunk(i,trunk);
        insert_leaf_trunk(i,lfs-i,trunk);//TODO TO CHECK
    }
}



template <class T>
void Best<T>::insert_key(const uint64_t key,uint level){
    leaf_filters[level]->insert_key(key);
}



template <class T>
void Best<T>::add_leaf(){
    leaf_filters.push_back(new Bloom(leaf_filters_size,number_hash_function,prefix+to_string(leaf_filters.size())));
}



template <class T>
void Best<T>::optimize(){
    for(uint64_t i = 0; i <leaf_filters.size();++i){
        leaf_filters[i]->optimize();
    }
}



template <class T>
void Best<T>::optimize(uint i){
    if(i<leaf_filters.size()){
        leaf_filters[i]->optimize();
    }else{
    }
}



template <class T>
void Best<T>::dump(uint i,bm::serializer<bm::bvector<> >& bvs){
    leaf_filters[i]->dump_disk(bvs);
    leaf_filters[i]->free_ram();
}



template <class T>
void Best<T>::free_ram(){
    for(uint i=0; i<leaf_filters.size(); ++i){
        leaf_filters[i]->free_ram();
    }
}



template <class T>
uint64_t Best<T>::rcb(uint64_t min)const {
  uint64_t res(0);
  uint64_t offset(1);
  offset<<=(2*K-2);
  for(uint i(0); i<K;++i){
    res+=(3-(min%4))*offset;
    min>>=2;
    offset>>=2;
  }
  return res;
}



template <class T>
void Best<T>::update_kmer(uint64_t& min, char nuc)const {
  min<<=2;
  min+=nuc2int(nuc);
  min%=offsetUpdatekmer;
}



template <class T>
void Best<T>::update_kmer_RC(uint64_t& min, char nuc)const {
  min>>=2;
  min+=(nuc2intrc(nuc)<<(2*K-2));
}



template <class T>
void Best<T>::insert_last_leaf_trunk(uint level,ExponentialBloom<T>* EB){
    Bloom* leon(leaf_filters[level]);
    bool free_necessary = false;
    uint64_t filter_size(EB->filter.size());
    if(not leon->available){
        free_necessary=true;
        leon->load_disk();
    }
    uint64_t value = leon->BV->get_first();
    uint64_t i=0;
    for(;i<filter_size and i<value;++i){
        if(EB->filter[i]!=0){
            EB->filter[i]++;
        }
    }
    if(value<filter_size){
        EB->filter[value]++;
    }else{
        cout<<"weird hash"<<endl;
    }
    number_bit_set++;
    i++;
    do{
        value = leon->BV->get_next(value);
        if (value!=0 and i<filter_size){
            for(;i<filter_size and i<value;++i){
                if(EB->filter[i]!=0){
                    EB->filter[i]++;
                }
            }
            if(value<filter_size){
                EB->filter[value]++;
            }else{
                cout<<"weird hash"<<endl;
            }
            number_bit_set++;
            i++;
        }else{
            break;
        }
    } while(true);
    for( ;i<filter_size;++i){
        if(EB->filter[i]!=0){
            EB->filter[i]++;
        }
    }
    if(free_necessary){
        leon->free_ram();
    }
}



template <class T>
void Best<T>::insert_leaf_trunk(uint level,uint indice,ExponentialBloom<T>* EB){
    Bloom* leon(leaf_filters[level]);
    bool free_necessary = false;
    if(not leon->available){
        free_necessary=true;
        leon->load_disk();
    }
    bm::bvector<>::enumerator en = leon->BV->first();
    bm::bvector<>::enumerator en_end = leon->BV->end();
    while (en < en_end)
    {
        if(EB->filter[*en]==0){
            EB->filter[*en]=indice;
        }
        ++en;  
        number_bit_set++;
    }
    if(free_necessary){
        leon->free_ram();
    }
}



template <class T>
vector<T> Best<T>::query_key(const uint64_t key){
    vector<T> result;
    int min_level=0;
    if(trunk!=NULL){
        min_level=leaf_filters.size()-trunk->check_key(key);
    }
    int max_level=leaf_filters.size();
    if(reverse_trunk!=NULL){
        max_level=reverse_trunk->check_key(key);
    }
    for(int i=min_level;i<max_level;i++){
        if(leaf_filters[i]->check_key(key)){
            result.push_back(i);
        }
    }
    return result;
}



template <class T>
void Best<T>::serialize(zstr::ostream* out,bool hot)const{
    void* point = &(trunk->filter[0]);
    out->write((char*)point, sizeof(T)*trunk_size);
    if(reverse_trunk!=NULL){
        void* point = &(reverse_trunk->filter[0]);
        out->write((char*)point, sizeof(T)*trunk_size);
    }
    if(hot){
        bm::serializer<bm::bvector<> > bvs;
        bvs.byte_order_serialization(false);
        bvs.gap_length_serialization(false);
        for(uint i = 0; i < leaf_filters.size(); ++i){
            leaf_filters[i]->dump_disk(bvs);
        }
    }
}



template <class T>
void Best<T>::load(zstr::ifstream* out,bool hot, uint64_t leaf_number, bool double_index ){
    trunk=new ExponentialBloom<T>(trunk_size,number_hash_function);
    void* point = &(trunk->filter[0]);
    out->read((char*)point, sizeof(T)*trunk_size);
    if(double_index){
        reverse_trunk=new ExponentialBloom<T>(trunk_size,number_hash_function);
        void* point = &(reverse_trunk->filter[0]);
        out->read((char*)point, sizeof(T)*trunk_size);
    }
   
    for(uint i = 0; i < leaf_number; ++i){
        leaf_filters.push_back(new Bloom(leaf_filters_size,number_hash_function,prefix+to_string(leaf_filters.size())));
        if(hot){
            leaf_filters[i]->load_disk();
        }
    }
}

