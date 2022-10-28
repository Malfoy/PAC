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



// template <class T>
// void Best<T>::construct_reverse_trunk(){
//     reverse_trunk=new ExponentialBloom<T>(size,number_hash_function);
//     uint64_t lfs(leaf_filters.size());
//     #pragma omp parallel for
//     for(int i = leaf_filters.size()-1; i>=0; i--){
//         insert_leaf_trunk(i,lfs-i,reverse_trunk);        
//     }
// }



// template <class T>
// uint64_t Best<T>::construct_trunk(){
//     trunk=new ExponentialBloom<T>(size,number_hash_function);
//     uint64_t lfs(leaf_filters.size());
//     #pragma omp parallel for
//     for(uint i = 0; i<lfs; i++){
//         insert_leaf_trunk(i,lfs-i,trunk);
//     }
//     return number_bit_set;
// }



template <class T>
void Best<T>::insert_key(const uint64_t key,uint level){
    leaf_filters[level]->insert_key(key);
}



template <class T>
void Best<T>::add_leaf(){
    leaf_filters.push_back(new Bloom<T>(this));
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
    disk_space_used+=leaf_filters[i]->dump_disk(bvs,out,i);
    leaf_filters[i]->free_ram();
}



template <class T>
void Best<T>::free_ram(){
    for(uint i=0; i<leaf_filters.size(); ++i){
        leaf_filters[i]->free_ram();
    }
    if(trunk!=NULL){
        delete trunk;
        trunk=NULL;
    }
    if(reverse_trunk!=NULL){
        delete reverse_trunk;
        reverse_trunk=NULL;
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



// template <class T>
// void Best<T>::insert_leaf_trunk(uint level,uint indice,ExponentialBloom<T>* EB){
//     Bloom<T>* leon(leaf_filters[level]);
//     bm::bvector<>::enumerator en = leon->BV->first();
//     bm::bvector<>::enumerator en_end = leon->BV->end();
//     while (en < en_end){
//         if((*(EB->filter))[*en]==0){
//             (*(EB->filter))[*en]=indice;
//             number_bit_set_abt++;
//         }
//         ++en;  
//         number_bit_set++;
//     }
// }

template <class T>
void Best<T>::insert_leaf_trunk_min(uint level,ExponentialBloom<T>* EB){
    Bloom<T>* leon(leaf_filters[level]);
    bm::bvector<>::enumerator en = leon->BV->first();
    bm::bvector<>::enumerator en_end = leon->BV->end();
    while (en < en_end){
        if((*(EB->filter))[*en]==0){
            (*(EB->filter))[*en]=level;
            number_bit_set_abt++;
        }else if((*(EB->filter))[*en]>level){
            (*(EB->filter))[*en]=level;
        }
        ++en;  
        number_bit_set++;
    }
}



template <class T>
void Best<T>::insert_leaf_trunk_max(uint level,ExponentialBloom<T>* EB){
    Bloom<T>* leon(leaf_filters[level]);
    bm::bvector<>::enumerator en = leon->BV->first();
    bm::bvector<>::enumerator en_end = leon->BV->end();
    while (en < en_end){
        if((*(EB->filter))[*en]==0){
            (*(EB->filter))[*en]=level;
            number_bit_set_abt++;
        }else if((*(EB->filter))[*en]<level){
            (*(EB->filter))[*en]=level;
        }
        ++en;  
        number_bit_set++;
    }
}



template <class T>
T Best<T>::query_key_min(const uint64_t key){
    T min_level=0;
    if(trunk!=NULL){
        min_level=leaf_number-trunk->check_key(key);
    }
    return min_level;
}



template <class T>
T Best<T>::query_key_max(const uint64_t key){
    T max_level=leaf_number;
    if(reverse_trunk!=NULL){
        max_level=leaf_number-reverse_trunk->check_key(key);
    }
    return max_level;
}



template <class T>
vector<T> Best<T>::query_key(const uint64_t key){
    vector<T> result;
    uint min_level=0;
    if(trunk!=NULL){
        min_level=leaf_number-trunk->check_key(key);
    }
    uint max_level=leaf_number;
    if(reverse_trunk!=NULL){
        max_level=reverse_trunk->check_key(key);
    }
    for(uint i=min_level;i<max_level;i++){
        if(leaf_filters[i]->check_key(key)){
            result.push_back(i);
        }
    }
    return result;
}



template <class T>
void Best<T>::serialize(){
    void* point = trunk->filter->data();
    out->write((char*)point, sizeof(T)*size);
    disk_space_used+=sizeof(T)*size;
    if(reverse_trunk!=NULL){
        void* point = reverse_trunk->filter->data();
        out->write((char*)point, sizeof(T)*size);
        disk_space_used+=sizeof(T)*size;
    }
    out->flush();
}



template <class T>
void Best<T>::load_bf(uint64_t leaf_number){
    zstr::ifstream in(prefix);
    leaf_filters.clear();
    for(uint i = 0; i < leaf_number; ++i){
        leaf_filters.push_back(new Bloom<T>(this));
    }
    for(uint i = 0; i < leaf_number; ++i){
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        leaf_filters[indiceBloom]->load_disk(&in);
    }
}



template <class T>
void Best<T>::load(uint64_t leaf_number, bool double_index ){
    cout<<"LOAD"<<endl;
    leaf_filters.clear();
    for(uint i = 0; i < leaf_number; ++i){
        leaf_filters.push_back(new Bloom<T>(this));
    }
    zstr::ifstream in(prefix);
    for(uint i = 0; i < leaf_number; ++i){
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        leaf_filters[indiceBloom]->load_disk(&in);
    }
    cout<<"LOAD BF OK "<<endl;
    trunk=new ExponentialBloom<T>(size,number_hash_function);
    void* point = (trunk->filter->data());
    in.read((char*)point, sizeof(T)*size);
    cout<<"LOAD ABF OK"<<endl;
    if(double_index){
        reverse_trunk=new ExponentialBloom<T>(size,number_hash_function);
        void* point = reverse_trunk->filter->data();
        in.read((char*)point, sizeof(T)*size);
    }
}



template <class T>
void Best<T>::load(uint64_t leaf_number, bool double_index, vector<bool>* BBV ){
    cout<<"LOAD"<<endl;
    cout<<leaf_number<<endl;
    cout<<leaf_filters.size()<<endl;
    leaf_filters.clear();
    cout<<"LOAD"<<endl;
    Bloom<T>* leon(new Bloom<T>(this));
    zstr::ifstream in(prefix);
    bm::bvector<>::enumerator en;
    bm::bvector<>::enumerator en_end;
    for(uint i = 0; i < leaf_number; ++i){
        uint32_t indiceBloom;
        in.read(reinterpret_cast<char*>(&indiceBloom), sizeof(indiceBloom));
        leon->load_disk(&in);
        en = leon->BV->first();
        en_end = leon->BV->end();
        while (en < en_end){
            (*BBV)[(*en)*leaf_number+indiceBloom]=true;
            ++en;  
        }
        leon->BV->clear();
    }
    delete leon;
    trunk=new ExponentialBloom<T>(size,number_hash_function);
    void* point = (trunk->filter->data());
    in.read((char*)point, sizeof(T)*size);
    if(double_index){
        reverse_trunk=new ExponentialBloom<T>(size,number_hash_function);
        void* point = reverse_trunk->filter->data();
        in.read((char*)point, sizeof(T)*size);
    }
}

