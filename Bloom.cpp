#include "Bloom.h"
#include "utils.h"
#include "BitMagic/src/bm.h"
#include "BitMagic/src/bmserial.h"
#include "BitMagic/src/bmundef.h"
#include <filesystem>



using namespace std;
using namespace filesystem;


template <class T>
 Bloom<T>::Bloom(Best<T>* Ifather){
    father=Ifather;
    BV=new bm::bvector<>(father->size+1,bm::BM_GAP);
    
}


template <class T>
void Bloom<T>::insert_key(uint64_t key){
    for(uint64_t i=0; i<father->number_hash_function;++i){
        uint64_t h=(hash_family(key,i))&(father->size);
        (*BV)[h]=true;
    }
}


template <class T>
bool Bloom<T>::check_key(uint64_t key){
    for(uint64_t i=0; i<father->number_hash_function;++i){
        uint64_t h=hash_family(key,i)&father->size;
        if((*BV)[h]==false){
            return false;
        }
    }
    return true;
}



template <class T>
uint64_t Bloom<T>::dump_disk(bm::serializer<bm::bvector<> >& bvs, zstr::ofstream* out,uint32_t i){
	bm::serializer<bm::bvector<> >::buffer sbuf;
    unsigned char* buf = 0;
    bvs.serialize(*(BV), sbuf);
    buf         = sbuf.data();
    uint64_t sz = sbuf.size();
    auto point2 = &buf[0];
    out->write(reinterpret_cast<const char*>(&i), sizeof(i));
    out->write(reinterpret_cast<const char*>(&sz), sizeof(sz));
    out->write((char*)point2, sz);
    out->flush();
    return sz;

}


template <class T>
void Bloom<T>::load_disk(zstr::ifstream* in){
    uint64_t sz;
    in->read(reinterpret_cast<char*>(&sz), sizeof(sz));
    uint8_t* buff = new uint8_t[sz];
    in->read((char*)buff, sz);
    bm::deserialize(*(BV), buff);
    delete[] buff;
}


template <class T>
void Bloom<T>::free_ram(){
    BV->reset();
}

