#include "Bloom.h"
#include "utils.h"
#include "BitMagic/src/bm.h"
#include "BitMagic/src/bmserial.h"
#include "BitMagic/src/bmundef.h"



void Bloom::insert_key(uint64_t key){
    for(uint64_t i=0; i<number_hash_functions;++i){
        uint64_t h=hash_family(key,i)%size;//TODO SIZE POWER OF TWO
        if((*BV)[h]==false){
            number_bit_set++;
            (*BV)[h]=true;
        }
    }
}



bool Bloom::check_key(uint64_t key){
    if(not available){
        load_disk();
    }
    for(uint64_t i=0; i<number_hash_functions;++i){
        uint64_t h=hash_family(key,i)%size;//TODO SIZE POWER OF TWO
        if((*BV)[h]==false){
            return false;
        }
    }
    return true;
}


void Bloom::dump_disk(){
    ofstream out(filename,iostream::trunc);
    bm::serializer<bm::bvector<> > bvs;
	bvs.byte_order_serialization(false);
	bvs.gap_length_serialization(false);
	bm::serializer<bm::bvector<> >::buffer sbuf;
    unsigned char* buf = 0;
    bvs.serialize(*(BV), sbuf);
    buf         = sbuf.data();
    uint64_t sz = sbuf.size();
    auto point2 = &buf[0];
    out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
    out.write((char*)point2, sz);

}



void Bloom::load_disk(){
    // cout<<filename<<endl;
    ifstream in(filename.c_str());
    uint64_t sz;
    in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    // cout<<sz<<endl;
    uint8_t* buff = new uint8_t[sz];
    in.read((char*)buff, sz);
    // cout<<"read ok"<<endl;
    bm::deserialize(*(BV), buff);
    // cout<<"deserialize_structure"<<endl;
    delete[] buff;
    available=true;
}

void Bloom::free_ram(){
    available=false;
    BV->clear();
}

