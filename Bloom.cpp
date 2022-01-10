#include "Bloom.h"
#include "utils.h"
#include "BitMagic/src/bm.h"
#include "BitMagic/src/bmserial.h"
#include "BitMagic/src/bmundef.h"
#include <filesystem>

using namespace std;

using namespace filesystem;




void Bloom::insert_key(uint64_t key){
    for(uint64_t i=0; i<number_hash_functions;++i){
        uint64_t h=hash_family(key,i)%size;//TODO SIZE POWER OF TWO
        if((*BV)[h]==false){
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




void Bloom::dump_disk(bm::serializer<bm::bvector<> >& bvs){
    filebuf fb;
    fb.open(filename, ios::out | ios::binary | ios::trunc);
	zstr::ostream out(&fb);

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
    zstr::ifstream in(filename.c_str());
    uint64_t sz;
    in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    uint8_t* buff = new uint8_t[sz];
    in.read((char*)buff, sz);
    bm::deserialize(*(BV), buff);
    delete[] buff;
    available=true;
}



void Bloom::free_ram(){
    available=false;
    BV->reset();
}

