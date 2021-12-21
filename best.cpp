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



using namespace std;


template <class T>
void Best<T>::construct_reverse_trunk(){
    reverse_trunk=new ExponentialBloom<T>(trunk_size,number_hash_function);
    for(int i = leaf_filters.size()-1; i>=0; i--){
        insert_last_leaf_trunk(i,reverse_trunk);
    }
}


template <class T>
void Best<T>::construct_trunk(){
    trunk=new ExponentialBloom<T>(trunk_size,number_hash_function);
    for(uint i = 0; i<leaf_filters.size(); i++){
        insert_last_leaf_trunk(i,trunk);
    }
}



template <class T>
void Best<T>::insert_key(const uint64_t key,uint level){
    if(level>=leaf_filters.size()){ 
        for(uint i(leaf_filters.size()-1);i<level;i++){
            leaf_filters.push_back(new Bloom(leaf_filters_size,number_hash_function));
        }
    }
    leaf_filters[level]->insert_key(key);
    nb_insertions++;
}



template <class T>
void Best<T>::optimize(){
    #pragma omp parallel for
    for(uint64_t i = 0; i <leaf_filters.size();++i){
        leaf_filters[i]->optimize();
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
    auto value = leon->filter.get_first();
    uint64_t i=0;
    for(;i<leon->size and i<value;++i){
        if(EB->filter[i]!=0){
            EB->filter[i]++;
        }
    }
    EB->filter[value]++;
    i++;
    do{
        value = leon->filter.get_next(value);
        if (value){
            for(;i<leon->size and i<value;++i){
                if(EB->filter[i]!=0){
                    EB->filter[i]++;
                }
            }
            EB->filter[value]++;
            i++;
        }else{
            break;
        }
    } while(true);
    for( ;i<leon->size;++i){
        if(EB->filter[i]!=0){
            EB->filter[i]++;
        }
    }
}






// template <class T>
// void Best<T>::insert_sequence(const string& reference) {
//   uint64_t S_kmer(str2numstrand(reference.substr(0,K-1)));//get the first kmer (k-1 bases)
//   uint64_t RC_kmer(rcb(S_kmer));//The reverse complement
//   for(uint i(0);i+K<reference.size();++i) {// all kmer in the genome
//     update_kmer(S_kmer,reference[i+K-1]);
//     update_kmer_RC(RC_kmer,reference[i+K-1]);
//     uint64_t canon(min(S_kmer,RC_kmer));//Kmer min, the one selected
//     insert_key(canon);
//   }
// }



template <class T>
void  Best<T>::insert_file(const string filename){
    char type=get_data_type(filename);
    zstr::ifstream in(filename);
    #pragma omp parallel
    {
        string ref;
        while(not in.eof()) {
            #pragma omp critical (inputfile)
            {
             Biogetline(&in,ref,type,K);
            }
            if(ref.size()>K) {
                insert_sequence(ref);
            }
        }
    }
    
}



template <class T>
void Best<T>::insert_file_of_file(const string filename){
    zstr::ifstream in(filename);
    string ref;
    while(not in.eof()) {
        getline(in,ref);
        if(exists_test(ref)) {
            cout<<"go insert_file"<<endl;
            insert_file(ref);
            change_level();
        }
    }
}



template <class T>
vector<T> Best<T>::query_key(const uint64_t key){
    // cout<<"QUERY KEY GO"<<endl;
    vector<T> result;
    int min_level=0;
    if(trunk!=NULL){
        min_level=leaf_filters.size()-trunk->check_key(key);
    }
    int max_level=leaf_filters.size();
    if(reverse_trunk!=NULL){
        max_level=reverse_trunk->check_key(key);
    }
    //  cout<<"CHEY KEY OK"<<endl;
    for(int i=min_level;i<max_level;i++){
        // cout<<leaf_filters.size()<<" "<<i<<endl;
        if(leaf_filters[i]->check_key(key)){
            result.push_back(i);
        }
    }
    //  cout<<"QUERY KEY END"<<endl;
    return result;
}



template <class T>
void Best<T>::serialize(zstr::ostream* out)const{
    void* point = &(trunk->filter[0]);
    // cout<<trunk_size<<endl;
    out->write((char*)point, sizeof(T)*trunk_size);
    bm::serializer<bm::bvector<> > bvs;
	bvs.byte_order_serialization(false);
	bvs.gap_length_serialization(false);
	bm::serializer<bm::bvector<> >::buffer sbuf;
    for(uint i = 0; i < leaf_filters.size(); ++i){
        unsigned char* buf = 0;
		bvs.serialize(leaf_filters[i]->filter, sbuf);
		buf         = sbuf.data();
		uint64_t sz = sbuf.size();
        // cout<<"size buff"<<sz<<endl;
		auto point2 = &buf[0];
		out->write(reinterpret_cast<const char*>(&sz), sizeof(sz));
		out->write((char*)point2, sz);
    }
}


template <class T>
void Best<T>::load(zstr::ifstream* out){
    void* point = &(trunk->filter[0]);
    // cout<<trunk_size<<endl;
    out->read((char*)point, sizeof(T)*trunk_size);
    // cout<<"read trunk"<<endl;
    for(uint i = 0; i < leaf_filters.size(); ++i){
        // cout<<"go leaf"<<endl;
        uint64_t sz;
        out->read(reinterpret_cast<char*>(&sz), sizeof(sz));
        // cout<<sz<<endl;
        uint8_t* buff = new uint8_t[sz];
        out->read((char*)buff, sz);
        //  cout<<"read ok"<<endl;
        bm::deserialize(leaf_filters[i]->filter, buff);
        delete[] buff;
    }
}






template <class T>
void Best<T>::get_stats()const{
    vector<uint> histograms(leaf_filters.size(),0);
    for(uint64_t i=0; i<trunk_size;++i){
        histograms[trunk->filter[i]]++;
    }
    for (size_t i = 0; i < histograms.size(); i++)
    {
        cout<<i<<" "<<intToString(histograms[i])<<endl;
    }
    
}