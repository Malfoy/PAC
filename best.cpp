#include "Bloom.h"
#include "ExponentialBloom.h"
#include "utils.h"
#include "best.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 



using namespace std;


template <class T>
void Best<T>::insert_key(const uint64_t key){
    leaf_filters[current_level]->insert_key(key);
}


template <class T>
void Best<T>::optimize(){
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
void Best<T>::insert_last_leaf_trunk(){
    Bloom* leon(leaf_filters[current_level]);
    for(uint64_t i=0; i<leon->size;++i){
        if(leon->filter[i] or trunk->filter[i]!=0){
            trunk->filter[i]++;
        }
    }
}


template <class T>
void Best<T>::change_level(){
    insert_last_leaf_trunk();
    leaf_filters[current_level]->optimize();
    current_level++;
    leaf_filters.push_back(new Bloom(leaf_filters_size,number_hash_function));
}



template <class T>
void Best<T>::insert_sequence(const string& reference) {
  uint64_t S_kmer(str2numstrand(reference.substr(0,K-1)));//get the first kmer (k-1 bases)
  uint64_t RC_kmer(rcb(S_kmer));//The reverse complement
  for(uint i(0);i+K<reference.size();++i) {// all kmer in the genome
    update_kmer(S_kmer,reference[i+K-1]);
    update_kmer_RC(RC_kmer,reference[i+K-1]);
    uint64_t canon(min(S_kmer,RC_kmer));//Kmer min, the one selected
    insert_key(canon);
  }
}


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
vector<uint> Best<T>::query_key(const uint64_t key){
    // cout<<"QUERY KEY GO"<<endl;
    vector<uint> result;
    uint level=trunk->check_key(key);
    //  cout<<"CHEY KEY OK"<<endl;
    for(int i=level;i>=0;i--){
        // cout<<leaf_filters.size()<<" "<<i<<endl;
        if(leaf_filters[i]->check_key(key)){
            result.push_back(i);
        }
    }
    //  cout<<"QUERY KEY END"<<endl;
    return result;
}



template <class T>
vector<uint> Best<T>::query_sequence(const string& reference) {
vector<uint> result(current_level,0);
vector<uint> colors;
  uint64_t S_kmer(str2numstrand(reference.substr(0,K-1)));//get the first kmer (k-1 bases)
  uint64_t RC_kmer(rcb(S_kmer));//The reverse complement
  for(uint i(0);i+K<reference.size();++i) {// all kmer in the genome
    update_kmer(S_kmer,reference[i+K-1]);
    update_kmer_RC(RC_kmer,reference[i+K-1]);
    uint64_t canon(min(S_kmer,RC_kmer));//Kmer min, the one selected
    colors=(query_key(canon));
    for(uint64_t i=0; i<colors.size();++i){
        result[colors[i]]++;
    }
  }
  return result;
}



template <class T>
void Best<T>::get_stats()const{
    vector<uint> histograms(current_level,0);
    for(uint64_t i=0; i<trunk_size;++i){
        histograms[trunk->filter[i]]++;
    }
    for (size_t i = 0; i < histograms.size(); i++)
    {
        cout<<i<<" "<<intToString(histograms[i])<<endl;
    }
    
}