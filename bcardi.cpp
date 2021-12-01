#include "Bloom.h"
#include "ExponentialBloom.h"
#include "utils.h"
#include "bcardi.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 



using namespace std;



void Bcardi::insert_key(const uint64_t key){
    if(current_level<big_filters.size()){
        if(current_level>0){
            if(big_filters[current_level-1]->check_key(key)){
                big_filters[current_level]->insert_key(key);
            }
        }else{
            big_filters[current_level]->insert_key(key);
        }
    }else{
        // cout<<current_level<<endl;
        uint exponent_level(current_level-big_filters.size());
        // cout<<exponent_level<<endl;
        if(exponent_level>0){
            if(little_filters->check_key(key,exponent_level)){
                little_filters->insert_key(key,exponent_level);
            }
        }else{
            if(big_filters[big_filters.size()-1]->check_key(key)){
                // cout<<"exp level    "<<exponent_level<<endl;
                little_filters->insert_key(key,exponent_level);
            }
        }
    }
}



vector<uint64_t> Bcardi::get_cardinalities()const{
    // cout<<"go get card"<<endl;
    uint64_t sum41(0);
    vector<uint64_t> result(current_level,0);
    for(int i=current_level;i>=0;--i){
        // cout<<"go level "<<i<<endl;
        uint64_t card_level;
        if(i>=(int)big_filters.size()){
            
            card_level=little_filters->get_cardinality(i-big_filters.size());
            // cout<<"C EXP    "<<intToString(card_level)<<endl;
        }else{
            card_level=big_filters[i]->get_cardinality();
            // cout<<"C BIG    "<<intToString(card_level)<<endl;
        }
        result[i]=card_level-sum41;
        sum41+=result[i];
        // cout<<"result "<<i<<"   "<<intToString(result[i])<<endl;
    }
    // cout<<"RECAP"<<endl;
    for(uint i(0);i<result.size();++i){
        cout<<i+1<<"  card:    "<<intToString(result[i])<<endl;
    }
    return result;
}





uint64_t Bcardi::rcb(uint64_t min)const {
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



void Bcardi::update_kmer(uint64_t& min, char nuc)const {
  min<<=2;
  min+=nuc2int(nuc);
  min%=offsetUpdatekmer;
}



void Bcardi::update_kmer_RC(uint64_t& min, char nuc)const {
  min>>=2;
  min+=(nuc2intrc(nuc)<<(2*K-2));
}



void Bcardi::change_level(){
    current_level++;
    // cout<<"LEVEL    "<<current_level<<endl;
}




void Bcardi::insert_sequence(const string& reference) {
  uint64_t S_kmer(str2numstrand(reference.substr(0,K-1)));//get the first kmer (k-1 bases)
  uint64_t RC_kmer(rcb(S_kmer));//The reverse complement
  for(uint i(0);i+K<reference.size();++i) {// all kmer in the genome
    update_kmer(S_kmer,reference[i+K-1]);
    update_kmer_RC(RC_kmer,reference[i+K-1]);
    uint64_t canon(min(S_kmer,RC_kmer));//Kmer min, the one selected
    // cout<<"go inesrt keu "<<i<<endl;
    insert_key(canon);
  }
}



void  Bcardi::insert_file(const string filename){
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
                // cout<<"go insert_seq"<<endl;
                insert_sequence(ref);
            }
        }
    }
}



void Bcardi::insert_file_of_file(const string filename){
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