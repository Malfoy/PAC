#include "Bloom.h"
#include "ExponentialBloom.h"
#include "utils.h"
#include "bestpart.h"
#include "best.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include <string>



using namespace std;



uint64_t str2num(const string& str) {
	uint64_t res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}


uint64_t revhash(uint64_t x) {
	// return hash64shift(x);
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x);
	return x;
}



uint64_t unrevhash(uint64_t x) {
	return hash64shift(x);
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}


template <class T>
void BestPart<T>::insert_keys( vector<uint64_t>& key,uint minimizer){
    omp_set_lock(&mutex_array[minimizer]);
    for(uint i=0;i<key.size(); ++i){
        buckets[minimizer]->insert_key(key[i]);
    }
    omp_unset_lock(&mutex_array[minimizer]);
    key.clear();
}


//TODO CAN BE IMPROVED?
template <class T>
uint64_t BestPart<T>::rcb(uint64_t min, uint64_t n)const {
  uint64_t res(0);
  uint64_t offset(1);
  offset<<=(2*n-2);
  for(uint i(0); i<n;++i){
    res+=(3-(min%4))*offset;
    min>>=2;
    offset>>=2;
  }
  return res;
}


template <class T>
uint64_t BestPart<T>::canonize(uint64_t x, uint64_t n) {
	return min(x, rcb(x, n));
}




template <class T>
uint64_t BestPart<T>::regular_minimizer_pos(uint64_t seq, uint64_t& position) {
	uint64_t mini, mmer;
	mmer = seq % large_minimizer_number;
	mini = mmer        = canonize(mmer, large_minimizer_size);
	uint64_t hash_mini = (unrevhash(mmer));
	position           = 0;
	for (uint64_t i(1); i <= K - large_minimizer_size; i++) {
		seq >>= 2;
		mmer          = seq % large_minimizer_number;
		mmer          = canonize(mmer, large_minimizer_size);
		uint64_t hash = (unrevhash(mmer));
		if (hash_mini > hash) {
			position  = K - large_minimizer_size - i;
			mini      = mmer;
			hash_mini = hash;
		}
	}
	return mini;
}



template <class T>
void BestPart<T>::updateK(uint64_t& min, char nuc)const {
  min<<=2;
  min+=nuc2int(nuc);
  min%=offsetUpdatekmer;
}


template <class T>
void BestPart<T>::updateM(uint64_t& min, char nuc)const {
  min<<=2;
  min+=nuc2int(nuc);
  min%=offsetUpdateminimizer;
}



template <class T>
void BestPart<T>::updateRCK(uint64_t& min, char nuc)const {
  min>>=2;
  min+=(nuc2intrc(nuc)<<(2*K-2));
}



template <class T>
void BestPart<T>::updateRCM(uint64_t& min, char nuc)const {
  min>>=2;
  min+=(nuc2intrc(nuc)<<(2*large_minimizer_size-2));
}



template <class T>
void BestPart<T>::change_level(){
    for(uint32_t i=0;i<bucket_number;++i){
        buckets[i]->change_level();
    }
    current_level++;
}



template <class T>
void BestPart<T>::insert_sequence(const string& ref) {
    uint64_t old_minimizer, minimizer;
    vector<uint64_t> superkmer;
    old_minimizer = minimizer = large_minimizer_number;
    uint64_t last_position(0);
    // FOREACH KMER
    uint64_t seq(str2num(ref.substr(0, K)));
    uint64_t rcseq=rcb(seq,K);
    uint64_t canon=min(seq,rcseq);
    uint64_t position_min;
    uint64_t min_seq = (str2num(ref.substr(K - large_minimizer_size, large_minimizer_size)));
    uint64_t min_rcseq(rcb(min_seq, large_minimizer_size));
    uint64_t min_canon(min(min_seq, min_rcseq));
    minimizer         = regular_minimizer_pos(seq, position_min);
    old_minimizer     = minimizer;
    uint64_t hash_min = unrevhash(minimizer);
    uint64_t i(0);
    for (; i + K < ref.size(); ++i) {
        updateK(seq, ref[i + K]);
        updateRCK(rcseq, ref[i + K]);
        canon=min(seq,rcseq);
        updateM(min_seq, ref[i + K]);
        updateRCM(min_rcseq, ref[i + K]);
        min_canon      = (min(min_seq, min_rcseq));
        uint64_t new_h = unrevhash(min_canon);
        // THE NEW mmer is a MINIMIZER
        if (new_h < hash_min) {
            minimizer    = (min_canon);
            hash_min     = new_h;
            position_min = i + K - large_minimizer_size + 1;
        } else {
            // the previous minimizer is outdated
            if (i >= position_min) {
                minimizer = regular_minimizer_pos(seq, position_min);
                hash_min  = unrevhash(minimizer);
                position_min += (i + 1);
            } else {
            }
        }
        // COMPUTE KMER MINIMIZER
        if (revhash(old_minimizer) % bucket_number != revhash(minimizer) % bucket_number) {
            old_minimizer = (revhash(old_minimizer) % bucket_number);
            insert_keys(superkmer,old_minimizer);
            last_position = i + 1;
            old_minimizer = minimizer;
        }
        superkmer.push_back(canon);
    }
    if (ref.size() - last_position > K - 1) {
        old_minimizer = (revhash(old_minimizer) % bucket_number);
        insert_keys(superkmer,old_minimizer);
    }
}



template <class T>
void  BestPart<T>::insert_file(const string filename){
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
void BestPart<T>::insert_file_of_file(const string filename){
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


