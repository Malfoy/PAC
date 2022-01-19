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
#include <chrono>
#include <ctime>
#include <filesystem>



using namespace std;
using namespace filesystem;



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
void BestPart<T>::insert_keys(const vector<uint64_t>& key,uint minimizer,uint level,Bloom<T>* unique_filter){
    if(filter){
        for(uint i=0;i<key.size(); ++i){
            if(unique_filter->check_key(key[i])){
                buckets[minimizer]->insert_key(key[i],level);
            }else{
                unique_filter->insert_key(key[i]);
            }
        }
    }else{
        for(uint i=0;i<key.size(); ++i){
            buckets[minimizer]->insert_key(key[i],level);
        }
    }
}



//TODO CAN BE IMPROVED
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
}



template <class T>
vector<pair<vector<uint64_t>,uint64_t> > BestPart<T>::get_super_kmers(const string& ref) {
    vector<pair<vector<uint64_t>,uint64_t> > result;
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
            result.push_back({superkmer,old_minimizer});
            superkmer.clear();
            last_position = i + 1;
            old_minimizer = minimizer;
        }
        superkmer.push_back(canon);
    }
    if (ref.size() - last_position > K - 1) {
        old_minimizer = (revhash(old_minimizer) % bucket_number);
        result.push_back({superkmer,old_minimizer});
        superkmer.clear();
    }
    return result;
}



template <class T>
void BestPart<T>::insert_sequence(const string& reference,uint level, Bloom<T>* unique_filter){
    auto V(get_super_kmers(reference));
    for(uint32_t i = 0; i < V.size();++i){
        insert_keys(V[i].first, V[i].second,level,unique_filter);
    }
}



template <class T>
void BestPart<T>::query_file(const string& filename, const string& output){
    cout<<"Query file "<<filename<<" in "<<output<<endl;
    auto  start = chrono::system_clock::now();
    zstr::ofstream out(output,ios::trunc);
    char type=get_data_type(filename);
    zstr::ifstream in(filename);
    path old(current_path());
    vector<uint32_t> colors_count(leaf_number,0);
    current_path(w_dir);
    vector< pair<string,vector<uint32_t> > > result;
    vector<vector<pair<uint64_t,uint32_t> > > colored_kmer_per_bucket(bucket_number);
    #pragma omp parallel num_threads(core_number)
    {
        string ref,header; 
        uint32_t query_id;
        while(not in.eof()) {
            #pragma omp critical (inputfile)
            {
                Biogetline(&in,ref,type,K,header);
                if(ref.size()>K){
                    result.push_back({header,colors_count});
                    query_id=result.size()-1;
                }
            }
            if(ref.size()>K) {
                load_super_kmer(colored_kmer_per_bucket,query_id,ref);
            }
        }
    }
    uint64_t sum_query=0;
    for(uint i=0;i<bucket_number;i++) {
        if(not colored_kmer_per_bucket[i].empty()){
            sum_query+=query_bucket(colored_kmer_per_bucket[i],result,i);
        }
    }
    for(uint i=0;i<result.size();i++) {
        {   
            out<<result[i].first<<"\n";
            for(uint j=0;j<result[i].second.size();j++){
                out<<result[i].second[j]<<' ';
            }
            out<<"\n";
        }
    }
    out.close();
    auto end = std::chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout <<  "Query time: " << elapsed_seconds.count() << "s for "<<intToString(sum_query)<<
    " kmer query or "<<intToString((double)sum_query/elapsed_seconds.count()) << " query per second  \n";

    current_path(old);
}



template <class T>
void BestPart<T>::load_super_kmer(vector<vector<pair<uint64_t,uint32_t> > >&colored_kmer_per_bucket,uint32_t query_id, const string& reference){
    auto V(get_super_kmers(reference));
    for(uint32_t i = 0; i < V.size();++i){
        omp_set_lock(&mutex_array[V[i].second]);
        for(uint32_t j = 0; j < V[i].first.size();++j){
            (colored_kmer_per_bucket[V[i].second]).push_back({V[i].first[j],query_id});
        }
        omp_unset_lock(&mutex_array[V[i].second]);
    }
}



template <class T>
uint64_t BestPart<T>::query_bucket(vector<pair<uint64_t,uint32_t> >& colored_kmer, vector< pair<string,vector<uint32_t> > >& result, uint bucket_id){
    cout<<'-'<<flush;
    uint64_t sum(0);
    buckets[bucket_id]->load(leaf_number,use_double_index);
    vector<T> minV(colored_kmer.size(),0);
    vector<T> maxV(colored_kmer.size(),leaf_number-1);
    #pragma omp parallel for num_threads(core_number)
    for(uint32_t i = 0; i < colored_kmer.size();++i){
        minV[i]=buckets[bucket_id]->query_key_min(colored_kmer[i].first);
    }
    if(use_double_index){
        maxV.resize(colored_kmer.size(),0);
        #pragma omp parallel for num_threads(core_number)
        for(uint32_t i = 0; i < colored_kmer.size();++i){
            maxV[i]=buckets[bucket_id]->query_key_max(colored_kmer[i].first);
        }
    }
    
    for(uint32_t l=(0);l<leaf_number;++l){
    #pragma omp parallel for num_threads(core_number)
        for(uint32_t i = 0; i < colored_kmer.size();++i){
            if(l>=minV[i] and l<=maxV[i]){
                #pragma omp atomic
                sum++;
                if(buckets[bucket_id]->leaf_filters[l]->check_key(colored_kmer[i].first)){
                    #pragma omp atomic
                    result[colored_kmer[i].second].second[l]++;
                }
            }
        }
    }
    
    delete buckets[bucket_id];
    return sum;
}





//~ template <class T>
//~ uint64_t BestPart<T>::query_bucket(vector<pair<uint64_t,uint32_t> >& colored_kmer, vector< pair<string,vector<uint32_t> > >& result, uint bucket_id){

    //~ uint64_t sum(0);
    //~ cout<<bucket_id<<endl;
    //~ buckets[bucket_id]->load(leaf_number,use_double_index);
    //~ #pragma omp parallel num_threads(core_number)
    //~ {
        //~ vector<T> colors;
        //~ #pragma omp for
        //~ for(uint32_t i = 0; i < colored_kmer.size();++i){
            //~ colors=buckets[bucket_id]->query_key(colored_kmer[i].first);
            //~ #pragma omp atomic
            //~ sum++;
            //~ for(uint32_t j = 0; j <colors.size();++j){
                //~ #pragma omp atomic
                //~ result[colored_kmer[i].second].second[colors[j]]++;
            //~ }
        //~ }
    //~ }
    //~ delete buckets[bucket_id];
    //~ return sum;
//~ }



template <class T>
vector<T> BestPart<T>::query_keys(const vector<uint64_t>& key,uint minimizer){
    vector<T> result;
    vector<T> colors;
    for(uint i=0;i<key.size(); ++i){
        colors=buckets[minimizer]->query_key(key[i]);
        result.insert(result.end(), colors.begin(), colors.end() );
        colors.clear();
    }
    return result;
}



template <class T>
void  BestPart<T>::insert_file(const string& filename, uint level, uint32_t indice_bloom){
    char type=get_data_type(filename);
    zstr::ifstream in(filename);
    Bloom<T>* unique_filter;
    if(filter){
        unique_filter=new Bloom<T>(buckets[0]);
    }
    {
        string ref;
        while(not in.eof()) {
            {
                Biogetline(&in,ref,type,K);
            }
            if(ref.size()>K) {
                insert_sequence(ref,level,unique_filter);
            }
        }
    }
    if(filter){
        delete unique_filter;
    }
    bm::serializer<bm::bvector<> > bvs;
	bvs.byte_order_serialization(false);
	bvs.gap_length_serialization(false);
    for(uint i=0;i<buckets.size();++i){
        buckets[i]->optimize(level);
        omp_set_lock(&mutex_array[i]);// TODO THE LOCK SHOULD NOT INCLUDE THE SERIALIZATION
        buckets[i]->out->write(reinterpret_cast<const char*>(&indice_bloom), sizeof(indice_bloom));
        buckets[i]->dump(level,bvs);
        omp_unset_lock(&mutex_array[i]);
    }
}



template <class T>
void BestPart<T>::insert_file_of_file(const string& filename){
    cout<<"I index "<<K<<"mers with Bloom filters of size " <<intToString(size)<<" with "<<number_hash_function<<" hash functions  using "<<intToString(1<<(2*small_minimizer_size))<<" partitions "<<endl;

    path old(current_path());
    path vodka(absolute(filename));
    zstr::ifstream in(vodka);
    cout<<"Insert file of file "<<filename<<endl;
    path initial_path=current_path();
    for(uint i(0);i<core_number;++i){
        add_leaf();
    }
    
    auto  start = chrono::system_clock::now();;
    #pragma omp parallel num_threads (core_number)
    {
        string ref;
        uint level;
        bool go;
        int idt=omp_get_thread_num ();
        //~ cout<<idt<<endl;
        while(not in.eof()) {
            #pragma omp critical (inputfile)
            {
                getline(in,ref);
                if(exists_test(ref)){
                    level=leaf_number;
                    leaf_number++;
                    //~ add_leaf();
                    go=true;
                }else{
                    go=false;
                }
            }
            if(go){
                insert_file(ref,idt,level);
                if(level%100==0){
                    cout<<level<<" files"<<endl;
                }
            }
        }
    }
    auto  middle = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = middle - start;
    cout <<  "Bloom construction time: " << elapsed_seconds.count() << "s\n";
    cout<<intToString(getMemorySelfMaxUsed())<<" KB total"<<endl;
    cout<<intToString(getMemorySelfMaxUsed()/(leaf_number))<<" KB per file ("<<leaf_number<<" files)"<<endl;
    index();
    nb_bit_set();
    auto  end = chrono::system_clock::now();
    elapsed_seconds = end - middle;
    cout <<  "Exponential Bloom construction time: " << elapsed_seconds.count() << "s\n";
    elapsed_seconds = end - start;
    cout <<  "Total Index time: " << elapsed_seconds.count() << "s\n";
    cout <<  "Index time per file: " << elapsed_seconds.count()/leaf_number << "s\n";
    cout<<intToString(getMemorySelfMaxUsed())<<" kB total"<<endl;
    cout<<intToString(getMemorySelfMaxUsed()/(leaf_number))<<" kB per file"<<endl;
    current_path(old);
}



template <class T>
void BestPart<T>::serialize()const{
    cout<<"Dump the index in "<<w_dir<<endl;
    path initial_path=current_path();
    path p {w_dir};
    if(exists(p)){
        
    }else{
        create_directory(p);
    }
    current_path(p);
    string filename("MainIndex");
    filebuf fb;
    fb.open(filename, ios::out | ios::binary | ios::trunc);
	zstr::ostream out(&fb);
    // VARIOUS INTEGERS
	out.write(reinterpret_cast<const char*>(&K), sizeof(K));
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    out.write(reinterpret_cast<const char*>(&number_hash_function), sizeof(number_hash_function));
    out.write(reinterpret_cast<const char*>(&small_minimizer_size), sizeof(small_minimizer_size));
    out.write(reinterpret_cast<const char*>(&large_minimizer_size), sizeof(large_minimizer_size));
    out.write(reinterpret_cast<const char*>(&leaf_number), sizeof(leaf_number));
    out.write(reinterpret_cast<const char*>(&use_double_index), sizeof(use_double_index));
    out<<flush;
    fb.close();
    current_path(initial_path);
    cout<<"Index use "<<intToString(disk_space_used)<<" Bytes on disk"<<endl;
}



template <class T>
void BestPart<T>::index(){
    #pragma omp parallel for num_threads(core_number)
    for(uint i=0;i<buckets.size();++i){
        buckets[i]->load_bf(leaf_number);
        uint64_t nbbs( buckets[i]->construct_trunk());
        if(use_double_index){
            buckets[i]->construct_reverse_trunk();
        }
        #pragma omp atomic
        number_bit_set+=nbbs;
        buckets[i]->serialize();
        nbbs=buckets[i]->disk_space_used;
        #pragma omp atomic
        disk_space_used+=nbbs;
        nbbs=buckets[i]->number_bit_set_abt;
        #pragma omp atomic
        number_bit_set_abt+=nbbs;
        delete buckets[i];
        buckets[i]=NULL;
    }
}



template <class T>
void BestPart<T>::add_leaf(){
    //~ #pragma omp parallel for
    for(uint i=0;i<buckets.size();++i){
        buckets[i]->add_leaf();
    }
}



template <class T>
uint64_t BestPart<T>::nb_bit_set()const{
    cout<<intToString(number_bit_set)<<" bits sets"<<endl;
    cout<<(double)100*number_bit_set/((double)leaf_number*size)<<"% of bits are set"<<endl;
    cout<<intToString(number_bit_set_abt)<<" non empty cell in ABT"<<endl;
    cout<<(double)100*number_bit_set_abt/(size)<<"% of cells are used"<<endl;
    return number_bit_set;
}


