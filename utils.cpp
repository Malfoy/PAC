#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include <sys/types.h>
#include <sys/stat.h>
#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "utils.h"
#include <sys/time.h>
#include <sys/resource.h>



uint64_t nuc2int(char c) {
  switch(c){
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return 0;
  }
  exit(0);
  return 0;
}



uint64_t nuc2intrc(char c) {
  switch(c) {
    case 'A': return 3;
    case 'C': return 2;
    case 'G': return 1;
    default: return 0;
  }
  exit(0);
  return 0;
}



uint64_t str2numstrand(const string& str) {
  uint64_t res(0);
  for(uint i(0);i<str.size();i++) {
    res<<=2;
    switch (str[i]){
      case 'A':res+=0;break;
      case 'C':res+=1;break;
      case 'G':res+=2;break;
      case 'T':res+=3;break;

      case 'a':res+=0;break;
      case 'c':res+=1;break;
      case 'g':res+=2;break;
      case 't':res+=3;break;
      default: return 0 ;break;
    }
  }
  return (res);
}



//pretty printing
string intToString(uint64_t n) {
  if (n < 1000) {
    return to_string(n);
  }
  string end(to_string(n % 1000));
  if (end.size() == 3) {
    return intToString(n / 1000) + "," + end;
  }
  if (end.size() == 2) {
    return intToString(n / 1000) + ",0" + end;
  }
  return intToString(n / 1000) + ",00" + end;
}



uint64_t xorshift64(uint64_t& state){
	uint64_t x = state;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	return state = x;
}



uint64_t hash64shift(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}



uint64_t getMemorySelfMaxUsed (){
	uint64_t result = 0;
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage)==0)  {  result = usage.ru_maxrss;  }
	return result;
}



uint64_t hash_family(uint64_t key,uint b){
    return hash64shift(key) + xorshift64(key)*b;
}



void Biogetline(zstr::ifstream* in,string& result,char type,uint K) {
  string discard;
  result.clear();
  switch(type){
    case 'Q':
      getline(*in,discard);
      getline(*in,result);
      getline(*in,discard);
      getline(*in,discard);
      break;
    case 'A':
      getline(*in,discard);
      char c=in->peek();
      while(c!='>' and c!=EOF){
        getline(*in,discard);
        result+=discard;
        c=in->peek();
      }
      break;
  }
  if(result.size()< K){
    result.clear();
  }
}



void Biogetline(zstr::ifstream* in,string& result,char type,uint K,string& header) {
  string discard;
  result.clear();
  switch(type){
    case 'Q':
      getline(*in,header);
      getline(*in,result);
      getline(*in,discard);
      getline(*in,discard);
      break;
    case 'A':
      getline(*in,header);
      char c=in->peek();
      while(c!='>' and c!=EOF){
        getline(*in,discard);
        result+=discard;
        c=in->peek();
      }
      break;
  }
  if(result.size()< K){
    result.clear();
    header.clear();
  }
}



uint64_t approx_power2(uint64_t n){
    uint64_t lower(1<<asm_log2(n));
    uint64_t upper((uint64_t)1<<(uint64_t)ceil(log2(n)));
    if(lower!=n){
        cout<<intToString(n)<<" is not a power of two, I will use "<<intToString(upper)<<" instead"<<endl;
    }
    return upper;
}



uint64_t asm_log2(const uint64_t x) {
	uint64_t y;
	asm("\tbsr %1, %0\n" : "=r"(y) : "r"(x));
	return y;
}



char get_data_type(const string& filename){
  if(filename.find(".fq")!=string::npos){
    return 'Q';
  }
  if(filename.find(".fastq")!=string::npos){
    return 'Q';
  }
  return 'A';
}



bool exists_test (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}



int directory_exists(string& path)
{
    struct stat info;

    int statRC = stat( path.c_str(), &info );
    if( statRC != 0 )
    {
        if (errno == ENOENT)  { return 0; } // something along the path does not exist
        if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
        return -1;
    }

    return ( info.st_mode & S_IFDIR ) ? 1 : 0;
}
