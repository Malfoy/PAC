#ifndef UTILH
#define UTILH

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath> 
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include "strict_fstream.hpp"
#include "zstr.hpp"



using namespace std;



uint64_t nuc2int(char c);

string intToString(uint64_t n);


uint64_t nuc2intrc(char c) ;



uint64_t str2numstrand(const string& str) ;



uint64_t xorshift64(uint64_t& state);



uint64_t hash64shift(uint64_t key) ;



uint64_t hash_family(uint64_t& key,uint b);



void Biogetline(zstr::ifstream* in,string& result,char type,uint K) ;



char get_data_type(const string& filename);



bool exists_test (const std::string& name);

#endif