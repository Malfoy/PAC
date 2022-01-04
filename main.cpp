#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <getopt.h>
#include "bcardi.h"
#include "best.h"
#include "bestpart.h"



using namespace std;
using namespace chrono;



string fof(""), dump_index(""), query_file(""),existing_index(""),query_output("output.gz");
uint16_t nb_hash_func(1), bit_encoding(16);
uint32_t kmer_size(31);
uint64_t bf_size(100000000); 
bool use_double_index(false);


void PrintHelp()
{
	cout <<
			"\n******************* ABC **************************************\n"

			"--help (-h)              :     Show help\n\n"

			"\n INDEX CONSTRUCTION\n"
            "--fof (-f)               :     Build index from file of files (FASTA/Q/GZ allowed in file of files)\n\n"
            "--load (-l)              :     Load index from file \n\n"
			"--dump (-d)              :     Write index in file\n"
			
			"\n INDEX QUERY\n"
			"--query (-q)             :     Query sequence file (FASTA/Q/GZ)\n"
            "--out (-o)               :     Write query output in file (default: output.gz)\n"
			
			"\n TWEAK PARAMETERS\n"
			"-k                       :     k-mer size (default: " << kmer_size << ")\n"
			"-b                       :     Bloom filter size (default " << intToString(bf_size) << ")\n"
			"-n                       :     Bloom filter's number of hash functions (default: " << intToString( nb_hash_func) << ")\n"
			"-e                       :     Bit encoding, possible values are 8 (max 256 files), 16 (max 16,384 files), 32 (Max 4,294,967,296 files) (default: " << bit_encoding << ")\n"
            "-i                       :     Build double index for faster queries \n";

	exit(1);
}




void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "k:d:q:b:f:e:l:n:hi";
	const option long_opts[] = 
	{
		{"index", no_argument, nullptr, 'i'},
		{"fof", required_argument, nullptr, 'f'},
		{"dump", required_argument, nullptr, 'd'},
        {"out", required_argument, nullptr, 'o'},
		{"help", no_argument, nullptr, 'h'},
		{"query", required_argument, nullptr, 'q'},
		{"load", required_argument, nullptr, 'l'},
		{nullptr, no_argument, nullptr, 0}
	};
	while (true)
	{
		const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

		if (-1 == opt)
			break;

		switch (opt)
		{
            case 'i':
				use_double_index=true;
				break;
			case 'k':
				kmer_size=stoi(optarg);
				break;
			case 'd':
				dump_index=optarg;
				break;
			case 'f':
				fof=optarg;
				break;
			case 'q':
				query_file=optarg;
				break;
            case 'o':
				query_output=optarg;
				break;
			case 'l':
				existing_index=optarg;
				break;
			case 'b':
				bf_size=stoi(optarg);
				break;
			case 'n':
				nb_hash_func=stoi(optarg);
				break;
			case 'e':
				bit_encoding=stoi(optarg);
				break;
			case 'h': // -h or --help
				PrintHelp();
				break;
			case '?': // Unrecognized option
				PrintHelp();
				break;
			default:
				PrintHelp();
				break;
		}
	}
}












int main(int argc, char **argv)
{
    omp_set_nested(1);
    // omp_set_num_threads(1);
	ProcessArgs(argc, argv);
	if (argc < 2)
	{
		PrintHelp();
	}

	cout << "************** SUPER TOOL GO ***************\n\n";

    // check bit encoding + build index
    if (!(bit_encoding == 8 or bit_encoding == 16 or bit_encoding == 32))
    {
        cerr << "[PARAMETER ERROR] Wrong bit encoding value\n";
        PrintHelp();
        return -1;
    }

    //WE BUILD THE INDEX
    if(existing_index != ""){
        cout << "LOADING index from file " << existing_index << endl;
        switch (bit_encoding)
        {
            case 8:
            {
                BestPart<uint8_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.load(existing_index);
                if(fof!=""){
                    ever.insert_file_of_file(fof);
                }
                 if(use_double_index){
                    ever.double_index();
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
            case 16:
            {
                BestPart<uint16_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.load(existing_index);
                if(fof!=""){
                    ever.insert_file_of_file(fof);
                }
                 if(use_double_index){
                    ever.double_index();
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
            case 32:
            {
                BestPart<uint32_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.load(existing_index);
                if(fof!=""){
                    ever.insert_file_of_file(fof);
                }
                 if(use_double_index){
                    ever.double_index();
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
        }
    }else if (fof!="")
	{
		

    

        cout << "Building index from file of file (" << fof << ")\n";
        switch (bit_encoding)
        {
            case 8:
            {
                BestPart<uint8_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.insert_file_of_file(fof);
                if(use_double_index){
                    ever.double_index();
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                if(dump_index!=""){
                    ever.serialize(dump_index);
                }
                break;
            }
            case 16:
            {
                BestPart<uint16_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.insert_file_of_file(fof);
                if(use_double_index){
                    ever.double_index();
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                if(dump_index!=""){
                    ever.serialize(dump_index);
                }
                break;
            }
            case 32:
            {
                BestPart<uint32_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.insert_file_of_file(fof);
                if(use_double_index){
                    ever.double_index();
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                if(dump_index!=""){
                    ever.serialize(dump_index);
                }
                break;
            }
        }
		return 0;
	}else 
	
	return 0;
}
