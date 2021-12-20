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



string fof(""), dump_index(""), query_file(""),existing_index("");
uint16_t nb_hash_func(1), bit_encoding(16);
uint32_t kmer_size(31);
uint64_t bf_size(100000000); //todo



void PrintHelp()
{
	cout <<
			"\n******************* Amazing tool **************************************\n"
			"******************* Does something ******************\n"

			"--help (-h)             :     Show help\n\n"
			"\n INDEX CONSTRUCION\n"
			"--index (-i)             :     build index\n"
            "--load (-l)              :     load index from disk \n\n"
			"--dump index (-d)         :     provide location to write index file\n"
			"--fof (-f)               :     provide file of files for index construction (FASTA/Q/GZ allowed in file of files)\n\n"
			"\n INDEX QUERY\n"
			"--query (-q)             :     provide sequence file (FASTA/Q/GZ) for the query\n"
			
			"\n OTHER OPTIONS\n"
			"-k                       :     k-mer size (default: " << kmer_size << ")\n"
			"-b                       :     Bloom filter size (default " << bf_size << ")\n"
			"-n                       :     Bloom filter's number of hash functions (default: " << nb_hash_func << ")\n"
			"-e                       :     bit encoding (possible values 8,16,32, default: " << bit_encoding << ")\n";

	exit(1);
}




void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "k:d:q:b:f:e:l:n:hi";
	const option long_opts[] = 
	{
		{"index", no_argument, nullptr, 'i'},
		{"fof", required_argument, nullptr, 'f'},
		{"index_dir", required_argument, nullptr, 'd'},
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
	if (fof!="")
	{
		

    

        cout << "Building index from file of file (" << fof << ")\n";
        switch (bit_encoding)
        {
            case 8:
            {
                BestPart<uint8_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.insert_file_of_file(fof);
                if(query_file!=""){
                    ever.query_file(query_file);
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
                if(query_file!=""){
                    ever.query_file(query_file);
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
                if(query_file!=""){
                    ever.query_file(query_file);
                }
                if(dump_index!=""){
                    ever.serialize(dump_index);
                }
                break;
            }
        }
		return 0;
	}else if(existing_index != ""){
        cout << "LOADING index from file " << existing_index << endl;
        switch (bit_encoding)
        {
            case 8:
            {
                BestPart<uint8_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.load(existing_index);
                cout<<"Loaded8 "<<endl;
                if(query_file!=""){
                    ever.query_file(query_file);
                }
                break;
            }
            case 16:
            {
                BestPart<uint16_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.load(existing_index);
                cout<<"Loaded16 "<<endl;
                if(query_file!=""){
                    ever.query_file(query_file);
                }
                break;
            }
            case 32:
            {
                BestPart<uint32_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
                ever.load(existing_index);
                if(query_file!=""){
                    ever.query_file(query_file);
                }
                break;
            }
        }
    }
	
	return 0;
}
