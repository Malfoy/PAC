#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <getopt.h>
#include "best.h"
#include "bestpart.h"
#include <filesystem>



using namespace std;
using namespace chrono;
using namespace filesystem;



string fof(""), w_dir("My_index"), query_file(""),existing_index(""),query_output("output.gz");
uint16_t nb_hash_func(1), bit_encoding(16);
uint32_t kmer_size(31);
uint32_t nb_partition(4);
uint64_t bf_size(134217728); 
uint64_t core_number(0);
bool use_double_index(false),filter_unique(false);



void PrintHelp()
{
	cout <<
			"\n******************* ABC **************************************\n"


			"\n INDEX CONSTRUCTION\n"
            "--fof (-f)               :     Build index from file of files (FASTA/Q/GZ allowed in file of files)\n\n"
            "--load (-l)              :     Load index from folder \n\n"
			"--dir (-d)               :     Write index in folder (default: My_index)\n"

			"\n INDEX QUERY\n"
			"--query (-q)             :     Query sequence file (FASTA/Q/GZ)\n"
            "--out (-o)               :     Write query output in file (default: output.gz)\n"
			
			"\n TWEAK PARAMETERS\n"
			"-k                       :     k-mer size (default: " << kmer_size << ")\n"
			"-b                       :     Bloom filter size (default " << intToString(bf_size) << ")\n"
			"-n                       :     Bloom filter's number of hash functions (default: " << intToString( nb_hash_func) << ")\n"
			"-e                       :     Bit encoding, possible values are 8 (max 256 files), 16 (max 65,636 files), 32 (Max 4,294,967,296 files) (default: " << bit_encoding << ")\n"
            "-u                       :     Filter unique Kmers \n"
            "-p                       :     Partitionning filters in 4^p files (default: "<<nb_partition<<" for "<<intToString(1<<(2*nb_partition))<<" file per filter)\n"
            "-c                       :     Number of core to use (default: all)\n"
            "-i                       :     Build double index for faster queries \n";

	exit(1);
}



void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "k:d:q:b:f:e:l:n:hiuo:p:c:";
	const option long_opts[] = 
	{
		{"index", no_argument, nullptr, 'i'},
		{"fof", required_argument, nullptr, 'f'},
		{"dir", required_argument, nullptr, 'd'},
        {"out", required_argument, nullptr, 'o'},
		{"query", required_argument, nullptr, 'q'},
		{"load", required_argument, nullptr, 'l'},
        {"uniqu", required_argument, nullptr, 'u'},
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
            case 'u':
				filter_unique=true;
				break;
			case 'k':
				kmer_size=stoi(optarg);
				break;
            case 'p':
				nb_partition=stoi(optarg);
				break;
            case 'c':
				core_number=stoi(optarg);
				break;
			case 'd':
				w_dir=(optarg);
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
				bf_size=(stoull(optarg));
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
    if(core_number!=0){
        omp_set_num_threads(1);
    }

    
	ProcessArgs(argc, argv);
	if (argc < 2)
	{
		PrintHelp();
	}

	cout << "************** PAC  ***************\n\n";

    // check bit encoding + build index
    if (!(bit_encoding == 8 or bit_encoding == 16 or bit_encoding == 32))
    {
        cerr << "[PARAMETER ERROR] Wrong bit encoding value\n";
        PrintHelp();
        return -1;
    }
    bf_size=approx_power2(bf_size);
    w_dir=filesystem::absolute(w_dir);
    cout<<w_dir<<endl;
    //WE BUILD THE INDEX
    if(existing_index != ""){
        existing_index=absolute(existing_index);
        switch (bit_encoding)
        {
            case 8:
            {
                BestPart<uint8_t> ever(existing_index);
                if(fof!=""){
                    ever.insert_file_of_file(fof);
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
            case 16:
            {
                BestPart<uint16_t> ever(existing_index);
                if(fof!=""){
                    ever.insert_file_of_file(fof);
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
            case 32:
            {
                BestPart<uint32_t> ever(existing_index);
                if(fof!=""){
                    ever.insert_file_of_file(fof);
                }
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
        }
    }else if (fof!="")
	{

        switch (bit_encoding)
        {
            case 8:
            {
                BestPart<uint8_t> ever(bf_size, nb_hash_func, kmer_size,filter_unique,w_dir,use_double_index,nb_partition);
                ever.insert_file_of_file(fof);
                ever.serialize();
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
            case 16:
            {
                BestPart<uint16_t> ever(bf_size, nb_hash_func, kmer_size,filter_unique,w_dir,use_double_index,nb_partition);
                ever.insert_file_of_file(fof);
                ever.serialize();
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
            case 32:
            {
                BestPart<uint32_t> ever(bf_size, nb_hash_func, kmer_size,filter_unique,w_dir,use_double_index,nb_partition);
                ever.insert_file_of_file(fof);
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                ever.serialize();
                if(query_file!=""){
                    ever.query_file(query_file,query_output);
                }
                break;
            }
        }
		return 0;
	}else 
	
	return 0;
}
