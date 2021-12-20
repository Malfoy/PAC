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



using namespace std;
using namespace chrono;

bool do_index(0), do_query(0);
string fof(""), index_files("./index_files"), query_file("");
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
			"--index_dir (-d)         :     provide location to write index files (default: " << index_files << ")\n\n"
			"--fof (-f)               :     provide file of files for index construction (FASTA/Q/GZ allowed in file of files)\n\n"
			"\n INDEX QUERY\n"
			"--query (-q)             :     provide sequence file (FASTA/Q/GZ) for the query\n"
			"--load (-l)              :     index files location (default:.)\n\n"
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
				index_files=optarg;
				break;
			case 'i':
				do_index=true;
				break;
			case 'f':
				fof=optarg;
				break;
			case 'q':
				query_file=optarg;
				break;
			case 'l':
				index_files=optarg;
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
	// BUILD INDEX
	cout << "************** SUPER TOOL GO ***************\n\n";
	if (do_index)
	{
		if (fof == "")
		{
			cerr << "[PARAMETER ERROR] Please provide a valid file of files\n";
			PrintHelp();
			return -1;
		}
		int systRet(0);
		//directory for index serialization
		if (not directory_exists(index_files)) 
		{
			cout << "Creating directory for index serialization in " << index_files << endl;
			systRet=system(("mkdir " + index_files).c_str());
		}
		else
		{
			cout << "Index file already exists in " << index_files << endl;
		}
		// check bit encoding + build index
		if (!(bit_encoding == 8 or bit_encoding == 16 or bit_encoding == 32))
		{
			cerr << "[PARAMETER ERROR] Wrong bit encoding value\n";
			PrintHelp();
			return -1;
		}
		else
		{
			cout << "Building index from file of file (" << fof << ") in " << index_files << "...\n";
			switch (bit_encoding)
			{
				case 8:
				{
					Best<uint8_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
					ever.Best<uint8_t>::insert_file_of_file(fof);
					break;
				}
				case 16:
				{
					Best<uint16_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
					ever.Best<uint16_t>::insert_file_of_file(fof);
					break;
				}
				case 32:
				{
					Best<uint32_t> ever(bf_size, bf_size, nb_hash_func, kmer_size);
					ever.Best<uint32_t>::insert_file_of_file(fof);
					break;
				}
			}
		}
		return 0;
	}
	// QUERY
	if (query_file != "")
	{
		// first load index //
		
		// then query //
		string sequence("");
		char data_type (get_data_type(query_file));
		// read query file
		zstr::ifstream * query_sequences = new zstr::ifstream(query_file); 
		while(not query_sequences->eof()) 
		{
			Biogetline(query_sequences, sequence, data_type, kmer_size);
			//~ ever.Best<uint8_t>::query_sequence(sequence);
		}
		delete query_sequences;
		return 0;
	}
	return 0;
}



//~ int main(int argc, char ** argv){
   
    //~ srand(time(NULL));
    //~ string  fof((argv[1]));
    //~ //TODO BENCH CONSTRUCTION AND QUERY TIME
    //~ Best<uint8_t> ever(1000*1000*1000, 1000*1000*1000*1,1,21);
    //~ ever.Best<uint8_t>::insert_file_of_file(fof);
    //~ ever.Best<uint8_t>::query_sequence("GACTGGTATGCGTTTTCACTGCCGAAAATGAAAGCCAGTAAAGAAGTTACAACGGACGATGAGTTACGTATCTGGAAATAAGGTTGAAAAATAAAAACGGCGCTAAAAAGCGCCGTTTTTTTTGACGGTGGTAAAGCCGATTAATCTTCCAGATCACCGCAGAAGCGATAACCTTCACCGTGAATGGTGGCGATGATTTCCGGCGTATCCGGCGTAGATTCGAAATGTTTACGAATACGGCGGATCGTCACGTCTACAGTACGGTCGTGCGGTTTCAGCTCACGGCCGGTCATTTTCTTCAGCAGTTCAGCACGGGACTGAATTTTGCCTGGGTTTTCACAGAAGTGAAGCATGGCGCGGAACTCGCTGCGCGGCAGCTTGTACTGCTCGCCATCAGGGCCGATCAACGAACGGCTGTTGATGTCCAGTTCCCAACCATTGAACTTGTAGCTTTCAACGCTACGACGTTCTTCGCTGACAGTACCCAGATTCATGGTACGGGACAGTAGGTTGCGTGCACGAATCGTCAGTTCACGCGGGTTGAACGGTTTGGTGATGTAGTCATCTGCACCGATATCGAGGCCGAGAATTTTATCGACTTCGTTGTCACGGCCAGTCAGGAACATCAACGCAACATTCGCCTGCTCGCGCAGTTCACGCGCTAACAGAAGACCGTTCTTACCCGGCAGATTGATATCCATGATCACCAGGTTGATGTCATATTCAGAGAGGATCTGATGCATTTCCGCGCCATCTGTCGCTTCGAAAACATCATAGCCTTCCGCTTCGAAAATACTTTTCAACGTGTTGCGTGTTACCAACTCGTCTTCAACGATAAGAATGTGCGGGGTCTGCATGTTTGCTACCTAAATTGCCAACTAAATCGAAACAGGAAGTACAAAAGTCCCTGACCTGCCTGATGCATGCTGCAAATTAACATGATCGGCGTAACATGACTAAAGTACGTAATTGCGTTCTTGATGCACTTTCCATCAACGTCAACAACATCATTAGCTTGGTCGTGGGTACTTTCCCTCAGGACCCGACAGTGTCAAAAACGGCTGTCATCCTAACCATTTTAACAGCAACATAACAGGCTAAGAGGGGCCGGACACCCAATAAAACTACGCTTCGTTGACATATATCAAGTTCAATTGTAGCACGTTAACAGTTTGATGAAATCATCGTATCTAAATGCTAGCTTTCGTCACATTATTTTAATAATCCAACTAGTTGCATCATACAACTAATAAACGTGGTGAATCCAATTGTCGAGATTTATTTTTTATAAAGTTATCCTAAGTAAACAGAAGGATATGTAGCATTTTTTAACAACTCAACCGTTAGTACAGTCAGGAAATAGTTTAGCCTTTTTTAAGCTAAGTAAAGGGCTTTTTCTGCGACTTACGTTAAGAATTTGTAAATTCGCACCGCGTAATAAGTTGACAGTGATCACCCGGTTCGCGGTTATTTGATCAAGAAGAGTGGCAATATGCGTATAACGATTATTCTGGTCGCACCCGCCAGAGCAGAAAATATTGGGGCAGCGGCGCGGGCAATGAAAACGATGGGGTTTAGCGATCTGCGGATTGTCGATAGTCAGGCACACCTGGAGCCAGCCACCCGCTGGGTCGCACATGGATCTGGTGATATTATTGATAATATTAAAGTTTTCCCGACATTGGCTGAATCGTTACACGATGTCGATTTCACTGTCGCCACCACTGCGCGCAGTCGGGCGAAATATCATTACTACGCCACGCCAGTTGCACTGGTGCCGCTGTTAGAGGAAAAATCTTCATGGATGAGCCATGCCGCGCTGGTGTTTGGTCGCGAAGATTCCGGGTTGACTAACGAAGAGTTAGCGTTGGCTGACGTTCTTACTGGTGTGCCGATGGTGGCGGATTATCCTTCGCTCAATCTGGGGCAGGCGGTGATGGTCTATTGCTATCAATTAGCAACATTAATACAACAACCGGCGAAAAGTGATGCAACGGCAGACCAACATCAACTGCAAGCTTTACGCGAACGAGCCATGACATTGCTGACGACTCTGGCAGTGGCAGATGACATAAAACTGGTCGACTGGTTACAACAACGCCTGGGGCTTTTAGAGCAACGAGACACGGCAATGTTGCACCGTTTGCTGCATGATATTGAAAAAAATATCACCAAATAAAAAACGCCTTAGTAAGTATTTTTC");
    //~ ever.Best<uint8_t>::get_stats();
    //~ // Bcardi Richelieu(1,(uint64_t)10*1000*1000*1000,1000*1000*1000,1,31);
    //~ // Richelieu.insert_file_of_file(fof);
    //~ // Richelieu.get_cardinalities();
//~ }
