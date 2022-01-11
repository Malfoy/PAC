#~ ******************* ABC **************************************

 #~ INDEX CONSTRUCTION
#~ --fof (-f)               :     Build index from file of files (FASTA/Q/GZ allowed in file of files)

#~ --load (-l)              :     Load index from folder 

#~ --dir (-d)               :     Write index in folder
#~ --hot (-h)               :     Keep all Bloom in RAM, (fast but expensive mode)

 #~ INDEX QUERY
#~ --query (-q)             :     Query sequence file (FASTA/Q/GZ)
#~ --out (-o)               :     Write query output in file (default: output.gz)

 #~ TWEAK PARAMETERS
#~ -k                       :     k-mer size (default: 31)
#~ -b                       :     Bloom filter size (default 134,217,728)
#~ -n                       :     Bloom filter's number of hash functions (default: 1)
#~ -e                       :     Bit encoding, possible values are 8 (max 256 files), 16 (max 65,636 files), 32 (Max 4,294,967,296 files) (default: 16)
#~ -u                       :     Filter unique Kmers 
#~ -p                       :     Partitionning filters in 4^p files (default: 3 for 64 file per filter)
#~ -c                       :     Number of core to use (default: all)
#~ -i                       :     Build double index for faster queries 



#~ Construct an index in a folder with default parameter
./best -f file_of_files.txt -d Index1;

#~ Load an index from a folder an perform a query
./best -l Index1 -q query_file.fa -o output_query.gz;

#~ Construct an index  with custom parameters (using 21mers, a larger bloom filter, four cores and filtering unique kmers)
./best -f file_of_files.txt -d Index2 -k 21 -b 2000000000 -c 4 -u;

#~ Construct an index  with low amount of partition and query it
./best -f file_of_files.txt -d Index3 -p 1 -q query_file.fa -o output_query_slow.gz;

#~ Construct an index  with high amount of partition and query it
./best -f file_of_files.txt -d Index4 -p 4 -q query_file.fa -o output_query_fast.gz;

#~ Construct an index  and build two ABF for faster queries
./best -f file_of_files.txt -d Index5 -p 4 -i -q query_file.fa -o output_query_faster.gz;

#~ Updating an index  by loading it and adding files to it 
./best -l Index1  -f file_of_files2.txt -d Index6  -q query_file.fa -o output_query_update.gz;

#~ To clean the mess
#~ rm -rf  Index1 Index2 Index3 Index4 Index5 Index6 output_query.gz output_query_slow.gz output_query_fast.gz output_query_faster.gz output_query_update.gz
