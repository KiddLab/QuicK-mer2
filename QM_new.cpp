#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>
#include <iterator>
#include <cmath>
#include <iomanip>
#include <inttypes.h>

#define Hash_size 0x2000000

struct QM_hash_struct
{
	uint64_t Kmer;
	uint64_t Next_loc_and_depth;
};

uint64_t Kmer_encode(char * kmer){
	uint64_t encoded = 0;
	//printf("%s\t",kmer);
	do {
		if (!*kmer) break;
		encoded <<= 2;
		/*
		switch (*kmer){
			case 'A': break;
			case 'T': encoded |= 1; break;
			case 'G': encoded |= 2; break;
			case 'C': encoded |= 3; break;
		}*/
		encoded |= (*kmer >> 1) & 3;
	}
	while (kmer++);
	//printf("%07X%08X\n",encoded >> 32, encoded);
	return encoded;
}

unsigned int DJBHash(char* str, unsigned int len)
{
	unsigned int hash = 5381;
	unsigned int i    = 0;
	for(i = 0; i < len; str++, i++)
	{   
		hash = ((hash << 5) + hash) + (*str);
	}
	return hash;
}

unsigned int DJBHash_encode(uint64_t kmer)
{
	unsigned int hash = 5381;
	unsigned int i    = 0;
	for(i = 0; i < 8; i++)
	{   
		hash = ((hash << 5) + hash) + ((kmer & 0xFF));
		kmer >>= 8;
	}
	return hash;
}

QM_hash_struct init_QM_hash;

int main(int argc, char** argv)
{
	std::ifstream control, kmer_list;
	//control.open(argv[2],std::ifstream::in);
	kmer_list.open(argv[1], std::ifstream::in);
	//output.open(argv[5],std::ofstream::out|std::ofstream::binary);
	unsigned int a,b,c,d, kmer_start, kmer_end;
	std::string chrom, randomstring, kmer;
	init_QM_hash.Kmer = 0;
	init_QM_hash.Next_loc_and_depth = 0;
	//std::vector<FAIDX> faidx;
	//std::map<std::string, unsigned int> Chrom_map;
	std::vector<QM_hash_struct> hash_array (Hash_size + (Hash_size >> 4), init_QM_hash);
	uint64_t last_index = 0;
	uint32_t worst_case = 0;
	puts("Memory allocated");
	uint32_t hist[8192] = {0};
	while (kmer_list >> chrom >> kmer_start >> kmer_end >> randomstring >> kmer)
	{
		
		char * kmer_cstr = (char * ) kmer.c_str();
		uint64_t kmer_value = Kmer_encode(kmer_cstr);
		uint32_t hash_index = DJBHash_encode(kmer_value); 
		//uint32_t hash_index = DJBHash(kmer_cstr, 30);
		hash_index &= (Hash_size-1);
		uint32_t worst = 0;
		while (hash_array[hash_index].Kmer != 0){
			hash_index++;
			worst++;
		}
		hash_array[hash_index].Kmer = kmer_value;
		hash_array[last_index].Next_loc_and_depth = hash_index << 16;
		last_index = hash_index;
		if (worst > worst_case) {
			worst_case = worst;
			printf("Worst %i\n", worst);
		}
		hist[worst]++;
	}
	uint32_t count = 0; 
	float average = 0;
	for (int k = 0; k <= worst_case; k++){
		count += hist[k];
		average += k * hist[k];
	}
	printf("Average %f, fill %f\%\n", average/count, ((float) count * 100)/ Hash_size);
}
