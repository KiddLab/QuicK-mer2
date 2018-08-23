#include "stdint.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "siphash24.h"

//#define Hash_size 0x100000000
#define Hash_size 0x2000000
#define buffer_size 1024*1024

const uint64_t Kmer_size = 30;
const uint8_t key = 42;

uint64_t Kmer_encode(char * kmer){
	uint64_t encoded = 0;
	uint64_t encoded_r = 0;
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
		uint8_t Letter = (*kmer >> 1) & 3;
		encoded |= Letter;
		Letter = (Letter - 2) & 3; //Very special conversion between A-T and G-C
		encoded_r |= (uint64_t)Letter << 60;
		encoded_r >>= 2;
	}
	while (kmer++);
	//printf("%07X%08X\n",encoded >> 32, encoded);
	if (encoded > encoded_r) return encoded_r;
	return encoded;
}

unsigned int DJBHash(char* str, unsigned int len)
{
	unsigned int hash = 5381;
	unsigned int i    = 0;
	for(i = 0; i < len; str++, i++)
	{
		hash = ((hash << 5) + hash) ^ (*str);
	}
	return hash;
}

uint64_t DJBHash_encode(uint64_t kmer)
{
	
	uint64_t hash = 5381;
	uint8_t i = 0;
	for(i = 0; i < 8; i++)
	{
		hash = ((hash << 5) + hash) + ((kmer & 0xFF));
		kmer >>= 8;
	}
	
	//uint64_t hash;
	//siphash(&hash,(uint8_t *) &kmer, 8, &key);
	
	return hash;
}

int main_hash(int argc, char ** argv)
{
	FILE *kmer_list = fopen(argv[0], "r");
	FILE * Hash_file = fopen(argv[1], "w");
	if (!Hash_file) {
		puts("File creation failed");
		return 1;
	}
	fseek(kmer_list, 0, SEEK_SET);
	if (!kmer_list) return 1;
	uint32_t kmer_start, kmer_end;
	char chrom[30];
	char kmer[60];
	char rs[60];
	//Malloc
	uint64_t * Kmer_hash = (uint64_t *) malloc((Hash_size + 16384) * sizeof(uint64_t));
	uint32_t * Kmer_next_index = (uint32_t *) malloc(sizeof(uint32_t) * (Hash_size + 16384));
	if (!Kmer_hash || !Kmer_next_index)
	{
		puts("Memory allocation failed");
		return 1;
	}
	uint64_t last_index = 0;
	uint32_t worst_case = 0;
	uint32_t first_index = 0;
	uint32_t hist[8192] = {0};
	while (fscanf(kmer_list, "%s\t%u\t%u\t%s\t%s", chrom, &kmer_start, &kmer_end, rs, kmer) == 5)
	{
		uint64_t kmer_value = Kmer_encode(kmer);
		uint32_t hash_index = DJBHash_encode(kmer_value);
		//uint32_t hash_index = DJBHash(kmer_cstr, 30);
		hash_index &= (Hash_size-1);
		uint32_t worst = 0;
		while (Kmer_hash[hash_index] != 0){
			hash_index++;
			worst++;
		}
		Kmer_hash[hash_index] = kmer_value;
		Kmer_next_index[last_index] = hash_index;
		if (!first_index) first_index = hash_index;
		last_index = hash_index;
		/*if (worst > worst_case) {
			worst_case = worst;
			printf("Worst %i\n", worst);
		}*/
		hist[worst]++;
	}
	Kmer_next_index[last_index] = 0xFFFFFFFF;
	uint32_t count = 0;
	float average = 0;
	for (uint32_t k = 0; k < 8192; k++){
		count += hist[k];
		average += k * hist[k];
	}
	printf("Average %f, fill %f\% \n", average/count, ((float) count * 100)/ Hash_size);
	rewind(Hash_file);
	const char Version[5] = "QM11";
	fwrite(Version, 1, 4, Hash_file);
	uint64_t Hashsize = Hash_size + 16384;
	fwrite(&Hashsize, 4, 1, Hash_file);
	fwrite(&first_index, 4, 1, Hash_file);
	fwrite((void *) Kmer_hash, sizeof(uint64_t), Hash_size + 16384, Hash_file);
	int count_kmer = 0;
	for (uint32_t k = 0; k < Hash_size +16384; k++){
		if (Kmer_hash[k]) count_kmer++;
	}
	printf("Total %i kmers\n", count_kmer);
	free(Kmer_hash);
	fwrite((void *) Kmer_next_index, sizeof(uint32_t), Hash_size + 16384, Hash_file);
	free(Kmer_next_index);
	fclose(Hash_file);
	return 0;
}

int main_count(int argc, char ** argv)
{
	FILE * Hash_file = fopen(argv[0], "r");
	uint32_t Hashsize;
	fseek(Hash_file, 4, 0);
	fread(&Hashsize, 1, 4, Hash_file);
	uint32_t first_idx;
	fread(&first_idx, 1, 4, Hash_file);
	printf("Hash Size: %i\nFirst location: %i\n", Hashsize, first_idx);
	uint64_t * Kmer_hash = (uint64_t *) malloc(Hashsize * sizeof(uint64_t));
	if (!Kmer_hash) {
		puts("Memory allocation failed");
		fclose(Hash_file);
		return 1;
	}
	printf("Read %i hash\n",fread(Kmer_hash, sizeof(uint64_t), Hashsize, Hash_file));
	uint16_t * Kmer_depth = (uint16_t *) calloc(Hashsize, sizeof(uint16_t));
	if (!Kmer_depth) {
		puts("Memory allocation failed");
		free(Kmer_hash);
		fclose(Hash_file);
		return 1;
	}
	FILE * fasta_input = fopen(argv[1], "r");
	if (!fasta_input) {
		puts("Input open fail");
	}
	char line[350];
	while (fgets(line, 350, fasta_input)){
		if (line[0] == '>') continue;
		char * seq_char_index = line;
		uint64_t encoded = 0;
		uint64_t encoded_r = 0;
		uint16_t cur_chars = 0;
		while (*seq_char_index != '\n') {
			if (*seq_char_index == 'N') {
				encoded = 0;
				encoded_r = 0;
				cur_chars = 0;
			}
			else {
				cur_chars++;
				uint8_t Letter = (*seq_char_index >> 1) & 3;
				encoded <<= 2;
				encoded |= Letter;
				Letter = (Letter - 2) & 3; //Very special conversion between A-T and G-C
				encoded_r |= (uint64_t)Letter << 60;
				encoded_r >>= 2;
				//Hash
				if (cur_chars >= Kmer_size) {
					uint64_t kmer = encoded & (((uint64_t)1 << (Kmer_size << 1)) - 1);
					uint32_t hash_index = DJBHash_encode(kmer) & (Hash_size - 1);
					if (kmer > encoded_r) kmer = encoded_r;
					while (Kmer_hash[hash_index] && Kmer_hash[hash_index] != kmer) hash_index++;
					if (Kmer_hash[hash_index]) {
						if (Kmer_depth[hash_index] != 65535)
							__sync_fetch_and_add(&Kmer_depth[hash_index], 1);
					}
				}
			}
			seq_char_index++;
		}
	}
	//Dump count file
	printf("Pileup finish\nRead chain file %i entries\n",fread(Kmer_hash, sizeof(uint32_t), Hashsize, Hash_file));
	uint32_t * Kmer_next_loc = (uint32_t *) Kmer_hash;
	uint32_t buf_count = 0;
	uint16_t buffer[buffer_size];
	FILE * OutFile = fopen(argv[2], "w");
	while(first_idx != 0xFFFFFFFF)
	{
		buffer[buf_count] = Kmer_depth[first_idx];
		first_idx = Kmer_next_loc[first_idx];
		buf_count++;
		if (buf_count == buffer_size) {
			fwrite(buffer, sizeof(uint16_t), buffer_size, OutFile);
			buf_count = 0;
		}
	}
	fwrite(buffer, sizeof(uint16_t), buf_count, OutFile);
	free(Kmer_hash);
	free(Kmer_depth);
	fclose(OutFile);
	return 0;
}

int main_search(int argc, char ** argv)
{
	//Malloc
	uint64_t * Kmer_hash = (uint64_t *) malloc((Hash_size + 16384) * sizeof(uint64_t));
	uint32_t * Kmer_next_index = (uint32_t *) malloc(sizeof(uint32_t) * (Hash_size + 16384));
	uint8_t * Kmer_occr = (uint8_t *) malloc(sizeof(uint8_t) * (Hash_size + 16384));
	uint8_t * Kmer_edit_depth = (uint8_t *) malloc(sizeof(uint8_t) * (Hash_size + 16384));
	if (!Kmer_hash || !Kmer_next_index)
	{
		puts("Memory allocation failed");
		return 1;
	}
	FILE * fasta = fopen(argv[0], "r");
	char fasta_buffer[200];
	uint8_t charge_size = 0;
	uint64_t encoded = 0;
	uint64_t encoded_r = 0;
	uint64_t processed = 0;
	uint32_t worst = 0;
	uint32_t hist[131071] = {0};
	//Loop through fasta lines
	while (fgets (fasta_buffer, 200, fasta) && fasta_buffer[0])
	{
		char * char_idx = fasta_buffer;
		if (*char_idx == '>')
		{
			charge_size = 0;
			encoded = 0;
			encoded_r = 0;
			printf("%s", fasta_buffer);
			continue;
		}
		while (*char_idx && *char_idx != '\n')
		{
			if (*char_idx == 'N'){
				charge_size = 0;
				encoded = 0;
				encoded_r = 0;
				char_idx++;
				continue;
			}
			uint8_t Letter = (*char_idx >> 1) & 3;
			char_idx++;
			encoded <<= 2;
			encoded |= Letter;
			Letter = (Letter - 2) & 3; //Very special conversion between A-T and G-C
			encoded_r |= (uint64_t)Letter << 60;
			encoded_r >>= 2;
			uint64_t kmer = encoded & (((uint64_t)1 << (Kmer_size << 1)) - 1);
			if (kmer > encoded_r) kmer = encoded_r;
			uint64_t hash_index = DJBHash_encode(kmer) & (Hash_size - 1);
			if (charge_size < Kmer_size) charge_size++;
			if (kmer && charge_size == Kmer_size)
			{
				//Add to hash memory
				uint32_t collision = 0;
				while (Kmer_hash[hash_index] && Kmer_hash[hash_index] != kmer)
				{
					hash_index++;
					collision++;
				}
				if (!Kmer_hash[hash_index])
				{
					if (collision > worst){
						worst = collision;
						printf("Worst %u\n", worst);
					}
					if (collision < 131071) hist[collision]++;
					else hist[131071]++;
					Kmer_hash[hash_index] = kmer;
				}
				Kmer_occr[hash_index]++;
			}
		}
		processed++;
		if (processed % 1666667 == 0){
			float average = 0;
			uint64_t count = 0;
			for (uint32_t k = 0; k < 131072; k++){
				count += hist[k];
				average += k * hist[k];
			}
			average /= count;
			printf("Processed %ubp, total %u Kmers, average collision %f\n", processed*60, count, average);
		}
	}
	float average = 0;
	uint64_t count = 0;
	for (uint32_t k = 0; k < 131072; k++){
		count += hist[k];
		average += k * hist[k];
	}
	printf("Average %f, fill %f\% \n", average/count, ((float) count * 100)/ Hash_size);
	uint64_t occr_idx = 0;
	uint64_t unique_count = 0;
	while (occr_idx < Hash_size+16384)
	{
		if (Kmer_occr[occr_idx] == 1) unique_count++;
		occr_idx++;
	}
	printf("Uniq count %u\n", unique_count);
	//Filter 
	
}

void printversion() {
	puts("QuicK-mer 2.0");
	puts("Operation modes: \n\tindex\tIndex a kmer list\n\tcount\tCNV estimate from library\n\tsearch\tSearch K-kmer in genome\n");
}

int main(int argc, char** argv)
{
	if (argc > 1) {
		if (strcmp(argv[1], "index") == 0)
			return main_hash(argc-2, argv+2);
		else if (strcmp(argv[1], "count") == 0)
			return main_count(argc-2, argv+2);
		else if (strcmp(argv[1], "search") == 0)
			return main_search(argc-2, argv+2);
		else {
			printversion();
			return 1;
		}
	}
	else {
		printversion();
		return 1;
	}
	return 0;
}
