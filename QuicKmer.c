#include "stdint.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pthread.h"
#include "unistd.h"
#include "time.h"
#include "semaphore.h"

#define buffer_size 1024*1024
#define FIFO_size 4096

uint64_t Hash_size = 0x2000000; //0x100000000
uint64_t Kmer_size = 30;
uint64_t Kmer_mask = ((uint64_t)1 << 60) - 1;

//Global Variables
uint64_t * Kmer_hash;
uint32_t * Kmer_next_index;
uint8_t * Kmer_occr;
uint8_t * Kmer_edit_depth;
uint16_t * Kmer_depth;
uint8_t edit_distance = 2;
uint8_t Edit_depth_thres = 100;
uint8_t thread_no_more_data = 0; //Set this flag when input stream finish

struct edit_dis_arg_struc {
	uint64_t start_idx;
	uint64_t end_idx;
	uint8_t thread_id;
};

struct FIFO_arg_struc {
	volatile uint64_t FIFO0[FIFO_size]; //FIFO
	volatile uint64_t FIFO1[FIFO_size]; //FIFO
	sem_t data_feed_sem;
	volatile uint8_t Write_count;
	volatile uint8_t Read_count;
	uint8_t thread_id;
};

void Set_Kmer_Size(uint8_t Size){
	Kmer_size = Size;
	Kmer_mask = ((uint64_t)1 << (Kmer_size * 2)) - 1;
}

uint64_t Kmer_encode(char * kmer){
	uint64_t encoded = 0;
	uint64_t encoded_r = 0;
	do {
		if (!*kmer) break;
		encoded <<= 2;
		uint8_t Letter = (*kmer >> 1) & 3;
		encoded |= Letter;
		Letter = (Letter - 2) & 3; //Very special conversion between A-T and G-C
		encoded_r |= (uint64_t)Letter << 60;
		encoded_r >>= 2;
	}
	while (kmer++);
	//printf("%llX\n",encoded);
	if (encoded > encoded_r) return encoded_r;
	return encoded;
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
	return hash;
}

void Permute_kmer(uint64_t *encoded_kmer, uint64_t *encoded_reverse, char position, char edit)
{
	//return the permuted_kmer
	uint64_t base = ((*encoded_kmer >> (position << 1)) & 3) + edit;
	base &= 3;
	(*encoded_kmer) &= Kmer_mask - (3 << (position << 1));
	(*encoded_kmer) |= base << (position << 1);
	base = (base - 2) & 3;
	(*encoded_reverse) &= Kmer_mask - (3 << (Kmer_size - 1 - position)*2);
	(*encoded_reverse) |= base << (Kmer_size - 1 - position)*2;
}

char Find_hash(uint64_t encoded_kmer, uint64_t * hash_index, uint64_t * Hash_dict)
{
	*hash_index = DJBHash_encode(encoded_kmer) & (Hash_size-1);
	while (Hash_dict[*hash_index] && Hash_dict[*hash_index] != encoded_kmer)
		(*hash_index)++;
	return Hash_dict[*hash_index] == encoded_kmer;
}

uint64_t Reverse_strand_encoded(uint64_t encoded_kmer)
{
	uint64_t encoded_reverse = 0;
	uint8_t index;
	for (index = 0; index < Kmer_size; index++){
		encoded_reverse <<= 2;
		encoded_reverse |= ((encoded_kmer & 3)-2) & 3;
		encoded_kmer >>= 2;
	}
	return encoded_reverse;
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
	char chrom[64];
	char kmer[33];
	//Malloc
	Kmer_hash = (uint64_t *) malloc((Hash_size + 16384) * sizeof(uint64_t));
	Kmer_next_index = (uint32_t *) malloc(sizeof(uint32_t) * (Hash_size + 16384));
	if (!Kmer_hash || !Kmer_next_index)
	{
		puts("Memory allocation failed");
		return 1;
	}
	uint64_t last_index = 0;
	uint32_t worst_case = 0;
	uint64_t first_index = 0;
	uint64_t Total_count = 0;
	uint32_t hist[8192] = {0};
	while (fscanf(kmer_list, "%s\t%u\t%u\t%s", chrom, &kmer_start, &kmer_end, kmer) == 4)
	{
		if (!Total_count) Set_Kmer_Size(strlen(kmer));
		uint64_t kmer_encoded = Kmer_encode(kmer);
		uint64_t hash_index = DJBHash_encode(kmer_encoded) & (Hash_size-1);
		uint32_t worst = 0;
		while (Kmer_hash[hash_index] != 0){
			hash_index++;
			worst++;
		}
		Kmer_hash[hash_index] = kmer_encoded;
		if (!Total_count) first_index = hash_index;
		else Kmer_next_index[last_index] = hash_index;
		last_index = hash_index;
		/*if (worst > worst_case) {
			worst_case = worst;
			printf("Worst %i\n", worst);
		}*/
		if (worst > 8191) worst = 8191;
		hist[worst]++;
		Total_count++;
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
	fwrite(&Kmer_size, 1, 1, Hash_file);
	fwrite(&edit_distance, 1, 1, Hash_file);
	fwrite(&Edit_depth_thres, 1, 1, Hash_file);
	fwrite(&Edit_depth_thres, 1, 1, Hash_file);
	uint64_t Hashsize = Hash_size + 16384;
	fwrite(&Hashsize, 8, 1, Hash_file);
	fwrite(&first_index, 8, 1, Hash_file);
	fwrite((void *) Kmer_hash, sizeof(uint64_t), Hash_size + 16384, Hash_file);
	int count_kmer = 0;
	for (uint64_t k = 0; k < Hash_size +16384; k++){
		if (Kmer_hash[k]) count_kmer++;
	}
	printf("Total %i kmers\n", count_kmer);
	free(Kmer_hash);
	fwrite((void *) Kmer_next_index, sizeof(uint32_t), Hash_size + 16384, Hash_file);
	free(Kmer_next_index);
	fclose(Hash_file);
	return 0;
}

void * Kmer_count_TSK(void *argvs)
{
	struct FIFO_arg_struc *argv = (struct FIFO_arg_struc *) argvs;
	//clock_t wait_before, thd_begin_time, total_time;
	//thd_begin_time = clock();
	//total_time = 0;
	//uint64_t total_processed = 0;
	//uint8_t last_wait = 0;
	while(1) {
		//wait_before = clock();
		sem_wait(&(argv -> data_feed_sem));
		if (thread_no_more_data) break;
		/*
		if ((argv -> Read_count) == (argv -> Write_count)) {
			if (last_wait == 0){
				last_wait = 1;
				wait_before = clock();
			}
			if (thread_no_more_data) break;
			else continue;
		}
		if (last_wait) {
			total_time += clock()-wait_before;
			last_wait = 0;
		}*/
		//total_time += clock()-wait_before;
		//Work on current batch
		volatile uint64_t * Kmers;
		if ((argv -> Read_count) & 1) Kmers = argv -> FIFO1;
		else Kmers = argv -> FIFO0;
		uint32_t i;
		for (i = 0; i < FIFO_size; i++)
		{
			uint64_t hash_index;
			if (Find_hash(Kmers[i], &hash_index, Kmer_hash))
				__sync_fetch_and_add(&Kmer_depth[hash_index], 1);
		}
		(argv -> Read_count)++;
	}
	//printf("T%i %f\n",argv -> thread_id, (float) total_time/(clock()-thd_begin_time));
}

int main_count(int argc, char ** argv)
{
	uint8_t thread_count = 0;
	char getopt_return;
	extern char *optarg;
	while ((getopt_return = getopt(argc, argv, "ht:")) != -1)
	{
		uint16_t len;
		switch (getopt_return)
		{
			case 'h':
				puts("Help information");
				return 1;
			case 't':
				thread_count = atoi(optarg);
				printf("[Option] Set %i threads\n",thread_count);
				break;
			case '?':
				puts("Option error, check help");
				return 1;
			default:
				return 1;
		}
	}
	FILE * Hash_file = fopen(argv[argc-3], "r");
	FILE * fasta_input = fopen(argv[argc-2], "r");
	if (!fasta_input) {
		puts("Input open fail");
	}
	FILE * OutFile = fopen(argv[argc-1], "w");
	fseek(Hash_file, 4, 0);
	fread(&Kmer_size, 1, 1, Hash_file);
	Set_Kmer_Size(Kmer_size);
	fseek(Hash_file, 8, 0);
	fread(&Hash_size, 8, 1, Hash_file);
	uint64_t first_idx;
	fread(&first_idx, 8, 1, Hash_file);
	printf("Hash Size: %i\nFirst location: %i\n", Hash_size, first_idx);
	Kmer_hash = (uint64_t *) malloc(Hash_size * sizeof(uint64_t));
	if (!Kmer_hash) {
		puts("Memory allocation failed");
		fclose(Hash_file);
		return 1;
	}
	printf("Read 0x%X hash\n",fread(Kmer_hash, sizeof(uint64_t), Hash_size, Hash_file));
	Kmer_depth = (uint16_t *) calloc(Hash_size, sizeof(uint16_t));
	if (!Kmer_depth) {
		puts("Memory allocation failed");
		free(Kmer_hash);
		fclose(Hash_file);
		return 1;
	}
	
	//Thread pool initialization
	uint8_t thd_idx;
	pthread_t *tid;
	struct FIFO_arg_struc *Thread_arg;
	if (thread_count)
	{
		tid = malloc(thread_count * sizeof(pthread_t));
		Thread_arg = malloc(thread_count * sizeof(struct FIFO_arg_struc));
		for (thd_idx = 0; thd_idx < thread_count; thd_idx++)
		{
			Thread_arg[thd_idx].Write_count = 0;
			Thread_arg[thd_idx].Read_count = 0;
			Thread_arg[thd_idx].thread_id = thd_idx;
			sem_init(&Thread_arg[thd_idx].data_feed_sem,0,0); //Initialize semaphore with first wait
			pthread_create(&tid[thd_idx], NULL, Kmer_count_TSK, &Thread_arg[thd_idx]);
		}
	}
	//End thread pool
	char line[350];
	thd_idx = 0;
	char pushed;
	uint32_t FIFO_write_idx = 0;
	uint64_t process_kmers = 0;
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
					if (kmer > encoded_r) kmer = encoded_r;
					if (thread_count)
					{
						//Thread modes
						if (Thread_arg[thd_idx].Write_count & 1)
							Thread_arg[thd_idx].FIFO1[FIFO_write_idx] = kmer;
						else Thread_arg[thd_idx].FIFO0[FIFO_write_idx] = kmer;
						FIFO_write_idx++;
						if (FIFO_write_idx == FIFO_size)
						{
							FIFO_write_idx = 0;
							Thread_arg[thd_idx].Write_count++;
							sem_post(&Thread_arg[thd_idx].data_feed_sem);
							do {
								thd_idx++;
								if (thd_idx == thread_count) thd_idx = 0;
							}
							while (Thread_arg[thd_idx].Write_count != Thread_arg[thd_idx].Read_count);
						}
					}
					else {
						uint64_t hash_index;
						if (Find_hash(kmer, &hash_index, Kmer_hash))
							Kmer_depth[hash_index]++;
					}
					process_kmers++;
					if ((process_kmers & 0x3FFFFFFF) == 0) printf("Read %iG kmers\n",process_kmers >> 30);
				}
			}
			seq_char_index++;
		}
	}
	//Join threads
	if (thread_count)
	{
		while (FIFO_write_idx != FIFO_size)
		{
			if (Thread_arg[thd_idx].Write_count & 1)
				Thread_arg[thd_idx].FIFO1[FIFO_write_idx] = 0;
			else Thread_arg[thd_idx].FIFO0[FIFO_write_idx] = 0;
			FIFO_write_idx++;
		}
		Thread_arg[thd_idx].Write_count++;
		sem_post(&Thread_arg[thd_idx].data_feed_sem);
		sleep(1);
		thread_no_more_data = 1;
		for (thd_idx = 0; thd_idx < thread_count; thd_idx++)
		{
			sem_post(&Thread_arg[thd_idx].data_feed_sem);
			pthread_join(tid[thd_idx], NULL);
			sem_destroy(&Thread_arg[thd_idx].data_feed_sem);
		}
		free(tid);
		free(Thread_arg);
	}
	
	//Dump count file
	printf("Pileup finish\nRead chain file %i entries\n",fread(Kmer_hash, sizeof(uint32_t), Hash_size, Hash_file));
	uint32_t * Kmer_next_loc = (uint32_t *) Kmer_hash;
	uint32_t buf_count = 0;
	uint16_t buffer[buffer_size];
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

void Recurse_edit(uint8_t Edits, uint64_t encoded, uint64_t encoded_r, uint8_t per_base, uint64_t kmer_idx)
{
	Edits--;
	for (char per_idx = 0; per_idx < per_base; per_idx++)
	{
		for (char per_value = 1; per_value < 4; per_value++)
		{
			uint64_t local_encoded = encoded;
			uint64_t local_encoded_r = encoded_r;
			Permute_kmer(&local_encoded, &local_encoded_r, per_idx, per_value);
			//Recursive edit the kmer
			if (Edits)
				Recurse_edit(Edits, local_encoded, local_encoded_r, per_idx, kmer_idx);
			//Find hash and finish
			uint64_t hash_index;
			if (local_encoded > local_encoded_r) local_encoded = local_encoded_r;
			if (Find_hash(local_encoded,&hash_index,Kmer_hash))
				Kmer_edit_depth[kmer_idx] += Kmer_occr[hash_index];
		}
	}
}

void * Kmer_filter_TSK(void *argvs)
{
	struct edit_dis_arg_struc *edit_dis_arg = (struct edit_dis_arg_struc *) argvs;
	uint64_t thread_cur_idx = edit_dis_arg -> start_idx;
	uint64_t thread_end_idx = edit_dis_arg -> end_idx;
	uint8_t thread_ID = edit_dis_arg -> thread_id;
	printf("Thread %d %u %u\n", thread_ID, thread_cur_idx, thread_end_idx);
	uint64_t count = 0;
	while (thread_cur_idx < thread_end_idx)
	{
		if (Kmer_occr[thread_cur_idx] == 1)
		{
			uint64_t encoded = Kmer_hash[thread_cur_idx];
			uint64_t encoded_r = Reverse_strand_encoded(encoded);
			Recurse_edit(edit_distance, encoded, encoded_r, Kmer_size, thread_cur_idx);
			count++;
			if ((count & 0xFFFFF) == 0) printf("Thread %u %uM\n", thread_ID, count >> 20);
		}
		thread_cur_idx++;
	}
	printf("Thread %i finished\n", thread_ID);
}

void hash_from_fasta(FILE * fasta){
	char fasta_buffer[200];
	uint8_t charge_size = 0;
	uint64_t encoded = 0;
	uint64_t encoded_r = 0;
	uint64_t processed = 0;
	uint32_t worst = 0;
	uint32_t hist[65536] = {0};
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
					if (collision < 65536) hist[collision]++;
					else hist[65535]++;
					Kmer_hash[hash_index] = kmer;
				}
				if (Kmer_occr[hash_index] < 255) Kmer_occr[hash_index]++;
			}
		}
		processed++;
		if (processed % 1666667 == 0){
			float average = 0;
			uint64_t count = 0;
			for (uint32_t k = 0; k < 65536; k++){
				count += hist[k];
				average += k * hist[k];
			}
			average /= count;
			printf("Processed %ubp, total %u Kmers, average collision %f\n", processed*60, count, average);
		}
	}
	float average = 0;
	uint64_t count = 0;
	for (uint32_t k = 0; k < 65536; k++){
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
	printf("Uniq count %u, total %u\n", unique_count, count);
}

uint64_t dump_kmer_list(char dump, char * kmerlist_fn, FILE * fasta){
	//Dump K-mer in second pass
	fseek(fasta, 0, SEEK_SET);
	FILE * Kmer_list;
	if (dump) Kmer_list = fopen(kmerlist_fn, "w");
	char fasta_buffer[200];
	uint8_t charge_size = 0;
	uint64_t encoded = 0;
	uint64_t encoded_r = 0;
	char chrom_name[128];
	char str_buf[256];
	uint64_t chr_pos;
	char Kmer_text[Kmer_size+1];
	uint64_t first_index;
	uint64_t last_index;
	uint64_t hash_index;
	uint64_t count = 0;
	while (fgets(fasta_buffer, 200, fasta) && fasta_buffer[0])
	{
		char * char_idx = fasta_buffer;
		if (*char_idx == '>')
		{
			charge_size = 0;
			encoded = 0;
			encoded_r = 0;
			chr_pos = 0;
			strncpy(chrom_name, char_idx+1, strlen(fasta_buffer)-2);
			continue;
		}
		while (*char_idx && *char_idx != '\n')
		{
			chr_pos++;
			memmove(Kmer_text,Kmer_text+1,Kmer_size-1);
			Kmer_text[Kmer_size-1] = *char_idx;
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
			if (charge_size < Kmer_size) charge_size++;
			if (kmer && charge_size == Kmer_size)
			{
				if (Find_hash(kmer, &hash_index, Kmer_hash))
				{
					if (!count) first_index = hash_index;
					else Kmer_next_index[last_index] = hash_index;
					last_index = hash_index;
					if (dump) {
						sprintf(str_buf, "%s\t%u\t%u\t%s\n", chrom_name, chr_pos-charge_size, chr_pos, Kmer_text);
						fputs(str_buf, Kmer_list);
					}
					count++;
				}
			}
		}
	}
	Kmer_next_index[last_index] = 0xFFFFFFFF;
	printf("Total output %i k-mers\n", count);
	if (dump) fclose(Kmer_list);
	return first_index;
}

int main_search(int argc, char ** argv)
{
	uint8_t thread_count = 1;
	char getopt_return;
	char * kmerlist_fn;
	extern char *optarg;
	char dump_kmer = 0;
	while ((getopt_return = getopt(argc, argv, "hk:t:s:e:d:f:")) != -1)
	{
		uint16_t len;
		switch (getopt_return)
		{
			case 'h':
				puts("Help information");
				return 1;
			case 'k':
				Set_Kmer_Size(atoi(optarg));
				printf("[Option] Set %i-mer\n",Kmer_size);
				break;
			case 't':
				thread_count = atoi(optarg);
				printf("[Option] Set %i threads\n",thread_count);
				break;
			case 's':
				len = strlen(optarg);
				switch (optarg[len-1]){
					case 'G':
						Hash_size = (uint64_t) atoi(optarg) << 30;
						break;
					case 'M':
						Hash_size = atoi(optarg) << 20;
						break;
					case 'K':
						Hash_size = atoi(optarg) << 10;
						break;
					default:
						Hash_size = atoi(optarg);
				}
				printf("[Option] Set hash space 0x%lX\n",Hash_size);
				break;
			case 'e':
				edit_distance = atoi(optarg);
				printf("[Option] Edit distance %i\n",edit_distance);
				break;
			case 'd':
				Edit_depth_thres = atoi(optarg);
				printf("[Option] Max repeat count with edit distance %i\n",Edit_depth_thres);
				break;
			case 'f':
				kmerlist_fn = optarg;
				dump_kmer = 1;
				printf("[Option] Final k-mer list enabled: %s\n",kmerlist_fn);
				break;
			case '?':
				puts("Option error, check help");
				return 1;
			default:
				return 1;
		}
	}
	
	//Malloc
	Kmer_hash = (uint64_t *) calloc(Hash_size + 16384, sizeof(uint64_t));
	Kmer_occr = (uint8_t *) calloc(Hash_size + 16384, sizeof(uint8_t));
	Kmer_edit_depth = (uint8_t *) calloc(Hash_size + 16384, sizeof(uint8_t));
	if (!Kmer_hash || !Kmer_edit_depth || !Kmer_occr)
	{
		puts("Memory allocation failed");
		return 1;
	}
	FILE * fasta = fopen(argv[argc-2], "r");
	hash_from_fasta(fasta);

	//Filter
	if (edit_distance)
	{
		uint64_t thread_seg_count = Hash_size / thread_count;
		uint64_t Thread_start_idx = 0;
		uint8_t thd_idx;
		pthread_t *tid = malloc(thread_count * sizeof(pthread_t));
		struct edit_dis_arg_struc *arg_arr = malloc(thread_count * sizeof(struct edit_dis_arg_struc));
		for (thd_idx = 1; thd_idx < thread_count; thd_idx++)
		{
			arg_arr[thd_idx].start_idx = Thread_start_idx;
			arg_arr[thd_idx].end_idx = Thread_start_idx + thread_seg_count;
			arg_arr[thd_idx].thread_id = thd_idx;
			Thread_start_idx += thread_seg_count;
			pthread_create(&tid[thd_idx], NULL, Kmer_filter_TSK, &arg_arr[thd_idx]);
		}
		arg_arr[0].start_idx = Thread_start_idx;
		arg_arr[0].end_idx = Hash_size + 16384;
		arg_arr[0].thread_id = 0;
		Kmer_filter_TSK(&arg_arr[0]); //Default thread
		//Join threads
		for (thd_idx = 1; thd_idx < thread_count; thd_idx++)
			pthread_join(tid[thd_idx], NULL);
		free(tid);
		free(arg_arr);
	}
	
	//High frequency mers stats and removal
	uint64_t occr_idx;
	uint64_t count = 0;
	for (occr_idx = 0; occr_idx < Hash_size+16384; occr_idx++)
	{
		if (Kmer_occr[occr_idx] == 0) continue;
		if (Kmer_occr[occr_idx] > 1 || Kmer_edit_depth[occr_idx] >= Edit_depth_thres)
		{
			//Kmers to get rid off
			Kmer_hash[occr_idx] = 0;
			count++;
		}
		else {
			uint64_t cur_kmer = Kmer_hash[occr_idx];
			Kmer_hash[occr_idx] = 0;
			uint64_t new_idx = DJBHash_encode(cur_kmer) & (Hash_size-1);
			while (Kmer_hash[new_idx] && Kmer_hash[new_idx] != cur_kmer) new_idx++;
			Kmer_hash[new_idx] = cur_kmer;
		}
	}
	free(Kmer_edit_depth);
	free(Kmer_occr);
	printf("%u kmers deleted\n", count);
	//Dump hash table into text format
	Kmer_next_index = (uint32_t *) calloc(Hash_size + 16384,sizeof(uint32_t));
	puts("Generate chaining table");
	uint64_t first_index = dump_kmer_list(dump_kmer, kmerlist_fn, fasta);
	//Save QM index binary
	FILE * Hash_file = fopen(argv[argc-1], "w");
	if (!Hash_file) {
		puts("File creation failed");
		return 1;
	}
	rewind(Hash_file);
	const char Version[5] = "QM11";
	fwrite(Version, 1, 4, Hash_file);
	fwrite(&Kmer_size, 1, 1, Hash_file);
	fwrite(&edit_distance, 1, 1, Hash_file);
	fwrite(&Edit_depth_thres, 1, 1, Hash_file);
	fwrite(&Edit_depth_thres, 1, 1, Hash_file);
	uint64_t Hashsize = Hash_size + 16384;
	fwrite(&Hashsize, 8, 1, Hash_file);
	fwrite(&first_index, 8, 1, Hash_file);
	fwrite((void *) Kmer_hash, sizeof(uint64_t), Hash_size + 16384, Hash_file);
	free(Kmer_hash);
	fwrite((void *) Kmer_next_index, sizeof(uint32_t), Hash_size + 16384, Hash_file);
	free(Kmer_next_index);
	fclose(Hash_file);
	puts("Kmer search finished!");
	return 0;
}

void printversion() {
	puts("QuicK-mer 2.0");
	puts("Operation modes: \n\tindex\tIndex a bed format kmer list\n\tcount\tCNV estimate from library\n\tsearch\tSearch K-kmer in genome\n");
}

int main(int argc, char** argv)
{
	if (argc > 1) {
		if (strcmp(argv[1], "index") == 0)
			return main_hash(argc-2, argv+2);
		else if (strcmp(argv[1], "count") == 0)
			return main_count(argc-1, argv+1);
		else if (strcmp(argv[1], "search") == 0)
			return main_search(argc-1, argv+1);
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
