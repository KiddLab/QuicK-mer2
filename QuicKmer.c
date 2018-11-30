#include "math.h"
#include "pthread.h"
#include "semaphore.h"
#include "stdint.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"
#include "unistd.h"

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
	int64_t scan_direction;
	if (*hash_index & (Hash_size >> 1)) scan_direction = -1;
	else scan_direction = 1;
	while (Hash_dict[*hash_index] && Hash_dict[*hash_index] != encoded_kmer)
		(*hash_index) += scan_direction;
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
	FILE *kmer_list = fopen(argv[argc - 2], "r");
	FILE * Hash_file = fopen(argv[argc - 1], "w");
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
	Kmer_hash = (uint64_t *) malloc(Hash_size * sizeof(uint64_t));
	Kmer_next_index = (uint32_t *) malloc(sizeof(uint32_t) * Hash_size);
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
		int64_t scan_direction;
		if (hash_index & (Hash_size >> 1)) scan_direction = -1;
		else scan_direction = 1;
		while (Kmer_hash[hash_index] != 0){
			hash_index+= scan_direction;
			worst++;
		}
		Kmer_hash[hash_index] = kmer_encoded;
		if (!Total_count) first_index = hash_index;
		else Kmer_next_index[last_index] = hash_index;
		last_index = hash_index;
		if (worst > worst_case) {
			worst_case = worst;
			printf("%s Worst %i\n", chrom, worst);
		}
		if (worst > 8191) worst = 8191;
		hist[worst]++;
		Total_count++;
	}
	Kmer_next_index[last_index] = first_index;
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
	uint64_t Hashsize = Hash_size;
	fwrite(&Hashsize, 8, 1, Hash_file);
	fwrite(&first_index, 8, 1, Hash_file);
	fwrite((void *) Kmer_hash, sizeof(uint64_t), Hash_size, Hash_file);
	int count_kmer = 0;
	for (uint64_t k = 0; k < Hash_size; k++){
		if (Kmer_hash[k]) count_kmer++;
	}
	printf("Total %i kmers\n", count_kmer);
	free(Kmer_hash);
	fwrite((void *) Kmer_next_index, sizeof(uint32_t), Hash_size, Hash_file);
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

void help_count(){
	puts("\nquicKmer count [Options] ref.fa sample.fast[a/q] Out_prefix\n\nOptions:");
	puts("-h\t\tShow this help information");
	puts("-t [num]\tNumber of threads");
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
				help_count();
				return 1;
			case 't':
				thread_count = atoi(optarg);
				printf("[Option] Set %u threads\n",thread_count);
				break;
			case '?':
				puts("Option error, check help");
				help_count();
				return 1;
			default:
				return 1;
		}
	}
	char path[65536];
	strcpy(path, argv[argc-3]);
	strcat(path, ".qm");
	FILE * Hash_file = fopen(path, "r");
	FILE * fasta_input = fopen(argv[argc-2], "r");
	if (!fasta_input) {
		puts("Input open fail");
	}
	strcpy(path, argv[argc-1]);
	strcat(path, ".bin");
	FILE * OutFile = fopen(path, "w");
	fseek(Hash_file, 4, 0);
	fread(&Kmer_size, 1, 1, Hash_file);
	Set_Kmer_Size(Kmer_size);
	fseek(Hash_file, 8, 0);
	fread(&Hash_size, 8, 1, Hash_file);
	uint64_t first_idx;
	fread(&first_idx, 8, 1, Hash_file);
	printf("Hash Size: 0x%lX\nFirst location: 0x%lX\n", Hash_size, first_idx);
	Kmer_hash = (uint64_t *) malloc(Hash_size * sizeof(uint64_t));
	if (!Kmer_hash) {
		puts("Memory allocation failed");
		fclose(Hash_file);
		return 1;
	}
	printf("Read 0x%lX hash\n",fread(Kmer_hash, sizeof(uint64_t), Hash_size, Hash_file));
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
	time_t start_time, end_time;
	time(&start_time);
	char line[350];
	thd_idx = 0;
	char pushed;
	uint32_t FIFO_write_idx = 0;
	uint64_t process_kmers = 0;
	fgets(line, 350, fasta_input);
	char set_fastq = 0;
	if (line[0] == '@') set_fastq = 1;
	else fseek(fasta_input, 0, SEEK_SET);
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
					if ((process_kmers & 0x3FFFFFFF) == 0) printf("Read %liG kmers\n",process_kmers >> 30);
				}
			}
			seq_char_index++;
		}
		if (set_fastq){
			fgets(line, 350, fasta_input);
			fgets(line, 350, fasta_input);
			fgets(line, 350, fasta_input);
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
	time(&end_time);
	printf("Counting elapse %u sec, total %lu kmers\n", end_time - start_time, process_kmers);
	//Dump count file
	printf("Pileup finish\nRead chain file %lu entries\n",fread(Kmer_hash, sizeof(uint32_t), Hash_size, Hash_file));
	strcpy(path, argv[argc-3]);
	strcat(path, ".qgc");
	FILE * GC_control = fopen(path, "r");
	if (GC_control == NULL)
		printf("GC control file %s absent. Continue without GC correction!\n", path);
	
	uint32_t * Kmer_next_loc = (uint32_t *) Kmer_hash;
	uint32_t buf_count = 0;
	uint16_t buffer[buffer_size];
	uint16_t gc_control[buffer_size];
	uint32_t chaining_idx = first_idx;
	double Control_curve[401] = {0};
	double Control_std[401] = {0};
	uint32_t Control_count[401] = {0};
	do {
		buffer[buf_count] = Kmer_depth[chaining_idx];
		chaining_idx = Kmer_next_loc[chaining_idx];
		if (GC_control != NULL) {
			if ((buf_count & (buffer_size -1)) == 0)
				fread(gc_control, sizeof(uint16_t), buffer_size, GC_control);
			if (gc_control[buf_count] & 0x8000) {
				Control_curve[gc_control[buf_count] & 0x1FF] += buffer[buf_count];
				Control_count[gc_control[buf_count] & 0x1FF]++;
				Control_std[gc_control[buf_count] & 0x1FF] += buffer[buf_count] * buffer[buf_count];
			}
		}
		buf_count++;
		if (buf_count == buffer_size) {
			fwrite(buffer, sizeof(uint16_t), buffer_size, OutFile);
			buf_count = 0;
		}
	}
	while(first_idx != chaining_idx);
	fwrite(buffer, sizeof(uint16_t), buf_count, OutFile);
	fclose(OutFile);
	free(Kmer_hash);
	free(Kmer_depth);
	//Save GC curve
	if (GC_control != NULL) {
		strcpy(path, argv[argc-1]);
		strcat(path, ".txt");
		FILE * GC_curve = fopen(path, "w");
		char str_buf[255];
		double total_depth = 0;
		uint64_t total_count = 0;
		for (uint16_t i = 0; i < 401; i++){
			total_count += Control_count[i];
			total_depth += Control_curve[i];
			if (Control_count[i]) {
				Control_curve[i] /= Control_count[i];
				Control_std[i] = Control_std[i] / Control_count[i] - Control_curve[i] * Control_curve[i];
			}
			sprintf(str_buf, "%.2f\t%f\t%i\t%f\n", i/4.0, Control_curve[i], Control_count[i], Control_std[i]);
			fputs(str_buf, GC_curve);
		}
		total_depth /= total_count;
		printf("Mean sequencing depth: %.2f\n",total_depth);
		fclose(GC_curve);
	}
	puts("Exit QuicK-mer count");
	return 0;
}

void help_est(){
	puts("quicKmer est ref.fa sample_prefix output.bed");
	puts("\tref.fa\t\tPrefix to genome reference. Program requires .qgc and .bed definition");
	puts("\tsample_prefix\tPrefix to sample.bin");
	puts("\toutput.bed\tOutput bedfile for copy number");
	puts("\nNo options available\n");
}

int main_estimate(int argc, char ** argv){
	FILE * GC_ctrl;
	FILE * Depth_Bin;
	FILE * Window_file;
	FILE * OutputFile;
	char path[65536];
	strcpy(path, argv[argc-3]);
	strcat(path, ".qgc");
	GC_ctrl = fopen(path, "r");
	if (!GC_ctrl) {
		puts("GC control file missing.");
		help_est();
		return 1;
	}
	strcpy(path, argv[argc-3]);
	strcat(path, ".bed");
	Window_file = fopen(path,"r");
	if (!Window_file) {
		puts("Window segmentation file missing.");
		help_est();
		return 1;
	}
	strcpy(path, argv[argc-2]);
	strcat(path, ".bin");
	Depth_Bin = fopen(path,"r");
	if (!Depth_Bin) {
		puts("Sample file missing (Use prefix).");
		help_est();
		return 1;
	}
	strcpy(path, argv[argc-2]);
	strcat(path, ".txt");
	FILE * GC_curve = fopen(path,"r");
	char str_buf[255];
	double total_depth = 0; //Later store average depth
	uint64_t total_count = 0;
	uint16_t gc_control[buffer_size];
	uint16_t depth[buffer_size];
	if (!GC_curve){
		puts("Depth control not found. Regenerating...");
		FILE * GC_curve = fopen(path, "w");
		uint32_t read_size;
		double Control_curve[401] = {0};
		double Control_std[401] = {0};
		uint32_t Control_count[401] = {0};
		while (read_size = fread(gc_control, 1, buffer_size, GC_ctrl)) {
			fread(depth, 1, read_size, Depth_Bin);
			read_size >>= 1;
			while (read_size){
				if (gc_control[read_size] & 0x8000){
					Control_curve[gc_control[read_size] & 0x1FF] += depth[read_size];
					Control_count[gc_control[read_size] & 0x1FF]++;
					Control_std[gc_control[read_size] & 0x1FF] += depth[read_size] * depth[read_size];
				}
				read_size -= 1;
			}
		}
		//Generate GC curve
		for (uint16_t i = 0; i < 401; i++){
			total_count += Control_count[i];
			total_depth += Control_curve[i];
			if (Control_count[i]) {
				Control_curve[i] /= Control_count[i];
				Control_std[i] = Control_std[i] / Control_count[i] - Control_curve[i] * Control_curve[i];
			}
			sprintf(str_buf, "%.2f\t%f\t%i\t%f\n", i/4.0, Control_curve[i], Control_count[i], Control_std[i]);
			fputs(str_buf, GC_curve);
		}
	}
	else {
		float percent, depth;
		uint32_t cur_count;
		while (fscanf(GC_curve, "%f\t%f\t%i\t%s\n", &percent, &depth, &cur_count, str_buf) == 4) {
			total_depth += depth * cur_count;
			total_count += cur_count;
		}
	}
	total_depth /= total_count;
	printf("Mean sequencing depth: %.2f\n",total_depth);
	fclose(GC_curve);
	strcpy(path, "smooth_GC_mrsfast.py ");
	strcat(path, argv[argc-2]);
	strcat(path, ".txt");
	//Open a pipe to call lowess smoothing
	FILE * pipe_in = popen(path, "r");
	OutputFile = fopen(argv[argc-1],"w");
	float correction_curve[401];
	fread(correction_curve,4,401,pipe_in);
	pclose(pipe_in);
	//Reset read pointer
	fseek(GC_ctrl, 0, SEEK_SET);
	fseek(Depth_Bin, 0, SEEK_SET);
	double cur_window_depth = 0.0;
	uint64_t kmer_idx = 0;
	char chrom[64];
	char win_begin[64];
	char win_end[64];
	uint32_t cur_window_right, cur_window_left;
	fscanf(Window_file, "%s\t%s\t%s\t%u\t%u\n", chrom, win_begin, win_end, &cur_window_left, &cur_window_right);
	uint32_t read_size;
	while (read_size = fread(gc_control, 1, buffer_size, GC_ctrl)) {
		fread(depth, 1, read_size, Depth_Bin);
		read_size >>= 1;
		uint32_t cur_idx = 0;
		while (cur_idx < read_size){
			if (kmer_idx > cur_window_right) {
				cur_window_depth /= cur_window_right - cur_window_left;
				cur_window_depth /= total_depth/2;
				sprintf(str_buf, "%s\t%s\t%s\t%f\n",chrom, win_begin, win_end,cur_window_depth);
				fputs(str_buf, OutputFile);
				if (fscanf(Window_file, "%s\t%s\t%s\t%u\t%u\n", chrom, win_begin, win_end, &cur_window_left, &cur_window_right)
					!= 5) break;
				cur_window_depth = 0.0;
			}
			if (kmer_idx < cur_window_right && kmer_idx >= cur_window_left){
				cur_window_depth += correction_curve[gc_control[cur_idx] & 0x1FF] * depth[cur_idx];
			}
			cur_idx++;
			kmer_idx++;
		}
	}
	fclose(OutputFile);
	return 0;
}

char Recurse_edit(uint8_t Edits, uint64_t encoded, uint64_t encoded_r, uint8_t per_base, uint64_t kmer_idx)
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
				if (Recurse_edit(Edits, local_encoded, local_encoded_r, per_idx, kmer_idx)) return 1;
			//Find hash and finish
			uint64_t hash_index;
			if (local_encoded > local_encoded_r) local_encoded = local_encoded_r;
			if (Find_hash(local_encoded,&hash_index,Kmer_hash)) {
				if (Kmer_edit_depth[kmer_idx] + Kmer_occr[hash_index] > Edit_depth_thres) {
					Kmer_edit_depth[kmer_idx] = Edit_depth_thres + 1;
					return 1;
				}
				else Kmer_edit_depth[kmer_idx] += Kmer_occr[hash_index];
			}
		}
	}
	return 0;
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

void Resize_hash_table(uint64_t old_size, uint64_t new_size){
	printf("Resize %f\n", (float)new_size/old_size);
	uint64_t index = 0;
	uint64_t hash_index;
	void * realloc_ptr;
	if (new_size > old_size) {
		if (realloc_ptr = realloc(Kmer_hash, new_size * sizeof(uint64_t))) Kmer_hash = realloc_ptr;
		else {
			puts("Realloc failed!\n");
			return;
		}
		if (Kmer_occr && (realloc_ptr = realloc(Kmer_occr, new_size))) Kmer_occr = realloc_ptr;
		else {
			puts("Realloc failed!\n");
			return;
		}
		puts("Allocation success");
		memset(Kmer_hash+old_size, 0, sizeof(uint64_t)*(new_size - old_size));
		Hash_size = new_size;
		index = old_size - 1;
		while (index >= (old_size >> 1)){
			if (Kmer_hash[index]){
				Find_hash(Kmer_hash[index], &hash_index, Kmer_hash);
				if (hash_index != index) {
					Kmer_hash[hash_index] = Kmer_hash[index];
					Kmer_hash[index] = 0;
					if (Kmer_occr) {
						Kmer_occr[hash_index] = Kmer_occr[index];
						Kmer_occr[index] = 0;
					}
				}
			}
			index--;
		}
		index = 0;
		while (index < (old_size >> 1)) {
			if (Kmer_hash[index]){
				Find_hash(Kmer_hash[index], &hash_index, Kmer_hash);
				if (hash_index != index) {
					Kmer_hash[hash_index] = Kmer_hash[index];
					Kmer_hash[index] = 0;
					if (Kmer_occr) {
						Kmer_occr[hash_index] = Kmer_occr[index];
						Kmer_occr[index] = 0;
					}
				}
			}
			index++;
		}
	}
	else {
		Hash_size = new_size;
		while (index < (old_size >> 1)) {
			if (Kmer_hash[index]){
				Find_hash(Kmer_hash[index], &hash_index, Kmer_hash);
				if (hash_index != index) {
					Kmer_hash[hash_index] = Kmer_hash[index];
					Kmer_hash[index] = 0;
					if (Kmer_occr) {
						Kmer_occr[hash_index] = Kmer_occr[index];
						Kmer_occr[index] = 0;
					}
				}
			}
			index++;
		}
		index = old_size - 1;
		while (index >= (old_size >> 1)){
			if (Kmer_hash[index]){
				Find_hash(Kmer_hash[index], &hash_index, Kmer_hash);
				if (hash_index != index) {
					Kmer_hash[hash_index] = Kmer_hash[index];
					Kmer_hash[index] = 0;
					if (Kmer_occr) {
						Kmer_occr[hash_index] = Kmer_occr[index];
						Kmer_occr[index] = 0;
					}
				}
			}
			index--;
		}
		if (realloc_ptr = realloc(Kmer_hash, new_size * sizeof(uint64_t))) Kmer_hash = realloc_ptr;
		if (realloc_ptr = realloc(Kmer_occr, new_size)) Kmer_occr = realloc_ptr;
	}
}

void hash_from_fasta(FILE * fasta){
	char fasta_buffer[200];
	uint8_t charge_size = 0;
	uint64_t encoded = 0;
	uint64_t encoded_r = 0;
	uint64_t processed = 0;
	uint64_t count = 0;
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
			if (charge_size < Kmer_size) charge_size++;
			if (kmer && charge_size == Kmer_size)
			{
				uint64_t hash_index = DJBHash_encode(kmer) & (Hash_size - 1);
				int64_t scan_direction;
				if (hash_index & (Hash_size >> 1)) scan_direction = -1;
				else scan_direction = 1;
				//Add to hash memory
				uint32_t collision = 0;
				while (Kmer_hash[hash_index] && Kmer_hash[hash_index] != kmer)
				{
					hash_index += scan_direction;
					collision++;
				}
				if (!Kmer_hash[hash_index])
				{
					count++;
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
		if (count > 0.8 * Hash_size) {
			Resize_hash_table(Hash_size, Hash_size << 1);
			memset(hist, 0, 4 * 65536);
			worst = 0;
		}
		processed++;
		if (processed % 1666667 == 0){
			float average = 0;
			uint64_t sum = 0;
			for (uint32_t k = 0; k < 65536; k++){
				sum += hist[k];
				average += k * hist[k];
			}
			average /= sum;
			printf("Processed %ubp, total %u Kmers, average collision %f\n", processed*60, count, average);
		}
	}
	float average = 0;
	uint64_t sum = 0;
	for (uint32_t k = 0; k < 65536; k++){
		sum += hist[k];
		average += k * hist[k];
	}
	printf("Average %f, fill %f\% \n", average/sum, ((float) count * 100)/ Hash_size);
	uint64_t occr_idx = 0;
	uint64_t unique_count = 0;
	while (occr_idx < Hash_size)
	{
		if (Kmer_occr[occr_idx] == 1) unique_count++;
		occr_idx++;
	}
	printf("Uniq count %lu, total %lu\n", unique_count, count);
}

uint64_t dump_kmer_list(FILE * Kmer_list, FILE * fasta, FILE * window_file, FILE * ctrl_file,
						FILE* GC_binout, uint16_t gc_win, uint16_t wsize){
	//Dump K-mer in second pass
	fseek(fasta, 0, SEEK_SET);
	char next_chr_name[200];
	char * chrom_name = next_chr_name+1;
	char str_buf[200];
	char ctrl_chr[200];
	uint32_t chr_pos;
	uint64_t win_start;
	uint32_t wstart;
	uint64_t count = 0;
	uint64_t first_index, last_index;
	char * fa_buf;
	uint16_t * GC_control_buf;
	uint32_t GC_control_buf_idx = 0;
	if (GC_binout != NULL) GC_control_buf = malloc(buffer_size*sizeof(uint16_t));
	fa_buf = malloc(256*1024*1024);
	while (fgets(next_chr_name, 200, fasta))
	{
		uint8_t charge_size = 0;
		uint64_t encoded = 0;
		uint64_t encoded_r = 0;
		if (next_chr_name[0] == '>')
		{
			charge_size = 0;
			wstart = count;
			encoded = 0;
			encoded_r = 0;
			chr_pos = 0;
			win_start = 0;
			while (!feof(fasta)) {
				uint64_t loc = ftell(fasta);
				fgets(&fa_buf[chr_pos], 200, fasta);
				uint8_t read_count = ftell(fasta)-loc;
				if (fa_buf[chr_pos] == '>') {
					fseek(fasta, -read_count, SEEK_CUR);
					fa_buf[chr_pos] = 0;
					break;
				}
				else {
					chr_pos += read_count - 1;
					if (fa_buf[chr_pos] != '\n') chr_pos++;
				}
			}
		}
		uint8_t ctrl_checked_same_chr = 0;
		uint8_t ctrl_this_chr_absent = 0; 
		uint16_t window_bp = 0;
		uint16_t GC_bp = 0;
		uint16_t N_bp = 0;
		uint32_t tl_idx = 0; //Trailing or leading index
		uint32_t chr_bp_count = chr_pos;
		uint32_t ctrl_w_s, ctrl_w_e;
		chrom_name[strlen(chrom_name)-1] = 0;
		printf("%s %i\n", chrom_name, chr_bp_count);
		//Charge kmer gc window
		for (chr_pos = 0; chr_pos < (gc_win - Kmer_size)/2; chr_pos++) {
			if (fa_buf[chr_pos] & 2) GC_bp++;
			else if (fa_buf[chr_pos] == 'N') N_bp++;
			window_bp++;
		}
		chr_pos = 0;
		char last_ctrl = 0;
		while (fa_buf[chr_pos])
		{
			if (chr_pos + (gc_win - Kmer_size) /2 < chr_bp_count){
				tl_idx = chr_pos + (gc_win - Kmer_size) /2;
				if (fa_buf[tl_idx] & 2) GC_bp++; //Is G or C
				else if (fa_buf[tl_idx] == 'N') N_bp++;
				window_bp++;
			}
			if (chr_pos >= (gc_win + Kmer_size)/2) {
				tl_idx = chr_pos - (gc_win + Kmer_size)/2;
				if (fa_buf[tl_idx] & 2) GC_bp--;
				else if (fa_buf[tl_idx] == 'N') N_bp--;
				window_bp--;
			}
			uint8_t Letter = (fa_buf[chr_pos] >> 1) & 3;
			if (fa_buf[chr_pos] == 'N'){
				charge_size = 0;
				encoded = 0;
				encoded_r = 0;
			}
			else {
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
					uint64_t hash_index;
					if (Find_hash(kmer, &hash_index, Kmer_hash))
					{
						if (GC_binout != NULL) {
							//GC calculation
							uint16_t GC_per = (N_bp + GC_bp * 2) * 200/ window_bp; //Use 401 bins
							//Control Window generation
							if (!ctrl_checked_same_chr) {
								fseek(ctrl_file, 0, SEEK_SET); //In case the bed file is not in the same chr sort order
								while (fscanf(ctrl_file, "%s\t%d\t%d", ctrl_chr, &ctrl_w_s, &ctrl_w_e) == 3 && 
									(ctrl_this_chr_absent = strcmp(ctrl_chr, chrom_name)));
								ctrl_checked_same_chr = 1;
							}
							else while (chr_pos > ctrl_w_e && !ctrl_this_chr_absent &&
								fscanf(ctrl_file, "%s\t%d\t%d", ctrl_chr, &ctrl_w_s, &ctrl_w_e) == 3)
								if (strcmp(ctrl_chr, chrom_name)) {
									ctrl_this_chr_absent = 1;
									break;
								}
							if (ctrl_checked_same_chr && !ctrl_this_chr_absent && chr_pos + 1 - Kmer_size > ctrl_w_s)
								GC_per |= 0x8000;
							GC_control_buf[GC_control_buf_idx] = GC_per;
							GC_control_buf_idx++;
							if ((GC_control_buf_idx & (buffer_size-1)) == 0) {
								GC_control_buf_idx = 0;
								fwrite(GC_control_buf, 2, buffer_size, GC_binout);
							}
						}
						//Chaining
						if (!count) first_index = hash_index;
						else Kmer_next_index[last_index] = hash_index;
						last_index = hash_index;
						count++;
						//Window generation
						if (window_file != NULL && count % wsize == 0){
							sprintf(str_buf, "%s\t%u\t%u\t%u\t%u\n", chrom_name, win_start, chr_pos, wstart, count);
							fputs(str_buf, window_file);
							win_start = chr_pos;
							wstart = count;
						}
					}
				}
			}
			chr_pos++;
		}
	}
	if (GC_binout != NULL) {
		fwrite(GC_control_buf, 2, GC_control_buf_idx, GC_binout);
		fclose(GC_binout);
	}
	Kmer_next_index[last_index] = first_index;
	printf("Total output %u k-mers\n", count);
	return first_index;
}

void help_search(){
	puts("\nquicKmer search [Options] ref.fa\n\nOptions:");
	puts("-h\t\tShow this help information");
	puts("-k [num]\tSize of K-mer. Must be between 3-32. Default 30");
	puts("-t [num]\tNumber of threads for edit distance search");
	puts("-s [num]\tSize of hash dictionary. Can use suffix G,M,K");
	puts("-e [num]\tEdit distance search. 0, 1 or 2");
	puts("-d [num]\tEdit distance depth threshold to keep. 1..255, Default 100");
	puts("-w [num]\tOutput window definition size. Default 1000");
	puts("-c [filename]\tInput bedfile for GC control regions");
	puts("");
}

int main_search(int argc, char ** argv)
{
	uint8_t thread_count = 1;
	char getopt_return;
	FILE * kmerdump;
	extern char *optarg;
	char dump_kmer = 0;
	FILE * window_file = NULL;
	FILE * GC_file = NULL;
	FILE * Control_bed = NULL;
	uint32_t window_size = 1000;
	if (argc < 2){
		help_search();
		return 1;
	}
	while ((getopt_return = getopt(argc, argv, "hk:t:s:e:d:f:w:c:")) != -1)
	{
		uint16_t len;
		switch (getopt_return)
		{
			case 'h':
				help_search();
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
				Hash_size = (uint64_t) 1 << (uint8_t) ceil(log2(Hash_size));
				printf("[Option] Set hash space 0x%lX\n",Hash_size);
				break;
			case 'c':
				Control_bed = fopen(optarg, "r");
				printf("[Option] Control bedfile %s\n",optarg);
				break;
			case 'w':
				window_size = atoi(optarg);
				printf("[Option] Write window size %s\n",optarg);
				break;
			case 'e':
				edit_distance = atoi(optarg);
				printf("[Option] Edit distance %i\n",edit_distance);
				break;
			case 'd':
				Edit_depth_thres = atoi(optarg);
				printf("[Option] Max repeat count with edit distance %i\n",Edit_depth_thres);
				break;
			case '?':
				puts("Option error, check help");
				help_search();
				return 1;
			default:
				return 1;
		}
	}
	
	//Malloc
	Kmer_hash = (uint64_t *) calloc(Hash_size, sizeof(uint64_t));
	Kmer_occr = (uint8_t *) calloc(Hash_size, sizeof(uint8_t));
	if (!Kmer_hash || !Kmer_occr)
	{
		puts("Memory allocation failed");
		return 1;
	}
	printf("ADDR %lX %lX %lX\n", Kmer_hash, Kmer_occr, Kmer_edit_depth);
	char path[65535];
	strcpy(path, argv[argc-1]);
	strcat(path, ".qm");
	FILE * Hash_file = fopen(path, "w");
	if (!Hash_file) {
		puts("File creation failed");
		return 1;
	}
	strcpy(path, argv[argc-1]);
	strcat(path, ".bed");
	window_file = fopen(path, "w");
	if (!window_file) {
		puts("Window creation failed");
		return 1;
	}
	FILE * fasta = fopen(argv[argc-1], "r");
	hash_from_fasta(fasta);

	//Filter
	if (edit_distance)
	{
		Kmer_edit_depth = (uint8_t *) calloc(Hash_size, sizeof(uint8_t));
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
		arg_arr[0].end_idx = Hash_size;
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
	for (occr_idx = 0; occr_idx < Hash_size; occr_idx++)
	{
		if (Kmer_occr[occr_idx] > 1 || (Kmer_edit_depth && Kmer_edit_depth[occr_idx] >= Edit_depth_thres))
		{
			//Kmers to get rid off
			Kmer_hash[occr_idx] = 0;
			count++;
		}
	}
	free(Kmer_edit_depth);
	free(Kmer_occr);
	printf("%u kmers deleted\n", count);
	uint64_t middle_buffer[1024] = {0};
	uint32_t middle_count = 0;
	for (occr_idx = 0; occr_idx < (Hash_size >> 1); occr_idx++)
	{
		if (Kmer_hash[occr_idx])
		{
			uint64_t cur_kmer = Kmer_hash[occr_idx];
			Kmer_hash[occr_idx] = 0;
			uint64_t new_idx = DJBHash_encode(cur_kmer) & (Hash_size-1);
			if (new_idx >= (Hash_size >> 1)){
				middle_buffer[middle_count] = cur_kmer;
				middle_count++;
			}
			else {
				Find_hash(cur_kmer, &new_idx, Kmer_hash);
				Kmer_hash[new_idx] = cur_kmer;
			}
		}
		uint64_t mirror_idx = Hash_size - 1 - occr_idx;
		if (Kmer_hash[mirror_idx])
		{
			uint64_t cur_kmer = Kmer_hash[mirror_idx];
			Kmer_hash[mirror_idx] = 0;
			uint64_t new_idx = DJBHash_encode(cur_kmer) & (Hash_size-1);
			if (new_idx < (Hash_size >> 1)){
				middle_buffer[middle_count] = cur_kmer;
				middle_count++;
			}
			else {
				Find_hash(cur_kmer, &new_idx, Kmer_hash);
				Kmer_hash[new_idx] = cur_kmer;
			}
		}
	}
	while (middle_count) {
		middle_count--;
		uint64_t new_idx = DJBHash_encode(middle_buffer[middle_count]) & (Hash_size-1);
		Find_hash(middle_buffer[middle_count], &new_idx, Kmer_hash);
		Kmer_hash[new_idx] = middle_buffer[middle_count];
	}
	//Dump hash table into text format
	Kmer_next_index = (uint32_t *) calloc(Hash_size,sizeof(uint32_t));
	puts("Generate chaining table");
	strcpy(path, argv[argc-1]);
	if (Control_bed != NULL) {
		strcat(path, ".qgc");
		GC_file = fopen(path, "w");
	}
	uint64_t first_index = dump_kmer_list(kmerdump, fasta, window_file, Control_bed, GC_file, 400, window_size);
	puts("Writing index file");
	//Save QM index binary
	
	rewind(Hash_file);
	const char Version[5] = "QM11";
	fwrite(Version, 1, 4, Hash_file);
	fwrite(&Kmer_size, 1, 1, Hash_file);
	fwrite(&edit_distance, 1, 1, Hash_file);
	fwrite(&Edit_depth_thres, 1, 1, Hash_file);
	fwrite(&Edit_depth_thres, 1, 1, Hash_file);
	fwrite(&Hash_size, 8, 1, Hash_file);
	fwrite(&first_index, 8, 1, Hash_file);
	printf("ADDR %lX %lX %lX\n", Kmer_hash, Kmer_occr, Kmer_edit_depth);
	fflush(stdout);
	fwrite((void *) Kmer_hash, sizeof(uint64_t), Hash_size, Hash_file);
	puts("Writing chain list");
	fflush(stdout);
	fwrite((void *) Kmer_next_index, sizeof(uint32_t), Hash_size, Hash_file);
	fclose(Hash_file);
	free(Kmer_hash);
	free(Kmer_next_index);
	puts("Kmer search finished!");
	return 0;
}

int main_sparse_kmer(int argc, char ** argv){
	char getopt_return;
	extern char *optarg;
	uint16_t window_size = 1000;
	FILE * Control_bed = NULL;
	while ((getopt_return = getopt(argc, argv, "hw:c:")) != -1)
	{
		uint16_t len;
		switch (getopt_return)
		{
			case 'h':
				puts("quicKmer sparse bp ref.fa");
				puts("\tbp\tReduce number of kmer to every 1/bp");
				return 1;
			case 'w':
				window_size = atoi(optarg);
				break;
			case 'c':
				Control_bed = fopen(optarg, "r");
				printf("[Option] Control bedfile %s\n",optarg);
				break;
			case '?':
				puts("Option error, check help");
				return 1;
			default:
				return 1;
		}
	}
	char path[65536];
	strcpy(path, argv[argc-1]);
	strcat(path, ".qm");
	FILE * Hash_file = fopen(path, "r");
	FILE * fasta = fopen(argv[argc-1], "r");
	uint32_t thin_frac = atoi(argv[argc-2]);
	fseek(Hash_file, 4, 0);
	fread(&Kmer_size, 1, 1, Hash_file);
	fread(&edit_distance, 1, 1, Hash_file);
	fread(&Edit_depth_thres, 1, 1, Hash_file);
	Set_Kmer_Size(Kmer_size);
	fseek(Hash_file, 8, 0);
	fread(&Hash_size, 8, 1, Hash_file);
	uint64_t first_idx;
	fread(&first_idx, 8, 1, Hash_file);
	printf("Hash Size: 0x%lX\nFirst location: 0x%X\n", Hash_size, first_idx);
	Kmer_hash = (uint64_t *) malloc(Hash_size * sizeof(uint64_t));
	if (!Kmer_hash) {
		puts("Memory allocation failed");
		fclose(Hash_file);
		return 1;
	}
	printf("Read 0x%lX hash\n",fread(Kmer_hash, sizeof(uint64_t), Hash_size, Hash_file));
	Kmer_next_index = (uint32_t *) malloc(Hash_size * sizeof(uint32_t));
	if (!Kmer_next_index) {
		puts("Memory allocation failed");
		free(Kmer_hash);
		fclose(Hash_file);
		return 1;
	}
	printf("Read 0x%lX chain file\n",fread(Kmer_next_index, sizeof(uint32_t), Hash_size, Hash_file));
	fclose(Hash_file);
	//Scan fasta
	char fasta_buffer[200];
	uint8_t charge_size = 0;
	uint64_t encoded = 0;
	uint64_t encoded_r = 0;
	uint64_t count = 0;
	uint64_t last_keep_index;
	uint32_t bp_count, last_kmer_bp;
	//Loop through fasta lines
	if (thin_frac > 1)
	{
		while (fgets (fasta_buffer, 200, fasta) && fasta_buffer[0])
		{
			char * char_idx = fasta_buffer;
			if (*char_idx == '>')
			{
				charge_size = 0;
				encoded = 0;
				encoded_r = 0;
				bp_count = 0;
				last_kmer_bp = 0;
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
				if (charge_size < Kmer_size) charge_size++;
				if (kmer && charge_size == Kmer_size)
				{	
					uint64_t hash_index;
					if (Find_hash(kmer, &hash_index, Kmer_hash))
					{
						if (bp_count - last_kmer_bp < thin_frac)
							Kmer_next_index[hash_index] = 0;
						else {
							last_kmer_bp = bp_count;
							if (count) Kmer_next_index[last_keep_index] = hash_index;
							else {
								first_idx = hash_index;
								Kmer_next_index[first_idx] = hash_index;
							}
							last_keep_index = hash_index;
							count++;
						}
					}
				}
				bp_count++;
			}
		}
		Kmer_next_index[last_keep_index] = first_idx;
		printf("%i kmers left\n", count);
		//Recalculate optimal hash size
		uint64_t Optimal_size = (uint64_t) 1 << (uint8_t) ceil(log2(count / 0.8));
		for (count = 0; count < Hash_size; count++)
			if (!Kmer_next_index[count]) {
				Kmer_hash[count] = 0;
				Kmer_next_index[count] = 0;
			}
		
		if (Optimal_size != Hash_size) Resize_hash_table(Hash_size, Optimal_size);
	}
	strcpy(path, argv[argc-1]);
	strcat(path, ".bed");
	FILE * window_file = fopen(path, "w");
	strcpy(path, argv[argc-1]);
	FILE * GC_file = NULL;
	if (Control_bed != NULL) {
		strcat(path, ".qgc");
		GC_file = fopen(path, "w");
	}
	fseek(fasta, 0, SEEK_CUR);
	first_idx = dump_kmer_list(NULL, fasta, window_file, Control_bed, GC_file, 400, window_size);

	//Dump new hash table
	strcpy(path, argv[argc-1]);
	strcat(path, ".rqm");
	Hash_file = fopen(path, "w");
	const char Version[5] = "QM11";
	fwrite(Version, 1, 4, Hash_file);
	fwrite(&Kmer_size, 1, 1, Hash_file);
	fwrite(&edit_distance, 1, 1, Hash_file);
	fwrite(&Edit_depth_thres, 1, 1, Hash_file);
	fwrite(&thin_frac, 1, 1, Hash_file);
	fwrite(&Hash_size, 8, 1, Hash_file);
	fwrite(&first_idx, 8, 1, Hash_file);
	fwrite((void *) Kmer_hash, sizeof(uint64_t), Hash_size, Hash_file);
	puts("Writing chain list");
	fwrite((void *) Kmer_next_index, sizeof(uint32_t), Hash_size, Hash_file);
	fclose(Hash_file);
	free(Kmer_hash);
	free(Kmer_next_index);
	puts("Sparse kmer finished");
	return 0;
}

void printversion() {
	puts("QuicK-mer 2.0");
	puts("Operation modes: \n\tindex\tIndex a bed format kmer list");
	puts("\tcount\tCNV estimate from library\n\tsearch\tSearch K-kmer in genome");
	puts("\test\tGC normalization into copy number");
	puts("\tsparse\tFractionate kmer in index for memory reduction or regenerate GC control/Window\n");
	puts("Simple operation:\n1. Construct a dictionary from fasta using \"search\"");
	puts("2. Count depth from sample fasta/fastq \"count\"");
	puts("3. Estimate copy number with \"est\"\n");
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
		else if (strcmp(argv[1], "est") == 0)
			return main_estimate(argc-1, argv+1);
		else if (strcmp(argv[1], "sparse") == 0)
			return main_sparse_kmer(argc-1, argv+1);
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
