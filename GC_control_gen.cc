/*
 * mrsFAST pipeline generating the control region
 * GC content value, control flag as minus sign
 * Masked control region
 */

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

//#define DEBUG
#define Use_strict 0

struct exclude
{
	unsigned int id;
	unsigned int start;
	unsigned int end;
};

struct FAIDX
{
	unsigned int length;
	unsigned int offset;
	unsigned short bp_per_line;
	unsigned short byte_per_line;
};
std::ifstream fasta, kmer_list;
char * seq_buffer;
unsigned int pointer = 0, pointer_t = 0;
unsigned short window_size;
unsigned short window_size_mem;
unsigned int readlength;
float out_buffer [1024*1024];
unsigned int buf_pointer = 0;
unsigned int kmer_pos;

unsigned short getbp(unsigned char bp_per_line)
{
	unsigned short base;
	if (pointer < readlength){
		base = seq_buffer[pointer] + (seq_buffer[pointer_t] << 8);
		pointer++;
		if (pointer % (bp_per_line + 1) == bp_per_line) pointer++;
	}
	else base = (seq_buffer[pointer_t] << 8);
	if (pointer > window_size_mem) pointer_t++;
	else base &= 0x00FF;
	if (pointer_t % (bp_per_line + 1) == bp_per_line) pointer_t++;
	return base;
}

void nextchrom(unsigned int offset, unsigned int length, unsigned char bp_per_line)
{
	readlength = length % bp_per_line + (length / bp_per_line) * (bp_per_line + 1);
	seq_buffer = new char [readlength];
	fasta.seekg(offset);
	fasta.read(seq_buffer, readlength);
	std::cout << "Read bytes: " << readlength << std::endl;
	window_size_mem = (window_size/bp_per_line) *(bp_per_line +1) + window_size % bp_per_line;
	pointer = 0;
	pointer_t = 0;
}



int main(int argc, char** argv)
{
	if (argc != 6)
	{
		std::cout << "Require 5 argumetns.\n Ref.fa Excluded.bed Kmer_list 400 Ref_GC_control.bin" << std::endl;
		std::cout << "Excluded.bed must be in same order as Ref.fa" << std::endl;
		return 0;
	}
	std::ifstream fa_idx;
	std::ifstream control;
	control.open(argv[2],std::ifstream::in);
	fasta.open(argv[1],std::ifstream::in|std::ifstream::binary);
	fa_idx.open(strcat(argv[1], ".fai"), std::ifstream::in);
	kmer_list.open(argv[3], std::ifstream::in);
	std::ofstream output, GCoutput;
	window_size = atoi(argv[4]);
	output.open(argv[5],std::ofstream::out|std::ofstream::binary);
	//GCoutput.open(strcat(argv[5],".txt"),std::ofstream::out);
	unsigned int a,b,c,d, kmer_start, kmer_end;
	std::string chrom, randomstring, kmer;
	std::vector<FAIDX> faidx;
	std::vector<std::string> chrom_name;
	std::map<std::string, unsigned int> Chrom_map;
	while (fa_idx >> chrom >> a >> b >> c >> d)
	{
		faidx.push_back({a, b, c, d});
		Chrom_map[chrom] = faidx.size() - 1;
		chrom_name.push_back(chrom);
	}
	std::string pre_chrom = "";
	std::vector<exclude> exclude_regions;
	int chr_id = -1;
	unsigned int LH_bp = window_size >> 1;
	unsigned int RH_bp = window_size - LH_bp;
	while (control >> chrom >> a >> b)
	{
		if (chrom != pre_chrom)
		{
			//Last chromosome window extend by half window
			//if (exclude_regions.size() && Use_strict) exclude_regions[exclude_regions.size() - 1].end += RH_bp;
			//New chromosome
			chr_id = Chrom_map[chrom];
			pre_chrom = chrom;
			exclude_regions.push_back({chr_id, a, b});
		}
		else {
			if (Use_strict)
			{
				//Use strict sense of control region
				//Boundary overlapping is prohibited
				if (a - exclude_regions[exclude_regions.size() - 1].end <= window_size)
					exclude_regions[exclude_regions.size() - 1].end = b;
				else {
					exclude_regions[exclude_regions.size() - 1].end += RH_bp;
					unsigned int next_start = a;
					if (a > LH_bp) next_start = a - LH_bp;
					exclude_regions.push_back({chr_id, next_start, b});
				}
			}
			else exclude_regions.push_back({chr_id, a, b});
		}
	}
//printf("%i\n", exclude_regions.size());
	unsigned int i = 0;
	unsigned int exc_this_chrom = 0;
	while (1) {
		exc_this_chrom += exclude_regions[i].end - exclude_regions[i].start;
		i++;
		if (i == exclude_regions.size()) break;
	}
	std::cout << "Total excluded bp: " << exc_this_chrom << std::endl;

#ifdef DEBUG
	for (int i = 0; i < exclude_regions.size(); i++){
		std::cout << "\t" << exclude_regions[i].start << "\t" << exclude_regions[i].end << std::endl;
	}
#endif

	control.close();

	kmer_list >> chrom >> kmer_start >> kmer_end >> randomstring >> kmer;
	kmer_pos = (kmer_start + kmer_end) >> 1;
	unsigned int total_ctrl_bp = 0;
	unsigned int total_bp = 0;
	unsigned int hist[401] = {0};
	unsigned int window_pointer = 0;
	unsigned int kmer_counter = 0;
	for (unsigned int faidx_i = 0; faidx_i < faidx.size(); faidx_i++)
	{
		unsigned int chr_len = faidx[faidx_i].length;
		unsigned int chr_offset = faidx[faidx_i].offset;
		unsigned short bp_per_line = faidx[faidx_i].bp_per_line;
		std::cout << faidx_i << '\t' << chr_len << '\t' << chr_offset <<std::endl;
		nextchrom(chr_offset, chr_len, bp_per_line);

		total_bp += chr_len;
		unsigned int AT = 0, GC = 0;
		//Precharge half buffer
		for (unsigned int chrom_idx = 0; chrom_idx < RH_bp; chrom_idx++)
		{
			unsigned short base = getbp(bp_per_line);
			char basenew = base & 0xFF;
			if (basenew == 'a' || basenew == 'A' || basenew =='t' || basenew == 'T') AT++;
			if (basenew == 'c' || basenew == 'C' || basenew =='g' || basenew == 'G') GC++;
		}
		//Actual GC sliding start
		for (unsigned int chrom_idx = 0; chrom_idx < chr_len; chrom_idx++)
		{
			unsigned short base = getbp(c);
			char basenew = base & 0xFF;
			char basediscard = base >> 8;
			if (chrom_idx < chr_len - RH_bp){
				if (basenew == 'a' || basenew == 'A' || basenew =='t' || basenew == 'T') AT++;
				if (basenew == 'c' || basenew == 'C' || basenew =='g' || basenew == 'G') GC++;
			}
			if (basediscard == 'a' || basediscard == 'A' || basediscard =='t' || basediscard == 'T') AT--;
			if (basediscard == 'c' || basediscard == 'C' || basediscard =='g' || basediscard == 'G') GC--;
			unsigned short total_meaningful_bp = AT + GC;
			float GC_content;
			{
				if (total_meaningful_bp != 0) GC_content = GC / (float)(AT + GC);
				else GC_content = 0.0;
			}
			//if (chrom_idx % 50 == 0) GCoutput << chrom_name[faidx_i] << '\t' << chrom_idx << '\t' << chrom_idx+1 << '\t' << GC_content << std::endl;
			/*Check control region*/
			char control_flag;
			if (faidx_i > exclude_regions[window_pointer].id || (faidx_i == exclude_regions[window_pointer].id && chrom_idx >= exclude_regions[window_pointer].end)) {
				if (window_pointer < exclude_regions.size()) window_pointer++;
			}
			if (chrom_idx >= exclude_regions[window_pointer].start && faidx_i == exclude_regions[window_pointer].id && window_pointer < exclude_regions.size()) {
				//within excluded region, GC is negative
				control_flag = 1;
			} else control_flag = 0;
			
			//Write to Buffer
			//Check position relative to kmer
			if (kmer_pos == chrom_idx) {
				if (control_flag) {
					total_ctrl_bp++;
					hist[(unsigned short) (GC_content * 400 + 0.5)]++;
				}
				else GC_content *= -1.0;
				out_buffer[buf_pointer] = GC_content;
				buf_pointer++;
				if (buf_pointer == 1024*1024)
				{
					buf_pointer = 0;
					output.write((char *)out_buffer, sizeof(out_buffer));
				}
				if (kmer_list >> chrom >> kmer_start >> kmer_end >> randomstring >> kmer)
				kmer_pos = (kmer_start + kmer_end) >> 1;
			}
		}
		delete seq_buffer;
	}

	output.write((char *)out_buffer, 4*buf_pointer);
	output.close();
	std::cout << total_ctrl_bp << " out of " << total_bp << " bp are control regions." << std::endl;
	for (int i = 0; i <= 400; i++) std::cout << i / 4.0 << '\t' << hist[i] << std::endl;
	fasta.close();
	fa_idx.close();
	return 0;
}
