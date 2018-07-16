
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


int main(int argc, char** argv)
{
	std::ifstream kmer;
	kmer.open(argv[1],std::ifstream::in);
	std::ofstream output;
	output.open(argv[2],std::ofstream::out);
	int skip = atoi(argv[3]);
	unsigned int a,b, next_kmer_start;
	std::string chrom, mer,c, chrom_region, random_text;
	std::string pre_chrom = "";
	next_kmer_start = 0;
	while (kmer >> chrom >> a >> b >> c >> mer)
	{
		if (chrom != pre_chrom)
		{
			pre_chrom = chrom;
			next_kmer_start = 0;
		}
		if (pre_chrom == "chrUn" || pre_chrom == "chrNovel")
			output << chrom << '\t' << a << '\t' << b << '\t' << c << '\t' << mer << std::endl;
		else {
			if (a >= next_kmer_start)
			{
				next_kmer_start = a + skip;
				output << chrom << '\t' << a << '\t' << b << '\t' << c << '\t' << mer << std::endl;
			}
		}
	}
	return 0;
}
