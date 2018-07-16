/*
 * GC correction for mrsFAST
 */

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <string.h>
#include <map>
#include <iterator>
#include <cmath>
#include <iomanip>
#include <stdint.h>

#define ReadLen 36
#define buffer_size 20 * 1024 * 1024

typedef uint16_t npy_uint16;
typedef uint32_t npy_uint32;

npy_uint16 floatbits_to_halfbits(npy_uint32 f);

unsigned short float2half(float value)
{
	return floatbits_to_halfbits(*(npy_uint32 *)(&value));
}

npy_uint16 floatbits_to_halfbits(npy_uint32 f)
{
    npy_uint32 f_exp, f_man;
    npy_uint16 h_sgn, h_exp, h_man;

    h_sgn = (npy_uint16) ((f&0x80000000u) >> 16);
    f_exp = (f&0x7f800000u);
    
    /* Exponent overflow/NaN converts to signed inf/NaN */
    if (f_exp >= 0x47800000u) {
        if (f_exp == 0x7f800000u) {
            /*
             * No need to generate FP_INVALID or FP_OVERFLOW here, as
             * the float/double routine should have done that.
             */
            f_man = (f&0x007fffffu);
            if (f_man != 0) {
                /* NaN - propagate the flag in the mantissa... */
                npy_uint16 ret = (npy_uint16) (0x7c00u + (f_man >> 13));
                /* ...but make sure it stays a NaN */
                if (ret == 0x7c00u) {
                    ret++;
                }
                return h_sgn + ret;
            } else {
                /* signed inf */
                return (npy_uint16) (h_sgn + 0x7c00u);
            }
        } else {
            /* overflow to signed inf */
#if HALF_GENERATE_OVERFLOW
            generate_overflow_error();
#endif
            return (npy_uint16) (h_sgn + 0x7c00u);
        }
    }
    
    /* Exponent underflow converts to denormalized half or signed zero */
    if (f_exp <= 0x38000000u) {
        /* 
         * Signed zeros, denormalized floats, and floats with small
         * exponents all convert to signed zero halfs.
         */
        if (f_exp < 0x33000000u) {
#if HALF_GENERATE_UNDERFLOW 
            /* If f != 0, we underflowed to 0 */
            if ((f&0x7fffffff) != 0) {
                generate_underflow_error();
            }
#endif
            return h_sgn;
        }
        /* It underflowed to a denormalized value */
#if HALF_GENERATE_UNDERFLOW 
        generate_underflow_error();
#endif
        /* Make the denormalized mantissa */
        f_exp >>= 23;
        f_man = (0x00800000u + (f&0x007fffffu)) >> (113 - f_exp);
        /* Handle rounding by adding 1 to the bit beyond half precision */
#if HALF_ROUND_TIES_TO_EVEN 
        /*
         * If the last bit in the half mantissa is 0 (already even), and
         * the remaining bit pattern is 1000...0, then we do not add one
         * to the bit after the half mantissa.  In all other cases, we do.
         */
        if ((f_man&0x00003fffu) != 0x00001000u) {
            f_man += 0x00001000u;
        }
#else
        f_man += 0x00001000u;
#endif
        h_man = (npy_uint16) (f_man >> 13);
        /*
         * If the rounding causes a bit to spill into h_exp, it will
         * increment h_exp from zero to one and h_man will be zero.
         * This is the correct result.
         */
        return (npy_uint16) (h_sgn + h_man);
    }

    /* Regular case with no overflow or underflow */
    h_exp = (npy_uint16) ((f_exp - 0x38000000u) >> 13);
    /* Handle rounding by adding 1 to the bit beyond half precision */
    f_man = (f&0x007fffffu);
#if HALF_ROUND_TIES_TO_EVEN 
    /*
     * If the last bit in the half mantissa is 0 (already even), and
     * the remaining bit pattern is 1000...0, then we do not add one
     * to the bit after the half mantissa.  In all other cases, we do.
     */
    if ((f_man&0x00003fffu) != 0x00001000u) {
        f_man += 0x00001000u;
    }
#else
    f_man += 0x00001000u;
#endif
    h_man = (npy_uint16) (f_man >> 13);
    /*
     * If the rounding causes a bit to spill into h_exp, it will
     * increment h_exp by one and h_man will be zero.  This is the
     * correct result.  h_exp may increment to 15, at greatest, in
     * which case the result overflows to a signed inf.
     */
#if HALF_GENERATE_OVERFLOW
    h_man += h_exp;
    if (h_man == 0x7c00u) {
        generate_overflow_error();
    }
    return h_sgn + h_man;
#else
    return h_sgn + h_exp + h_man;
#endif
}

int main(int argc, char** argv)
{
	/*
	 * reference fasta index
	 * GC_control_region
	 * input sam stream
	 * output file prefix
	 *  mrsfast-3.3.8
	 */

	
	unsigned int cumulate;
	std::ifstream GC_region;
	std::ifstream jfquery;
	GC_region.open(argv[1], std::ifstream::in | std::ifstream::binary);
	GC_region.seekg(0, std::ios::end);
	cumulate = GC_region.tellg() >> 2;
	GC_region.seekg(0);
	jfquery.open(argv[2], std::ifstream::in);
	std::string line;
	unsigned short * Depth = new unsigned short [cumulate];
	unsigned int site = 0;
	while (!jfquery.eof())
	{
		std::getline(jfquery, line);
		Depth[site] = atoi(line.c_str());
		site++;
	}
	jfquery.close();

	//GC control curve
	std::vector<unsigned int> Bin_count (401,0);
	std::vector<float> Average_depth (401,0.0);
	float * GC_buffer = new float[buffer_size];
	unsigned int bp_pos = 0;
	while (bp_pos < cumulate)
	{
		unsigned int readsize;
		if ((cumulate - bp_pos) >= buffer_size) readsize = buffer_size;
		else readsize = cumulate - bp_pos;
		GC_region.read((char *) GC_buffer, sizeof(float)* buffer_size);
		for (unsigned int i = 0; i < readsize; i++){
			if (GC_buffer[i] > 0)
			{
				//Addition of depth to control region bins
				unsigned short bin_idx = (unsigned short) (GC_buffer[i]*400+0.5);
				Average_depth[bin_idx] += Depth[bp_pos];
				Bin_count[bin_idx]++;
			}
			bp_pos++;
		}
	}
	puts("Finish GC control scan\n");

	//Average Depth
	float average_depth = 0;
	unsigned int total_ctrl_window = 0;
	std::ofstream GC_curve_file;
	std::string argv4 = argv[3];
	std::string filename;
	filename = argv4 + ".txt";
	GC_curve_file.open(filename.c_str(), std::ofstream::out);
	for (int i = 0; i <= 400; i++)
	{
		average_depth += Average_depth[i];
		if (Bin_count[i]) Average_depth[i] /= Bin_count[i];
		total_ctrl_window += Bin_count[i];
		GC_curve_file << i/4.0 << '\t' << Average_depth[i] << '\t' << Bin_count[i] << "\t0" << std::endl;
	}
	average_depth /= total_ctrl_window;
	GC_curve_file.close();
	puts("Finish GC curve pile up\n");
	//Lowess
	float * lowess_factor = new float[401];

	FILE *lowess_pipe;

	filename = "smooth_GC_mrsfast.py " + argv4 + ".txt";
	lowess_pipe = popen(filename.c_str(), "r");

	if (!lowess_pipe) return 1;
	fread((char *) lowess_factor, sizeof(float), 401, lowess_pipe);
	pclose(lowess_pipe);

	//Apply correction
	unsigned short * out_buffer = new unsigned short[buffer_size];
	std::ofstream Output_file;
	filename = argv4 + ".bin";
	Output_file.open(filename.c_str(), std::ofstream::out | std::ofstream::binary);
	//Reset seek pointer
	GC_region.clear();
	GC_region.seekg(0);
	bp_pos = 0;

	while (bp_pos < cumulate)
	{
		unsigned int readsize;
		if ((cumulate - bp_pos) >= buffer_size) readsize = buffer_size;
		else readsize = cumulate - bp_pos;
		GC_region.read((char *) GC_buffer, sizeof(float)* readsize);
		for (unsigned int i = 0; i < readsize; i++)
		{
			unsigned short bin_idx = (unsigned short) (abs(GC_buffer[i]*400) + 0.5);
			float corr_depth = Depth[bp_pos] * lowess_factor[bin_idx];
			if (GC_buffer[i] >= 0) out_buffer[i] = float2half(corr_depth);
			else out_buffer[i] = float2half(-corr_depth);
			bp_pos++;
		}
		Output_file.write((char *)out_buffer, sizeof(uint16_t)*readsize);
	}
	Output_file.close();
	GC_region.close();
	delete Depth;
	delete lowess_factor;
	delete GC_buffer;
	delete out_buffer;
	return 0;
}
