#!/usr/bin/env python3
# Jeff Kidd
# read fastq from stdin, write out fasta to stdout folded
# to make input for quicKmer2 from long pacbio reads

import sys


n = 50  # 50 bp per line
while True:    
    r = sys.stdin.readline()
    if r == '': # end of file
        break
    if r[0] != '@':
        sys.stdout.write('ERROR! record does not start with @!\n')
        sys.stdout.write(r[0:100])
        sys.exit()

    #print out name, read name does not matter
    sys.stdout.write('>read\n')    
    seq = sys.stdin.readline().rstrip() 
    seqLen = len(seq)
    for i in range(0,seqLen,n):
        sys.stdout.write(seq[i:i+n])
        sys.stdout.write('\n')

    
    
    r = sys.stdin.readline() # the + line
    r = sys.stdin.readline() # the qual line
            