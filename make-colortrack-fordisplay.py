#!/usr/bin/env python3
import sys
import os
import subprocess

from optparse import  OptionParser


###############################################################################
USAGE = """
make-colortrack-fordisplay.py  --cn <depth file>  --name <name for track/sample>

takes file of copy number in windows and makes colored bed file, as well
as merged colored bed file suitable for loading in UCSC browser

"""

parser = OptionParser(USAGE)

parser.add_option('--cn',dest='copyNumberFile',help='copy number file')
parser.add_option('--name',dest='sampleName',help='name of sample for track')


(options,args)=parser.parse_args()

if options.copyNumberFile is None:
    parser.error('copyNumberFile not given')
if options.sampleName is None:
    parser.error('sampleName not given')

###############################################################################



###############################################################################
CNtoColor = {}
CNtoColor[0]=   '224,224,224'
CNtoColor[1]=   '160,160,160'
CNtoColor[2]=   '0,0,0'
CNtoColor[3]=   '0,0,153'
CNtoColor[4]=   '51,51,255'
CNtoColor[5]=   '0,255,255'
CNtoColor[6]=   '0,153,0'
CNtoColor[7]=   '255,255,0'
CNtoColor[8]=   '255,153,51'
CNtoColor[9]=   '153,76,0'
CNtoColor[10]=  '204,0,0'
###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print('command failed')
        print(cmd)
        sys.exit(1)
###############################################################################
def make_bed_color(trackName,CNfile,outFileName):
    inFile = open(CNfile,'r')
    outFile = open(outFileName,'w')
    for line in inFile:
        line = line.strip()
        line = line.split()
        c = line[0]
        b = line[1]
        e = line[2]
        cn = int(round(float(line[3]) ))
        if cn > 10:
            cn = 10
        
        # with the many short contigs, we sometimes get no unmasked positions
        # this gives a negative depth, which we will code as 0 for color tracks
        
        if cn < 0:
            cn = 0

            
        col = CNtoColor[cn]
        nl = [c,b,e,trackName,'0','.',b,e,col]
        nl = '\t'.join(nl) + '\n'
        outFile.write(nl)
    inFile.close()
    outFile.close()
###############################################################################
###############################################################################
def check_merge(l1,l2):
    if l1[0] != l2[0]:
        return False
    if l1[8] != l2[8]:
        return False    
    if l2[1] == l1[2]: # adjacent:
        return True
    return False
###############################################################################
def merge_bed9(inFileName,outFileName):
	prevLine = ''
	inFile = open(inFileName,'r')
	outFile = open(outFileName,'w')
	for line in inFile:
		line = line.rstrip()
		line = line.split()
		line[1] = int(line[1])
		line[2] = int(line[2])
	
	
		# if there is no prev line
		if prevLine == '':
			prevLine = line
			continue
	
		# compare with prev line    
		to_merge = check_merge(prevLine,line)
		if to_merge is True:
			prevLine[2] = line[2]
			prevLine[7] = line[7]
		else:
			# print out prev line and rest
			nl = [str(j) for j in prevLine]
			nl = '\t'.join(nl) + '\n'
			outFile.write(nl)
			prevLine = line
	inFile.close()

	# print out prevLine at end
	nl = [str(j) for j in prevLine]
	nl = '\t'.join(nl) + '\n'
	outFile.write(nl)
	inFile.close()
	outFile.close()
###############################################################################



outBedDisplay = options.copyNumberFile
sampleName = options.sampleName
print('initial copy number file is',outBedDisplay)
tmpCol = outBedDisplay + '.tmp'
outBedColor = outBedDisplay + '.bedColor'


make_bed_color(sampleName,outBedDisplay,tmpCol)
# do merge bed
merge_bed9(tmpCol,outBedColor)
# rm tmp
cmd = 'rm %s ' % tmpCol
print(cmd)
runCMD(cmd)

