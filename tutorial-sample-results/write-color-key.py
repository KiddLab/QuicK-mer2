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


outFile = open('color-track.bed','w')

for i in range(0,11):
    b = '0'
    e = '1000'
    c = 'chr1'
    trackName = str(i)
    if i == 10:
        trackName = '10+'
    col = CNtoColor[i] 
    nl = [c,b,e,trackName,'0','.',b,e,col]
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)


outFile.close()