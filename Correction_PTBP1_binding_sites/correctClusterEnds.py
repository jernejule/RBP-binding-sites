'''
Created on Aug 20, 2013

@author: Nejc

The script will correct cluster ends after bedtools closest intersect.
'''


import sys

def setClusterEnds(fin_fname, fout_fname):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    fout_old = open(fout_fname + "-NOT.bed", "w")
    line = fin.readline()

    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        clusterStart = col[1]
        clusterEnd = col[2]
        NewClusterStart  = col[7]
        NewClusterEnd  = col[8]
        distance = int(col[12])
        strand = col[5]

        if distance <= 40 and distance > 0: #maximum read length is 40 amd -1 means there were none in that direction
            if strand == '+':   #we write corrected positions and original one in 3th column
                fout.write(chr + '\t' + str(clusterStart) + '\t' + str(NewClusterEnd) + '\t' + chr+':'+clusterStart+':'+clusterEnd + '\t' + '' + '\t' + strand + '\n')
            else:
                fout.write(chr + '\t' + str(NewClusterStart) + '\t' + str(clusterEnd) + '\t' + chr+':'+clusterStart+':'+clusterEnd + '\t' + '' + '\t' + strand + '\n')
        else: #we didn't extand the cluster ends so we write original cluster to a new fil
            fout_old.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[4] + '\t' + col[5] + '\n')
            
        line = fin.readline()
    fout.close()
    fin.close()
    fout_old.close()
    


if sys.argv.__len__() == 3:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    setClusterEnds(fin_fname, fout_fname)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python correctClusterEnds.py input_fname output_fname"

