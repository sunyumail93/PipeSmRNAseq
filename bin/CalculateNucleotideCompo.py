#!/usr/bin/python
#CalculateNucleotideCompo.py
#This script is part of PipeSmRNAseq pipeline to calculate nucleotide composition, especially for piRNA
#Output: Nt composition table
#Nt order, 0, 1, 2, 3 : A, C, G, T
#Yu Sun, 2020-11-25

import sys

def NucleotideCompo():
    fi=open(sys.argv[1],'r')
    fo=open(sys.argv[1]+".NtCompo",'w')

    Matrix=[[0]*23, [0]*23, [0]*23, [0]*23]
    for line in fi:
        line=line.strip()
        for nt in range(23):
            if line[nt] == 'A' or line[nt] == 'a':
                Matrix[0][nt]+=1
            elif line[nt] == 'C' or line[nt] == 'c':
                Matrix[1][nt]+=1
            elif line[nt] == 'G' or line[nt] == 'g':
                Matrix[2][nt]+=1
            elif line[nt] == 'T' or line[nt] == 't':
                Matrix[3][nt]+=1
            else:
                pass
    fo.write("A"+"\t"+'\t'.join(map(str,Matrix[0]))+"\n")
    fo.write("C"+"\t"+'\t'.join(map(str,Matrix[1]))+"\n")
    fo.write("G"+"\t"+'\t'.join(map(str,Matrix[2]))+"\n")
    fo.write("T"+"\t"+'\t'.join(map(str,Matrix[3]))+"\n")

    fi.close()
    fo.close()

if len(sys.argv) != 2:
    print("This script is part of PipeSmRNAseq pipeline to calculate nucleotide composition, especially for piRNA")
    print("Usage: [ExtractSAMRawReads.py] [23mer]")
else:
    NucleotideCompo()
