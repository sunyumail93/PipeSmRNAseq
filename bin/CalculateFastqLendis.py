#!/usr/bin/python
#CalculateFastqLendis.py
#This script is part of PipeSmRNAseq pipeline to calculate unmapped fastq lendis
#Output: Data.Lendis
#Yu Sun, 2020-11-25

import sys

def CalcLendis():
    fi=open(sys.argv[1],'r')
    fo=open(sys.argv[1]+".Lendis",'w')

    Lendis={}
    for i in range(18,36):
        Lendis[i]=0

    while True:
        FirstLine=fi.readline().strip()
        SecondLine=fi.readline().strip()
        ThirdLine=fi.readline().strip()
        ForthLine=fi.readline().strip()
        if not FirstLine:
            break
        else:
            Length=len(SecondLine)
            if (Length > 17) and (Length < 36):
                Lendis[Length]+=1

    for i in range(18,36):
        fo.write(str(i)+"\t"+str(Lendis[i])+"\n")
        
    fi.close()
    fo.close()

if len(sys.argv) != 2:
    print("This script is part of PipeSmRNAseq pipeline to calculate unmapped fastq lendis")
    print("Output: Data.Lendis")
    print("Usage: [CalculateFastqLendis.py] [Fastq]")
else:
    CalcLendis()
