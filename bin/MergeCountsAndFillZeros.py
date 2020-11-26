#!/usr/bin/python
#MergeCountsAndFillZeros.py
#This script is part of PipeSmRNAseq.sh to fill zeros to undetected miRNAs, and sort the output miRNA list based on a desired order
#Yu Sun, 2020-11-17
                                                                          
import sys

def MergeCountsAndFillZeros():
    fi=open(sys.argv[1],'r')
    forder=open(sys.argv[2],'r')
    fo=open(sys.argv[3],'w')
    
    Counts={}
    miRNAList=[]
    for miRNA in forder:
        miRNA=miRNA.strip()
        miRNAList.append(miRNA)
        Counts[miRNA]=0

    for line in fi:
        Name=line.strip().split()[0]
        count=line.strip().split()[1]
        Counts[Name]=count

    for i in miRNAList:
        fo.write(i+"\t"+str(Counts[i])+"\n")

    fi.close()
    forder.close()
    fo.close()

if len(sys.argv) != 4:
    print("This script is part of PipeSmRNAseq.sh to fill zeros to undetected miRNAs, and sort the output miRNA list based on a desired order")
    print("Usage: [MergeCountsAndFillZeros.py] [CountsTable|Two columns: Name, Counts] [OrderList] [OutputFile]")
else:
    MergeCountsAndFillZeros()
