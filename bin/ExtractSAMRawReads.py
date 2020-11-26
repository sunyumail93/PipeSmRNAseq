#!/usr/bin/python
#ExtractSAMRawReads.py
#This script is part of PipeSmRNAseq pipeline to extract raw reads from SAM
#For the reads mapped to the reverse strand (i.e. column 2 == 16)
#Output: Data.Lendis
#Yu Sun, 2020-11-25

import sys

def Reverse(String):
    StringList=list(String)
    for i in range(len(StringList)):
        if StringList[i]=="A": StringList[i]="T"
        elif StringList[i]=="a": StringList[i]="t"
        elif StringList[i]=="T": StringList[i]="A"
        elif StringList[i]=="t": StringList[i]="a"
        elif StringList[i]=="G": StringList[i]="C"
        elif StringList[i]=="g": StringList[i]="c"
        elif StringList[i]=="C": StringList[i]="G"
        elif StringList[i]=="c": StringList[i]="g"
        else: StringList[i]="N"
    return "".join(map(str,StringList))

def ExtractSAMRawReads():
    fi=open(sys.argv[1],'r')
    fo=open(sys.argv[1]+".Reads",'w')
    
    for line in fi:
        if not line.startswith("@"):
            currline=line.strip().split()
            if currline[1] == "0":
                Seq=currline[9]
                fo.write(Seq+"\n")
            elif currline[1] == "16":
                Seq=Reverse(currline[9])
                Seq=Seq[::-1]
                fo.write(Seq+"\n")
            else:
                pass

    fi.close()
    fo.close()

if len(sys.argv) != 2:
    print("This script is part of PipeSmRNAseq pipeline to extract raw reads from SAM")
    print("For the reads mapped to the reverse strand (i.e. column 2 == 16)")
    print("Usage: [ExtractSAMRawReads.py] [SAM]")
else:
    ExtractSAMRawReads()
