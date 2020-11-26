#!/usr/bin/python
#SeparateFastq_siRNA_piRNA.py
#This script is part of PipeSmRNAseq pipeline to separate fastq by length
#21-23 nt reads will be siRNAs, and >23 nt will be piRNAs. The rest <21 nt will be the rest miRNAs
#Output: Prefix.Less21.fastq, Prefix.21-23.fastq, Prefix.Over23.fastq
#Yu Sun, 2020-11-22

import sys

def LengthSeparator():
    fi=open(sys.argv[1],'r')
    fosmall=open(sys.argv[2]+".Less21.fastq",'w')
    fosiRNA=open(sys.argv[2]+".21-23.fastq",'w')
    fopiRNA=open(sys.argv[2]+".Over23.fastq",'w')

    Counter=0
    while True:
        FirstLine=fi.readline().strip()
        SecondLine=fi.readline().strip()
        ThirdLine=fi.readline().strip()
        ForthLine=fi.readline().strip()
        if not FirstLine:
            break
        else:
            Counter+=1
            Length=len(SecondLine)
            if Length < 21:
                fosmall.write(FirstLine+"\n")
                fosmall.write(SecondLine+"\n")
                fosmall.write(ThirdLine+"\n")
                fosmall.write(ForthLine+"\n")
            elif Length > 23:
                fopiRNA.write(FirstLine+"\n")
                fopiRNA.write(SecondLine+"\n")
                fopiRNA.write(ThirdLine+"\n")
                fopiRNA.write(ForthLine+"\n")
            else:
                fosiRNA.write(FirstLine+"\n")
                fosiRNA.write(SecondLine+"\n")
                fosiRNA.write(ThirdLine+"\n")
                fosiRNA.write(ForthLine+"\n")
        
    fi.close()
    fosmall.close()
    fosiRNA.close()
    fopiRNA.close()
    print(Counter)

if len(sys.argv) != 3:
    print("This script is part of PipeSmRNAseq pipeline to separate fastq by length")
    print("21-23 nt reads will be siRNAs, and >23 nt will be piRNAs. The rest <21 nt will be the rest miRNAs")
    print("Usage: [SeparateFastq_siRNA_piRNA.py] [Fastq] [Prefix]")
else:
    LengthSeparator()
