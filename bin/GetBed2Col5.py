#!/usr/bin/python
#This script is part of PipeSmRNAseq.sh
#This script take one column value, and repeat the non-one value that value times
#for example: [1 2 3 4 2] -> [1 2 2 3 3 3 4 4 4 4 2 2]
#Version: 2020-12-03, Yu Sun
import sys

def GetBed2Col5():
    fi=open(sys.argv[1],'r')
    fo=open(sys.argv[2],'w')
    counting=fi.readlines()
    
    for i in counting:
        i=int(i.strip())
        if i==1:
            fo.write(str(i)+'\n')
        else:
            for j in range(0,i):
                fo.write(str(i)+'\n')

    fi.close()
    fo.close()

if len(sys.argv) != 3:
    print("This script is part of PipeSmRNAseq.sh")
    print("This script take one column value, and repeat the non-one value that value times")
    print("Usage: [GetBed2Col5.py] [input] [output]")
else:
    GetBed2Col5()
