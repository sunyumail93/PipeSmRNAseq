#!/usr/bin/env bash
#DefineMatureMiRNACoordinates.sh
#This pipeline maps mature miRNAs to precursors
#It generats a bed6 index for each mature miRNA, mapping to the precursors.
#The result can server as an reference of mature miRNA coordinates on precursors. This pipeline will only report 1 location for multi-mappings.
#After generating the coordinates is will define the 5' or 3' arms de novo, and add the info to column 5
#Version: Yu Sun, 2020-11-15

if [ ! -n "$3" ]
then
  echo "    This pipeline maps mature miRNAs to precursors"
  echo "    It generats a bed6 index for each mature miRNA, mapping to the precursors."
  echo "    The result can server as an reference of mature miRNA annotation. This pipeline will only report 1 location for multi-mappings."
  echo "    After generating the matching is will define the 5' or 3' arms de novo (mapped to left half or right half of precursors), and add the info to column 5"
  echo "    Usage: `basename $0` [hairpin.fa] [mature.fa] [output.bed6]"
  echo "    Output: BED6 output coordinates of mature miRNAs"
else
  calculating(){ awk "BEGIN { print "$*" }"; }
  Hairpin=$1
  Mature=$2
  OutputFinal=$3
  Output=$OutputFinal.t
  OutputFinalFull=$OutputFinal.full
  rm -rf $Output
  rm -rf $OutputFinalFull $OutputFinal.t
  module load bowtie
  module load samtools
  echo "1. Indexing"
  bowtie-build --quiet $Hairpin hairpin
  awk '{print $1}' $Hairpin|grep ">" |sed 's/>//' > temp.hairpin.names
  awk '{print $1}' $Hairpin|grep -v ">" |awk '{print length($0)}' > temp.hairpin.len
  paste temp.hairpin.names temp.hairpin.len > hairpin.sizes && rm -rf temp.hairpin.names temp.hairpin.len

  echo "2. Mapping"
#  bowtie -S -f -v 0 -k 1 hairpin $Mature 1> /dev/stdout 2> /dev/null| samtools view -uS -|bedtools_piPipes bamtobed -i -|awk '$6=="+"' > $Output
  bowtie -S -f -v 0 -k 1 --no-unal hairpin $Mature > ${Mature}.mapped.sam
  samtools view -b ${Mature}.mapped.sam > ${Mature}.mapped.bam 
  bedtools bamtobed -i ${Mature}.mapped.bam | awk '$6=="+"' > $Output

  echo "3. Defining arms"
  #Preparing
  Start=1
  Lines=`wc -l $Output|awk '{print $1}'`
  for i in $(eval echo {$Start..$Lines})
  do
    whole=`head -${i} $Output | tail -1`
    Precursor=`head -${i} $Output | tail -1|awk '{print $1}'`
    PosStart=`head -${i} $Output | tail -1|awk '{print $2}'`
    PosEnd=`head -${i} $Output | tail -1|awk '{print $3}'`
    MatureName=`head -${i} $Output | tail -1|awk '{print $4}'`
    LastTwo=${MatureName: -2}
    if [[ $LastTwo == "3p" ]];then
      PreArm="3p"
    elif [[ $LastTwo == "5p" ]];then
      PreArm="5p"
    else
      PreArm="xx"
    fi
    echo $whole
    echo $PreArm
    HalfLengthp=`awk -v t=$Precursor '{if ($1==t) print $2/2}' hairpin.sizes`
    HalfLength=${HalfLengthp%.*}
    echo $HalfLength
    if [[ "$PosEnd" -lt "$HalfLength" ]];then
      DefineArm="5p"
    elif [[ "$HalfLength" -lt "$PosStart" ]];then
      DefineArm="3p"
    else
      #Rescue some miRNAs across the middle base
      DefineArm="none"
      leftover=`calculating $HalfLength-$PosStart`
      rightover=`calculating $PosEnd-$HalfLength`
      if [[ "$leftover" -lt "$rightover" ]];then
        DefineArm="3p"
      else
        DefineArm="5p"
      fi
    fi
    echo "DefineArm: "$DefineArm

    if [[ $PreArm == $DefineArm ]];then
      Type="consistent"
      echo "consistent"
      echo "Final arm: "$DefineArm
    elif [[ $PreArm != $DefineArm ]];then
      if [[ $PreArm == "xx" ]];then
        echo "Novel defined: "$DefineArm
        Type="novel"
      else
        #This final case is very rare, for mouse, only has one case (mmu-mir-466i), when the miRNA across the middle of the hairpin
        echo -e "Inconsistent: "$PreArm"\t"$DefineArm
        DefineArm=$PreArm  #Assign the annotation
        Type="inconsistent"
        echo "Final arm: "$DefineArm
      fi
    fi
    echo $whole|awk -v newarm=$DefineArm 'BEGIN{OFS="\t"}{$5=newarm;print $0}'
    echo $whole|awk -v newarm=$DefineArm 'BEGIN{OFS="\t"}{$5=newarm;print $0}' >> $OutputFinal
    echo $whole|awk -v newarm=$DefineArm -v type=$Type 'BEGIN{OFS="\t"}{$5=newarm;$7=type;print $0}'
    echo $whole|awk -v newarm=$DefineArm -v type=$Type 'BEGIN{OFS="\t"}{$5=newarm;$7=type;print $0}' >> $OutputFinalFull

  done

  echo "4. Done"
  rm -rf $Output
  rm -rf hairpin.*ebwt
  rm -rf hairpin.sizes ${Mature}.mapped.sam ${Mature}.mapped.bam
fi
