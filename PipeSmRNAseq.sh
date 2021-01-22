#!/bin/bash
#PipeSmRNAseq.sh
#This is a pipeline for small RNAseq FASTQ analysis
#Inputs: Trimmed (e.g use FastqAdapterTimmer) small RNAseq single-end raw reads in FASTQ format
#This pipeline will separate rRNA, tRNA, snRNA, snoRNA, reads before charactering miRNAs, and finally piRNAs
#    For miRNA analysis, we map reads to known annotation, all other species miRNA annotation, and also predict novel miRNAs
#	 For piRNA analysis, we perform genome mapping and TE mapping.
#Increase the number of genome mapping (or using -reportAll) times will cause the pipeline and piRNA cluster defining step very slow.
#Example: 
#PipeSmRNAseq.sh -i Embryo.SmRNAseq.fastq.gz -g mm10 -noqc
#Version1 basic write up: 2020-11-15, Y.S

array=( "$@" )

#Usage
if [ ! -n "$1" ]
then
  echo "********************************************************************************************"
  echo "*                   PipeSmRNAseq: pipeline for Small RNAseq analysis.                      *"
  echo "*                            Version 1, 2020-11-15, Y.S                                    *"
  echo "* Usage: `basename $0`                                                                   *"
  echo "*        Required (single end data):                                                       *"
  echo "*                  -i [Data.fastq.gz]                                                      *"
  echo "*                  -g [mm10/susScr11/Ppup]                                                 *"
  echo "*        Optional: -p [Number of CPUs, default=1]                                          *"
  echo "*                  -m [Genome mapping mismatch, default=1]                                 *"
  echo "*                  -k [int] [Genome mapping times, default=1]                              *"
  echo "*                  -reportAll [Report all genome mapping locations, default off]           *"
  echo "*                  -pre Sequences.fa [Run pre-mapping to Sequences.fa]                     *"
  echo "*                  -noqc [Suppress fastqc]                                                 *"
  echo "*                  -normFactor [Num] [User defined normalization factor, to divide]        *"
  echo "*                  -piCluster  [Predict piRNA cluster, default off]                        *"
  echo "* Inputs: The fastq file needs to be trimmed, and has .fastq (or.fastq.gz) as suffix       *"
  echo "* Run: Default to run fastqc, rRNA, tRNA, snRNA, snoRNA, miRNA & genomic mapping           *"
  echo "*      Figures will be generated in /plots folder, and bigWig files in /tracks folder      *"
  echo "* Indexes: prepare genome.rRNA.fa, as weel as tRNA, snRNA, snoRNA, hairpin, mature and TE  *"
  echo "*      under genome/Sequence path, and this pipeline will generate bowtie indexes for you  *"
  echo "* Normalization: use per million genomic mapping reads: RPM                                *"
  echo "* Outputs: All output files will be generated in the same folder as the pipeline submitted *"
  echo "* This pipeline requires additional scripts in PipelineHomeDir/bin folder                  *"
  echo "********************************************************************************************"
  exit 1
fi

echo "*****************************************************************************"
echo "*            PipeSmRNAseq: pipeline for Small RNAseq analysis.              *"
echo "*                       Version 1, 2020-11-15, Y.S                          *"
echo "*****************************************************************************"
echo "0. Loading softwares:"
#Get current time
St=`date`

#Load softwares: This part may need to be changed when using different cluster
#Main tools:STAR, cufflinks, salmon, featureCounts (from Subread package)
#Other tools: samtools

#Module load for Bluehive:
#module load fastqc
#module load bowtie
#module load bedtools
#module load samtools
#module load fastqc         #Optional when -noqc specified.
#module load rnastar        #Use STAR to rescue some piRNA junction reads.
#module load r
#module load weblogo

#Get pipeline directory
HomeDir=$(dirname `readlink -f $0`)
#For mac local check only:
#HomeDir="/Users/yusun/Downloads/PipelineHomeDir"
echo "   Home Directory:"
echo "   "$HomeDir

echo "1. Resolving inputs:"

#Get parameters
Data="unassigned"
datastats=0
genome="unassigned"
genomestats=0
CPU=1
GenomeMM=1               #Default Genome mapping MM
LibraryType="fr-firststrand"
StrandType="unassigned"
runfastqc=1              #You can change this default to allow/suppress fastqc
premap=0
premapData="unassigned"
gMappingTimes=1
normDepth=1
reportAll=0
normFactor="unassigned"
piCluster=0

for arg in "$@"
do
 if [[ $arg == "-i" ]]
  then
    Data=${array[$counter+1]}
    echo '   Single end data: '$Data
 elif [[ $arg == "-g" ]]
  then
    genome=${array[$counter+1]}
    echo '   Genome: '$genome
 elif [[ $arg == "-p" ]]
  then
    CPU=${array[$counter+1]}
    echo '   CPU: '$CPU
 elif [[ $arg == "-m" ]]
  then
    GenomeMM=${array[$counter+1]}
    echo '   Genome mapping mismatch: '$GenomeMM
 elif [[ $arg == "-k" ]]
  then
    gMappingTimes=${array[$counter+1]}
    echo '   Pre mapping file: '$gMappingTimes
 elif [[ $arg == "-pre" ]]
  then
    premapData=${array[$counter+1]}
    premap=1
    echo '   Pre mapping file: '$premapData
 elif [[ $arg == "-reportAll" ]]
  then
    reportAll=${array[$counter+1]}
    reportAll=1
    echo '   Report all mapping positions'
 elif [[ $arg == "-noqc" ]]
  then
    runfastqc=0
    echo '   Suppress fastqc quality control'
 elif [[ $arg == "-normDepth" ]]
  then
    normCDS=1
    echo '   Use genome mapping reads for normalization'
 elif [[ $arg == "-normFactor" ]]
  then
    normFactor=${array[$counter+1]}
    echo '   Use user defined norm factor: '$normFactor
 elif [[ $arg == "-piCluster" ]]
  then
    piCluster=1
    echo '   Run piRNA cluster prediction: '
 fi
  let counter=$counter+1
done

#Get current directory and create folders or files
[ $runfastqc == "1" ] && FastqcDir=fastqc && mkdir -p $FastqcDir
[ $premap == "1" ] && PremapDir=pre_mapping && mkdir -p $PremapDir
OtherDir=smRNA_other && mkdir -p $OtherDir
miRNADir=smRNA_miRNA && mkdir -p $miRNADir
piRNADir=smRNA_piRNA && mkdir -p $piRNADir
GenomeMappingDir=genome_mapping && mkdir -p $GenomeMappingDir
FiguresDir=plots && mkdir -p $FiguresDir
TracksDir=tracks && mkdir -p $TracksDir

#Check data
if [ $Data == "unassigned" ];then
  echo "     >>> [Error]: Please input data!"
else
  datastats=1
fi

#Getting Output Suffix and default strand types
OutputSuffix=$(echo $Data|sed 's/.fastq.*//g')
TABLE=${OutputSuffix}.summary

#Check genome
if [ $genome == "unassigned" ];then
  echo "     >>> [Error]: Please assign genome file!"
else
  if [ -d $HomeDir/$genome ];  then
    echo "     >>> This genome is supported."
    genomestats=1
  else
    echo "     >>> [Error]: Genome unsupported!"
  fi
fi

#Export checking status
if [ $datastats == "1" -a $genomestats ==  "1" ];then
    echo "   Output Suffix: "$OutputSuffix
  echo "   Status check: pass"
else
  echo "   Status check: failed, stop the pipeline"
  exit 1
fi

calculating(){ awk "BEGIN { print "$*" }"; }

#In case we don't have miRNA coordinates bed6 generated before, this function will generate for you.
GenerateMiRNACoordinates () {
  Hairpin=$1
  Mature=$2
  OutputFinal=$3
  Output=$OutputFinal.t
  OutputFinalFull=$OutputFinal.full
  rm -rf $Output $OutputFinal.t $OutputFinalFull
  calculating(){ awk "BEGIN { print "$*" }"; }
  
  echo "     1. Indexing"
  bowtie-build --quiet $Hairpin hairpin
  awk '{print $1}' $Hairpin|grep ">" |sed 's/>//' > temp.hairpin.names
  awk '{print $1}' $Hairpin|grep -v ">" |awk '{print length($0)}' > temp.hairpin.len
  paste temp.hairpin.names temp.hairpin.len > hairpin.sizes && rm -rf temp.hairpin.names temp.hairpin.len
  
  echo "     2. Mapping"
  bowtie -S -f -v 0 -k 1 --no-unal hairpin $Mature > ${Mature}.mapped.sam
  samtools view -b ${Mature}.mapped.sam > ${Mature}.mapped.bam 
  bedtools bamtobed -i ${Mature}.mapped.bam | awk '$6=="+"' > $Output

  echo "     3. Defining 5p and 3p"
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
    HalfLengthp=`awk -v t=$Precursor '{if ($1==t) print $2/2}' hairpin.sizes`
    HalfLength=${HalfLengthp%.*}
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

    if [[ $PreArm == $DefineArm ]];then
      Type="consistent"
    elif [[ $PreArm != $DefineArm ]];then
      if [[ $PreArm == "xx" ]];then
        Type="novel"
      else
        #This final case is very rare, for mouse, only has one case (mmu-mir-466i), when the miRNA across the middle of the hairpin
        DefineArm=$PreArm  #Assign the annotation
        Type="inconsistent"
      fi
    fi
    echo $whole|awk -v newarm=$DefineArm 'BEGIN{OFS="\t"}{$5=newarm;print $0}' >> $OutputFinal
    echo $whole|awk -v newarm=$DefineArm -v type=$Type 'BEGIN{OFS="\t"}{$5=newarm;$7=type;print $0}' >> $OutputFinalFull
  done

  echo "     4. Done"
  rm -rf $Output
  rm -rf hairpin.*ebwt
  rm -rf hairpin.sizes ${Mature}.mapped.sam ${Mature}.mapped.bam
}

#Calculate length distribution for unique mapping sam files
GetUniqueLendis() {
	Data=$1
	RefLength=$2
	grep -v "@" $Data |awk '{print length($10)}' |sort|uniq -c|awk '{OFS="\t";print $2,$1}' > ${Data%.*}.Lendis.t
	MergeCountsAndFillZeros.py ${Data%.*}.Lendis.t $RefLength ${Data%.*}.Lendis
	rm -rf ${Data%.*}.Lendis.t
}

#Calculate length distribution for multi mapping sam files
GetMultiLendis() {
	Data=$1
	RefLength=$2
	grep -v "@" $Data |awk '{print $1,length($10)}' |uniq|awk '{print $2}' |sort|uniq -c|awk '{OFS="\t";print $2,$1}' > ${Data%.*}.Lendis.t
	MergeCountsAndFillZeros.py ${Data%.*}.Lendis.t $RefLength ${Data%.*}.Lendis
	rm -rf ${Data%.*}.Lendis.t
}

echo "*****************************************************************************"
echo "2. Checking dependencies:"
echo "Annotations:"
#Genome files and index check
#This check includes folders: /Annotation, /Sequence, /Index
#                    files: ${genome}.RefSeq.gtf, ${genome}.RefSeq.bed12
#                    files: $genome.fa, $genome.RefSeq.fa
#If indexes not generated, then check these files under /Sequence:
#${genome}.rRNA.fa, ${genome}.tRNA.fa, ${genome}.snRNA.fa, ${genome}.snoRNA.fa
#${genome}.hairpin.fa, ${genome}.mature.fa, ${genome}.TE.fa
dependenciescount=0

#Check FASTA files
if [ -s $HomeDir/$genome/Sequence/${genome}.fa ];then
  echo "   Genome: "${genome}.fa
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] Genome lost "
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.hairpin.fa ];then
  echo "   annotated miRNA precursors: "${genome}.hairpin.fa
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] annotated miRNA precursor lost "
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.mature.fa ];then
  echo "   annotated mature miRNAs: "${genome}.mature.fa
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] annotated mature miRNA lost "
fi

###Check Index paths
echo "Index files:"
if [ -d $HomeDir/$genome/Index/Sm_rRNAIndex ];then
  echo "   rRNA index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] rRNA index lost"
  echo "   Try to build one..."
  mkdir $HomeDir/$genome/Index/Sm_rRNAIndex
  bowtie-build --quiet $HomeDir/$genome/Sequence/${genome}.rRNA.fa $HomeDir/$genome/Index/Sm_rRNAIndex/Sm_rRNAIndex && \
      let dependenciescount=$dependenciescount+1 && echo "   Index done"
fi
if [ -d $HomeDir/$genome/Index/Sm_tRNAIndex ];then
  echo "   tRNA index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] tRNA index lost"
  echo "   Try to build one..."
  mkdir $HomeDir/$genome/Index/Sm_tRNAIndex
  bowtie-build --quiet $HomeDir/$genome/Sequence/${genome}.tRNA.fa $HomeDir/$genome/Index/Sm_tRNAIndex/Sm_tRNAIndex && \
      let dependenciescount=$dependenciescount+1 && echo "   Index done"
fi
if [ -d $HomeDir/$genome/Index/Sm_snRNAIndex ];then
  echo "   snRNA index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] snRNA index lost "
  echo "   Try to build one..."
  mkdir $HomeDir/$genome/Index/Sm_snRNAIndex
  bowtie-build --quiet $HomeDir/$genome/Sequence/${genome}.snRNA.fa $HomeDir/$genome/Index/Sm_snRNAIndex/Sm_snRNAIndex && \
      let dependenciescount=$dependenciescount+1 && echo "   Index done"
fi
if [ -d $HomeDir/$genome/Index/Sm_snoRNAIndex ];then
  echo "   snoRNA index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] snoRNA index lost"
  echo "   Try to build one..."
  mkdir $HomeDir/$genome/Index/Sm_snoRNAIndex
  bowtie-build --quiet $HomeDir/$genome/Sequence/${genome}.snoRNA.fa $HomeDir/$genome/Index/Sm_snoRNAIndex/Sm_snoRNAIndex && \
      let dependenciescount=$dependenciescount+1 && echo "   Index done"
fi
if [ -d $HomeDir/$genome/Index/Sm_miRNAIndex ];then
  echo "   miRNA precursor index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] miRNA precursor index lost "
  echo "   Try to build one..."
  mkdir $HomeDir/$genome/Index/Sm_miRNAIndex
  bowtie-build --quiet $HomeDir/$genome/Sequence/${genome}.hairpin.fa $HomeDir/$genome/Index/Sm_miRNAIndex/Sm_miRNAIndex && \
      let dependenciescount=$dependenciescount+1 && echo "   Index done"
fi
if [ -d $HomeDir/$genome/Index/Sm_miRNAMatureIndex ];then
  echo "   mature miRNA index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] mature miRNA index lost "
  echo "   Try to build one..."
  mkdir $HomeDir/$genome/Index/Sm_miRNAMatureIndex
  bowtie-build --quiet $HomeDir/$genome/Sequence/${genome}.mature.fa $HomeDir/$genome/Index/Sm_miRNAMatureIndex/Sm_miRNAMatureIndex && \
      let dependenciescount=$dependenciescount+1 && echo "   Index done"
fi
if [ -d $HomeDir/$genome/Index/Sm_miRNAAllSpeciesIndex ];then
  echo "   miRNA all species precursor index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] miRNA all species precursor index lost "
  echo "   Try to build one..."
  mkdir $HomeDir/$genome/Index/Sm_miRNAAllSpeciesIndex
  bowtie-build --quiet $HomeDir/bin/AllSpecies.hairpin.fa $HomeDir/$genome/Index/Sm_miRNAAllSpeciesIndex/Sm_miRNAAllSpeciesIndex && \
      let dependenciescount=$dependenciescount+1 && echo "   Index done"
fi
if [ -d $HomeDir/$genome/Index/Sm_TE ];then
  echo "   TE index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] TE index lost "
  echo "   Try to build one..."
  mkdir $HomeDir/$genome/Index/Sm_TE
  bowtie-build --quiet $HomeDir/$genome/Sequence/${genome}.TE.fa $HomeDir/$genome/Index/Sm_TE/Sm_TE && \
      let dependenciescount=$dependenciescount+1 && echo "   Index done"
fi
if [ -d $HomeDir/$genome/Index/STARIndex ];then
  echo "   STAR index exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] STAR index lost"
fi

#Check other dependencies
if [ -s $HomeDir/$genome/Annotation/${genome}.miRNA.Matching.bed6 ];then
  echo "   annotated miRNA coordinates: "${genome}.miRNA.Matching.bed6
  let dependenciescount=$dependenciescount+1
else
  echo "   [Warning] annotated miRNA coordinates lost, try to build one:"
  echo "   Running: PipelineHomeDir/bin/DefineMatureMiRNACoordinates.sh"
  GenerateMiRNACoordinates $HomeDir/${genome}/Sequence/${genome}.hairpin.fa $HomeDir/${genome}/Sequence/${genome}.mature.fa ${genome}.miRNA.Matching.bed6
  mv ${genome}.miRNA.Matching.bed6* $HomeDir/${genome}/Annotation && echo "   Successfully built ${genome}.miRNA.Matching.bed6" && let dependenciescount=$dependenciescount+1
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ];then
  echo "   ChromInfo: "${genome}.ChromInfo.txt
  let dependenciescount=$dependenciescount+1
else
  $HomeDir/bin/faSize -tab -detailed $HomeDir/$genome/Sequence/${genome}.fa > $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt && \
  echo "   [Error] ChromInfo lost, generate a new one.." && let dependenciescount=$dependenciescount+1
fi

#echo $dependenciescount
if [ $dependenciescount == "14" ];then
  echo "   Dependencies check: pass"
else
  echo "   [Error] Some dependencies lost "
fi

calculating(){ awk "BEGIN { print "$*" }"; }

echo "*****************************************************************************"
echo "3. Run fastqc quality control:"
if [ $runfastqc == "1" ];then
echo "   Running fastqc single-end mode"
fastqc \
    -f fastq \
    -o $FastqcDir \
    $Data \
    2> $FastqcDir/${OutputSuffix}.fastqc.log && \
    echo "   Done fastqc"
elif [ $runfastqc == "0" ];then
echo "   Skipping fastqc..."
fi

echo "*****************************************************************************"
echo "4. rRNA mapping and removal:"
echo "   Running bowtie single-end mode"
bowtie \
	-S \
	-k 1 \
	-v 3 \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/Sm_rRNAIndex/Sm_rRNAIndex \
    -q $Data \
    --un $OtherDir/${OutputSuffix}.${genome}.No_rRNA.fastq \
    1> $OtherDir/${OutputSuffix}.${genome}.rRNA.sam \
    2> $OtherDir/${OutputSuffix}.${genome}.rRNA.log && \
    samtools view -b $OtherDir/${OutputSuffix}.${genome}.rRNA.sam > $OtherDir/${OutputSuffix}.${genome}.rRNA.bam && \
    echo "   Done rRNA mapping"
	rRNAReads=`head -2 $OtherDir/${OutputSuffix}.${genome}.rRNA.log | tail -1 | awk '{print $9}'`
	for i in {18..35};do echo $i;done > ${OutputSuffix}.RefLendis
	echo "   Calculating length distribution"
	GetUniqueLendis $OtherDir/${OutputSuffix}.${genome}.rRNA.sam ${OutputSuffix}.RefLendis
	rm -rf $OtherDir/${OutputSuffix}.rRNA.${genome}.sam

echo "*****************************************************************************"
echo "5. tRNA mapping and removal:"
echo "   Running bowtie single-end mode"
bowtie \
	-S \
	-k 1 \
 	-v 3 \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/Sm_tRNAIndex/Sm_tRNAIndex \
    -q $OtherDir/${OutputSuffix}.${genome}.No_rRNA.fastq  \
    --un $OtherDir/${OutputSuffix}.${genome}.No_rRNA.No_tRNA.fastq \
    1> $OtherDir/${OutputSuffix}.${genome}.tRNA.sam \
    2> $OtherDir/${OutputSuffix}.${genome}.tRNA.log && \
    samtools view -b $OtherDir/${OutputSuffix}.${genome}.tRNA.sam > $OtherDir/${OutputSuffix}.${genome}.tRNA.bam && \
    echo "   Done tRNA mapping"
	tRNAReads=`head -2 $OtherDir/${OutputSuffix}.${genome}.tRNA.log | tail -1 | awk '{print $9}'`
	echo "   Calculating length distribution"
	GetUniqueLendis $OtherDir/${OutputSuffix}.${genome}.tRNA.sam ${OutputSuffix}.RefLendis
	rm -rf $OtherDir/${OutputSuffix}.${genome}.tRNA.sam

echo "*****************************************************************************"
echo "6. snRNA mapping and removal:"
echo "   Running bowtie single-end mode"
bowtie \
	-S \
	-k 1 \
 	-v 3 \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/Sm_snRNAIndex/Sm_snRNAIndex \
    -q $OtherDir/${OutputSuffix}.${genome}.No_rRNA.No_tRNA.fastq  \
    --un $OtherDir/${OutputSuffix}.${genome}.No_rRNA.No_tRNA.No_snRNA.fastq \
    1> $OtherDir/${OutputSuffix}.${genome}.snRNA.sam \
    2> $OtherDir/${OutputSuffix}.${genome}.snRNA.log && \
    samtools view -b $OtherDir/${OutputSuffix}.${genome}.snRNA.sam > $OtherDir/${OutputSuffix}.${genome}.snRNA.bam && \
    echo "   Done snRNA mapping"
	snRNAReads=`head -2 $OtherDir/${OutputSuffix}.${genome}.snRNA.log | tail -1 | awk '{print $9}'`
	echo "   Calculating length distribution"
	GetUniqueLendis $OtherDir/${OutputSuffix}.${genome}.snRNA.sam ${OutputSuffix}.RefLendis
	rm -rf $OtherDir/${OutputSuffix}.${genome}.snRNA.sam

echo "*****************************************************************************"
echo "7. snoRNA mapping and removal:"
echo "   Running bowtie single-end mode"
bowtie \
	-S \
	-k 1 \
 	-v 3 \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/Sm_snoRNAIndex/Sm_snoRNAIndex \
    -q $OtherDir/${OutputSuffix}.${genome}.No_rRNA.No_tRNA.No_snRNA.fastq  \
    --un $OtherDir/${OutputSuffix}.${genome}.No_rRNA.No_tRNA.No_snRNA.No_snoRNA.fastq \
    1> $OtherDir/${OutputSuffix}.${genome}.snoRNA.sam \
    2> $OtherDir/${OutputSuffix}.${genome}.snoRNA.log && \
    samtools view -b $OtherDir/${OutputSuffix}.${genome}.snoRNA.sam > $OtherDir/${OutputSuffix}.${genome}.snoRNA.bam && \
    echo "   Done snoRNA mapping"
	snoRNAReads=`head -2 $OtherDir/${OutputSuffix}.${genome}.snoRNA.log | tail -1 | awk '{print $9}'`
	echo "   Calculating length distribution"
	GetUniqueLendis $OtherDir/${OutputSuffix}.${genome}.snoRNA.sam ${OutputSuffix}.RefLendis
	rm -rf $OtherDir/${OutputSuffix}.${genome}.snoRNA.sam

echo "*****************************************************************************"
echo "8. Map cleaned reads to annotated miRNAs:"
echo "   a, Mapping to miRNA precursors for quantification"
echo "      Multi-mapping miRNAs will be reported to only one location"
bowtie \
	-S \
	-k 1 \
 	-v 1 \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/Sm_miRNAIndex/Sm_miRNAIndex \
    -q $OtherDir/${OutputSuffix}.${genome}.No_rRNA.No_tRNA.No_snRNA.No_snoRNA.fastq  \
    --un $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.fastq \
    1> $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.sam \
    2> $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.log && \
    echo "      Done miRNA precursor mapping"
    echo "      Calculating length distribution"
    GetUniqueLendis $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.sam ${OutputSuffix}.RefLendis
    echo "      Summarizing miRNA counts"
	miRNAAnnoReads=`head -2 $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.log | tail -1 | awk '{print $9}'`
	samtools view -b $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.sam > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bam
	bedtools bamtobed -i $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bam > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bed6.t
	grep -v "@" $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.sam | awk '{print $10}' > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.sam.seq
	paste $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bed6.t $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.sam.seq > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bed7 && \
		rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bed6.t
	
	Matching=$HomeDir/$genome/Annotation/${genome}.miRNA.Matching.bed6
	awk '{print $4}' $Matching > $miRNADir/${OutputSuffix}.miRNAAnno.list
	rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.summary.txt && touch $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.summary.txt
	bedtools intersect -a $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bed7 -b $Matching -wa -wb -s -f 0.8 > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.overlap
	miRNAAnnoReads_sense=`wc -l $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.overlap | awk '{print $1}'`
	awk '{print $11}' $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.overlap |sort|uniq -c|awk '{OFS="\t";print $2,$1}' > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.overlap.t
	python $HomeDir/bin/MergeCountsAndFillZeros.py $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.overlap.t $miRNADir/${OutputSuffix}.miRNAAnno.list $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.summary.txt
	rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.sam* $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bed7 $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.overlap.t
	mv $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.overlap $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno.bed13

echo "   b, Mapping to mature miRNAs for quantification (strict mode)"
echo "      Multi-mapping miRNAs will be reported to only one location"
bowtie \
	-S \
	-k 1 \
 	-v 1 \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/Sm_miRNAMatureIndex/Sm_miRNAMatureIndex \
    -q $OtherDir/${OutputSuffix}.${genome}.No_rRNA.No_tRNA.No_snRNA.No_snoRNA.fastq  \
    1> $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.sam \
    2> $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.log && \
    echo "      Done mature miRNA mapping"
    echo "      Calculating length distribution"
    GetUniqueLendis $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.sam ${OutputSuffix}.RefLendis
    echo "      Summarizing miRNA counts"
	miRNAAnnoReads_mature=`head -2 $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.log | tail -1 | awk '{print $9}'`
	samtools view -F 0x10 -b $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.sam > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bam   #Only keep sense mapping
	bedtools bamtobed -i $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bam > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bed6.t
	samtools view $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bam | awk '{print $10}' > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bam.seq
	paste $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bed6.t $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bam.seq > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bed7 && \
		rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bed6.t $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bam.seq
	miRNAAnnoReads_mature_sense=`wc -l $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bed7 | awk '{print $1}'`

	Matching=$HomeDir/$genome/Annotation/${genome}.miRNA.Matching.bed6
	rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.summary.txt && touch $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.summary.txt
	awk '{print $1}' $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.bed7 | sort|uniq -c|awk '{OFS="\t";print $2,$1}' > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.summary.t
	python $HomeDir/bin/MergeCountsAndFillZeros.py $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.summary.t $miRNADir/${OutputSuffix}.miRNAAnno.list $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.summary.txt
	rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.summary.t $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAnno_mature.sam

echo "*****************************************************************************"
echo "9. Map the rest reads to miRNAs from all known species in miRBase:"
echo "   Multi-mapping miRNAs will be reported to only one location"
bowtie \
	-S \
	-k 1 \
 	-v 1 \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/Sm_miRNAAllSpeciesIndex/Sm_miRNAAllSpeciesIndex \
    -q $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.fastq  \
    --al $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.miRNAAllSpecies.fastq \
    --un $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.fastq \
    1> $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.sam \
    2> $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.log && \
    echo "   Done all species miRNA precursor mapping"
    echo "   Calculating length distribution"
    GetUniqueLendis $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.sam ${OutputSuffix}.RefLendis
    echo "   Summarizing miRNA counts"
	miRNAAnnoReads_AllSpecies=`head -2 $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.log | tail -1 | awk '{print $9}'`
	samtools view -F 0x10 -b $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.sam > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bam   #Only keep sense mapping
	FilteredMapped=`samtools view $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bam | wc -l`
	bedtools bamtobed -i $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bam > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bed6.t
	samtools view $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bam | awk '{print $10}' > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bam.seq
	paste $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bed6.t $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bam.seq > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bed7 && \
		rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bed6.t $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bam.seq
	MatchingAll=$HomeDir/bin/AllSpecies.miRNA.Matching.bed6
	awk '{print $4}' $MatchingAll > $miRNADir/${OutputSuffix}.miRNAAllSpecies.list
	
	rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.summary.txt && touch $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.summary.txt
	bedtools intersect -a $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bed7 -b $MatchingAll -wa -wb -s -f 0.8 > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.overlap
	miRNAAnnoReads_AllSpecies_sense=`wc -l $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.overlap | awk '{print $1}'`
	awk '{print $11}' $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.overlap |sort|uniq -c|awk '{OFS="\t";print $2,$1}' > $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.overlap.t
	python $HomeDir/bin/MergeCountsAndFillZeros.py $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.overlap.t $miRNADir/${OutputSuffix}.miRNAAllSpecies.list $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.summary.txt
	rm -rf $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.summary.t $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.overlap.t $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.sam
	mv $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.overlap $miRNADir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.bed13

	if [[ "$miRNAAnnoReads_AllSpecies" == 0 ]]; then
		echo "   No reads mapped to other species"
	else
		echo "   Some reads mapped to other species"
		echo "   Perform genome mapping for new miRNA discovery, reporting at most 5 positions"
		bowtie \
		-S \
		-k 5 \
		-v 1 \
		-p $CPU \
		--no-unal \
		$HomeDir/$genome/Index/bowtieIndex/bowtieIndex \
		-q $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.miRNAAllSpecies.fastq  \
		1> $GenomeMappingDir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.${genome}.sam \
		2> $GenomeMappingDir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.${genome}.log && \
		samtools view -b $GenomeMappingDir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.${genome}.sam > $GenomeMappingDir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.${genome}.bam && \
		rm -rf $GenomeMappingDir/${OutputSuffix}.clean.${genome}.miRNAAllSpecies.${genome}.sam && \
		echo "   Done miRNA genome mapping using miRNAs mapped to all species"
	fi

echo "*****************************************************************************"
echo "10. Separateing small miRNAs (<21), siRNA (21-23 nt) and piRNA length (>23 nt)"
echo "   Generating: "${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.Less21.fastq
echo "   Generating: "${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.21-23.fastq
echo "   Generating: "${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.Over23.fastq
RestReads=`python $HomeDir/bin/SeparateFastq_siRNA_piRNA.py $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.fastq $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies`

echo "*****************************************************************************"
echo "11. Genome mapping with "$GenomeMM" mismatch(es)"
echo "   Report "$gMappingTimes" time(s) for multimapping reads, choosing the best one"
echo "   a, Aligning small RNAs: length < 21 nt"
if [ -s "$miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.Less21.fastq" ];then
bowtie \
	-S \
	--best \
	-k $gMappingTimes \
 	-v $GenomeMM \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/bowtieIndex/bowtieIndex \
    -q $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.Less21.fastq  \
    --un $GenomeMappingDir/${OutputSuffix}.final.Less21.unmapped.fastq \
    1> $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.sam \
    2> $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.log
    echo "      Calculating length distribution"
    if [ $gMappingTimes == "1" ];then
    	GenomeMappedLess21=`head -2 $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.log | tail -1 | awk '{print $9}'`
    	GetUniqueLendis $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.sam ${OutputSuffix}.RefLendis
    else
    	GenomeMappedLess21=`grep -v "@" $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.sam | awk '{print $1}'|sort|uniq|wc -l`
    	GetMultiLendis $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.sam ${OutputSuffix}.RefLendis
    fi
    samtools view -b $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.sam > $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.bam && \
    rm -rf $GenomeMappingDir/${OutputSuffix}.final.Less21.${genome}.sam
    echo "      Done small RNA mapping: length < 21 nt"
    
else
	echo "      No small RNAs to map, length < 21 nt"
fi

echo "   b, Aligning small RNAs: 21 nt < length < 23 nt"
if [ -s "$miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.21-23.fastq" ];then
bowtie \
	-S \
	--best \
	-k $gMappingTimes \
 	-v $GenomeMM \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/bowtieIndex/bowtieIndex \
    -q $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.21-23.fastq  \
    --un $GenomeMappingDir/${OutputSuffix}.final.21-23.unmapped.fastq \
    1> $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.sam \
    2> $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.log
    echo "      Calculating length distribution"
    if [ $gMappingTimes == "1" ];then
    	GenomeMappedLess2123=`head -2 $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.log | tail -1 | awk '{print $9}'`
    	GetUniqueLendis $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.sam ${OutputSuffix}.RefLendis
    else
    	GenomeMappedLess21=`grep -v "@" $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.sam | awk '{print $1}'|sort|uniq|wc -l`
    	GetMultiLendis $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.sam ${OutputSuffix}.RefLendis
    fi
    samtools view -b $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.sam > $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.bam && \
    rm -rf $GenomeMappingDir/${OutputSuffix}.final.21-23.${genome}.sam
    echo "      Done small RNA mapping: 21 nt < length < 23 nt"
    
else
	echo "      No small RNAs to map, length 21 nt < length < 23 nt"
fi

echo "   c, Aligning small RNAs: length > 23 nt"
if [ -s "$miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.Over23.fastq" ];then
bowtie \
	-S \
	--best \
	-k $gMappingTimes \
 	-v $GenomeMM \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/bowtieIndex/bowtieIndex \
    -q $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.Over23.fastq \
    --un $GenomeMappingDir/${OutputSuffix}.final.Over23.unmapped.fastq \
    1> $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.sam \
    2> $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.log
    echo "      Calculating length distribution"
    if [ $gMappingTimes == "1" ];then
    	GenomeMappedOver23=`head -2 $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.log | tail -1 | awk '{print $9}'`
    	GetUniqueLendis $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.sam  ${OutputSuffix}.RefLendis
    else
    	GenomeMappedOver23=`grep -v "@" $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.sam | awk '{print $1}'|sort|uniq|wc -l`
    	GetMultiLendis $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.sam ${OutputSuffix}.RefLendis
    fi
    samtools view -b $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.sam > $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.bam
    echo "      Extracting sequences"
    python $HomeDir/bin/ExtractSAMRawReads.py $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.sam
    echo "      Done small RNA mapping: length > 23 nt"
	cd $piRNADir && ln -s ../genome_mapping/${OutputSuffix}.final.Over23.${genome}.bam ${OutputSuffix}.final.Over23.${genome}.bam && cd ../
	rm -rf $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.sam
else
	echo "      No small RNAs to map, length > 23 nt"
fi

echo "*****************************************************************************"
echo "12. Direct transposon mapping using reads > 23 nt, reporting all alignments"
if [ -s "$miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.Over23.fastq" ];then
bowtie \
	-S \
	--best \
	-a \
 	-v 2 \
	-p $CPU \
	--no-unal \
    $HomeDir/$genome/Index/Sm_TE/Sm_TE \
    -q $miRNADir/${OutputSuffix}.clean.${genome}.No_miRNAAnno.No_miRNAAllSpecies.Over23.fastq  \
    1> $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam \
    2> $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.log
	TEMapped=`grep -v "@" $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam | awk '{print $1}'|sort|uniq|wc -l`
    GetMultiLendis $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam ${OutputSuffix}.RefLendis
    samtools view -b $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam > $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.TE.bam && \
    rm -rf $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.TE.sam
    echo "   Done small RNA mapping: length > 23 nt"
else
	echo "   No small RNAs to map, length > 23 nt"
fi

echo "*****************************************************************************"
echo "13. Merging length distributions of each step"
rm -rf $TABLE
InputReads=`grep 'reads processed' $OtherDir/${OutputSuffix}.${genome}.rRNA.log |head -1| awk '{print $4}'`
CleanningReads=$((rRNAReads+tRNAReads+snRNAReads+snoRNAReads))
GenomeUnmappedReads=`calculating $RestReads-$GenomeMappedLess21-$GenomeMappedLess2123-$GenomeMappedOver23`
TotalMapped=`calculating $CleanningReads+$miRNAAnnoReads+$miRNAAnnoReads_AllSpecies+$GenomeMappedLess21+$GenomeMappedLess2123+$GenomeMappedOver23`
if [ $normFactor == "unassigned" ];then
  echo "   Use Reads Per Million (RPM) as normalization factor"
  Norm=`calculating $TotalMapped/1000000`
else
  Norm=$normFactor
fi

echo "   Calculating unmapped reads length distribution"
python $HomeDir/bin/CalculateFastqLendis.py $GenomeMappingDir/${OutputSuffix}.final.Less21.unmapped.fastq
python $HomeDir/bin/CalculateFastqLendis.py $GenomeMappingDir/${OutputSuffix}.final.21-23.unmapped.fastq
python $HomeDir/bin/CalculateFastqLendis.py $GenomeMappingDir/${OutputSuffix}.final.Over23.unmapped.fastq
paste $GenomeMappingDir/${OutputSuffix}.final.Less21.unmapped.fastq.Lendis $GenomeMappingDir/${OutputSuffix}.final.21-23.unmapped.fastq.Lendis \
	$GenomeMappingDir/${OutputSuffix}.final.Over23.unmapped.fastq.Lendis | awk '{OFS="\t";print $1,$2+$4+$6}' > $GenomeMappingDir/${OutputSuffix}.final.unmapped.Lendis
echo "   Generating Lendis figures"
Rscript --vanilla $HomeDir/bin/PlotLendis_cmd.R ${OutputSuffix} ${genome} ${Norm} > /dev/null && mv ${OutputSuffix}*pdf $FiguresDir && \
	echo "   8 figures generated"

echo "*****************************************************************************"
echo "14. Getting statistics:"

echo -e "   total input reads:\t${InputReads}"
echo -e "   rRNA_reads:\t${rRNAReads}"
echo -e "   tRNA_reads:\t${tRNAReads}"
echo -e "   snRNA_reads:\t${snRNAReads}"
echo -e "   snoRNA_reads:\t${snoRNAReads}"
echo -e "   4_above_types_sum:\t${CleanningReads}"
echo -e "   miRNA_reads_annotated:\t${miRNAAnnoReads}"
echo -e "   miRNA_reads_annotated_sense_counted:\t${miRNAAnnoReads_sense}"
echo -e "   miRNA_reads_annotated_strict:\t${miRNAAnnoReads_mature}"
echo -e "   miRNA_reads_annotated_strict_sense_counted:\t${miRNAAnnoReads_mature_sense}"
echo -e "   miRNA_reads_all_species:\t${miRNAAnnoReads_AllSpecies}"
echo -e "   miRNA_reads_all_species_sense_counted:\t${miRNAAnnoReads_AllSpecies_sense}"
echo -e "   After removing above reads:"
echo -e "   clean_genome_mapping_Less21nt:\t${GenomeMappedLess21}"
echo -e "   clean_genome_mapping_21-23nt:\t${GenomeMappedLess2123}"
echo -e "   clean_genome_mapping_Over23nt:\t${GenomeMappedOver23}"
echo -e "   TE_mapping:\t${TEMapped}"
echo -e "   Total_mapped:\t${TotalMapped}"
echo -e "   Normalization:\t${Norm}"
echo -e "   Unmapped:\t${GenomeUnmappedReads}"

echo -e "total input reads:\t${InputReads}" > $TABLE
echo -e "rRNA_reads:\t${rRNAReads}" >> $TABLE
echo -e "tRNA_reads:\t${tRNAReads}" >> $TABLE
echo -e "snRNA_reads:\t${snRNAReads}" >> $TABLE
echo -e "snoRNA_reads:\t${snoRNAReads}" >> $TABLE
echo -e "4_above_types_sum:\t${CleanningReads}" >> $TABLE
echo -e "miRNA_reads_annotated:\t${miRNAAnnoReads}" >> $TABLE
echo -e "miRNA_reads_annotated_sense_counted:\t${miRNAAnnoReads_sense}" >> $TABLE
echo -e "miRNA_reads_annotated_strict:\t${miRNAAnnoReads_mature}" >> $TABLE
echo -e "miRNA_reads_annotated_strict_sense_counted:\t${miRNAAnnoReads_mature_sense}" >> $TABLE
echo -e "miRNA_reads_all_species:\t${miRNAAnnoReads_AllSpecies}" >> $TABLE
echo -e "miRNA_reads_all_species_sense_counted:\t${miRNAAnnoReads_AllSpecies_sense}" >> $TABLE
echo -e "After removing above reads:" >> $TABLE
echo -e "clean_genome_mapping_Less21nt:\t${GenomeMappedLess21}" >> $TABLE
echo -e "clean_genome_mapping_21-23nt:\t${GenomeMappedLess2123}" >> $TABLE
echo -e "clean_genome_mapping_Over23nt:\t${GenomeMappedOver23}" >> $TABLE
echo -e "TE_mapping:\t${TEMapped}" >> $TABLE
echo -e "Total_mapped:\t${TotalMapped}" >> $TABLE
echo -e "Normalization:\t${Norm}" >> $TABLE
echo -e "Unmapped:\t${GenomeUnmappedReads}" >> $TABLE

echo "*****************************************************************************"
echo "15. piRNA analysis"
echo "   Formatting mapped piRNA reads (> 23 nt)"
bedtools bamtobed -i $piRNADir/${OutputSuffix}.final.Over23.${genome}.bam > $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed
paste $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed $GenomeMappingDir/${OutputSuffix}.final.Over23.${genome}.sam.Reads | awk '{OFS="\t";print $1,$2,$3,$7,0,$6}' > \
	$piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6
sort $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6 > $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted

echo "   Generating piRNA summary file for PILFER: "${OutputSuffix}.final.Over23.${genome}.bed6.sorted.pilfer
uniq -c $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted | awk '{OFS="\t";print $2,$3,$4,$5"::PI",$1,$7}' |sort -V -k1,1 -k2,2 > \
	$piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.pilfer

echo "   piRNA nucleotide composition calculation (use the first 23nt)"
cut -f4 $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted |cut -c1-23 > $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer
python $HomeDir/bin/CalculateNucleotideCompo.py $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer
cat $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer.NtCompo|awk '{print "   ",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}'

echo "   generating piRNA composition figure"
Rscript --vanilla $HomeDir/bin/CompositionPlot_cmd.R $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer.NtCompo "genome mapping" > /dev/null
mv $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer.NtCompo.pdf $FiguresDir/

echo "   Analyzing TE mapping piRNAs"
python $HomeDir/bin/ExtractSAMRawReads.py $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam
samtools view -b $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam > $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.bam
bedtools bamtobed -i $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.bam > $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.bed
paste $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.bed $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.Reads |awk '{print $4,$7}'|uniq|awk '{print $2}'|cut -c1-23 > \
	$piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer
python $HomeDir/bin/CalculateNucleotideCompo.py $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer
echo "   piRNA nucleotide composition calculation (use the first 23nt)"
cat $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer.NtCompo|awk '{print "   ",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}'

echo "   generating piRNA composition figure for TE mapping piRNAs"
Rscript --vanilla $HomeDir/bin/CompositionPlot_cmd.R $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer.NtCompo "TE mapping" > /dev/null
mv $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer.NtCompo.pdf $FiguresDir/

Weblogo=`weblogo --help 2> /dev/null|grep Usage|wc -l`
if [ $Weblogo == "1" ];then
	echo "   Checking... Weblogo is installed. "
	echo "   Plot SeqLogo figures for genome mapping piRNAs"
	sed 's/T/U/g' $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer > $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer.RNA
	weblogo --composition none -F PDF -A rna -W 15 -n 62 -c classic -S 2 < $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer.RNA > \
		$piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer.RNA.pdf
	
	echo "   Plot SeqLogo figures for TE mapping piRNAs"	
	sed 's/T/U/g' $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer > $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer.RNA
	weblogo --composition none -F PDF -A rna -W 15 -n 62 -c classic -S 2 < $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer.RNA > \
		$piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer.RNA.pdf
		
	mv $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.23mer.RNA.pdf $FiguresDir/
	mv $piRNADir/${OutputSuffix}.final.Over23.${genome}.TE.sam.23mer.RNA.pdf $FiguresDir/
fi

if [ $piCluster == "1" ];then
	echo "   Run piRNA cluster prediction using PILFER... This can be time consuming"
	python $HomeDir/bin/pilfer.py -i $piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.pilfer > \
		$piRNADir/${OutputSuffix}.final.Over23.${genome}.bed6.sorted.pilfer.piCluster.bed
fi

#Cleaning
rm -rf ${OutputSuffix}.RefLendis


echo "Time Started: "$St
Ed=`date`
echo "Time Ended:   "$Ed
echo "*                           End of the pipeline                             *"
echo "*****************************************************************************"
