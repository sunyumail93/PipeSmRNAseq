# PipeSmRNAseq
A comprehensive pipeline for small RNAseq data analysis.

Some path/files are shared with other Pipe* pipelines, such as the PipelineHomeDir folder, and scripts in /bin folder.

## Software prerequisites
This pipeline is designed to run on Linux servers, and requires the following softwares:
```
R
Python2
STAR
bowtie
bedtools
samtools
fastqc (optional)
weblogo (optional)
PILFER (The main script from PILFER has been included in the ./bin, and no need to install again)
```
PipeSmRNAseq.sh is the pipeline main script. Other dependencies are in ./bin folder

One UCSC tool (from http://hgdownload.soe.ucsc.edu/admin/exe/) is used: faSize.

For more information about the PILFER piRNA cluster prediction tool: https://github.com/rishavray/PILFER/blob/master/pipeline.sh

## Pipeline setup

Here is an example of mm10 genome setup.

1, Download scripts from github to Linux server:

```
git clone https://github.com/sunyumail93/PipeSmRNAseq.git
mv PipeSmRNAseq PipelineHomeDir
```

2, Set up index files for genome mapping
Use Pteromalus puparum (Ppup) genoem as an example

```
cd PipelineHomeDir/Ppup/Sequence

```
3, Set up index files:
```
mkdir PipelineHomeDir/Ppup/Index
cd PipelineHomeDir/Ppup/Index

#STAR index:
mkdir STARIndex
STAR --runMode genomeGenerate --genomeDir STARIndex --genomeFastaFiles ../Sequence/Ppup.fa --sjdbGTFfile ../Annotation/Ppup.RefSeq.reduced.bed12.geneid.gtf --sjdbOverhang 100

#bowtie index:
mkdir bowtieIndex
bowtie-build ../Sequence/Ppup.fa bowtieIndex/bowtieIndex

#Other index files can be automatically generated when running the pipline.
#Make sure the sequence files for rRNA, tRNA, snRNA, snoRNA, TE are existed in the ./Ppup/Sequence folder
```
4, Compile the C++ program if necessary

The FastqAdapterTimmer binary file was compiled on Linux system. It may need to be compiled again:

```
g++ FastqAdapterTimmer.cpp -o FastqAdapterTimmer
```

5, Add executable permissions

```
chmod +x PipeSmRNAseq.sh
chmod +x ./bin/faSize
```

## Pipeline components
```
PipelineHomeDir/
    ├── PipeSmRNAseq.sh
    ├── bin/
    └── Ppup/
      └── Annotation/
        ├── Ppup.RefSeq.reduced.bed12
        └── Ppup.RefSeq.reduced.bed12.geneid.gtf
      └── Index/
        ├── bowtieIndex
        ├── Sm_miRNAAllSpeciesIndex/
        ├── Sm_miRNAIndex/
        ├── Sm_miRNAMatureIndex/
        ├── Sm_rRNAIndex/
        ├── Sm_snoRNAIndex/
        ├── Sm_snRNAIndex/
        ├── Sm_TE/
        ├── Sm_tRNAIndex/
        └── STARIndex/
      └── Sequence/
        ├── Ppup.ChromInfo.txt
        ├── Ppup.fa
        ├── Ppup.fai
        ├── Ppup.hairpin.fa
        ├── Ppup.mature.fa
        ├── Ppup.rRNA.fa
        ├── Ppup.snoRNA.fa
        ├── Ppup.snRNA.fa
        ├── Ppup.TE.fa
        └── mPpup.tRNA.fa
    └── susScr11/
       ...
```

Notes: 

1, For Annotation folder, the GTF file is required for STAR index generation

2, For Index folder, STAR and bowtie whole genome indexes are big, and it would be better to generate them separately. Other directly bowtie mapping indexes are small and can be generated automatically when running the pipeline.

3, Whole genome sequence and other direct mapping sequences (rRNA, tRNA, snRNA, snoRNA, TE) are in the Sequence folder.

## Usage

Type the pipeline name, then you will see the manual page:

```
PipeSmRNAseq.sh
```

Manual page:

![](images/Usages.png)

The trimming script also has manual page:

```
FastqAdapterTimmer
```

![](images/TrimmerUsages.png)

Sometimes we need more trimming at the 5'-end, for example, trim the first 7nt from the left:

```
FastqSeqClipLeft.py
```



## Examples

A regular run using mostly default parameters:

```
PipeSmRNAseq.sh -i Data.fastq.gz -g Ppup
```

More parameters used. The piRNA cluster prediction step may be time consuming:

```
PipeSmRNAseq.sh -i Data.fastq.gz -g Ppup -p 8 -piCluster
```

## Run a real data to test the pipeline

1, Download data

Use a public dataset: [GEO SRA: SRR9851822](https://www.ncbi.nlm.nih.gov/sra/SRX6606439[accn])

`fastq-dump` is part of [NCBI SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software):

```
fastq-dump --split-3 SRR9851822
```

2, Trim adaptor using the trimming script from /bin repository: FastqAdapterTimmer

```
FastqAdapterTimmer -i SRR9851822.fastq -a AGATCG -o Data.trimmed.fastq
gzip Data.trimmed.fastq
```

3, Run the small RNAseq pipeline:

```
PipeSmRNAseq.sh -i Data.trimmed.fastq.gz -g Ppup -p 8 -piCluster
```

## Outputs

1, Length distribution of small RNA all reads, separated by categories:

2, Nucleotide composition analysis
