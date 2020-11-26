#PlotLendis_cmd.R
#This is part of the PipeSmRNA.sh to plot length distribution based on different small RNA categories
#Output will be Contaminants, miRNA, all smRNA, with and without normalization plots. Total 6 plots.
#Yu Sun, 2020-11-25

# CurrPath <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(CurrPath)

args <- commandArgs(TRUE)

Prefix <- args[1]
genome <- args[2]
Norm <- as.numeric(args[3])

rRNA <- read.table(paste0("smRNA_other/",Prefix,".",genome,".rRNA.Lendis"), row.names = 1, col.names = c("Length","Count"))
tRNA <- read.table(paste0("smRNA_other/",Prefix,".",genome,".tRNA.Lendis"), row.names = 1, col.names = c("Length","Count"))
snRNA <- read.table(paste0("smRNA_other/",Prefix,".",genome,".snRNA.Lendis"), row.names = 1, col.names = c("Length","Count"))
snoRNA <- read.table(paste0("smRNA_other/",Prefix,".",genome,".snoRNA.Lendis"), row.names = 1, col.names = c("Length","Count"))

AnnoMiRNA <- read.table(paste0("smRNA_miRNA/",Prefix,".clean.",genome,".miRNAAnno.Lendis"), row.names = 1, col.names = c("Length","Count"))
OtherSpeciesMiRNA <- read.table(paste0("smRNA_miRNA/",Prefix,".clean.",genome,".miRNAAllSpecies.Lendis"), row.names = 1, col.names = c("Length","Count"))

GenomeMapping21 <- read.table(paste0("genome_mapping/",Prefix,".final.Less21.",genome,".Lendis"), row.names = 1, col.names = c("Length","Count"))
GenomeMapping2123 <- read.table(paste0("genome_mapping/",Prefix,".final.21-23.",genome,".Lendis"), row.names = 1, col.names = c("Length","Count"))
GenomeMapping23 <- read.table(paste0("genome_mapping/",Prefix,".final.Over23.",genome,".Lendis"), row.names = 1, col.names = c("Length","Count"))

Unmapped <- read.table(paste0("genome_mapping/",Prefix,".final.unmapped.Lendis"), row.names = 1, col.names = c("Length","Count"))

FourTypes <- t(cbind(rRNA,tRNA,snRNA,snoRNA))
MiRNAs <- t(cbind(AnnoMiRNA,OtherSpeciesMiRNA))
FinalTable <- t(cbind(Unmapped,rRNA,tRNA,snRNA,snoRNA,AnnoMiRNA,OtherSpeciesMiRNA,
                      GenomeMapping21,GenomeMapping2123,GenomeMapping23))
Mapped <- t(cbind(rRNA,tRNA,snRNA,snoRNA,AnnoMiRNA,OtherSpeciesMiRNA,
                      GenomeMapping21,GenomeMapping2123,GenomeMapping23))

clustercol <- c("gray65","tomato","slateblue1","tan1","#984EA3","steelblue1","gold","peru","orchid","darkseagreen","dimgrey","firebrick","deepskyblue1","darkgoldenrod1","darkolivegreen1","bisque4","coral","chartreuse1","deeppink","aquamarine1","burlywood1")
FourTypes_col <- clustercol[2:5]
MiRNAs_col <- clustercol[6:7]
FinalTable_col <- clustercol[1:10]
Mapped_col <- clustercol[2:10]

#Figures using raw reads
pdf(paste0(Prefix,".Lendis.raw.4contaminants.pdf"),width = 12,height = 8)
par(mar = c(4.5,6,3,3))
barplot(FourTypes,ylab="",xlab="Length",col=FourTypes_col,
        names.arg=18:35,space=0.4,las=1,main=paste0("Length distribution of ",Prefix," 4 contaminants"))
legend("topright", y=max(FourTypes)*0.85,legend = rev(c("rRNA","tRNA","snRNA","snoRNA")),fill =rev(FourTypes_col), cex = 0.75)
mtext("Reads",side = 2, line = 4.5)
dev.off()

pdf(paste0(Prefix,".Lendis.raw.miRNAs.pdf"),width = 12,height = 8)
par(mar = c(4.5,6,3,3))
barplot(MiRNAs,ylab="",xlab="Length",col=MiRNAs_col,
        names.arg=18:35,space=0.4,las=1,main=paste0("Length distribution of ",Prefix," miRNAs"))
legend("topright", y=max(MiRNAs)*0.85,legend = rev(c("miRNA annotated","miRNA other species")),fill =rev(MiRNAs_col), cex = 0.75)
mtext("Reads",side = 2, line = 4.5)
dev.off()

pdf(paste0(Prefix,".Lendis.raw.all.pdf"),width = 12,height = 8)
par(mar = c(4.5,6,3,3))
barplot(FinalTable,ylab="",xlab="Length",col=FinalTable_col,
        names.arg=18:35,space=0.4,las=1,main=paste0("Length distribution of ",Prefix," all small RNAs"))
legend("topright", y=max(FinalTable)*0.85,
       legend = c("Unmapped","rRNA","tRNA","snRNA","snoRNA","miRNA annotated","miRNA other species",
                                                     "Other smRNA","siRNA","piRNA"),
       fill =FinalTable_col, cex = 0.75)
mtext("Reads",side = 2, line = 4.5)
dev.off()

pdf(paste0(Prefix,".Lendis.raw.mapped.pdf"),width = 12,height = 8)
par(mar = c(4.5,6,3,3))
barplot(Mapped,ylab="",xlab="Length",col=Mapped_col,
        names.arg=18:35,space=0.4,las=1,main=paste0("Length distribution of ",Prefix," mapped small RNAs"))
legend("topright", y=max(Mapped)*0.85,
       legend = c("rRNA","tRNA","snRNA","snoRNA","miRNA annotated","miRNA other species",
                  "Other smRNA","siRNA","piRNA"),
       fill =Mapped_col, cex = 0.75)
mtext("Reads",side = 2, line = 4.5)
dev.off()


#Figures using normalized reads
pdf(paste0(Prefix,".Lendis.normalized.4contaminants.pdf"),width = 12,height = 8)
par(mar = c(4.5,6,3,3))
barplot(FourTypes/Norm,ylab="",xlab="Length",col=FourTypes_col,
        names.arg=18:35,space=0.4,las=1,main=paste0("Length distribution of ",Prefix," 4 contaminants"))
legend("topright", y=max(FourTypes/Norm)*0.85,legend = rev(c("rRNA","tRNA","snRNA","snoRNA")),fill =rev(FourTypes_col), cex = 0.75)
mtext("Normalized reads (RPM, reads per million)",side = 2, line = 4.5)
dev.off()

pdf(paste0(Prefix,".Lendis.normalized.miRNAs.pdf"),width = 12,height = 8)
par(mar = c(4.5,6,3,3))
barplot(MiRNAs/Norm,ylab="",xlab="Length",col=MiRNAs_col,
        names.arg=18:35,space=0.4,las=1,main=paste0("Length distribution of ",Prefix," miRNAs"))
legend("topright", y=max(MiRNAs/Norm)*0.85,legend = rev(c("miRNA annotated","miRNA other species")),fill =rev(MiRNAs_col), cex = 0.75)
mtext("Normalized reads (RPM, reads per million)",side = 2, line = 4.5)
dev.off()

pdf(paste0(Prefix,".Lendis.normalized.all.pdf"),width = 12,height = 8)
par(mar = c(4.5,6,3,3))
barplot(FinalTable/Norm,ylab="",xlab="Length",col=FinalTable_col,
        names.arg=18:35,space=0.4,las=1,main=paste0("Length distribution of ",Prefix," all small RNAs"))
legend("topright", y=max(FinalTable)*0.85,
       legend = c("Unmapped","rRNA","tRNA","snRNA","snoRNA","miRNA annotated","miRNA other species",
                  "Other smRNA","siRNA","piRNA"),
       fill =FinalTable_col, cex = 0.75)
mtext("Normalized reads (RPM, reads per million)",side = 2, line = 4.5)
dev.off()

pdf(paste0(Prefix,".Lendis.normalized.mapped.pdf"),width = 12,height = 8)
par(mar = c(4.5,6,3,3))
barplot(Mapped/Norm,ylab="",xlab="Length",col=Mapped_col,
        names.arg=18:35,space=0.4,las=1,main=paste0("Length distribution of ",Prefix," mapped small RNAs"))
legend("topright", y=max(Mapped)*0.85,
       legend = c("rRNA","tRNA","snRNA","snoRNA","miRNA annotated","miRNA other species",
                  "Other smRNA","siRNA","piRNA"),
       fill =Mapped_col, cex = 0.75)
mtext("Normalized reads (RPM, reads per million)",side = 2, line = 4.5)
dev.off()
