#CompositionPlot_}cmd.R
#This is part of the PipeSmRNA.sh to plot length distribution based on different small RNA categories
#Yu Sun, 2020-11-25

# CurrPath <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(CurrPath)

args <- commandArgs(TRUE)

FileName <- args[1]
Type <- args[2]

Data <- read.table(FileName, header = F, row.names = 1)
colnames(Data) <- 1:23
Sum <- apply(Data, 2, FUN=sum)

FinalPerc <- Data/Sum*100
Colors <- c("lightskyblue", "lemonchiffon1", "seagreen2", "lightsalmon")

pdf(paste0(FileName, ".pdf"), height = 9, width = 12)
barplot(as.matrix(FinalPerc),ylab="Percentage (%)",xlab="Nucleotide",col=Colors,
        names.arg=1:23,space=0.1,las=1,main=paste0("Nucleotide composition of ",Type," piRNAs"))
mtext("Top to bottom: T, G, C, A")
dev.off()