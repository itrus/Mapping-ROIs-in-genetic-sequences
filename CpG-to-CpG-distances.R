# Initialization
start.time <- Sys.time()
cat("\014")
library("seqinr")
library("xlsx")

cat("Reading sequences\n")
seq <- read.alignment(file="PR8.fas", format="fasta")

for (CurrentSeq in 1:seq$nb) {
  data2<-words.pos("cg", seq$seq[[CurrentSeq]])
  data3<-NULL
  for (i in 2:length(data2)) data3<-c(data3,data2[i]-data2[i-1]-2)
  hist(data3, las=1, breaks=10, xlim=c(1,140), col="red",
       main=paste(seq$nam[[CurrentSeq]],",",length(data2),"CpGs per",round(length(s2c(seq$seq[[CurrentSeq]]))/1000,1),"kb, M±SD=",
       round(mean(data3),1),"±", round(sd(data3),1)),
       xlab="CpG to CpG distance, nt", ylim=c(0,42))
  write.xlsx(data3, "distances.xlsx", sheetName=seq$nam[CurrentSeq], append = T, col.names = F, row.names = F)
}

# Cleaning of workaround
cat("\n")
print(Sys.time()-start.time)
