# Initialization
start.time <- Sys.time()
cat("\014")
library("seqinr")
library("xlsx")

# Reading samples.
# (They should be aligned, no gaps, same length.)
seq <- read.alignment(file="PR8.fas", format="fasta")
seq.length <- nchar(seq$seq[[1]]) # End position
cpg_matrix<-matrix(0, seq$nb, seq.length,dimnames=list(as.list(seq$nam)))

for (i in 1:seq$nb) {
    for (k in 1:(seq.length-1)) {
      if (substr(seq$seq[i], k, k)=="c")  
        if (substr(seq$seq[i], k+1, k+1)=="g") cpg_matrix[i,k]<-1
  }
}

cpg_matrix2<-matrix(0, seq$nb, seq.length)
for (i in 1:seq$nb) {
  frame<-500
  for (k in 1:(seq.length-frame)) {
    if ((cpg_matrix[i,k])==1) cpg_matrix2[i,k:(k+frame)]<-cpg_matrix2[i,k:(k+frame)]+1
  }
}
matplot(t(cpg_matrix2/10), type='o', ylim=c(0,2), xlab='Position', ylab=paste("CpG sites per",frame/10,"nt", sep=" "), pch=15:19)
legend("topright", inset=0.01, legend=seq$nam, col=c(1:6), pch=15:19, bg= ("white"))

# exporting to excel file
write.xlsx(t(cpg_matrix), "positions.xlsx", append = T, col.names = F, row.names = F)

# Cleaning of workaround
cat("\n")
print(Sys.time()-start.time)