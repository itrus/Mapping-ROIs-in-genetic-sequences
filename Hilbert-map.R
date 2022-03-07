# Initialization
start.time <- Sys.time()
cat("\014")
library(seqinr)
library(HilbertCurve)
library(circlize)
library(shape)
library(ade4)
library(IRanges)

# Reading samples.
# (They should be aligned, no gaps, same length, if possible)
seq <- read.alignment(file="PR8.fas", format="fasta")
seq.length <- nchar(seq$seq[[1]]) # End position
cpg_matrix<-matrix(0, seq$nb, seq.length,dimnames=list(as.list(seq$nam)))

for (i in 1:seq$nb) {
  seq.length <- nchar(seq$seq[[i]]) # End position
  for (k in 1:(seq.length-1)) {
    if (substr(seq$seq[i], k, k)=="c")  
      if (substr(seq$seq[i], k+1, k+1)=="g") cpg_matrix[i,k]<-1
  }
  hc = HilbertCurve(1, seq.length, level = 5, title=row.names(cpg_matrix)[i])
  s<-0
  e<-0
  for (k in 1:seq.length){
    if (cpg_matrix[i,k]=="1") {
      s<-c(s,k)
      e<-c(e,k)
    }
  }
  ir = IRanges(s[2:length(s)], e[2:length(e)])
  hc_segments(hc, IRanges(1, seq.length))
  hc_points(hc, ir, np = NULL, size = unit(4,"mm"), pch = 17)
}

# Cleaning of workaround
cat("\n")
print(Sys.time()-start.time)

