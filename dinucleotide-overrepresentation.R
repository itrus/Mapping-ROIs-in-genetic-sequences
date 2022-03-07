library("seqinr")
library("xlsx")
seq<-read.alignment(file="PR8.fas", format="fasta")
underrepresentation<-c()
estfreq<-function(nuc){
  k<-count(s2c(nuc), word=1, freq=T)
  return(c(k[1]*k[1],k[1]*k[2],k[1]*k[3],k[1]*k[4],
           k[2]*k[1],k[2]*k[2],k[2]*k[3],k[2]*k[4],
           k[3]*k[1],k[3]*k[2],k[3]*k[3],k[3]*k[4],
           k[4]*k[1],k[4]*k[2],k[4]*k[3],k[4]*k[4]))
}
for (i in 1:seq$nb) {
  underrepresentation<-rbind(underrepresentation,(count(s2c(seq$seq[[i]]), word=2, freq=T)/estfreq((seq$seq[[i]]))))
}
boxplot(underrepresentation)
write.xlsx(cbind(seq$nam,underrepresentation),"results.xlsx")
