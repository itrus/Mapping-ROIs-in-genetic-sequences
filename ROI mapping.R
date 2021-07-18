# Initialization
library("seqinr")
library("magrittr")
library("stringr")

#Defining the reference genome and its length
ReferenceGenome <- read.alignment(file="PR8.fas", format="fasta")
seq.length <- nchar(ReferenceGenome$seq)

#Defining sliding frame size and ROIs (Regions of Interest) and current ROI
FrameSize <- 70
#ROI <- c("a", "t", "c", "g")
ROI <- c("aa", "at", "ag", "ac",
         "ta", "tt", "tg", "tc",
         "ga", "gt", "gg", "gc",
         "ca", "ct", "cg", "cc")
#ROI <- c(
#'aaa', 'caa', 'taa', 'gaa', 
#'aac', 'cac', 'tac', 'gac', 
#'aat', 'cat', 'tat', 'gat', 
#'aag', 'cag', 'tag', 'gag', 
#'aca', 'cca', 'tca', 'gca', 
#'acc', 'ccc', 'tcc', 'gcc', 
#'act', 'cct', 'tct', 'gct', 
#'acg', 'ccg', 'tcg', 'gcg', 
#'ata', 'cta', 'tta', 'gta', 
#'atc', 'ctc', 'ttc', 'gtc', 
#'att', 'ctt', 'ttt', 'gtt', 
#'atg', 'ctg', 'ttg', 'gtg', 
#'aga', 'cga', 'tga', 'gga', 
#'agc', 'cgc', 'tgc', 'ggc', 
#'agt', 'cgt', 'tgt', 'ggt', 
#'agg', 'cgg', 'tgg', 'ggg')
CurrentROI <- ROI[13]
print(CurrentROI)

#Let's find ROIs
ROIperSegment <- list()
for (i in 1:ReferenceGenome$nb) {
  a <- rep(0, seq.length[i])
  ROIperSegment[[i]] <- a
}
names(ROIperSegment) <- ReferenceGenome$nam
for (i in 1:ReferenceGenome$nb) {
  a <- str_locate_all(ReferenceGenome$seq[[i]], CurrentROI)[[1]][ ,1]
  b <- as.vector(ROIperSegment[i])[[1]]
  b[a] <- 1
  ROIperSegment[i][[1]] <- b
}

#Let's smooth our results
ROIperSegmentSmoothed <- list()
for (i in 1:ReferenceGenome$nb) {
  a <- as.vector(ROIperSegment[i])[[1]]
  b <- array(0, seq.length[i])
  for (z in 1:(seq.length[i] - FrameSize)) {
    if (a[z] == 1)
      b[z:(z + FrameSize)] <- b[z:(z + FrameSize)] + 1
  }
  ROIperSegmentSmoothed[i][[1]] <- b
}
names(ROIperSegmentSmoothed) <- ReferenceGenome$nam

#Let's export our results
results <- do.call(rbind, lapply(ROIperSegmentSmoothed, as.data.frame)) %>%
  c %>% unlist %>% as.numeric 
results.matrix <- matrix(NA, nrow = max(seq.length), ncol = ReferenceGenome$nb)
results.cumsum <- cumsum(seq.length)
a <- 1
for(i in 1:ReferenceGenome$nb) {
  results.matrix[1:seq.length[i], i] <- results[a:results.cumsum[i]]
  a <- results.cumsum[i]+1
}
results.df <- as.data.frame(results.matrix)
colnames(results.df) <- ReferenceGenome$nam
results.df <- cbind(1:max(seq.length), results.df)
colnames(results.df)[1] <- "Position"
results.df %>% write.table(file = paste0(CurrentROI, ".csv"), sep = ";",
            row.names = F, na = "")
