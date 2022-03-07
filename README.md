# Mapping-ROIs-in-genetic-sequences

## ROI-mapping.R
Simple R script that maps density of simple motifs (e.g. CpG) in genetic sequences.

This is a simple script that provides the density of specific motifs (or ROIs; Region Of Interest) across the provided genetic sequences.
* The input FASTA file (RNA, DNA or AA) should be not-aligned, without gaps and special symbols. A sample of an input file containt the H1N1 IAV PR8 strain sequence is provided ("PR8.fas").
* Output CSV-file contains numbers (frequencies) of the specific motif (ROI) across the genetic sequence. A sliding window with a default frame size of 70 nt/aa is implemented to smooth the output. A sample of an output file is provided ("ca.csv").

## dinucleotide-overrepresentation.R

Simple R script that makes a chart and Excel file with dinucleotide over- or underrepresentation.

* The input FASTA file (RNA, DNA or AA) should be not-aligned, without gaps and special symbols. A sample of an input file containt the H1N1 IAV PR8 strain sequence is provided ("PR8.fas").

## Hilbert-map.R

Simple R script that makes a chart with Hilbert curve for the CpG dinucleotide.

* The input FASTA file (DNA) should be not-aligned, without gaps and special symbols. A sample of an input file containt the H1N1 IAV PR8 strain sequence is provided ("PR8.fas").
* You need to have the required packages:
```
install.packages("seqinr")
install.packages("circlize ")
install.packages("ade4")
install.packages("BiocManager")
BiocManager::install("HilbertCurve")
BiocManager::install("IRanges")
```
## CpG-to-CpG-distances.R

Simple R script that makes a chart with distances calculated between CpG dinucleotides.
* The input FASTA file (DNA) should be not-aligned, without gaps and special symbols. A sample of an input file containt the H1N1 IAV PR8 strain sequence is provided ("PR8.fas").
