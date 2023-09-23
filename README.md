# Mapping-ROIs-in-genetic-sequences

## ROI-mapping.R
Simple R script that maps density of simple motifs (e.g. CpG) in genetic sequences.

This is a simple script that provides the density of specific motifs (or ROIs; Region Of Interest) across the provided genetic sequences.
* The input FASTA file (RNA, DNA or AA) should be not-aligned, without gaps and special symbols. A sample of an input file containing the H1N1 IAV PR8 strain sequence is provided ("PR8.fas").
* Output CSV-file contains numbers (frequencies) of the specific motif (ROI) across the genetic sequence. A sliding window with a default frame size of 70 nt/aa is implemented to smooth the output. A sample of an output file for the **CA** dinucleotide is provided ("ca.csv").

## dinucleotide-overrepresentation.R

Simple R script that makes a chart and Excel file with dinucleotide over- or underrepresentation.

* The input FASTA file (RNA, DNA or AA) should be not-aligned, without gaps and special symbols. A sample of an input file containt the H1N1 IAV PR8 strain sequence is provided ("PR8.fas").

Typical output chart:

![Rplot01](https://user-images.githubusercontent.com/9166776/157060514-25d45f6d-028e-4938-a191-74a973b5d829.png)


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

Typical output chart:

![Rplot02](https://user-images.githubusercontent.com/9166776/157060563-fc23db8d-3b47-495f-ae2f-a7fe79fd55d5.png)

## CpG-to-CpG-distances.R

Simple R script that makes a chart with distances calculated between CpG dinucleotides.

* The input FASTA file (DNA) should be not-aligned, without gaps and special symbols. A sample of an input file containt the H1N1 IAV PR8 strain sequence is provided ("PR8.fas").

Typical output chart:

![Rplot](https://user-images.githubusercontent.com/9166776/157060843-02efa4e8-43b1-4081-b519-10ccabafc47a.png)

## CpG-map-linear.R

Simple R script that makes a chart with the density of CpG dinucleotides.

* The input FASTA file (DNA) should be not-aligned, without gaps and special symbols. A sample of an input file containt the H1N1 IAV PR8 strain sequence is provided ("PR8.fas").

Typical output chart:

![Rplot](https://user-images.githubusercontent.com/9166776/157075731-0debf368-f4fa-491b-b097-a37f11c9cce6.png)

