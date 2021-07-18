# Mapping-ROIs-in-genetic-sequences
Simple R script that maps the density of simple motifs (e.g. CpG) in genetic sequences

This is a simple script that provides the density of specific motifs (or ROIs; Region Of Interest) across the provided genetic sequences.
The input FASTA file (RNA, DNA or AA) should be not aligned and without gaps and special symbols. A sample of an input file is provided.
Output CSV-file contains numbers (frequencies) of the specific motif (ROI) across the genetic sequence. A sliding window with a default frame size of 70 nt/aa is implemented to smooth the output. A sample of an output file is provided.
