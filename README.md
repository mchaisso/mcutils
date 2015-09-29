This is a set of tools I find helpful for working with whole genome
alignments stored in the SAM format.  Blasr (from this repository 
only: https://github.com/mchaisso/blasr ), and bwa-mem may both be
used to generate whole-genome alignments.

To build:
    cd src && make && make install

The blasr command is:
		sawriter target_genome.fasta
    blasr query_genome.fasta target_genome.fasta -alignContigs -maxAnchorGap 30000 -sam -minMapQV 30


The programs are:
+ samToBed - This version of a sam <-> bed converter allows the calculation of identity
+ samToDot - Create a file that may be rendered into a dotplot using RenderPlot.R
+ samLiftover - Transform coordinates from the query genome to target using the SAM file as a mapping.




