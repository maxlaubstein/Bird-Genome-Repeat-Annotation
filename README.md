# Bird-Genome-Repeat-Annotation
A script automating repeat annotation (using RepeatModeler &amp; RepeatMasker) for avian genomes, following the protocol of Daren Card: https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
Assumes you have working installations of RepeatModeler, RepeatMasker, seqkit, ProcessRepeats, rmOutToGFF3custom, and bedtools in your path.

Usage:
./bird_repeat_annotation.sh /path/to/genome/assembly species_shorthand_name
e.g.
./bird_repeat_annotation.sh /ncbi_dataset/data/GCF_048771995.1/GCF_048771995.1_bTaeGut7.mat_genomic.fna taeGut
