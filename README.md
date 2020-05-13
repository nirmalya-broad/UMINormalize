

## Objective
Dual transcriptional profiling of host and bacteria during infection is challenging due to the low abundance of bacterial mRNA. We report Pathogen Hybrid Capture (PatH-Cap), a method to enrich for bacterial mRNA and deplete bacterial rRNA simultaneously from dual RNA-seq libraries using transcriptome-specific probes. 

PatH-Cap generates host and pathogen reads using the [scDual-Seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1340-x) protocol. To circumvent the problem of duplicate reads due to PCR amplifcation, the library preparation step uses <b>U</b>nique <b>M</b>olecular <b>I</b>dentifier (<b>UMI</b>). 

The existing methods for UMI normalization (or UMI collapse) are applicable to reads from eukaryotes and not from prokaryotes. However, due to the presence of operons the transcripts from bacteria can span multiple-gene boundary. Hence, the existing UMI methods for eukaryotes are not able to estimate the transcript count from scDual-Seq accurately. 

## Algorithm
Within each UMI, we assume that reads from each transcript are densely distributed. So, we sort all the reads from an aligned bam with respect to the UMI string and inside UMI we sort them by alignment coordinates of the reads on the genome. If two consecutive reads are spaced by at least N bases (e.g. 500 bases), the algorithm considers that those reads are coming from different transcripts. Otherwise, if the gap is less than M bases, we decide that those two transcripts are coming from the same transcript.

We also have an implementation of UMI normalization for eukaryotes that mostly follow published algorithms.

## Implementation
Our algorithm uses a specialized <b>out of memory k-way merge sort</b> in the order of (1) UMI barcode,  (2) strand specificity of the alignments and finally (3) coordinates of the aligned reads (mainly bacterial reads), boundary of the features/genes (mainly for eukaryotic host reads).

After sorting is done, UMINormalize algorithm segments and collapses the reads based on a heuristic (the gap between the cluster described earlier in the <b>Algorithm</b> section) to obtain reads of mRNA transcripts without any PCR duplicates. <b>Out of memory sort</b> implemented in C++ allows processing of unlimited number of reads.

## Statistical modeling
In this version of the project we segment the reads based on a fix gap (500 bases). However, in another project we developed a statistical model of the gaps based on a mixture of distributions. That part of the project is still in progress. But for concept and code please see the [NB_EM](https://github.com/nirmalya-broad/NB_EM) repository.

## Running UMINormalize
After compiling the way to run UMINormalize is
```
umi_norm -i <infile> -o <outdir> -p <prefix> -c <collapse_type>
```
where, 

<b>infile</b> contains the input sam/bam file containing aligned reads.<br>
<b>outdir</b> points to the path of the output directory.<br>
<b>prefix</b> is a string used as a prefix of output files.<br>
<b>collapse_type</b> is used to specify if the umi collapse is based on coordinates (for bacterial reads) or feature boundaries (used for eukaryotic host reads).

## Cite the project
Betin, V., Penaranda, C., Bandyopadhyay, N. et al. Hybridization-based capture of pathogen mRNA enables paired host-pathogen transcriptional analysis. Sci Rep 9, 19244 (2019). https://doi.org/10.1038/s41598-019-55633-6



 

