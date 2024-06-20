# polySV
Code used in:<br/>
Hämälä T, Moore C, Cowan L, Carlile M, Gopaulchan D, Brandrud MK, Birkeland S, Loose M, Kolář F, Koch MA & Yant L (2024). Impact of whole-genome duplications on structural variant evolution in _Cochlearia_. Nature Communications. https://doi.org/10.1038/s41467-024-49679-y<br>
<br>
prune_ld.c: A program for conducting LD-pruning on mixed ploidy VCF files.<br>
poly_sfs.c: A program for estimating SFS from mixed ploidy VCF files.<br>
poly_fst.c: A program for estimating pairwise Fst and Dxy from mixed ploidy VCF files.<br>
poly_freq.c: A program for estimating allele frequencies from mixed ploidy VCF files.<br>
est_sfs_updog.r: An R script for estimating SFS and Tajima's D from genotype probabilities.<br>
est_cov_pca.r: An R script for conducting PCA on mixed ploidy VCF files.<br>
est_adapt_dist.r: An R script for estimating and plotting the distance between SV and SNP-based climatic landscapes.<br>
<br>
Instructions for compiling and running the C software in Unix-like operating systems are provided in the individual source code files. The software and scripts were tested on macOS 14 (C compiler: Clang 15.0.0) and R v4.3.0.  
