## ASM-Clust

#### Overview:
The three scripts in this repository comprise ASM-Clust, an approach to de novo classify complex protein superfamilies.
ASM-Clust is implemented in bash with two helper scripts in perl, and will take a protein fasta file as the sole input. 

Fasta files are processed with ASM_clust.sh, which then:

- randomly selects a subset of n sequences (default 1000)
- aligns the entire dataset to the subset of n sequences
- combines all scores into a matrix (inserting 0 for query-database pairs that did not produce an alignment)
- reduces the matrix to 2 dimensions using t-SNE (Van der Maaten and Hinton 2008; Van der Maaten 2014)


#### Usage:
ASM_clust.sh [options] fasta_file   
Options:   
-s INTEGER (subset of sequences for matrix, default 1000)   
-p INTEGER (t-SNE perplexity value, default 1000)   
-m INTEGER (max iterations of t-SNE, default 5000)   
-a mmseqs2/diamond/blast (aligner, one of three options, default mmseqs2)   
-t INTEGER (threads for aligner)   
-f FILE (fasta file of sequences to use as references)   


Although the clustering is generally similar with multiple randomly chosen subsets, the subset can be provided as a fasta file with the -f flag for reproducibility.
The output of ASM_clust.sh can be visualized as a scatterplot where each dot represents a sequence, and clusters are readily apparent. This format allows additional annotation with sequence features, such as taxonomy, length, or composition.


#### External dependencies:
1) The python wrapper for the Barnes-Hut implementation of t-SNE available here: https://github.com/lvdmaaten/bhtsne 

2) Alignment software. For flexible usage, ASM-Clust supports alignment using:
- MMSeqs2 (Steinegger and Söding 2017) (default aligner)
- DIAMOND (Buchfink, Xie, and Huson 2015)
- BLAST (Altschul et al. 1990)


#### Recommended setup:
For reproducibility I recommend running ASM-Clust in a conda environment with the dependencies installed.

1) Make a directory for ASM-Clust and enter the directory.   
```mkdir ASM_Clust```   
```cd ASM_Clust```   

2) clone this repository    
```git clone https://github.com/dspeth/ASM_clust.git```   

3) clone the bh-tSNE repository, and compile the bh_tsne executable   
```git clone https://github.com/lvdmaaten/bhtsne.git```   
```cd bhtsne```   
```g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2```   

4) set up an anaconda environment with the alignment software installed   
The alignment software versions specified below were used when testing ASM-Clust.   
```conda create -n ASM_clust```   
```conda activate ASM_clust```   
```conda install -c bioconda -c conda-forge mmseqs2=11.e1a1c```   
```conda install -c bioconda -c conda-forge blast=2.9.0```   
```conda install -c bioconda -c conda-forge diamond=0.9.24```   

5) add links to the relevant executables to the anaconda environment /bin directory to place them in you PATH when the environment is loaded.   
```ln -s /absolute/path/to/ASM_clust.sh /absolute/path/to/ASM_clust/conda/env/bin```   
```ln -s /absolute/path/to/merge_tab_files.pl /absolute/path/to/ASM_clust/conda/env/bin```   
```ln -s /absolute/path/to/tab_seq_lookup.pl /absolute/path/to/ASM_clust/conda/env/bin```   
```ln -s /absolute/path/to/bhtsne.py /absolute/path/to/ASM_clust/conda/env/bin```   
```ln -s /absolute/path/to/bh_tsne /absolute/path/to/ASM_clust/conda/env/bin```   

© 2020 GitHub, Inc.
