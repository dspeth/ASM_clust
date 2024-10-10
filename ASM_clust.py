#!/usr/bin/env python

# import required libraries
import argparse
import sys
from pathlib import Path
import random
from Bio import SeqIO
import subprocess
import pandas as pd

# set up command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("fasta_file", help="fastA file to cluster")

parser.add_argument("-V", "--version", action="store_true", help="show script version and exit")
parser.add_argument("-a", "--aligner", default="diamond", type=str, help="choose aligner, can be \"mmseqs\", \"diamond\", or \"blastp\"")
parser.add_argument("-ss", "--subset_size", type=int, help="Define number of sequences to be picked as reference subset, not compatible with defined subset file")
parser.add_argument("-sf", "--subset_file", help="Path to fastA file containing reference subset, not compatible with defined subset size")
parser.add_argument("-p", "--threads", type=int, default=1, help="Set number of threads to use for aligner, default 1")

parser.add_argument("-t", "--tsne", action="store_true", help="Use tSNE for dimentsionality reduction, can be used in addition to or instead of UMAP")
parser.add_argument("-tp", "--tsne_perplexity", default=500, type=int, help="Set upper tSNE perplexity parameter for multiscale tsne. Default 500, change for smaller datasets")
parser.add_argument("-tm", "--tsne_iterations", default=500, type=int, help="Set iterations for tSNE, applies to both exaggerated and non-exaggerated phase of tSNEembedding, default 500")
parser.add_argument("-tx", "--tsne_exaggeration", default=6, type=int, help="Set early exaggeration parameter for tSNE, important for clustering. Default 6")

parser.add_argument("-u", "--umap", action="store_true", help="Use UMAP for dimensionality reduction, can be used in addition to or instead of tSNE")
parser.add_argument("-un", "--umap_n_neighbors", help="Set n_neighbor parameter of umap. Default ???")
parser.add_argument("-ud", "--umap_min_dist", help="Set min_dist parameter of umap. Range 0-0.99, default ???")

parser.add_argument("-mp", "--metadata_protein", help="path to tab delimited metadata file referencing protein sequence IDs, header must be 'prot_ID'")
parser.add_argument("-mg", "--metadata_genome", help="path to tab delimited metadata file referencing genome accession IDs, header must be 'genome_ID'")

# parse args
args = parser.parse_args()

fasta = Path(args.fasta_file)

# check existence of input files
if args.version:
    print("python3 implementation of ASM-clust, version 1")
    sys.exit(0)

if not fasta.is_file():
    print("provided input fastA file does not exist")
    sys.exit(0)

if args.subset_file:
    subset = Path(args.subset_file)
    if not subset.is_file():
        print("subset specified, but subset fastA file does not exist")
        sys.exit(0)

if args.metadata_genome:
    genome_metadata = Path(args.metadata_genome)
    if not genome_metadata.is_file():
        print("genome metadata file specified, but metadata file does not exist")
        sys.exit(0)

if args.metadata_protein:
    protein_metadata = Path(args.metadata_protein)
    if not protein_metadata.is_file():
        print("protein metadata file specified, but metadata file does not exist")
        sys.exit(0)

# check for used aligner in path
def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    from shutil import which
    return which(name) is not None

if not is_tool(args.aligner):
    print("You selected {} to align sequences, but it is not available in PATH".format(args.aligner))
    sys.exit(0)

# check first line of fasta file and ref file for ">"
def fasta_check(name):
    """Check whether 'name' is a fasta file"""
    with open(name) as f:
        firstline = f.readline().rstrip()
        if not firstline.startswith(">"):
            print("{} does not seem to be a fastA file".format(name))
            sys.exit(0)

fasta_check(fasta)

if args.subset_file:
    subset = Path(args.subset_file)
    fasta_check(subset)

### possible other issues:
# illegal characters in seq ids? both files to <- don't know yet
# - perp can not be larger than data/3 <- automated in openTSNE


# get number of sequences
num_seqs = len([1 for line in open(fasta) if line.startswith(">")])
print("There are {} sequences in the dataset, proceeding to subset selection".format(num_seqs))


# pick subset

basename = fasta.stem

def fasta_subsample(fasta_file, subset_file, subset_size):
    """
    randomly sample 'subset_size' number of seqs from a fasta file
    and write a new file containing the subset in single line fasta
    """

    seq_index = SeqIO.index(fasta_file, "fasta")
    subset_keys = random.sample(seq_index.keys(), subset_size)
    subset_dict = {k: seq_index[k] for k in subset_keys}
    SeqIO.write(subset_dict.values(), subset_file, "fasta-2line")

    print("written {} randomly selected sequences from {} to {}, proceeding to alignment".format(subset_size, fasta_file, subset_file))

if not args.subset_file:
    if args.subset_size is None:
        args.subset_size = 1000
    subset_size = args.subset_size
    if subset_size > num_seqs:
        print("WARNING: subset size ({}) is larger than number of sequences in dataset ({}) setting subset size to {}".format(subset_size, num_seqs, num_seqs))
        subset_size = num_seqs
    random_subset = basename+"_subset_"+str(subset_size)+".faa"
    fasta_subsample(str(fasta), random_subset, subset_size)
else:
    print("using specified file {} as subset, proceeding to alignment".format(subset))


# align
if args.subset_file:
    align_subset = str(subset)
else:
    align_subset = random_subset

if args.aligner == "diamond":
    dbname = basename+"_subset_dmnd"
    run_dmnd_makedb = ["diamond", "makedb", "--in", align_subset, "-d", dbname, "-p", str(args.threads)]
    run_dmnd_blastp = ["diamond", "blastp", "-q", str(fasta), "-d", dbname, "-o", basename+"_align", "-p", str(args.threads), "-k", str(subset_size), "--sensitive", "--masking", "0", "--outfmt", "6", "qseqid", "sseqid", "score", "--comp-based-stats", "0"]
    subprocess.run(run_dmnd_makedb, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run(run_dmnd_blastp, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    print("aligned total sequence set to subset using DIAMOND, proceeding to matrix construction")

elif args.aligner == "mmseqs":
    run_mmseqs = ["mmseqs", "easy-search", str(fasta), align_subset, basename+"_align", "tmp/", "--threads", str(args.threads), "--max-seqs", str(subset_size), "--mask", "0", "--format-output", "query,target,raw", "--comp-bias-corr", "0"]
    subprocess.run(run_mmseqs, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    print("aligned total sequence set to subset using MMSeqs2, proceeding to matrix construction")

elif args.aligner == "blastp":
    run_blast_makedb = ["makeblastdb", "-in", align_subset, "-dbtype", "prot"]
    run_blastp = ["blastp", "-query", str(fasta), "-db", align_subset, "-num_threads", str(args.threads), "-evalue", "0.001", "-outfmt", "6", "qseqid", "sseqid", "score", "-max_target_seqs", subset_size, "-max_hsps", "1", "-out", basename+"_align"]
    subprocess.run(run_blast_makedb, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run(run_blastp, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    print("aligned total sequence set to subset using blastp, proceeding to matrix construction")

# generate matrix
align_df = pd.read_csv(basename+"_align", sep="\t", header=None)
ASM_df = align_df.pivot(index=0, columns=1, values=2).fillna(0).astype(int)
ASM_df.to_csv(basename+"_ASM.tsv", sep="\t")
print("written the Alignment score matrix to {}, proceeding to embedding".format(basename+"_ASM.tsv"))
ASM_matrix = ASM_df.values

# tsne
if args.tsne:
    import openTSNE
    from sklearn.cluster import DBSCAN

    # key tSNE parameters:
    # perplexity - https://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf
    # kernel degree of freedom (dof) - https://link.springer.com/chapter/10.1007/978-3-030-46150-8_8
    # exaggeration - https://arxiv.org/abs/2007.08902
    # learning rate - https://www.nature.com/articles/s41467-019-13056-x
    ASM_affinities = openTSNE.affinity.Multiscale(ASM_matrix, perplexities=[20,args.tsne_perplexity], metric="cosine", n_jobs=args.threads)
    ASM_PCA_init = openTSNE.initialization.pca(ASM_matrix)
    ASM_tsne_embedding = openTSNE.TSNEEmbedding(ASM_PCA_init, ASM_affinities, n_jobs=args.threads, dof=1, learning_rate="auto")

    # first 500 iterations of tSNE embedding with exaggeration
    ASM_tsne_embedding.optimize(n_iter=args.tsne_iterations, exaggeration=args.tsne_exaggeration, momentum=0.5, inplace=True, n_jobs=args.threads)

    # cluster on early embedding and save plot for inspection
    ASM_tsne_clustering = DBSCAN(eps=1, min_samples=10).fit(ASM_tsne_embedding)

    tsne_embed_early = pd.DataFrame(data=ASM_tsne_embedding)
    tsne_embed_early.columns = ['tsne1', 'tsne2']
    tsne_embed_early['prot_ID'] = ASM_df.index
    tsne_embed_early['locus_nr'] = tsne_embed_early['prot_ID'].astype("string").str.rsplit("_", n=1).str[1]
    tsne_embed_early['genome_ID'] = tsne_embed_early['prot_ID'].astype("string").str.rsplit("_", n=1).str[0]
    tsne_embed_early['cluster'] = ASM_tsne_clustering.labels_
    tsne_embed_early = tsne_embed_early[['prot_ID', 'locus_nr', 'genome_ID', 'tsne1', 'tsne2', 'cluster']]

    if args.metadata_protein:
        protein_metadata_df = pd.read_csv(protein_metadata, sep="\t")
        tsne_embed_early = pd.merge(tsne_embed_early, protein_metadata_df, on="prot_ID")
    if args.metadata_genome:
        genome_metadata_df = pd.read_csv(genome_metadata, sep="\t")
        tsne_embed_early = pd.merge(tsne_embed_early, genome_metadata_df, on="genome_ID")

    tsne_embed_early.to_csv(basename+"_tsne_early_clust.tsv", sep="\t", index=False)

    # finalize embedding for visualization
    ASM_tsne_embedding.optimize(n_iter=args.tsne_iterations, momentum=0.8, inplace=True, n_jobs=args.threads)

    # convert embedding to dataframe, add clustering and utility columns
    tsne_embed_df = pd.DataFrame(data=ASM_tsne_embedding)
    tsne_embed_df.columns = ['tsne1', 'tsne2']
    tsne_embed_df['prot_ID'] = ASM_df.index
    tsne_embed_df['cluster'] = ASM_tsne_clustering.labels_
    tsne_embed_df['locus_nr'] = tsne_embed_df['prot_ID'].astype("string").str.rsplit("_", n=1).str[1]
    tsne_embed_df['genome_ID'] = tsne_embed_df['prot_ID'].astype("string").str.rsplit("_", n=1).str[0]
    tsne_embed_df = tsne_embed_df[['prot_ID', 'locus_nr', 'genome_ID', 'tsne1', 'tsne2', 'cluster']]

    if args.metadata_protein:
        protein_metadata_df = pd.read_csv(protein_metadata, sep="\t")
        tsne_embed_df = pd.merge(tsne_embed_df, protein_metadata_df, on="prot_ID")
    if args.metadata_genome:
        genome_metadata_df = pd.read_csv(genome_metadata, sep="\t")
        tsne_embed_df = pd.merge(tsne_embed_df, genome_metadata_df, on="genome_ID")

    tsne_embed_df.to_csv(basename+"_tsne_embed.tsv", sep="\t", index=False)

# umap enmbedding and clustering
if args.umap:
    NUMBA_NUM_THREADS = args.threads # this has to be set before importing umap/numba
    import umap
    from sklearn.cluster import DBSCAN

    ASM_umap_embedding = umap.UMAP(min_dist=0.5,n_neighbors=200).fit_transform(ASM_matrix)

    ASM_tsne_clustering = DBSCAN(eps=1, min_samples=10).fit(ASM_umap_embedding)

    umap_embed_df = pd.DataFrame(data=ASM_umap_embedding)
    umap_embed_df.columns = ['umap1', 'umap2']
    umap_embed_df['prot_ID'] = ASM_df.index
    umap_embed_df['cluster'] = ASM_umap_clustering.labels_
    umap_embed_df = umap_embed_df[['prot_ID', 'tsne1', 'tsne2', 'cluster']]
    umap_embed_df.to_csv(basename+"_umap_embed.tsv", sep="\t", index=False)

