# Anchorage

The Anchorage algorithm is a method to assemble synthetic long reads (SLR) with anchors. Anchors are known short sequences ligated to both ends of a sequenced molecule. It is applicable to LoopSeq, LoopSeq Solo, and other sequencing technologies that have known sequences at the endpoints of a molecule.

### Brief description of algorithm

Modern sequencing technologies often ligate short synthetic adapters with known sequences to both ends of a captured molecule. Those adapters often play functional roles in sequencing, for example, the capture of target molecules, template switch, and amplification in PCR. Those adaptors are often sequenced together with the target molecules (and then trimmed away by tailored algorithms, e.g. trimmomatic). Those adaptors are known and often dissimilar from the sequenced target, so they can serve to indicate the termini of molecules in sequence assembly. Additional such oligo-nucleotides, "anchors", can be synthesized and ligated to the target molecule, for example, in LoopSeq. They are extremely useful for assembly as they mark the two endpoints of molecules.

Anchorage starts with a kmer-based approach for precise estimation of molecule lengths. It then formulates the assembly problem as finding an optimal path that connects the two nodes determined by anchors in the underlying compact de Bruijn graph. The optimality is defined as maximizing the weight of the smallest node while matching the estimated sequence length. Anchorage uses a modified dynamic programming algorithm to efficiently find the optimal path. Anchorage is robust to sequencing artifacts and errors, particularly under very high sequencing depths.

More detailed description of the algorithm is available at https://doi.org/10.4230/LIPIcs.WABI.2024.22.



# Installation

Anchorage can be installed by cloning this repository.

```sh
git clone https://github.com/Shao-Group/anchorage.git
```

Anchorage is implemented in python3. Anchoarge requires [kmc3](https://github.com/refresh-bio/KMC), [seqkit](https://github.com/shenwei356/seqkit), and [spades](https://ablab.github.io/spades/index.html), for kmer counting, sequence statistics, and graph construction. It is also dependent on python package [networkx](https://networkx.org/) and [numpy](https://numpy.org/). If preferred, users may install the dependencies separately and manually, but all executables and libraries should be callable from `$PATH`.

Otherwise, an easier way to install and use those dependencies is via conda with the `.yml` file.
```sh
conda env create -f src/python/environment.yml
conda activate anchorage
```

# Usage

The input of Anchorage should be de-multiplexed, anchor-preserved reads. In other words, the input reads should (1) be from the same cell/sample and molecular target, i.e. having the same cell barcode (CB) and unique molecular identifier (UMI), if any; and (2) have anchor sequences in (some of) the reads. Quality control prior to assembly is permitted as long as the anchor sequences are not trimmed away. Prior knowledge sequence of start and end anchor sequences is required. Users can consult [loop-core](https://github.com/Elembio/loop-core) for LoopSeq's anchor sequences. For some other sequencing technologies, TSO and oligo-dT/polyA adaptors (or other tech-specific sequences) might be able to serve as anchors, since those adaptors are also known and ligated to the termini of a sequenced target. Anchorage was extensively evaluated on LoopSeq Solo assemblies with internal gold standards.

### Quick start:

For MacOS, run `ulimit -n 2048` before running anchorage. This is needed by kmc. See https://github.com/refresh-bio/KMC/issues/140. Linux system does not need this step.

Run `anchorage`  to assemble a synthetic long read:

```bash
python anchorage.py -s1 <start-anchor> -s2 <end-anchor> -r1 <read1.fa> -r2 <read2.fa> [options]
```

The default output is in `anchorage_contig.fa`. If a prefix is specified by argument `-o`/`--output_prefix`, the output will be `prefix.fa` 



### Required:

`-s1`, `--anchor_start`: str  
	&nbsp;&nbsp;&nbsp;&nbsp;Sequence of start anchor.

`-s2`, `--anchor_end`: str  
	&nbsp;&nbsp;&nbsp;&nbsp;Sequence of end anchor.

`-r1`, `--read1_fq`: str  
	&nbsp;&nbsp;&nbsp;&nbsp;Read1 file of fastq/fasta.

`-r2`, `--read2_fq`: str  
	&nbsp;&nbsp;&nbsp;&nbsp;Read2 file of fastq/fasta.



### Options:

`-o`, `--output_prefix`: str  
	&nbsp;&nbsp;&nbsp;&nbsp;Output file prefix, default: `anchorage_contig`.    
	&nbsp;&nbsp;&nbsp;&nbsp;The assembled contig will be in `prefix.fa`.

`-k`: comma-separated-int  
	&nbsp;&nbsp;&nbsp;&nbsp;A series of k-mer sizes, default: 21,33,55,77,99.  
	&nbsp;&nbsp;&nbsp;&nbsp;This is used for the construction of the de Bruijn Graph.

`-t`, `--threads`: int  
	&nbsp;&nbsp;&nbsp;&nbsp;Number of threads for spades to construct the de Bruijn Graph, default: 8.    
	&nbsp;&nbsp;&nbsp;&nbsp;Please be advised that the assembly step of Anchorage runs in a single thread.

`--contig_barcode_len`: int  
	&nbsp;&nbsp;&nbsp;&nbsp;Length of contig barcode for LoopSeq Solo, default: 0.    
	&nbsp;&nbsp;&nbsp;&nbsp;LoopSeq Solo adds "contig barcode" to the 5' of target molecule as an identifier. This barcode will be used as the name of the contig in the output `.fa` file.

`--no-trim_barcode`: no arguments  
	&nbsp;&nbsp;&nbsp;&nbsp;NOT trimming contig barcode from 5' of the assembled contig, default: trim.  
	&nbsp;&nbsp;&nbsp;&nbsp;The "contig barcode" is ligated to the target molecule during library preparation/sequencing. By default, it will be trimmed away after assembly since it's not part of the actual target molecule. When the trim is undesired or such nucleotides are absent, use `--no-trim_barcode` to prevent the trim.

`--max_nm_anchors`: int  
	&nbsp;&nbsp;&nbsp;&nbsp;Maximum number of mismatches/indels in anchors permitted, default: 2.  
	&nbsp;&nbsp;&nbsp;&nbsp;Since sequencing errors may occur in the anchor, Anchorage iteratively searches for anchor nodes by tolerating `x` mismatches+indels in the anchor nodes sequences, when aligning the node sequence with the anchor sequence. Anchorage iteratively increments `x` from 0 to the max tolerance `max_nm_anchors`, until such anchors nodes and whereafter a satisfying contig is found (or teminate with a failure if the max tolerance is reached but no satisfying contig is found). The value of mismatches + indels is calculated using Smith-Waterman algorithm.

`--verbose`: no arguments  
	&nbsp;&nbsp;&nbsp;&nbsp;Verbosely print messages, default: false.
	
`-h`, `--help`: no arguments  
	&nbsp;&nbsp;&nbsp;&nbsp;Show help message and exit.

# Example
We provided an example dataset and its ground truth sequence in the `example` directory.
The reads are simulated using NG_042068.1.fa as a reference and added CGCAGAGTACAT/TTGGAGTTAAAG as the start/end anchors.

Users can assemble this SLR by using:
```sh
ulimit -n 2048  # for MacOS. Not needed for Linux

python anchorage.py -o anchorage-example -s1 CGCAGAGTACAT -s2 TTGGAGTTAAAG -r1 example/sample_01_1.fasta -r2 example/sample_01_2.fasta  --no-trim_barcode --contig_barcode_len 0 
```
The assembled fasta sequence will be in `anchorage-example.fa`

# License

Anchorage is freely available under a [BSD-3-Clause license](./LICENSE). This repository provides only the anchorage assembly algorithm and parts of its code are forked and modified from [loop-core](https://github.com/Elembio/loop-core) ([LICENSE](https://github.com/Elembio/loop-core/blob/main/LICENSE)). For the complete LoopSeq data analysis pipeline, including trimming, demultiplex, read grouping by UMI, assembly, and aggregation of results, please refer to [loop-core](https://github.com/Elembio/loop-core). 

# Citation

If you find this repository useful, please consider citing:

>Xiaofei Carl Zang, Xiang Li, Kyle Metcalfe, Tuval Ben-Yehezkel, Ryan Kelley, and Mingfu Shao. Anchorage Accurately Assembles Anchor-Flanked Synthetic Long Reads. In 24th International Workshop on Algorithms in Bioinformatics (WABI 2024). Leibniz International Proceedings in Informatics (LIPIcs), Volume 312, pp. 22:1-22:17, Schloss Dagstuhl – Leibniz-Zentrum für Informatik (2024)
>https://doi.org/10.4230/LIPIcs.WABI.2024.22
