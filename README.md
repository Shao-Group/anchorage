# Anchorage

The Anchorage algorithm is a method to assemble synthetic long reads (SLR) with anchors, that are known short sequences ligated to both ends ligated of a sequenced molecule. It is applicable to LoopSeq and LoopSeq Solo. 

Compared to other state-of-art assemblers, Anchorage efficiently utilize the following information to perform assembly:

1. strand-specificity and monoclonality of LoopSeq data
2. contig anchors to mark start and end of the contig
3. accurate estimation of barcoded contig length and sequencing coverage from raw reads

# License

Anchorage is available under a [BSD-3-Clause license](./LICENSE). This repository provides only the anchorage assembly algorithm and parts of its code are forked from [loop-core](https://github.com/Elembio/loop-core) (see [loop-core's LICENSE](https://github.com/Elembio/loop-core/blob/main/LICENSE)). For the complete LoopSeq data analysis pipeline, including trimming, demultiplex, read grouping by UMI, assembly, aggregation of results, please refer to [loop-core](https://github.com/Elembio/loop-core). 

# Installation


## Install manually
Anchoarge is implemented in python3. Anchoarge requires [kmc3](https://github.com/refresh-bio/KMC) and [seqkit](https://github.com/shenwei356/seqkit) and is dependent on python package [networkx](https://networkx.org/) and [numpy](https://numpy.org/). 
If desired, users may install the dependencies separately and manually, but all executables and libraries should be callable from `$PATH`.

Otherwise, an easier way to install and use those dependencies is via conda with the `.yml` file .
```
conda env create -f src/python/environment.yml
conda activate anchorage
```

# Usage
```
usage: anchorage.py -s1 -s2 -r1 -r2 [options]

Required arguments:
  -s1, --anchor_start
  -s2, --anchor_end
  -r1, --read1_fq
  -r2, --read2_fq

optional arguments:
  -h, --help                 show this help message and exit
  -o, --output_prefix        output file prefix, default: anchor_guide_contig
  -k                         a series of k-mer sizes, default 21,33,55,77,99.
  -t, --threads              number of threads for spades
  --contig_barcode_len       length of contig barcode (for LoopSeq Solo)
  --no-trim_barcode          NOT trimming contig barcode from contig
  --max_nm_anchors           maximum number of mismatches/indels in anchors permitted
  --verbose                  
```

## Additional information for MacOS
run `ulimit -n 2048` before running anchorage. This is needed by kmc. See https://github.com/refresh-bio/KMC/issues/140.
