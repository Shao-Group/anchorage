#!/usr/bin/env python3

"""
Description:
    The main pipeline of anchorage. It performs contig assembly for LoopSeq SLR.

BSD 3-Clause License

Copyright (c) 2024, Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University
Copyright (c) 2024, Element Biosciences 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
import networkx as nx

from anchorGuidedAssembler import AnchorGuidedAssembler
from contigStatistician import ContigStatistician
from aligner import PseudoAligner
import aligner
from anchorageUtil import RunShellCommand

from sys import argv
import sys
import argparse
import os

def anchorage(args):
    anchorage_assembly(args)
    print("Anchorage completed running!")
    return 0

def anchorage_assembly(args):
    spades_gfa = os.path.join(args.output_prefix, "assembly_graph_with_scaffolds.gfa")
    if not os.path.isfile(spades_gfa):
        execute_spades(args)

    try:
        execute_AGA(args, spades_gfa)
        return 0
    except nx.NetworkXNoPath:
        pass

    try:
        args.max_nm_anchors += 1
        execute_AGA(args, spades_gfa)
        return 0
    except nx.NetworkXNoPath:
        pass

    try:
        gfa = execute_spades_gbuilder(args)
        execute_AGA(args, gfa)
        return 0
    except Exception as e:
        raise e


def execute_spades(args):
    cmd_spades = []
    cmd_spades.append("spades.py")
    cmd_spades.append("-t  {}".format(args.threads))
    cmd_spades.append("--pe-1  0 {}".format(args.read1_fq))
    cmd_spades.append("--pe-2  0 {}".format(args.read2_fq))
    cmd_spades.append("-k  {}".format(args.k))
    cmd_spades.append("-o {}".format(args.output_prefix))
    if(str(args.read1_fq).endswith(".fq") or 
       str(args.read1_fq).endswith(".fq.gz") or 
       str(args.read1_fq).endswith(".fastq") or
       str(args.read1_fq).endswith(".fastq.gz")):
        cmd_spades.append("--careful")
        cmd_spades.append("--phred-offset  33")
    else:
        cmd_spades.append("--only-assembler")
    cmd_spades.append("--sc -t 8")       
    cmd_spades.append("--disable-gzip-output")
    cmd_line = " ".join(cmd_spades)
    msg = ''
    if args.verbose:
        msg = 'Calling spades for gfa construction'
    RunShellCommand(cmd_line, msg)
    return 0

def execute_spades_gbuilder(args):
    #TODO: Graph simplification
    gfa_output = "{}.gbuild.gfa".format(args.output_prefix)
    if os.path.isfile(gfa_output):
        return gfa_output

    gbuild_kmer_size = 99
    data_yml_file = "{}_spades_guild_k{}.yaml".format(args.output_prefix, gbuild_kmer_size)
    data_yml = """
    [
      {{
        orientation: "fr",
        type: "paired-end",
        right reads: [
          "{}"
        ],
        left reads: [
          "{}"
        ]
      }}
    ]
    """.format(os.path.abspath(args.read1_fq), os.path.abspath(args.read2_fq))
    with open(data_yml_file, 'w') as f:
        f.write(data_yml)
    
    cmdgbuild = []
    cmdgbuild.append("spades-gbuilder")
    cmdgbuild.append(data_yml_file)
    cmdgbuild.append(gfa_output) #TODO:
    cmdgbuild.append("-c --gfa")
    cmdgbuild.append("-t  {}".format(args.threads))
    cmdgbuild.append("-k {}".format(gbuild_kmer_size))
    cmdgbuild.append("-tmp-dir {}_tmp_dir_gbuild_k{}".format(args.output_prefix, gbuild_kmer_size))

    cmd_line = " ".join(cmdgbuild)
    msg = ''
    if args.verbose:
        msg = 'Calling spades-gbuilder'
    RunShellCommand(cmd_line, msg)
    return gfa_output



def execute_AGA(args, gfa):
    """
        output file: args.output_prefix + ".fa", default "anchorage_contig.fa"
    """
    AGA = AnchorGuidedAssembler(
                        # args.gfa,
                        gfa,
                        args.anchor_start,
                        args.anchor_end,
                        left_fq=args.read1_fq,
                        right_fq=args.read2_fq,
                        output_contig_file = args.output_prefix,
                        ht_index_1=args.ht_index_1,
                        ht_index_2=args.ht_index_2,
                        barcode_trim = args.trim_barcode,
                        max_nm_anchors = args.max_nm_anchors,
                        barcode_len = args.contig_barcode_len,   
                        do_strategy_low_cov= False,
                        # algorithm=args.algorithm,
                        algorithm='dp',
                        verbose=args.verbose)
    return 0

def execute_argparse():
    parser = argparse.ArgumentParser()
    requiredarg = parser.add_argument_group('Required arguments')
    recommendarg = parser.add_argument_group('Recommended arguments')
    configarg = parser.add_argument_group('Algorithm configuration arguments')
    
    requiredarg.add_argument('-s1', '--anchor_start', required=True, metavar='\b')
    requiredarg.add_argument('-s2', '--anchor_end', required=True, metavar='\b')
    requiredarg.add_argument('-r1', '--read1_fq', required=True, metavar='\b')
    requiredarg.add_argument('-r2', '--read2_fq', required=True, metavar='\b')

    recommendarg.add_argument('-h1', '--ht_index_1', default="", metavar='\b',
                        help="A variable start anchor following anchor_start immediately in molecule, e.g. HT_index_1 or equivalent")
    recommendarg.add_argument('-h2', '--ht_index_2', default="", metavar='\b',
                        help="A variable end anchor following anchor_end immediately in molecule, e.g. HT_index_2 or equivalent")   
    recommendarg.add_argument('-o', '--output_prefix', default="anchorage_contig", required=False, metavar='\b',
                        help='output file prefix, default: anchor_guide_contig')
    recommendarg.add_argument('-k', default="21,33,55,77,99", required=False, metavar='\b',
                        help='a series of k-mer sizes, default 21,33,55,77,99.')
    # recommendarg.add_argument('-n', default="100000", required=False, metavar='\b',
    #                     help='Default number of reads to be used')
    
    # configarg.add_argument('-a', '--algorithm', default='dp', type=str, required=False, metavar='\b',
    #                     help=r'algorithm for assembly, must be one of ["dp", "all-paths"]')
    
    configarg.add_argument('-t', '--threads', default=8, type=int, required=False, metavar='\b',
                        help='number of threads for spades')
    
    configarg.add_argument('--contig_barcode_len', default=0, type=int, metavar='\b',
                        help='length of contig barcode')
    
    configarg.add_argument('--no-trim_barcode', 
                           dest='trim_barcode', action='store_false', 
                           help='NOT trimming contig barcode from contig ')
    configarg.set_defaults(trim_barcode=True)

    configarg.add_argument('--max_nm_anchors', type=int, default=2, metavar='\b',
                        help='maximum number of anchors permitted')

    # configarg.add_argument('--no-tackle_low_cov',
    #                        dest='tackle_low_cov', action='store_false', 
    #                        help='NOT tackle low coverage samples by selecting long contig with descent coverage')
    # configarg.set_defaults(tackle_low_cov=True)
    
    configarg.add_argument('--verbose',
                           dest='verbose', action='store_true')
    configarg.set_defaults(verbose=False)

    args = parser.parse_args()
    # if args.algorithm not in ["dp", "all-paths"]:
    #     raise ValueError(str(r'--algorithm must be one of ["dp", "all-paths"] but got ')+ str(args.algorithm))
    return args


if __name__ == "__main__":
    args = execute_argparse()
    anchorage(args)
