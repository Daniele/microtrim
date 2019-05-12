#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

from microlib.aligner import bowtie, htseq

parser = argparse.ArgumentParser()
parser.add_argument(
    "--aligner",
    help="select the aligner [bowtie, bowtie_htseq]",
    choices=["bowtie", "bowtie_htseq"],
    default="bowtie_htseq",
)
parser.add_argument(
    "--index",
    help="specify the index file",
    default="data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index",
)
parser.add_argument("--gff", help="specify the gff file", default="data/hsa.gff3")
parser.add_argument("-o", "--out-name", help="output name", required=True)
# parser.add_argument("--sam", help="specify the sam file", default="data/eg2.sam")
# parser.add_argument(
#     "--count", help="specify the htseq count file", default="data/count.tsv"
# )
parser.add_argument(
    "-i", "--in-file", help="input file", default="data/SRR8311267.trimmed.fastq"
)
parser.add_argument("--workers", type=int, help="number of parallel workers", default=4)
# TODO: implement quiet and verbose
# parser.add_argument('-q', '--quiet', help='suppress output',
#                     action='store_true')
# parser.add_argument('-v', '--verbose',
#                     help='print additional information', action='store_true')


def main():
    args = parser.parse_args()

    print()
    print(f"Input file name: {args.in_file}")
    print(f"Used {args.workers} workers")
    print()

    # Align results
    if args.aligner in ("bowtie", "bowtie_htseq"):
        args.sam = f"{args.out_name}.sam"
        bowtie(args)

    if args.aligner == "bowtie_htseq":
        args.count = f"{args.out_name}.tsv"
        htseq(args)


if __name__ == "__main__":
    main()
