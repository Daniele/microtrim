#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import time
from queue import Empty
from multiprocessing import Process, Queue

from fastqandfurious import fastqandfurious as ff
from fastqandfurious._fastqandfurious import entrypos as entrypos_c

from microlib.matcher import adaptergen, adaptergen_faster, leven, ndleven, ssw
from microlib.aligner import bowtie, htseq

MATCHER_BUILDER = {
    "adagen": adaptergen.build,
    "adagen-fast": adaptergen_faster.build,
    "leven": leven.build,
    "ndleven": ndleven.build,
    "ssw": ssw.build,
}
EOF = "EOF"

parser = argparse.ArgumentParser()
parser.add_argument(
    "-m",
    "--matcher",
    help="the matcher to use [adagen, adagen-fast, leven, ndleven, ssw]",
    choices=["adagen", "adagen-fast", "leven", "ndleven", "ssw"],
    required=True,
)
parser.add_argument(
    "--aligner",
    help="select the aligner [bowtie, bowtie_htseq]",
    choices=["bowtie", "bowtie_htseq"],
)
parser.add_argument(
    "--index",
    help="specify the index file",
    default="data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index",
)
parser.add_argument("--sam", help="specify the sam file", default="data/eg2.sam")
parser.add_argument("--gff", help="specify the gff file", default="data/hsa.gff3")
parser.add_argument(
    "--count", help="specify the htseq count file", default="data/count.tsv"
)
parser.add_argument(
    "-a", "--adapter", help="adapter to remove", default="TGGAATTCTCGGGTGCCAAGG"
)
parser.add_argument(
    "--trim-first", type=int, help="number of initial bases to trim", default=0
)
parser.add_argument(
    "--trim-last", type=int, help="number of final bases to trim", default=0
)
parser.add_argument(
    "--trim-to", type=int, help="number of final bases to trim", default=23
)
parser.add_argument(
    "--match-only", type=int, help="match only the first X bases of adapter", default=15
)
parser.add_argument(
    "-i", "--in-file", help="input file", default="data/SRR8311267.fastq"
)
parser.add_argument(
    "-o", "--out-file", help="output file", default="data/SRR8311267.trimmed.fastq"
)
# TODO: implement quiet and verbose
# parser.add_argument('-q', '--quiet', help='suppress output',
#                     action='store_true')
# parser.add_argument('-v', '--verbose',
#                     help='print additional information', action='store_true')
parser.add_argument(
    "--max-distance",
    type=float,
    help="maximum string distance (used only in ndleven)",
    default=0.1,
)
parser.add_argument(
    "--stop-after",
    type=int,
    help="stop after 1/X of the string (used only in leven and ndleven)",
    default=2,
)
parser.add_argument("--workers", type=int, help="number of parallel workers", default=4)
parser.add_argument(
    "--chunk", type=int, help="number of chunks send to the workers", default=500
)
parser.add_argument(
    "--debug-limit",
    type=int,
    help="read only the first N read from the input file",
    default=-1,
)


def trim_partition(partition, trimFirst, trimLast, trimTo, match_fun):
    """
    Trim a partition of lines with the passed match fun
    """
    for i, seq in enumerate(partition):
        comment = seq[0].decode("utf-8")
        line = seq[1].decode("utf-8")
        quality = seq[2].decode("utf-8")
        match = match_fun(line)
        tFirst = 0
        tLast = 0
        count = 0

        if trimTo and match:
            lineLen = len(line[:match])
            while lineLen > trimTo:
                if count % 2:
                    tFirst += 1
                else:
                    tLast += 1
                count += 1
                lineLen -= 1

        tFirst = max(tFirst, trimFirst)
        tLast = max(tLast, trimLast)

        if match:
            line = line[tFirst : match - tLast]
            quality = quality[tFirst : match - tLast]
        else:
            line = line[tFirst : tFirst + trimTo]
            quality = quality[tFirst : tFirst + trimTo]

        partition[i] = f"@{comment}\n{line}\n+\n{quality}\n"
    return partition


def worker_fun(q1, q2, trimFirst, trimLast, trimTo, match_fun):
    for p in iter(q1.get, EOF):
        q2.put(trim_partition(p, trimFirst, trimLast, trimTo, match_fun))
    q2.put(EOF)


def collector_fun(outFilePath, queues):
    with open(outFilePath, "w") as outFile:
        count = 0
        eof_count = 0
        while eof_count < len(queues):
            for q in queues:
                try:
                    partition = q.get_nowait()
                    if partition == EOF:
                        eof_count += 1
                        continue
                    for p in partition:
                        outFile.write(p)
                        count += 1
                except Empty:
                    pass
    print(f"Received {count} elements")


def main():
    args = parser.parse_args()
    matcher_name = args.matcher
    adapter = args.adapter
    trimFirst = args.trim_first
    trimLast = args.trim_last
    trimTo = args.trim_to
    inFilePath = args.in_file
    outFilePath = args.out_file
    maxThread = args.workers
    chunk = args.chunk
    debugLimit = args.debug_limit

    print()
    if trimFirst == 0:
        print(f"Not trimming any initial bases")
    else:
        print(f"Trimming the first {trimFirst} bases")
    print(f"Trimming adapter: {adapter}")
    # if version == 2:
    print(f"The matcher '{matcher_name}' is used to find the adapter")
    # print(f'Considering only first {matchOnly} bases of adapter: {adapter[:matchOnly]}')
    # if version < 3:
    #     print(f'Considering {len(adapters)} possible variants of the adapter')
    # else:
    #     print(f'Using Levenshtein-Damerau distance to find adapter variants')
    print(f"Trimming all bases after the adapter (if present)")
    if trimLast == 0:
        print(f"Not trimming any other bases after adapter removal")
    else:
        print(f"Trimming the last {trimLast} bases after adapter removal")
    print(f"Saving to file: {outFilePath}")
    print("Used", f"{maxThread} workers" if maxThread > 0 else "sequential version")
    print()

    # get the matcher function
    matcher_builder = MATCHER_BUILDER[matcher_name]
    matcher = matcher_builder(adapter, args)

    if maxThread > 0:
        # build the parallel topology
        process = [None] * maxThread
        queues1 = [None] * maxThread
        queues2 = [None] * maxThread
        for i in range(maxThread):
            queues1[i] = Queue()
            queues2[i] = Queue()
            out_queue = queues2[i]
            process[i] = Process(
                target=worker_fun,
                args=(queues1[i], out_queue, trimFirst, trimLast, trimTo, matcher),
            )
            process[i].start()
        collector = Process(target=collector_fun, args=(outFilePath, queues2))
        collector.start()

        # start file read
        t_start = time.perf_counter() * 1000
        with open(inFilePath, "r+b") as infile:
            t = 0
            sequence = ff.readfastq_iter(infile, fbufsize=50000, _entrypos=entrypos_c)
            for i, seq in enumerate(sequence):
                p = i % chunk
                if p == 0:
                    partition = [None] * chunk

                if i == debugLimit:
                    break

                partition[p] = seq

                if p == chunk - 1:
                    queues1[t].put(partition)
                    t = (t + 1) % maxThread
            if p < chunk - 1:
                partition = partition[:p]
                queues1[t].put(partition)

        print(f"Sent {i} elements to the workers")
        for q in queues1:
            q.put(EOF)

        print("Wait process")
        for p in process:
            p.join()
        collector.join()
        t_end = time.perf_counter() * 1000
        time_match = math.floor(t_end - t_start)

        print(f"Matching time: {time_match}")
    else:
        # Sequential version
        t_start = time.perf_counter() * 1000
        with open(inFilePath, "r+b") as infile:
            sequence = ff.readfastq_iter(infile, fbufsize=50000, _entrypos=entrypos_c)
            with open(outFilePath, "w") as outFile:
                for i, seq in enumerate(sequence):
                    comment = seq[0].decode("utf-8")
                    line = seq[1].decode("utf-8")
                    quality = seq[2].decode("utf-8")
                    match = matcher(line)
                    tFirst = 0
                    tLast = 0
                    count = 0

                    if trimTo and match:
                        lineLen = len(line[:match])
                        while lineLen > trimTo:
                            if count % 2:
                                tFirst += 1
                            else:
                                tLast += 1
                            count += 1
                            lineLen -= 1

                    tFirst = max(tFirst, trimFirst)
                    tLast = max(tLast, trimLast)

                    if match:
                        line = line[tFirst : match - tLast]
                        quality = quality[tFirst : match - tLast]
                    else:
                        line = line[tFirst : tFirst + trimTo]
                        quality = quality[tFirst : tFirst + trimTo]

                    outFile.write(f"@{comment}\n{line}\n+\n{quality}\n")

        t_end = time.perf_counter() * 1000
        time_match = math.floor(t_end - t_start)
        print(f"Processed {i} elements")
        print(f"Matching time: {time_match}")

    # Align results
    if args.aligner:
        print("Start alignment")

    # Align results
    if args.aligner in ("bowtie", "bowtie_htseq"):
        bowtie(args)

    if args.aligner == "bowtie_htseq":
        htseq(args)


if __name__ == "__main__":
    main()
