#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import re
import time
from subprocess import PIPE, run
from multiprocessing import Process, Queue

from fastqandfurious import fastqandfurious as ff

from matcher import adaptergen, adaptergen_faster, leven, ndleven, ssw

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
    "--aligner", help="select the aligner [bowtie, htsec]", choices=["bowtie", "htsec"]
)
parser.add_argument(
    "--index",
    help="Specified the index file",
    default="data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index",
)
parser.add_argument(
    "-a", "--adapter", help="adapter to remove", default="TGGAATTCTCGGGTGCCAAGG"
)
parser.add_argument(
    "--trim-first", type=int, help="number of initial bases to trim", default=4
)
parser.add_argument(
    "--trim-last", type=int, help="number of final bases to trim", default=4
)
# TODO: implement smart trimming
# parser.add_argument('--trim-to', type=int,
#                     help='trim to a specific number of bases', default=28)
parser.add_argument(
    "--match-only", help="match only the first X bases of adapter", default=15
)
parser.add_argument(
    "-i", "--in-file", help="input file", default="./data/SRR8311267.fastq"
)
parser.add_argument(
    "-o", "--out-file", help="output file", default="./data/SRR8311267.trimmed.fastq"
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
parser.add_argument(
    "--single-queue",
    help="use a single output queue instead on one per worker queues",
    action="store_true",
    default=False,
)


def process_queue(q, f):
    """
    Utility function to iterate over the a queue
    """
    partition = q.get()
    while partition != EOF:
        f(partition)
        partition = q.get()


def trim_partiton(partition, trimFirst, trimLast, match_fun):
    """
    Trim a partition of lines with the passed match fun
    """
    for i, seq in enumerate(partition):
        comment = seq[0].decode("utf-8")
        line = seq[1].decode("utf-8")
        quality = seq[2].decode("utf-8")
        match = match_fun(line)

        if match:
            line = line[trimFirst : match - trimLast]
            quality = quality[trimFirst : match - trimLast]

        partition[i] = f"@{comment}\n{line}\n+\n{quality}\n"
    return partition


def worker_fun(q1, q2, trimFirst, trimLast, match_fun):
    for p in iter(q1.get, EOF):
        q2.put(trim_partiton(p, trimFirst, trimLast, match_fun))
    q2.put(EOF)


def parse_bowtie2(out):
    ignored = 0
    for line in out.stderr.decode("utf-8").split("\n"):
        m1 = re.search(r"\(([0-9.]*)%\) aligned 0 times", line)
        if m1:
            no_aligned = float(m1.group(1))
        m2 = re.search(r"([0-9.]*)% overall alignment rate", line)
        if m2:
            aligned = float(m2.group(1))
        if "Warning: skipping read" in line:
            ignored += 1
    return aligned, no_aligned, ignored


def main():
    args = parser.parse_args()
    matcher_name = args.matcher
    adapter = args.adapter
    trimFirst = args.trim_first
    trimLast = args.trim_last
    inFilePath = args.in_file
    outFilePath = args.out_file
    maxThread = args.workers
    chunk = args.chunk
    debugLimit = args.debug_limit
    singleQueue = args.single_queue

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
    # TODO is maxThreads is 0 use the sequential version of the code
    print(f"Used {maxThread} workers")
    print()

    # get the matcher function
    matcher_builder = MATCHER_BUILDER[matcher_name]
    matcher = matcher_builder(adapter, args)

    # build the parallel topology
    process = [None] * maxThread
    queues1 = [None] * maxThread
    if singleQueue:
        out_queue = Queue()
    else:
        queues2 = [None] * maxThread
    for i in range(maxThread):
        queues1[i] = Queue()
        if not singleQueue:
            queues2[i] = Queue()
            out_queue = queues2[i]
        process[i] = Process(
            target=worker_fun,
            args=(queues1[i], out_queue, trimFirst, trimLast, matcher),
        )
        process[i].start()

    # start file read
    t_start = time.perf_counter() * 1000
    with open(inFilePath, "r+b") as infile:
        t = 0
        # TODO: find the optimal size of the buffer "fbufsize" and "chunk"
        for i, seq in enumerate(ff.readfastq_iter(infile, fbufsize=20_000_000)):
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

    print("Wait results")
    with open(outFilePath, "w") as outFile:
        count = 0

        def write_partition(partition):
            nonlocal count
            for p in partition:
                outFile.write(p)
                count += 1

        if singleQueue:
            for p in range(maxThread):
                process_queue(out_queue, write_partition)
        else:
            for q in queues2:
                process_queue(q, write_partition)

    print(f"Received {count} elements")
    print("Wait process")
    for p in process:
        p.join()
    t_end = time.perf_counter() * 1000
    time_match = math.floor(t_end - t_start)

    # Align results
    if args.aligner:
        print("Start alignment")

    if args.aligner == "bowtie":
        cmd = "bowtie2 -x {} {} -S /tmp/eg2.sam -p {}" "".format(
            args.index, args.out_file, args.workers
        )
        t_start = time.perf_counter() * 1000
        out = run(cmd, check=True, stdout=PIPE, stderr=PIPE, shell=True)
        t_end = time.perf_counter() * 1000

        time_align = math.floor(t_end - t_start)
        alined, not_aligned, ignored = parse_bowtie2(out)
        ignored = ignored / count * 100.0

    if args.aligner == "htsec":
        time_align = 0
        alined = 0
        ignored = 0
        not_aligned = 0

    # Print results
    print()
    print(f"time_match:{time_match}")
    if args.aligner:
        print(f"time_align:{time_align}")
        print(f"alined:{alined:2.2f}")
        print(f"ignored:{ignored:2.2f}")
        print(f"not_aligned:{not_aligned:2.2f}")


if __name__ == "__main__":
    main()
