#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import re
import csv
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
    "--aligner", help="select the aligner [bowtie, bowtie_htseq]", choices=["bowtie", "bowtie_htseq"]
)
parser.add_argument(
    "--index",
    help="specify the index file",
    default="data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index",
)
parser.add_argument(
    "--sam",
    help="specify the sam file",
    default="data/eg2.sam",
)
parser.add_argument(
    "--gff",
    help="specify the gff file",
    default="data/hsa.gff3",
)
parser.add_argument(
    "--count",
    help="specify the htseq count file",
    default="data/count.tsv",
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


def parse_bowtie2(out):
    ignored = 0
    for line in out.stderr.decode("utf-8").split("\n"):
        m1 = re.search(r"\(([0-9.]*)%\) aligned 0 times", line)
        if m1:
            not_aligned = float(m1.group(1))
        m2 = re.search(r"([0-9.]*)% overall alignment rate", line)
        if m2:
            aligned = float(m2.group(1))
        if "Warning: skipping read" in line:
            ignored += 1
    return aligned, not_aligned, ignored


def parse_htseq(count_file):
    ignored = 0
    with open(count_file, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        values = [int(row[1]) for row in reader]
        aligned = sum(values[:-5])
        bad_aligned = sum(values[-5:-2])
        not_aligned = values[-2]
        not_unique = values[-1] # not needed for total
        total = aligned + bad_aligned + not_aligned
    return 100*aligned/total, 100*bad_aligned/total, 100*not_aligned/total


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
            args=(queues1[i], out_queue, trimFirst, trimLast, trimTo, matcher),
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
    
    print(f"Matching time: {time_match}")

    # Align results
    if args.aligner:
        print("Start alignment")

    if args.aligner in ("bowtie", "bowtie_htseq"):
        cmd = "bowtie2 -x {} {} -S {} -p {}" "".format(
            args.index, args.out_file, args.sam, args.workers
        )
        t_start = time.perf_counter() * 1000
        out = run(cmd, check=True, stdout=PIPE, stderr=PIPE, shell=True)
        t_end = time.perf_counter() * 1000

        time_align = math.floor(t_end - t_start)
        aligned, not_aligned, ignored = parse_bowtie2(out)
        ignored = ignored / count * 100.0
        
        # Print results
        print()
        print("Bowtie alignment")
        print(f"- Alignment time: {time_align}")
        print(f"- Aligned: {aligned:2.2f}%")
        print(f"- Ignored: {ignored:2.2f}%")
        print(f"- Not aligned: {not_aligned:2.2f}%")

    if args.aligner == "bowtie_htseq":
        cmd = "python3 -m HTSeq.scripts.count -t miRNA -i Name {} {} > {}".format(
            args.sam, args.gff, args.count
        )
        t_start = time.perf_counter() * 1000
        out = run(cmd, check=True, stdout=PIPE, stderr=PIPE, shell=True)
        t_end = time.perf_counter() * 1000

        time_align = math.floor(t_end - t_start)
        aligned, bad_aligned, not_aligned = parse_htseq(args.count)
        ignored = ignored / count * 100.0

        # Print results
        print()
        print("HTSeq alignment")
        print(f"- Alignment time: {time_align}")
        print(f"- Aligned: {aligned:2.2f}%")
        print(f"- Badly aligned: {bad_aligned:2.2f}%")
        print(f"- Not aligned: {not_aligned:2.2f}%")
        print()

if __name__ == "__main__":
    main()
