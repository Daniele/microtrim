# -*- coding: utf-8 -*-
import time
import math
import re
import csv
from subprocess import PIPE, run


def parse_bowtie2(out):
    out_str = out.stderr.decode("utf-8")
    m = re.search(r"([0-9]+) reads;", out_str)
    if m:
        total = int(m.group(1))
    m = re.search(r"\(([0-9.]+)%\) aligned 0 times", out_str)
    if m:
        not_aligned = float(m.group(1))
    m = re.search(r"([0-9.]+)% overall alignment rate", out_str)
    if m:
        aligned = float(m.group(1))
    ignored = 0
    for line in out_str.split("\n"):
        if "Warning: skipping read" in line:
            ignored += 1
    return total, aligned, not_aligned, ignored


def parse_htseq(count_file):
    with open(count_file, "r") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        values = [int(row[1]) for row in reader]
        aligned = sum(values[:-5])
        bad_aligned = sum(values[-5:-2])
        not_aligned = values[-2]
        not_unique = values[-1]  # not needed for total
        total = aligned + bad_aligned + not_aligned
    return 100 * aligned / total, 100 * bad_aligned / total, 100 * not_aligned / total


def bowtie(args):
    cmd = "bowtie2 -x {} {} -S {} -p {}" "".format(
        args.index, args.in_file, args.sam, args.workers
    )
    t_start = time.perf_counter() * 1000
    out = run(cmd, check=True, stdout=PIPE, stderr=PIPE, shell=True)
    t_end = time.perf_counter() * 1000

    time_align = math.floor(t_end - t_start)
    count, aligned, not_aligned, ignored = parse_bowtie2(out)
    ignored = ignored / count * 100.0

    # Print results
    print()
    print("Bowtie alignment")
    print(f"- Alignment time: {time_align}")
    print(f"- Aligned: {aligned:2.2f}%")
    print(f"- Ignored: {ignored:2.2f}%")
    print(f"- Not aligned: {not_aligned:2.2f}%")


def htseq(args):
    cmd = "python3 -m HTSeq.scripts.count -t miRNA -i Name {} {} > {}".format(
        args.sam, args.gff, args.count
    )
    t_start = time.perf_counter() * 1000
    run(cmd, check=True, stdout=PIPE, stderr=PIPE, shell=True)
    t_end = time.perf_counter() * 1000

    time_align = math.floor(t_end - t_start)
    aligned, bad_aligned, not_aligned = parse_htseq(args.count)

    # Print results
    print()
    print("HTSeq alignment")
    print(f"- Alignment time: {time_align}")
    print(f"- Aligned: {aligned:2.2f}%")
    print(f"- Badly aligned: {bad_aligned:2.2f}%")
    print(f"- Not aligned: {not_aligned:2.2f}%")
    print()
