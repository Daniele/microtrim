#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import os
import sys
import mmap
import json
import time
import math
import argparse
import Levenshtein
from multiprocessing import Process, Queue
from pyxdameraulevenshtein import damerau_levenshtein_distance as dleven
from pyxdameraulevenshtein import normalized_damerau_levenshtein_distance as ndleven
from fastqandfurious import fastqandfurious as ff

parser = argparse.ArgumentParser()
parser.add_argument('-V', '--version', help='algorithm version', default = 3)
parser.add_argument('-a', '--adapter', help='adapter to remove', default = 'TGGAATTCTCGGGTGCCAAGG')
parser.add_argument('-f', '--trim_first', help='number of initial bases to trim', default = 4)
parser.add_argument('-l', '--trim_last', help='number of final bases to trim', default = 4)
parser.add_argument('-T', '--trim_to', help='trim to a specific number of bases', default = 28)
parser.add_argument('-m', '--match_only', help='match only the first X bases of adapter', default = 15)
parser.add_argument('-i', '--in_file', help='input file', default = './data/SRR8311267.fastq')
parser.add_argument('-o', '--out_file', help='output file', default = './trimmed.fq')
parser.add_argument('-q', '--quiet', help='suppress output', action = 'store_true')
parser.add_argument('-v', '--verbose', help='print additional information', action = 'store_true')
parser.add_argument('-d', '--max_distance', help='maximum string distance', default = .1)
parser.add_argument('-s', '--stop_after', help='stop after 1/X of the string', default = 2)
parser.add_argument('-t', '--max_threads', help='maximum threads', default = 4)
parser.add_argument('-c', '--chunks', help='number of chunks', default = 100)
args = parser.parse_args()

version = args.version
adapter = args.adapter
trimFirst = int(args.trim_first)
trimLast = int(args.trim_last)
trimTo = int(args.trim_to) # TODO implement smart trimming
matchOnly = args.match_only
inFilePath = args.in_file
outFilePath = args.out_file
maxThread = int(args.max_threads)
chunk = int(args.chunks)
maxDistance = float(args.max_distance)
stopAfter = int(args.stop_after)
verbose = args.verbose
quiet = args.quiet

abc = ('A', 'C', 'G', 'T')
adapters = set()
adLength = len(adapter)
cutAdapter = adapter[:matchOnly][::-1]
verbose = False


def addNewAdapterToSet(ad, adSet):
    adSet.add(ad)
    if verbose:
        print(f'Adding {ad}')
        time.sleep(.1)
    return adSet

# V2 (~1000 adapter variants - slower)


def makeAdaptersV2(adapters):
    adapters.add(cutAdapter)
    if verbose:
        print(f'Adding {cutAdapter}')
        time.sleep(.1)

    if verbose:
        print('\nSubstitutions\n')
    for i, x in enumerate(cutAdapter):
        # Adapter with wrong substitution of 1 base
        for j in abc:
            ad = cutAdapter[:i] + j + cutAdapter[i+1:]
            adapters = addNewAdapterToSet(ad[:matchOnly], adapters)

    if verbose:
        print('\nAdditions (1)\n')
    newAdapters1 = set()
    for adapter in adapters:
        for i, x in enumerate(adapter):
            # Adapter with wrong addition of 1 base
            for l in abc:
                ad = adapter[:i] + l + adapter[i:]
                newAdapters1 = addNewAdapterToSet(
                    ad[:matchOnly], newAdapters1)

    if verbose:
        print('\nRemovals (1)\n')
    newAdapters2 = set()
    for adapter in adapters:
        for i, x in enumerate(adapter):
            # Adapter with wrong removal of 1 base
            ad = adapter[:i] + adapter[i+1:]
            newAdapters2 = addNewAdapterToSet(ad[:matchOnly], newAdapters2)

    adapters = adapters.union(newAdapters1).union(newAdapters2)

    return adapters

# V1 (~100 adapter variants - faster)


def makeAdaptersV1(adapters):
    adapters.add(adapter[:-8])
    for i, x in enumerate(adapter):
        for j in abc:
            ad = adapter[:i] + j + adapter[i+1:]
            adapters.add(ad[:-8])
            for l in abc:
                ad = adapter[:i] + l + adapter[i:]
                adapters.add(ad[:-8])
            ad = adapter[:i] + adapter[i+1:]
            adapters.add(ad[:-8])
            ad = adapter[:i] + adapter[i+2:]
            adapters.add(ad[:-8])
    return adapters

def leven(s1, s2):
    return Levenshtein.ratio(s1, s2)

def matchAdapters(line, adapters):
    for adapter in adapters:
        if adapter in line:
            return line.find(adapter)
    return None

def matchLeven(line, adapter):
    for j, char in enumerate(line[:math.floor(len(line)/stopAfter)]):
        possibleMatch = line[j:j+len(adapter)]
        #l = leven(adapter, possibleMatch)
        #dl = dleven(adapter, possibleMatch)
        ndl = ndleven(adapter, possibleMatch)
        if ndl <= maxDistance:
            return -(j+len(adapter))
    return None

# def group_lines(it):
#     for i, l in enumerate(it):
#         if i % 4 == 0:
#             lines = []
#         lines.append(l.decode('utf-8'))
#         if i % 4 == 3:
#             yield lines

if version == 1:
    adapters = sorted(makeAdaptersV1(set()))
elif version == 2:
    adapters = sorted(makeAdaptersV2(set()))

def analysis(partition):
    for i, seq in enumerate(partition):
        comment = seq[0].decode('utf-8')
        line = seq[1].decode('utf-8')
        quality = seq[2].decode('utf-8')
        match = matchLeven(line[::-1], cutAdapter)
        
        if match:
            line = line[trimFirst:match-trimLast]
            quality = quality[trimFirst:match-trimLast]

        partition[i] = f'@{comment}\n{line}\n+\n{quality}\n'
    return partition

def worker_fun2(q1, q2):
    print("process start!", os.getpid())
    partition = q1.get()
    while partition != "EOF":
        res = 1
        for i in range(1000000):
            res += math.sin(res)
        # print (f"{os.getpid()} after sin", res)
        partition = q1.get()
    q2.put(res)
    q2.put("EOF")
    print("quit")

def worker_fun(q1, q2):
    print("process start!", os.getpid())
    partition = q1.get()
    while partition != "EOF":
        partition = analysis(partition)
        q2.put(partition)
        partition = q1.get()
    q2.put("EOF")
    print("quit")

def main():
    print()
    if trimFirst == 0:
        print(f'Not trimming any initial bases')
    else:
        print(f'Trimming the first {trimFirst} bases')
    print(f'Trimming adapter: {adapter}')
    if version == 2:
        print(
            f'Considering only first {matchOnly} bases of adapter: {adapter[:matchOnly]}')
    if version < 3:
        print(f'Considering {len(adapters)} possible variants of the adapter')
    else:
        print(f'Using Levenshtein-Damerau distance to find adapter variants')
    print(f'Trimming all bases after the adapter (if present)')
    if trimLast == 0:
        print(f'Not trimming any other bases after adapter removal')
    else:
        print(f'Trimming the last {trimLast} bases after adapter removal')
    print(f'Saving to file: {outFilePath}')
    print(f'Used {maxThread} worker threads')
    print()

    # count = 0
    # matches = 0

    # def file_size(infile):
    #     m = mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ)
    #     n=0
    #     for i,_ in enumerate(iter(m.readline, b'')):
    #         n+=1
    #     return n

    process = [None] * maxThread
    queues1 = [None] * maxThread
    queues2 = [None] * maxThread
    for i in range(maxThread):
        queues1[i] = Queue()
        queues2[i] = Queue()
        process[i] = Process(target=worker_fun, args=(queues1[i], queues2[i]))
        process[i].start()

    with open(inFilePath, 'r+b') as infile:
        #m = mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ)
        # m = infile

        n_lines = 4
        size = chunk * n_lines
        count = 0
        t = 0

        for i, seq in enumerate(ff.readfastq_iter(infile, 20000000, ff.entryfunc)):
            # if i == 1000*4:
            #     break
            p = count % size
            if p == 0:
                partition = [None] * size

            partition[p] = seq

            if p == size - 1:
                # print("put element")
                # print(f"send to {t}")
                queues1[t].put(partition)
                t = (t + 1) % maxThread
            count+=1

    for q in queues1:
        q.put("EOF")

    print("wait results")
    with open(outFilePath, 'w') as outFile:
        p_end = maxThread
        for q in queues2:
            lines = q.get()
            while lines != 'EOF':
                outFile.write(''.join(lines))
                lines = q.get()
        print("received", lines)

    print("wait process")
    for p in process:
        p.join()


if __name__ == "__main__":
    main()
