#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import os
import sys
import mmap
import json
import time
import math
from multiprocessing import Process, Queue
import Levenshtein
from pyxdameraulevenshtein import damerau_levenshtein_distance as dleven
from pyxdameraulevenshtein import normalized_damerau_levenshtein_distance as ndleven
from itertools import product

# TODO - Make options configurable through argparse
version = 1
inDirPath = './data/'
adapter = 'TGGAATTCTCGGGTGCCAAGG'
matches = 0
trimFirst = 4
ignoreAfter = 10
outFilePath = './data/SRR8311267.trimmed-parallel.fq'
maxThread = int(sys.argv[1]) if len(sys.argv) > 1 else 4
chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 100

abc = ('A', 'C', 'G', 'T')
adapters = set()
adLength = len(adapter)
cutAdapter = adapter[:ignoreAfter]
verbose = False


def addNewAdapterToSet(ad, adSet):
    adSet.add(ad)
    if verbose:
        print(f'Adding {ad}')
        time.sleep(.05)
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
            adapters = addNewAdapterToSet(ad[:ignoreAfter], adapters)

    if verbose:
        print('\nAdditions (1)\n')
    newAdapters1 = set()
    for adapter in adapters:
        for i, x in enumerate(adapter):
            # Adapter with wrong addition of 1 base
            for l in abc:
                ad = adapter[:i] + l + adapter[i:]
                newAdapters1 = addNewAdapterToSet(
                    ad[:ignoreAfter], newAdapters1)

    if verbose:
        print('\nRemovals (1)\n')
    newAdapters2 = set()
    for adapter in adapters:
        for i, x in enumerate(adapter):
            # Adapter with wrong removal of 1 base
            ad = adapter[:i] + adapter[i+1:]
            newAdapters2 = addNewAdapterToSet(ad[:ignoreAfter], newAdapters2)

    adapters = adapters.union(newAdapters1).union(newAdapters2)

    print(len(adapters))
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

# def group_lines(it):
#     for i, l in enumerate(it):
#         if i % 4 == 0:
#             lines = []
#         lines.append(l.decode('utf-8'))
#         if i % 4 == 3:
#             yield lines

if version == 1:
    adapters = sorted(makeAdaptersV1(set()))
else:
    adapters = sorted(makeAdaptersV2(set()))
def analysis(partition):
    for i, line in enumerate(partition):
        line = line.decode('utf-8')
        # print(line)
        if i % 4 == 1:
            match = matchAdapters(line, adapters)
            if match:
                line = line[trimFirst:match] + "\n"
        if i % 4 == 3 and match:
            line = line[trimFirst:match] + "\n"
            match = False
        partition[i] = line
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
            f'Considering only first {ignoreAfter} bases of adapter: {adapter[:ignoreAfter]}')
    print(f'Considering {len(adapters)} possible variants of the adapter')
    print(f'Trimming all bases after the adapter (if present)')
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

    with open(inDirPath + 'SRR8311267.fastq', 'r+b') as infile:
        m = mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ)
        # m = infile

        n_lines = 4
        size = chunk * n_lines
        count = 0
        t = 0
        for line in iter(m.readline, b''):
            # if i == 1000*4:
            #     break
            p = count % size
            if p == 0:
                partition = [None] * size

            partition[p] = line
            
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
