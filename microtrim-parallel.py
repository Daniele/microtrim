#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import os
import sys
import mmap
import json
import time
from concurrent.futures import ProcessPoolExecutor
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
maxThread = int(sys.argv[1]) if len(sys.argv) > 1 else None
chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 50000

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
def analyse(partition):
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

    with open(inDirPath + 'SRR8311267.fastq', 'r+b') as infile:
        # m = mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ)
        m = infile
        
        group_lines = 4
        p_size = chunk * group_lines
        f_results= []
        with ProcessPoolExecutor(max_workers=maxThread) as executor:
            data_iter = iter(m.readline, b'')
            # for i in range(0, maxThread):
            stop = False
            while True:
                print('before create')
                partition = [None]*p_size
                try:
                    for i in range(p_size):
                        el = next(data_iter)
                        partition[i] = el
                except StopIteration:
                    stop = True
                    partition = partition[:i]
                print('after create')

                # print(f"send partition of {len(partition)}")
                # print(f'first element is {partition[0]}')
                # print(f'last element is {partition[-1]}')
                print(f'send chunk N {len(f_results)} of size {len(partition)}', end='\n')
                f_results.append(executor.submit(analyse, partition))
                print(f'after send')
                if stop:
                    break
                # p = end
    
    print("send all results")

    with open(outFilePath, 'w') as outFile:
        for f in f_results:
            for r in f.result():
                outFile.write(r)
        # m = infile
        # def group_lines(it):
        #     for i, l in enumerate(it):
        #         if i % 4 == 0:
        #             lines = []
        #         lines.append(l.decode('utf-8'))
        #         if i % 4 == 3:
        #             yield lines
        
        # def limit(it):
        #     for i, el in enumerate(it):
        #         sys.stdout.write(f'\r{i}')
        #         # if i == 170000000: # DEBUG
        #         #     break
        #         yield el
        
        # def chunk(it):
        #     size = len(m)//4
        #     p = 0
        #     for i in range(0, maxThread):
        #         yield [next(it) for i in range(p, p+size)]
        #         p += size
        # for lines in executor.map(analyse, limit(group_lines(iter(m.readline, b""))), chunksize=1000):
        # with ProcessPoolExecutor(max_workers=maxThread) as executor:
        #     for lines in executor.map(analyse, chunk(group_lines(iter(m.readline, b""))), chunksize=1):
        #         count += 1
        #         if lines:
        #             matches += 1
        #             for line in lines:
        #                 result += line

        # print(f'\nTrimmed {count} lines, found {matches} matches\n')

        # TODO - Fix Levenshtein distance matching
        """
        else:
            for j, char in enumerate(line[10:-len(adapter)-1]):
                possibleMatch = line[j:j+len(adapter)-1]
                #l = leven(adapter, possibleMatch)
                #dl = dleven(adapter, possibleMatch)
                ndl = ndleven(adapter, possibleMatch)
                if ndl < .1:
                    matches += 1
                    continue
                    print(f'{i}. {adapter}')
                    print(f'{i}. {possibleMatch} • {l} • {dl} • {ndl}')
                    print()
        sys.stdout.write(f'\r{matches}\r')
        """
        # Used to print matches
        # if matches >= 5000:
        # sys.exit()
        # time.sleep(1)
        #matches = re.findall('', line)
        # for match in matches:
        #    adapterMatch = line[index:index+len(adapter)]
        #    print(adapterMatch)
        #    print(leven(adapter, adapterMatch))
        #    print()
        #    time.sleep(1)

    # with open(outFilePath, 'w') as outFile:
    #     outFile.write(result)


if __name__ == "__main__":
    main()
