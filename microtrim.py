#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import os
import sys
import mmap
import json
import time
import argparse
import Levenshtein
from pyxdameraulevenshtein import damerau_levenshtein_distance as dleven
from pyxdameraulevenshtein import normalized_damerau_levenshtein_distance as ndleven
from itertools import product

parser = argparse.ArgumentParser()
parser.add_argument('-V', '--version', help='algorithm version', default = 1)
parser.add_argument('-a', '--adapter', help='adapter to remove', default = 'TGGAATTCTCGGGTGCCAAGG')
parser.add_argument('-f', '--trim_first', help='number of initial bases to trim', default = 4)
parser.add_argument('-l', '--trim_last', help='number of final bases to trim', default = 4)
parser.add_argument('-T', '--trim_to', help='trim to a specific number of bases', default = 28)
parser.add_argument('-m', '--match_only', help='match only the first X bases of adapter', default = 15)
parser.add_argument('-i', '--in_file', help='input file', default = './data/SRR8311267.fastq')
parser.add_argument('-o', '--out_file', help='output file', default = './trimmed.fq')
parser.add_argument('-q', '--quiet', help='suppress output', action = 'store_true')
parser.add_argument('-v', '--verbose', help='print additional information', action = 'store_true')
parser.add_argument('-d', '--max_distance', help='maximum string distance', default = .05)
args = parser.parse_args()

version = args.version
adapter = args.adapter
trimFirst = args.trim_first
trimLast = args.trim_last
trimTo = args.trim_to # TODO implement smart trimming
matchOnly = args.match_only
inFilePath = args.in_file
outFilePath = args.out_file
maxDistance = args.max_distance
verbose = args.verbose
quiet = args.quiet

matches = 0
abc = ('A','C','G','T')
adapters = set()
adLength = len(adapter)
#cutAdapter = adapter[:matchOnly]
cutAdapter = adapter[::-1]

def addNewAdapterToSet(ad, adSet):
    adSet.add(ad)
    if verbose:
        print(f'Adding {ad}')
        time.sleep(maxDistance)
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
                newAdapters1 = addNewAdapterToSet(ad[:matchOnly], newAdapters1)

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
    adapters.add(adapter[:matchOnly])
    for i, x in enumerate(adapter):
        for j in abc:
            ad = adapter[:i] + j + adapter[i+1:]
            adapters.add(ad[:matchOnly])
            for l in abc:
                ad = adapter[:i] + l + adapter[i:]
                adapters.add(ad[:matchOnly])
            ad = adapter[:i] + adapter[i+1:]
            adapters.add(ad[:matchOnly])
            ad = adapter[:i] + adapter[i+2:]
            adapters.add(ad[:matchOnly])
    return adapters

def leven(s1, s2):
    return Levenshtein.ratio(s1, s2)

def matchAdapters(line, adapters):
    for adapter in adapters:
        if adapter in line:
            return line.find(adapter)
    return None

def matchLeven(line, adapter):
    for j, char in enumerate(line[:int(len(line)/2)]):
        possibleMatch = line[j:j+len(adapter)]
        #l = leven(adapter, possibleMatch)
        #dl = dleven(adapter, possibleMatch)
        ndl = ndleven(adapter, possibleMatch)
        if ndl < .05:
            #if verbose:
            #   print(f'{i}. {adapter}')
            #   print(f'{i}. {possibleMatch} • {ndl}')
            #   print()
            return -(j+len(adapter))
    return None

def getSubstrings(string,length):
    subs = set()
    for j in string[10:-length-1]:
        subs.add(line[j:j+length-1])
    return subs

if version == 1:
    adapters = sorted(makeAdaptersV1(set()))
else:
    adapters = sorted(makeAdaptersV2(set()))

if not quiet:
    print()
    if trimFirst == 0:
        print(f'Not trimming any initial bases')
    else:
        print(f'Trimming the first {trimFirst} bases')
    print(f'Trimming adapter: {adapter}')
    if version == 2:
        print(f'Considering only first {matchOnly} bases of adapter: {adapter[:matchOnly]}')
    print(f'Considering {len(adapters)} possible variants of the adapter')
    print(f'Trimming all bases after the adapter (if present)')
    if trimLast == 0:
        print(f'Not trimming any other bases after adapter removal')
    else:
        print(f'Trimming the last {trimLast} bases after adapter removal')
    print(f'Saving to file: {outFilePath}')
    print()

result = ""

with open(inFilePath, 'r+b') as infile:
    m = mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ)
    length = len(adapter)
    finalLength = str(51-length) # should be set based on first line
    match = None
    isRead = False
    count = 0
    for i, line in enumerate(iter(m.readline, b"")):
        if (i == 13):
            print(line[::-1], cutAdapter)
        # if i == 5*4:  # DEBUG
        #     break
        line = line.decode("utf-8")
        #for i in range(0,len(m)):
        #line = str(m.readline())
        #print(line)
        #if adapters.intersection(getSubstrings(line,length)):
        #    matches += 1
        #    sys.stdout.write(f'\r{i} • {matches}\r')
        #continue
        if isRead and (i+3) % 4 == 0:
            count += 1
            #match = matchAdapters(line, adapters)
            match = matchLeven(line[::-1], cutAdapter)
            if match:
                matches += 1
                result += f'{prevLine}\n{line[trimFirst:match-trimLast]}\n'
                if not quiet:
                    sys.stdout.write(f'\r{count} • {matches}\r')
                #if len(line[:match]) < 10:
                #    print(adapter)
                #    print(match)
                #    print(line)
                #    time.sleep(.3)
            elif verbose:
                # Print unmatched sequences
                print(f'{i:06}. {line[:-1]}')
                #time.sleep(.1)
        elif i % 2 != 0 and match:
            result += f'{prevLine}\n{line[trimFirst:match-trimLast]}\n'
            match = None
        else:
            isRead = line[0] == '@'
            prevLine = line[:-1]# + finalLength

    if not quiet:
        print(f'Trimmed {count} lines, found {matches} matches\n')

with open(outFilePath, 'w') as outFile:
    outFile.write(result)
