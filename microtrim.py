#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import os
import sys
import mmap
import json
import time
import Levenshtein
from pyxdameraulevenshtein import damerau_levenshtein_distance as dleven
from pyxdameraulevenshtein import normalized_damerau_levenshtein_distance as ndleven
from itertools import product

# TODO - Make options configurable through argparse
version = 1
inDirPath = './'
adapter = 'TGGAATTCTCGGGTGCCAAGG'
matches = 0
trimFirst = 4
ignoreAfter = 10
outFilePath = './SRR8311267.trimmed.fq'

abc = ('A','C','G','T')
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
                newAdapters1 = addNewAdapterToSet(ad[:ignoreAfter], newAdapters1)

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

def getSubstrings(string,length):
    subs = set()
    for j in string[10:-length-1]:
        subs.add(line[j:j+length-1])
    return subs

if version == 1:
    adapters = makeAdaptersV1(set())
else:
    adapters = makeAdaptersV2(set())

print()
if trimFirst == 0:
    print(f'Not trimming any initial bases')
else:    
    print(f'Trimming the first {trimFirst} bases')
print(f'Trimming adapter: {adapter}')
if version == 2:
    print(f'Considering only first {ignoreAfter} bases of adapter: {adapter[:ignoreAfter]}')
print(f'Considering {len(adapters)} possible variants of the adapter')
print(f'Trimming all bases after the adapter (if present)')
print(f'Saving to file: {outFilePath}')
print()

result = ""

with open(inDirPath + 'SRR8311267.fastq', 'r+b') as infile:
    m = mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ)
    length = len(adapter)
    finalLength = str(51-length)
    match = None
    isRead = False
    count = 0
    for i, line in enumerate(iter(m.readline, b"")):
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
            match = matchAdapters(line, adapters)
            if match:
                matches += 1
                sys.stdout.write(f'\r{count} • {matches}\r')
                result += prevLine + '\n' + line[trimFirst:match] + '\n'
                #if len(line[:match]) < 10:
                #    print(adapter)
                #    print(match)
                #    print(line)
                #    time.sleep(.3)
            #else:
                # Print unmatched sequences
                #print(f'{i:06}. {line[:-1]}')
                #time.sleep(.3)
        elif i % 2 != 0 and match:
            result += prevLine + '\n' + line[trimFirst:match] + '\n'
            match = None
        else:
            isRead = line[0] == '@'
            prevLine = line[:-3] + finalLength
    
    print(f'Trimmed {count} lines, found {matches} matches\n')
    
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
    #if matches >= 5000:
        #sys.exit()
        #time.sleep(1) 
        #matches = re.findall('', line)
        #for match in matches:
        #    adapterMatch = line[index:index+len(adapter)]
        #    print(adapterMatch)
        #    print(leven(adapter, adapterMatch))
        #    print()
        #    time.sleep(1)

with open(outFilePath, 'w') as outFile:
    outFile.write(result)
