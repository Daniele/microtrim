'''
Adapter matcher
'''
import time

VERBOSE = False
BASES = ('A', 'C', 'G', 'T')


def addNewAdapterToSet(ad, adSet):
    adSet.add(ad)
    if VERBOSE:
        print(f'Adding {ad}')
        time.sleep(.1)
    return adSet


def makeAdapters(adapter, matchOnly):
    adapters = set()
    adapters.add(adapter)
    if VERBOSE:
        print(f'Adding {adapter}')
        time.sleep(.1)

    if VERBOSE:
        print('\nSubstitutions\n')
    for i, x in enumerate(adapter):
        # Adapter with wrong substitution of 1 base
        for j in BASES:
            ad = adapter[:i] + j + adapter[i+1:]
            adapters = addNewAdapterToSet(ad[:matchOnly], adapters)

    if VERBOSE:
        print('\nAdditions (1)\n')
    newAdapters1 = set()
    for adapter in adapters:
        for i, x in enumerate(adapter):
            # Adapter with wrong addition of 1 base
            for l in BASES:
                ad = adapter[:i] + l + adapter[i:]
                newAdapters1 = addNewAdapterToSet(
                    ad[:matchOnly], newAdapters1)

    if VERBOSE:
        print('\nRemovals (1)\n')
    newAdapters2 = set()
    for adapter in adapters:
        for i, x in enumerate(adapter):
            # Adapter with wrong removal of 1 base
            ad = adapter[:i] + adapter[i+1:]
            newAdapters2 = addNewAdapterToSet(ad[:matchOnly], newAdapters2)

    adapters = adapters.union(newAdapters1).union(newAdapters2)
    return adapters


def build(adapter, args):
    adapter = adapter[:args.match_only]
    adapters = sorted(makeAdapters(adapter, args.match_only))

    def match(line):
        for adapter in adapters:
            if adapter in line:
                return line.find(adapter)
    return match
