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

def makeAdapters(adapter):
    adapters = set()
    adapters.add(adapter[:-8])
    for i, x in enumerate(adapter):
        for j in BASES:
            ad = adapter[:i] + j + adapter[i+1:]
            adapters.add(ad[:-8])
            for l in BASES:
                ad = adapter[:i] + l + adapter[i:]
                adapters.add(ad[:-8])
            ad = adapter[:i] + adapter[i+1:]
            adapters.add(ad[:-8])
            ad = adapter[:i] + adapter[i+2:]
            adapters.add(ad[:-8])
    return adapters


def build(adapter, args):
    adapters = sorted(makeAdapters(adapter))

    def match(line):
        for adapter in adapters:
            if adapter in line:
                return line.find(adapter)
    return match
