'''
NDLeven matcher
'''

import math

from pyxdameraulevenshtein import normalized_damerau_levenshtein_distance as ndleven


def build(adapter, args):
    adapter = adapter[:args.match_only][::-1]

    def match(line):
        rline = line[::-1]
        for j, char in enumerate(rline[:math.floor(len(rline)/args.stop_after)]):
            possibleMatch = rline[j:j+len(adapter)]
            if ndleven(adapter, possibleMatch) <= args.max_distance:
                return -(j+len(adapter))

    return match
