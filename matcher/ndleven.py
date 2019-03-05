'''
NDLeven matcher
'''

import math

from pyxdameraulevenshtein import normalized_damerau_levenshtein_distance as ndleven


def build(adapter, args):
    adapter = adapter[:args.match_only][::-1]

    def match(line):
        for j, char in enumerate(line[:math.floor(len(line)/args.stop_after)]):
            possibleMatch = line[j:j+len(adapter)]
            if ndleven(adapter, possibleMatch) <= args.max_distance:
                return -(j+len(adapter))

    return match
