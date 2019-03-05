'''
NDLeven matcher
'''

import math

import Levenshtein
from pyxdameraulevenshtein import normalized_damerau_levenshtein_distance as ndleven


def build(adapter, args):
    max_distance = args.max_distance
    stop_after = args.stop_after

    def match(line):
        for j, char in enumerate(line[:math.floor(len(line)/stop_after)]):
            possibleMatch = line[j:j+len(adapter)]
            if ndleven(adapter, possibleMatch) <= max_distance:
                return -(j+len(adapter))

    return match
