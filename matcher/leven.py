'''
Leven matcher
'''
import math

import Levenshtein


def leven(s1, s2):
    return Levenshtein.ratio(s1, s2)


def build(adapter, args):
    stop_after = args.stop_after

    def match(line):
        for j, char in enumerate(line[:math.floor(len(line)/stop_after)]):
            possibleMatch = line[j:j+len(adapter)]
            if leven(adapter, possibleMatch) >= .9:
                return -(j+len(adapter))

    return match
