'''
Leven matcher
'''
import math

import Levenshtein


def build(adapter, args):
    adapter = adapter[:args.match_only][::-1]

    def match(line):
        rline = line[::-1]
        for j, char in enumerate(rline[:math.floor(len(rline)/args.stop_after)]):
            possibleMatch = rline[j:j+len(adapter)]
            if Levenshtein.ratio(adapter, possibleMatch) >= .9:
                return -(j+len(adapter))

    return match
