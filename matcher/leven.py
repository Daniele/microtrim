'''
Leven matcher
'''
import math

import Levenshtein


def build(adapter, args):
    adapter = adapter[:args.match_only][::-1]

    def match(line):
        for j, char in enumerate(line[:math.floor(len(line)/args.stop_after)]):
            possibleMatch = line[j:j+len(adapter)]
            if Levenshtein.ratio(adapter, possibleMatch) >= .9:
                return -(j+len(adapter))

    return match
