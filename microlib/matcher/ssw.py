'''
Striped Smith-Waterman matcher
'''

import math

from ssw import SSW


def build(adapter, args):
    '''
    Build a Striped Smith–Waterman matcher with parameters:
      - match_only
    '''
    adapter = adapter[:args.match_only][::-1]

    def match(line):
        rline = line[::-1]
        matcher = SSW(1)
        matcher.setRead(adapter)
        matcher.setReference(rline)
        align = matcher.align()
                
        if align.optimal_score/len(adapter) >= 1 - args.max_distance:
            return -(align.reference_start+len(adapter)) + 1
        return None

    return match
