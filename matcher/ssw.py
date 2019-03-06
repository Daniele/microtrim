'''
Striped Smith-Waterman matcher
'''

import math

from skbio.alignment import StripedSmithWaterman as ssw


def build(adapter, args):
    '''
    Build a Striped Smith–Waterman matcher with parameters:
      - match_only
    '''
    adapter = adapter[:args.match_only][::-1]

    def match(line):
        rline = line[::-1]
        target = ssw(rline, match_score=1)
        match = target(adapter)
        align = match['optimal_alignment_score']/len(adapter)
        if align >= .8:
            return -(match['query_begin']+len(adapter))
        return None

    return match
