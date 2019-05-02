from __future__ import print_function

import re

from lib import run_test


def parse_time(string):
    str_re = re.findall(r"Matching time: (.*)", string)
    return int(str_re[0])


# Damerau-Levenshtein
run_test(
    (
        "taskset -c 0 python microtrim.py "
        + "-m leven --out-file leven_seq.trimmed.fastq --max-distance .1 "
        + "--match-only 18 --trim-to 23 --trim-first 0 --trim-last 0 "
        + "--workers 0 "
        + "2> /dev/null"
    ),
    {},
    n=3,
    simulate=False,
    time_unit="ms",
    time_parser=parse_time,
)
