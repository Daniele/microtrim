from __future__ import print_function

import re

from lib import run_test


def parse_time(string):
    str_re = re.findall(r"Matching time: (.*)", string)
    return int(str_re[0])


# Adagen Fast
run_test(
    (
        "taskset -c 0 python microtrim.py "
        + "-m adagen-fast --out-file adagen_fast_seq.trimmed.fastq "
        + "--trim-to 23 --trim-first 0 --trim-last 0 "
        + "--workers 0 "
        + "2> /dev/null"
    ),
    {},
    n=3,
    simulate=False,
    time_unit="ms",
    time_parser=parse_time,
)
