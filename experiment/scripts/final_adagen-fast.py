from __future__ import print_function

import re

from lib import run_test


def parse_time(string):
    str_re = re.findall(r"Matching time: (.*)", string)
    return int(str_re[0])


MAX_CORE = 48
WORKERS = [1] + list(range(0, MAX_CORE + 1, 4))[1:]
CORES = map(lambda x: x - 1, WORKERS)

# Adagen Fast
run_test(
    (
        "taskset -c 0-{0[worker][1]} python microtrim.py "
        + "-m adagen-fast --out-file adagen_fast.trimmed.fastq "
        + "--trim-to 23 --trim-first 0 --trim-last 0 "
        + "--workers {0[worker][0]} --chunk {0[chunk]} "
        + "2> /dev/null"
    ),
    {"worker": list(zip(WORKERS, CORES)), "chunk": [1000]},
    n=3,
    simulate=False,
    time_unit="ms",
    time_parser=parse_time,
)
