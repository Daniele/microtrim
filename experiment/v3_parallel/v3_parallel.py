#!/usr/bin/python
from __future__ import print_function

import subprocess
import math
import re
import sys
from itertools import product
import time


def get_time(command, times=5):
    data = []
    for _ in range(0, times):
        t = time.perf_counter()
        subprocess.check_output(command, shell=True)
        elapsed_s = time.perf_counter() - t
        elapsed = math.floor(elapsed_s * 1000)
        data.append(elapsed)
    avg = sum(data) / len(data)
    std = int(math.sqrt(sum([(e - avg)**2 for e in data]) / len(data)))
    vmin, vmax = min(data), max(data)
    return avg, std, vmin, vmax


MAX_THREADS = 24
# CMD = '''taskset -c 0-{0[worker]} python microtrim-parallel4.py \\
#          --max-threads={0[worker]} --chunks={0[chunks]} -V={0[version]} \\
#          > /dev/null
# '''
CMD = '''python microtrim.py --matcher={0[matcher]} \\
         --worker={0[worker]} --chunk={0[chunks]} \\
         > /dev/null
'''
WORKER = range(2, MAX_THREADS + 1)
CHUNKS = [1000]
MATCHER = ['ndleven']
# MATCHER = ['adagen', 'adagen-fast', 'leven', 'ndleven', 'ssw']


def get_cmd(values):
    return CMD.format(values)


print('matcher', 'worker', 'chunk', 'time (ms)', 'time (s)')

for c, w, m in product(CHUNKS, WORKER, MATCHER):
    di = {'worker': w,
          'chunks': c,
          'matcher': m}
    cmd = get_cmd(di)
    # print(cmd)
    avg, std, vmin, vmax = get_time(cmd, times=3)
    print(di['matcher'], di['worker'], di['chunks'], avg, avg/1000.0)
