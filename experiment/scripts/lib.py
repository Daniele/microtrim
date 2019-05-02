from __future__ import print_function
from __future__ import division

import subprocess
import math
import sys
import os
from itertools import product


def time_parser_wrapper(fun, string):
    assert isinstance(string, (str, bytes))
    if isinstance(string, bytes):
        string = string.decode("UTF-8", "ignore")
    return fun(string)


def default_time_parser(string):
    return int(string)


def get_time(command, time_parser, times=5):
    data = []
    for _ in range(0, times):
        cmd_out = subprocess.check_output(command, shell=True)
        ctime = time_parser_wrapper(time_parser, cmd_out)
        data.append(ctime)
    avg = sum(data) / len(data)
    std = math.sqrt(sum((e - avg) ** 2 for e in data) / len(data))
    vmin, vmax = min(data), max(data)
    return avg, std, vmin, vmax


def tsend(s, simulate):
    if simulate:
        print(s)
    else:
        subprocess.call(["tsend", s])


def print_percentage(i, tot_inst, my_name, every, simulate):
    slot = math.ceil(tot_inst / every)
    if i % slot == 0:
        tsend("{} {}% complete".format(my_name, i / slot / every * 100), simulate)


def print_header(f, par_name, time_unit):
    f.write(
        "\t".join(par_name)
        + "\t"
        + "\t".join(
            map(
                lambda x: x.format(time_unit),
                ["tc({})", "std({})", "min({})", "max({})"],
            )
        )
        + "\n"
    )


def run_test(
    cmd,
    par,
    name=None,
    time_unit="s",
    n=5,
    simulate=False,
    path=".",
    every=4,
    time_parser=None,
):
    time_parser = default_time_parser if time_parser is None else time_parser
    try:
        par_name = par.keys()
        par_value = par.values()
        my_name = os.path.splitext(sys.argv[0])[0] if name is None else name
        my_name += "-sim" if simulate else ""
        file_name = "{}.log".format(my_name)
        if not simulate:
            i = 1
            while os.path.isfile(file_name):
                file_name = "{}-{}.log".format(my_name, i)
                i += 1
        # TODO: fix path
        file_path = os.path.join(path, file_name)
        with open(file_path, "w") as f:
            print_header(f, par_name, time_unit)
            prod_list = list(product(*par_value))
            tot_inst = len(prod_list)
            for i, inst in enumerate(product(*par_value)):
                p = dict(zip(par_name, inst))
                cmd_str = cmd.format(p)
                if simulate:
                    print("  " + cmd_str)
                    t_res = 1.2, 0.4, 1.0, 2.3
                else:
                    t_res = get_time(cmd_str, time_parser, times=n)
                f.write(
                    "\t".join(
                        map(
                            lambda x: str(x[0]) if isinstance(x, tuple) else str(x),
                            inst + t_res,
                        )
                    )
                    + "\n"
                )
                print_percentage(i, tot_inst, my_name, every, simulate)
        tsend("{} ok".format(my_name), simulate)
    except Exception as e:
        tsend("{} error {}".format(my_name, e), simulate)
        print(e)
