#!/usr/bin/env python3

import os
import ast
import re
import json
import types
import shutil
import fileinput
import argparse
import sys


def natural_sort(l):
    """
    Natural sorting function.

    Args:
        l: list of strings where each string contains a number, either on 
            format of 1,2,3... or 00001, 00002 ect.

    Returns:
        A sorted list.
    """

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('(\d+)', key)]
    return sorted(l, key=alphanum_key)



def build_cmd(run_config, field_cfg, args, run_type=None):
    """Builds a single argument list"""

    cmd = ["python2", "createJobs"]

    if args.dryrun:
        cmd.append("--dryrun")

    cmd.append("load")
    cmd.append(run_config)
    cmd.append("-s")
    cmd.append("{0:s}".format(args.system))

    if run_type == "io":
        num_cfgs_args = ["-NCf", "10"]
    else:
        num_cfgs_args = []

    if not isinstance(field_cfg, type(None)):
        # Checks if we are loading a folder or single config.
        if os.isdir(field_cfg):
            cmd.append("-lcfg")
        else:
            cmd.append("-lcfgr")
            # At least 1 config must be generated with lcfgr
            num_cfgs_args = ["-NCf", "1"]
        cmd.append(field_cfg)

    cmd.append("--ignore_tasks_per_node")

    return cmd


def main():
    desc_str = \
        """Small program for initiating multiple jobs on slurm or torque."""
    parser = argparse.ArgumentParser(prog="MultiJobSetup",
                                     description=desc_str)

    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    parser.add_argument("-dr", "--dryrun", default=False, action="store_true",
                        help="Dryrun to not perform any permanent actions.")
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="A more verbose output.")

    parser.add_argument("-s", "--system", type=str, default="abel",
                        help="System we are running on.")

    parser.add_argument("folder", default=False,
                        help="Folder of configs to initiate job for.")
    parser.add_argument("-lscfg", "--load_single_field_config", default=False,
                        type=str,
                        help=("Loads a the same single field configuration "
                              "for each of the run configurations."))
    parser.add_argument("-lmcfgs", "--load_multiple_field_configs",
                        type=str, default=False,
                        help=("Loads a unique single field "
                              "configuration for each of the "
                              "run configurations,"
                              " based on the run name"))

    args = parser.parse_args()

    if args.dryrun:
        print("*"*10 + " DRYRUN " + "*"*10)

    # # Assures that we have at least lscfg or lmcfgs arguemnt
    # if ((not args.load_single_field_config) and
    #         (not args.load_multiple_field_configs)):
    #     exit("At least one argument of either lscfg or lmcfgs allowed.")

    if args.load_single_field_config and args.load_multiple_field_configs:
        exit("Only one argument of either lscfg or lmcfgs allowed.")

    # Gets base path name and run type
    try:
        base_path, scaling_type, run_type = args.folder.split("/")
    except ValueError:
        paths = args.folder.split("/")
        run_type = paths[-1]

    # Gets number of run configs
    run_cfgs = natural_sort(os.listdir(args.folder))
    N_configs = len(run_cfgs)

    # Gets field configurations of each of the run configs
    field_cfgs = []
    if args.load_single_field_config != False:
        _cfg = args.load_single_field_config
        field_cfgs = [_cfg for i in range(N_configs)]

        # if run_type == "io":

        # elif run_type == "cfg_gen":
        # elif run_type == "flow":
        # else:
        #     raise NameError("Config folder not recognized: {0}".format(run_type))

    if args.load_multiple_field_configs != False:
        _cfgs = os.listdir(args.load_multiple_field_configs)

        def _filter_func(_c):
            """Function for filtering out configs
             not being of correct run type(e.g. flow, io, cfg_gen)."""
            if _c in run_type:
                return _c

        _tmp_cfgs = natural_sort(filter(_filter_func, _cfgs))

        print(_tmp_cfgs)

        # for run_cfg in run_cfgs:

        #     if "weak" in run_type:
        #         field_cfgs.append()

    if ((not args.load_single_field_config) and
            (not args.load_multiple_field_configs)):
        field_cfgs = [None for i in range(N_configs)]

    cmds = []
    for run_cfg, field_cfg in zip(run_cfgs, field_cfgs):
        cmds.append(build_cmd(os.path.join(args.folder, run_cfg),
                              field_cfg, args, run_type=run_type))

    for cmd in cmds:
        print("> " + " ".join(cmd))
        if not args.dryrun:
            raise NotImplementedError("running createJobs not implemented yet")
            # proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            # print(proc.stdout.read())
    print("Started {0:d} jobs.".format(len(cmds)))


# python2 createJobs.py $DRYRUN load $filename - s $SYSTEM $GAUGECONFIGARGUMENT - -ignore_tasks_per_node $ADDITIONAL_ARGS
if __name__ == '__main__':
    main()
