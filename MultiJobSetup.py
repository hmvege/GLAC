#!/usr/bin/env python3

import os
import re
import argparse
import subprocess


def natural_sort(l):
    """
    Natural sorting function.

    Args:
        l: list of strings where each string contains a number, either on 
            format of 1,2,3... or 00001, 00002 ect.

    Returns:
        A sorted list.
    """

    def convert(text): return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): return [convert(c) for c in re.split(r'(\d+)', key)]
    return sorted(l, key=alphanum_key)


def build_cmd(run_config, field_cfg, args, run_type=None, add_Ncf=None):
    """Builds a single argument list"""

    cmd = ["python2", "createJobs.py"]

    if args.dryrun:
        cmd.append("--dryrun")

    cmd.append("load")
    cmd.append(run_config)
    cmd.append("-s")
    cmd.append("{0:s}".format(args.system))

    if not isinstance(add_Ncf, type(None)):
        num_cfgs_args = ["-NCf", "%d" % add_Ncf]
    else:
        num_cfgs_args = []

    if not isinstance(field_cfg, type(None)):
        # Checks if we are loading a folder or single config.
        if os.path.isdir(field_cfg):
            cmd.append("-lcfg")
        else:
            cmd.append("-lcfgr")
            # At least 1 config must be generated with lcfgr
            # num_cfgs_args = ["-NCf", "1"]
        cmd.append(field_cfg)

    cmd += num_cfgs_args

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
    folder = os.path.normpath(args.folder)
    try:
        base_path, scaling_type, run_type = folder.split("/")
        scaling_type = scaling_type.split("_")[0]
    except ValueError:
        paths = folder.split("/")
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
        output_folder = args.load_multiple_field_configs
        _cfgs = os.listdir(output_folder)

        # Because we were silly, and didnt make gen the default.
        # if run_type == "cfg_gen":
        #     _run_type = "gen"
        # else:
        #     _run_type = run_type

        # Configurations to use in lmcfgs is always in gen folder
        _run_type = "gen" 

        def _filter_func(_c):
            """Function for filtering out configs
             not being of correct run type(e.g. flow, io, cfg_gen)."""
            c_contents = _c.split("_")
            if ((_run_type == c_contents[-1]) and 
                ("np" in c_contents[-2]) and
                (c_contents[0] == scaling_type)):
                return _c

        # Temporary sorts the different runs we are going to look for
        # configs in.
        _tmp_cfgs = natural_sort(filter(_filter_func, _cfgs))

        # Gets the processor size (REDUNDANT)
        _tmp_proc_size = [re.findall(r"(\w*np\d+\w*)", _i)[0]
                          for _i in _tmp_cfgs]

        # Assures that we actually picks up a config to use
        if run_type != "cfg_gen":
            _tmp_cfgs = [_c.replace("flow", "gen")
                         for _c in _tmp_cfgs]

        # Builds up the field config folders
        _field_cfg_folders = [os.path.join(output_folder, _c,
                                           "field_configurations")
                              for _c in _tmp_cfgs]

        # Makes sure we have at least one config in the folder to load.
        assert_msg = (
            "Config folders for {0:s} does not contain any configs.".format(
                " ".join(_field_cfg_folders)))
        for elem in [len(os.listdir(_c)) for _c in _field_cfg_folders]:
            assert elem > 0, assert_msg

        # Builds up the field config folders with full field .bin paths
        if run_type == "flow":
            field_cfgs = _field_cfg_folders
        else:
            field_cfgs = [os.path.join(_c, natural_sort(os.listdir(_c))[0])
                          for _c in _field_cfg_folders]

    run_cfgs = [os.path.join(args.folder, _c) for _c in run_cfgs]

    if ((not args.load_single_field_config) and
            (not args.load_multiple_field_configs)):
        field_cfgs = [None for i in range(N_configs)]
    elif (args.load_multiple_field_configs and
          not args.load_single_field_config):

        for run_cfg, field_cfg in zip(run_cfgs, field_cfgs):
            # Checks that we have use the correct number of processors
            numprocs_run_cfg = re.findall(r"np(\d+)", run_cfg)[0]
            numprocs_field_cfg = re.findall(r"np(\d+)", field_cfg)[0]
            assert numprocs_run_cfg == numprocs_field_cfg, (
                "Number of processors used do not "
                "match:\n    '{0}'\nand\n    '{1}'".format(
                    run_cfg, field_cfg))

    if ((not args.load_single_field_config) and
            (not args.load_multiple_field_configs)):
        add_Ncf = None
    else:
        if run_type == "io":
            add_Ncf = 10
        elif run_type == "cfg_gen":
            add_Ncf = 1
        elif run_type == "flow":
            add_Ncf = None
        else:
            raise ValueError("{0:s} not recognized.".format(run_type))

    cmds = []
    for run_cfg, field_cfg in zip(run_cfgs, field_cfgs):
        cmds.append(build_cmd(run_cfg, field_cfg, args,
                              run_type=run_type, add_Ncf=add_Ncf))

    for cmd in cmds:
        print("> " + " ".join(cmd))
        # if not args.dryrun:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        print(proc.stdout.read())

    print("Started {0:d} jobs.".format(len(cmds)))


if __name__ == '__main__':
    main()
