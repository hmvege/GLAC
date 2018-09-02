#!/usr/bin/env python3

import os
import re
import argparse


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

    def alphanum_key(key): return [convert(c) for c in re.split('(\d+)', key)]
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
        if os.path.isdir(field_cfg):
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
        output_folder = args.load_multiple_field_configs
        _cfgs = os.listdir(output_folder)

        print(_cfgs)

        # Because we were silly, and didnt make gen the default.
        if run_type == "cfg_gen":
            _run_type = "gen"
        else:
            _run_type = run_type

        def _filter_func(_c):
            """Function for filtering out configs
             not being of correct run type(e.g. flow, io, cfg_gen)."""
            if _run_type in _c and "np" in _c:
                return _c

        # Temporary sorts the different runs we are going to look for
        # configs in.
        _tmp_cfgs = natural_sort(filter(_filter_func, _cfgs))

        # Gets the processor size (REDUNDANT)
        _tmp_proc_size = [re.findall(r"(\w*np\d+\w*)", _i)[0]
                          for _i in _tmp_cfgs]

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
