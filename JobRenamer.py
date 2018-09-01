#!/usr/bin/env python3

import sys
import os
import shutil
import re
import argparse
import json
import numpy as np

__NAME_STR = "Batch name:"
__TIME_STR = "Program complete. Time used:"
__TOT_UP_TIME = "Total update time for 600 updates: "
__LAT_DIMS = "Lattice dimensions(spatial, temporal): "
__SUB_DIMS = "Sub lattice dimensions:                "

_regex_name_str = re.compile(r"%s[ ]*(\w)" % __NAME_STR)
_regex_time_str = re.compile(r"(%s[ ]\w+)" % __TIME_STR)


def JobRenamer(folder, get_run_times=False, verbose=False, dryrun=False):
    """Function for renaming jobs."""

    if dryrun:
        print("*"*10 + " DRYRUN " + "*"*10)

    def filter_func(f):
        if len(re.findall(r"(.o\d{8})", f)) == 1:
            return f, "torque"
        if os.path.splitext(f)[-1] == ".out":
            return f, "slurm"

    folder_list = filter(filter_func, os.listdir(folder))

    run_statistics = []

    if get_run_times:
        json_dict = {"runs": []}

    for f in folder_list:
        fpath = os.path.join(folder, f)
        if not __check_if_complete(fpath):
            continue

        runname, job_dict = __get_job_content(fpath)

        if get_run_times:
            json_dict["runs"].append(job_dict)

        if not runname in f:
            fhead, fext = os.path.splitext(f)
            new_path = os.path.join(folder, fhead + "-" + runname + fext)

            if verbose:
                print("File: {0:20s} JobName: {1:40s} TimeUsed: {2:15f} NewName: {3:s}".format(
                    f, runname, job_dict["time"], new_path))

            print("> mv {} {}".format(fpath, new_path))

            if not dryrun:
                shutil.move(fpath, new_path)

    if get_run_times:
        json_fpath = os.path.join(folder, "run_times.json")
        print("Writing json job run times file at location {0:<s}".format(
            json_fpath))

        if verbose:
            print(json.dumps(json_dict, indent=4, separators=(", ", ": ")), "\n")

        if not dryrun:
            with file(json_fpath, "w+") as json_file:
                json.dump(json_dict, json_file, indent=4)


def __check_if_complete(fpath):
    """Checks if job is complete."""
    is_complete = False
    with open(fpath, "r") as f:
        for i, l in enumerate(f):
            if "Program complete" in l:
                is_complete = True
                break
        else:
            is_complete = False
    return is_complete


def __get_job_content(fpath):
    found_first_elem = False
    sub_dims_found = False
    lat_dims_found = False
    tot_up_time = False

    runname = ""
    seconds_used = 0
    update_time = 0
    size_N = 0
    size_NT = 0
    sub_dims = []

    with open(fpath, "r") as f:
        for i, l in enumerate(f):

            if not found_first_elem and __NAME_STR in l:
                runname = l.split(" ")[-1].strip("\n")
                found_first_elem = True

            if not lat_dims_found and __LAT_DIMS in l:
                __temp_sizes = l.split(__LAT_DIMS)[-1].strip("\n").split(" ")
                size_N, size_NT = map(int, __temp_sizes)
                lat_dims_found = True

            if not sub_dims_found and __SUB_DIMS in l:
                sub_dims = l.split(__SUB_DIMS)[-1].strip("\n ").split(" ")
                sub_dims = map(int, sub_dims)
                sub_dims_found = True

            if not tot_up_time and __TOT_UP_TIME in l:
                update_time = l.split(__TOT_UP_TIME)
                update_time = update_time[-1].strip("\n ")
                update_time = float(update_time.split(" ")[0])
                tot_up_time = True

            if found_first_elem and __TIME_STR in l:
                seconds_used = re.findall(r"\(([\.\w]+){1} seconds\)", l)[0]
                seconds_used = float(seconds_used)
                break

    if runname == "" or seconds_used == 0:
        raise ValueError(
            ("missing runname {0:s} or seconds used {1:d}: is job complete?".format(
                runname, seconds_used)))

    results_dictionary = {
        "runname": runname,
        "time": seconds_used,
        "update_time": update_time,
        "subdims": sub_dims,
        "subdimsize": np.prod(sub_dims),
        "N": size_N,
        "NT": size_NT,
        "totsize": size_N*size_N*size_N*size_NT,
    }

    return runname, results_dictionary


def main():
    desc_str = \
        """Quick, small program for appending createJobs.py "runname" 
        to .out files from slurm or torque."""
    parser = argparse.ArgumentParser(prog="JobRenamer",
                                     description=desc_str)

    # Base arguments
    parser.add_argument(
        "--version", action="version", version="%(prog)s 1.0")
    parser.add_argument(
        "-dr", "--dryrun", default=False, action="store_true",
        help="Dryrun to not perform any permanent actions.")
    parser.add_argument(
        "-v", "--verbose", default=False, action="store_true",
        help="A more verbose output.")

    # Main arguments for program
    parser.add_argument("folder", type=str,
                        help="Folder to modify file names in.")
    parser.add_argument("-extms", "--extract_times",
                        default=False, action="store_true",
                        help=("Retrieves time spend in program(s) and "
                              "saves in a text file."))

    if len(sys.argv) == 1:
        arguments = ["test_out", "--extract_times"]
        # arguments = ["testoutfolder", "--dryrun",
        #              "--verbose", "--extract_times"]
        args = parser.parse_args(arguments)
    else:
        args = parser.parse_args()

    JobRenamer(args.folder, get_run_times=args.extract_times,
               verbose=args.verbose, dryrun=args.dryrun)


if __name__ == '__main__':
    main()
