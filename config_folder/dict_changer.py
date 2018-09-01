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
from collections import OrderedDict


def change_config_folder(folder, replace_values, append_to_values=None,
                         modify_values=None, dryrun=False, verbose=False):

    if dryrun:
        print("**** DRYRUN ****")

    assert_err_msg = ("Please provide dict with "
                      "values to replace in config folder %s" % folder)
    assert isinstance(replace_values, dict), assert_err_msg

    _temp = [f for f in os.listdir(folder) if f.endswith(".py")]
    files = sorted(_temp, key=lambda s: int(re.findall(r"\d+", s)[0]))

    for i, f in enumerate(files):
        if verbose:
            print("Modifying file %s" % f)

        fpath = os.path.join(folder, f)

        list_name_appended = False

        ordered_keys = []
        with open(fpath, "r") as raw_file:
            for l in raw_file:
                found_words = re.findall(r"\w+", l)
                if len(found_words) < 2:

                    # In cases where we have json formatting
                    if (len(found_words) != 0) and \
                        (isinstance(found_words[0], str)) and \
                            (not list_name_appended):

                        ordered_keys.append(found_words[0])
                        list_name_appended = True

                    elif len(found_words) == 0:
                        list_name_appended = False

                    else:
                        continue

                else:
                    ordered_keys.append(found_words[0])

        data_dict = ast.literal_eval(open(fpath, "r").read())

        # Modifies dictionary values
        for item in replace_values.items():
            key, val = item

            if verbose:
                if isinstance(val, list):
                    if append_to_values and isinstance(val[0], str):
                        if val[0] in data_dict[key]:
                            print_val = data_dict[key]
                        else:
                            print_val = data_dict[key] + val[0]

                        if key == "runName":
                            if f.split(".")[0] + str(val[i]) == data_dict[key]:
                                print_val = data_dict[key]
                            else:
                                print_val = f.split(".")[0] + str(val[i])
                        else:
                            if val[i] in data_dict[key]:
                                print_val = data_dict[key]
                            else:
                                print_val += str(val[i])

                    else:
                        print_val = val[i]
                else:
                    print_val = val

                if data_dict[key] != print_val:
                    print("    %25s: %30s -> %s" % (key, data_dict[key], print_val))

            # If we have a list, we first check if we are to append
            # Then we check if we are to modify the value somehow through a function
            # Then we simply replace value in file
            if isinstance(val, list):

                if len(val) == len(files):
                    if not isinstance(append_to_values, type(None)) and \
                            key in append_to_values:

                        # Checks if we have already appended string on previous
                        # run.
                        # if val[i] in data_dict[key]:
                        if key == "runName":
                            # print f.split(".")[0]+str(val[i]) == data_dict[key], data_dict[key], f.split(".")[0]+str(val[i])
                            if f.split(".")[0] + str(val[i]) == data_dict[key]:
                                continue
                            else:
                                data_dict[key] = f.split(".")[0] + str(val[i])
                        else:
                            if val[i] in data_dict[key]:
                                continue
                            else:
                                data_dict[key] += str(val[i])


                    elif not isinstance(modify_values, type(None)):
                        data_dict[key] += modify_values(val[i])

                    else:
                        data_dict[key] = val[i]
            else:
                data_dict[key] = val

        # Reorganized dictionary values
        ord_data_dict = OrderedDict(
            [(key, data_dict[key]) for key in ordered_keys])

        json_fpath = os.path.join(folder, os.path.splitext(f)[0] + ".py")

        # Writes json to file
        if verbose:
            print(("Writing python configuration file at "
                   "location {0:<s}".format(json_fpath))
                  )
        if not dryrun:
            with open(json_fpath, "w+") as json_file:
                json.dump(ord_data_dict, json_file, indent=4)

            # Re-names false->False, true->True
            f1 = open(json_fpath, "r")
            f2_path = os.path.join(
                folder, os.path.splitext(f)[0] + "_temp2.py")
            f2 = open(f2_path, "w")
            replacement_values = [("false", "False"), ("true", "True")]
            for l1 in f1:
                # Replaces correct true/false values for python
                for old, new in replacement_values:
                    if old in l1:
                        f2.write(l1.replace(old, new))
                        break
                else:
                    f2.write(l1)
            f1.close()
            f2.close()
            shutil.copy(f2_path, json_fpath)
            os.remove(f2_path)

        if verbose:
            print("")
        # print(json.dumps(ord_data_dict, indent=4, separators=(", ", ": ")), "\n")
        # exit("Exits after changing file %s" % f)

    print("Changed config files in folder %s." % folder)

    # exit("Exits after folder %s" % folder)


def main():
    desc_str = \
        """Program for changer folder of configurations for GluonAction."""
    parser = argparse.ArgumentParser(prog="Job configuration changer",
                                     description=desc_str)

    # Base arguments
    parser.add_argument(
        "--version", action="version", version="%(prog)s 0.1")
    parser.add_argument(
        "--dryrun", default=False, action="store_true",
        help="Dryrun to not perform any permanent actions.")
    parser.add_argument(
        "-v", "--verbose", default=False, action="store_true",
        help="A more verbose output.")

    # subparser = parser.add_subparsers(dest="subparser")

    # preset_choices = [
    #     "weak_flow", "weak_cfg_gen", "weak_io", "strong_flow",
    #     "strong_cfg_gen", "strong_io"]
    # preset_parser = subparser.add_parser(
    #     "preset", help="runs a specific preset arguments.")
    # preset_parser.add_argument(
    #     "folder_type", choices=preset_choices,
    #     help="Type of folder to implement preset values of")

    # old_method_parser = subparser.add_parser(
    #     "old_method", help="Runs old methods of file parsing.")

    # folder_parser = parser.ArgumentParser("folder")

    # Folder to change stuff in
    parser.add_argument(
        "folder", type=str,
        help="Folder of .py files with dictionaries to change.")

    # Arguments to change
    parser.add_argument(
        "-rn_app", "--runname_append", default=None,
        type=str,
        help="filename extension to add to config files in given folder.")
    parser.add_argument(
        "-rn_rep", "--runname_replace", default=None,
        type=str,
        help=("filename to replace in 'runname' section in config files "
              "in given folder."))
    parser.add_argument(
        "-NCf", "-NCfg", "-NCfgs", "--NConfigs", type=int,
        help="Number of configurations")
    parser.add_argument(
        "-NCor", "-NCorr", "--NCor", default=None,
        type=int, help="Number of correlation updates")
    parser.add_argument(
        "-NUp", "--NUpdates", default=None,
        type=int, help="Number of single link updates")
    parser.add_argument(
        "-NFlows", "--NFlows", default=None, type=int,
        help="Number of number of flow steps")
    parser.add_argument(
        "-NTh", "--NTherm", default=None, type=int,
        help="number of thermalization steps")
    parser.add_argument(
        "-sc", "--storeCfgs", default=False,
        action="store_true", help="Store configurations(Default=false)")
    parser.add_argument(
        "-chr", "--cpu_approx_runtime_hr", default=None, nargs="*",
        type=int, help="Number of approximate cpu-hours for *each* config.")
    parser.add_argument(
        "-cmin", "--cpu_approx_runtime_min", default=None, nargs="*",
        type=int, help="Number of approximate cpu-min for *each* config.")

    if len(sys.argv) == 1:
        arguments = ["--dryrun", "--verbose",
                     "weak_scaling/flow", "-NCf", "10"]
        # arguments.remove("--dryrun")

        weak_flow_args = [
            "-v", "weak_scaling/flow", "-rn_app", "_flow", "-NCf", "0",
            "-NCor", "0", "-NUp", "0", "-NFlows", "1000", "-chr",
            "2", "2", "2", "3", "3", "4", "4", "6", "6", "8", "8"]

        weak_cfg_gen_args = [
            "-v", "weak_scaling/cfg_gen", "-rn_app", "_gen", "-NCf", "1",
            "-NCor", "600", "-NUp", "30", "-NFlows", "0", "-NTh", "1", "-sc", "-chr",
            "2", "2", "2", "3", "3", "4", "4", "6", "6", "8", "8"]

        weak_io_args = [
            "-v", "weak_scaling/io", "-rn_app", "_io", "-NCf", "10", "-NTh", 
            "0", "-NCor", "1", "-NUp", "1", "-NFlows", "0", "-sc", "-chr",
            "2", "2", "2", "3", "3", "4", "4", "6", "6", "8", "8"]

        strong_flow_args = [
            "-v", "strong_scaling/flow", "-rn_app", "_flow", "-NCf", "0",
            "-NCor", "0", "-NUp", "0", "-NFlows", "1000", "-chr",
            "64", "64", "32", "16", "8", "4", "2", "1", "1"]

        strong_cfg_gen_args = [
            "-v", "strong_scaling/cfg_gen", "-rn_app", "_gen", "-NCf", "1",
            "-NCor", "600", "-NUp", "30", "-NTh", "0", "-chr",
            "64", "64", "32", "16", "8", "4", "2", "1", "1"]

        strong_io_args = [
            "-v", "strong_scaling/io", "-rn_app", "_io", "-NCf", "10",
            "-NCor", "1", "-NUp", "1", "-NFlows", "0", "-sc", "-chr",
            "64", "64", "32", "16", "8", "4", "2", "1", "1"]

        arguments = weak_flow_args
        arguments = weak_cfg_gen_args
        # arguments = weak_io_args
        # arguments = strong_flow_args
        # arguments = strong_cfg_gen_args
        # arguments = strong_io_args

        # arguments.insert(0, "--dryrun")
        # arguments.remove("-v")

        args = parser.parse_args(arguments)
    else:
        args = parser.parse_args()

    dryrun = args.dryrun
    verbose = args.verbose

    dict_folder = args.folder
    files = [f for f in os.listdir(dict_folder) if f.endswith(".py")]
    values_to_replace = {}
    append_to_values = None
    modify_values = None

    if not isinstance(args.runname_append, type(None)):
        append_to_values = "runName"
        new_runnames = ["%s" %
                        args.runname_append for i in range(len(files))]
        values_to_replace["runName"] = new_runnames
    if not isinstance(args.runname_replace, type(None)):
        new_runnames = ["%s" %
                        args.runname_replace for i in range(len(files))]
    if not isinstance(args.NConfigs, type(None)):
        values_to_replace["NCf"] = args.NConfigs
    if not isinstance(args.NCor, type(None)):
        values_to_replace["NCor"] = args.NCor
    if not isinstance(args.NUpdates, type(None)):
        values_to_replace["NUpdates"] = args.NUpdates
    if not isinstance(args.NTherm, type(None)):
        values_to_replace["NTherm"] = args.NTherm
    if not isinstance(args.NFlows, type(None)):
        values_to_replace["NFlows"] = args.NFlows
    if not isinstance(args.storeCfgs, type(None)):
        values_to_replace["storeCfgs"] = args.storeCfgs
    if not isinstance(args.cpu_approx_runtime_hr, type(None)):
        values_to_replace["cpu_approx_runtime_hr"] = \
            args.cpu_approx_runtime_hr
    if not isinstance(args.cpu_approx_runtime_min, type(None)):
        values_to_replace["cpu_approx_runtime_min"] = \
            args.cpu_approx_runtime_min


    change_config_folder(
        dict_folder, values_to_replace,
        append_to_values=append_to_values,
        modify_values=modify_values, dryrun=dryrun,
        verbose=verbose)

    # scaling_types = ["weak", "strong"]

    # weak_flow_args = [
    #     "-v", "weak_scaling/flow", "-rn_app", "_flow", "-NCf", "0",
    #     "-NCor", "0", "-NUp", "0", "-NFlows", "0"]

    # for scaling in scaling_types:

    #     dict_folder = "%s_scaling/flow" % scaling
    #     files = [f for f in os.listdir(dict_folder) if f.endswith(".py")]
    #     append_to_values = ["runName"]
    #     modify_values = None
    #     value_to_replace = {
    #         "runName": ["_flow" for i in range(len(files))],
    #         "NCf": 0,
    #         "NCor": 0,
    #         "NUp": 0,
    #         "NFlows": 1000,
    #         "storeCfgs": False,
    #         # "cpu_approx_runtime_hr": [],
    #     }
    #     change_config_folder(dict_folder, value_to_replace,
    #         append_to_values=append_to_values,
    #         modify_values=modify_values,
    #         dryrun=dryrun)

    # weak_cfg_gen_args = [
    #     "-v", "weak_scaling/cfg_gen", "-rn_app", "_flow", "-NCf", "1",
    #     "-NCor", "600", "-NUp", "30", "-NFlows", "1000", "-chr",
    #     "20", "12", "6", "4", "3", "2", "2", "2", "2", "2", "2"]

    #     dict_folder = "%s_scaling/cfg_gen" % scaling
    #     files = [f for f in os.listdir(dict_folder) if f.endswith(".py")]
    #     append_to_values = ["runName"]
    #     modify_values = None
    #     value_to_replace = {
    #         "runName": ["_cfg_gen" for i in range(len(files))],
    #         "NCf": 1,
    #         "NCor": 600,
    #         "NUp": 30,
    #         "NFlows": 1000,
    #         "storeCfgs": True,
    #         "cpu_approx_runtime_hr": [20, 12, 6, 4, 3, 2, 2, 2, 2, 2, 2],
    #     }
    #     change_config_folder(dict_folder, value_to_replace,
    #         append_to_values=append_to_values,
    #         modify_values=modify_values,
    #         dryrun=dryrun)

    # weak_io_args = [
    #     "-v", "weak_scaling/io", "-rn_app", "_io", "-NCf", "100",
    #     "-NCor", "1", "-NUp", "1", "-NFlows", "0", "-chr", "-sc",
    #     "20", "12", "6", "4", "3", "2", "2", "2", "2", "2", "2"]

    # exit("Exiting after run for %s scaling." % scaling)


if __name__ == '__main__':
    main()
