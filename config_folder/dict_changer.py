import os
import ast
import re
import json
import types
import shutil
import fileinput
from collections import OrderedDict

def change_config_folder(folder, replace_values, append_to_values=None, 
    modify_values=None, dryrun=False):
    
    if dryrun:
        print("**** DRYRUN ****")

    assert isinstance(replace_values, dict), ("Please provide dict with "
        "values to replace in config folder %s" % folder)

    _temp = [f for f in os.listdir(folder) if f.endswith(".py")]
    files = sorted(_temp, key=lambda s: int(re.findall(r"\d+", s)[0]))

    for i, f in enumerate(files):
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

        # Modifies dicitonary values
        for item in replace_values.items():
            key, val = item


            if isinstance(val, list):

                if len(val) == len(files):
                    if not isinstance(append_to_values, type(None)) and \
                        key in append_to_values:

                        # Checks if we have already appended string on previous
                        # run.
                        if val[i] in data_dict[key]:
                            continue
                        else:
                            data_dict[key] += str(val[i])

                    elif not isinstance(modify_values, type(None)):
                        data_dict[key] += modify_values(val[i])

                    else:
                        data_dict[key] += val[i]
            else:
                data_dict[key] = val


        # Reorganized dictionary values
        ord_data_dict = OrderedDict([(key, data_dict[key]) for key in ordered_keys])

        json_fpath = os.path.join(folder, os.path.splitext(f)[0] + ".py")

        # Writes json to file        
        print("Writing python configuration file at location {0:<s}".format(json_fpath))
        if not dryrun:
            with open(json_fpath, "w+") as json_file:
                json.dump(ord_data_dict, json_file, indent=4)

            # Re-names false->False, true->True
            f1 = open(json_fpath, "r")
            f2_path = os.path.join(folder, os.path.splitext(f)[0] + "_temp2.py")
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

        # print(json.dumps(ord_data_dict, indent=4, separators=(", ", ": ")), "\n")
        # exit("Exits after changing file %s" % f)

    print("Changed config files in folder %s." % folder)

    # exit("Exits after folder %s" % folder)

def main():
    dryrun = False

    scaling_types = ["weak", "strong"]

    for scaling in scaling_types:

        dict_folder = "%s_scaling/flow" % scaling
        files = [f for f in os.listdir(dict_folder) if f.endswith(".py")]
        append_to_values = ["runName"]
        modify_values = None
        value_to_replace = {
            "runName": ["_flow" for i in range(len(files))],
            "NCf": 0,
            "NCor": 0,
            "NUp": 0,
            "NFlows": 1000,
            "storeCfgs": False,
            # "cpu_approx_runtime_hr": [],
        }
        change_config_folder(dict_folder, value_to_replace, 
            append_to_values=append_to_values, 
            modify_values=modify_values,
            dryrun=dryrun)


        dict_folder = "%s_scaling/cfg_gen" % scaling
        files = [f for f in os.listdir(dict_folder) if f.endswith(".py")]
        append_to_values = ["runName"]
        modify_values = None
        value_to_replace = {
            "runName": ["_cfg_gen" for i in range(len(files))],
            "NCf": 1,
            "NCor": 600,
            "NUp": 30,
            "NFlows": 1000,
            "storeCfgs": True,
            "cpu_approx_runtime_hr": [20, 12, 6, 4, 3, 2, 2, 2, 2, 2, 2],
        }
        change_config_folder(dict_folder, value_to_replace, 
            append_to_values=append_to_values, 
            modify_values=modify_values,
            dryrun=dryrun)


        # dict_folder = "%s_scaling/full" % scaling
        # files = [f for f in os.listdir(dict_folder) if f.endswith(".py")]
        # append_to_values = ["runName"]
        # modify_values = None
        # value_to_replace = {
        #     "runName": ["_full" for i in range(len(files))],
        #     "NCf": 5,
        #     "NCor": 600,
        #     "NUp": 30,
        #     "NFlows": 1000,
        #     "storeCfgs": True,
        #     "cpu_approx_runtime_hr": [20, 12, 6, 4, 3, 2, 2, 2, 2, 2, 2],
        # }
        # change_config_folder(dict_folder, value_to_replace, 
        #     append_to_values=append_to_values, 
        #     modify_values=modify_values,
        #     dryrun=dryrun)

        exit("Exiting after run for %s scaling." % scaling)

if __name__ == '__main__':
    main()