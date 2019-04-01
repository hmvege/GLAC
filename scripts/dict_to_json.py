import os

dryrun = True
if dryrun:
    print("**** DRYRUN ****")

path_configs_to_convert = "../config_folder"

count = 0

# Walks through config_folder and converts every python file to a json file
for w in os.walk(path_configs_to_convert):
    
    for f in w[-1]:
        if not f.startswith("."):
            folder_path = w[0]
            file_path = os.path.join(w[0], f)

            # Splits file path and .py extension
            file_body, file_ext = os.path.splitext(file_path)
            if file_ext != ".py":
                continue
            new_file_path = file_body + ".json"

            if not dryrun:
                count += 1
                os.rename(file_path, new_file_path)
            print("> mv {} {}".format(file_path, new_file_path))

print("{} files renamed to .json format.".format(count))