import os

dryrun = False
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

                # First open file
                # change False -> false
                # change True -> true
                # change # -> //
                with open(file_path, "r") as _f:
                    file_contents = _f.read()

                file_contents = file_contents.replace("True","true")
                file_contents = file_contents.replace("False","false")
                file_contents = file_contents.replace("#","//")

                with open(file_path, "w") as _f:
                    _f.write(file_contents)

                os.rename(file_path, new_file_path)

            print("> mv {} {}".format(file_path, new_file_path))

print("{} files renamed to .json format.".format(count))