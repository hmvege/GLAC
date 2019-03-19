# GLAC - Gluonic Action

Lattice Quantum Chromo Dynamics program for generating pure gauge field configurations and flowing them. Created by Hans Mathias Mamen Vege at UiO(University of Oslo).

![action density of the field](docs/field_density_b62_b6.200000_N32_NT64_np512_config00600.gif)
*The action density of a gauge field.*

#### Program folder structure
* `GLAC` contains main program. 
* `config_folder` contains parameter-files for different runs.
* 
```
GLAC
config_folder
docs
python_scripts
JobRenamer.py 
createJobs.py
multiJobSetup.py
```

#### Installation
Compile with following libraries:
```
GNU/4.9
OpenMPI/1.10.0
Qt/5.6.2
```
To compile, run the following inside the `build` folder:
```
qmake ../GLAC/GLAC.pro
make
```

#### A guide to GLAC
To run GLAC through its Python interface, you will for the moment being have to use python 2.7.

Example:
```
> python createJobs.py load weinberg_config.py -s slurm -bf {base directory}
```
This will automatically create a job script and submit it to the cluster the script is running on.

`-s {system}` specifies the run system. Either slurm, torque or local.

`-bf {base_directory}` is the full path leading up to where to output directory is situated. This would be the scratch space or wherever you want to store the output data. `createJobs.py` will automatically setup all relevant folders in a following structure,
```
{base directory}/
|--- output/
    |--- {run-name}/
        |--- field_configurations/
            |-- {field configurations generated}.bin
        |--- flow_observables/
            |--- {observables}/
                |--- {observable output files for each flow time}.dat
        |--- observables/
            |--- {non-flowed observable data}.dat
        |--- scalar_fields/
            |--- {observables}/
                |--- {the observables for the entire lattice}.bin
|--- input/
    |--- {run-name}.json
```

If you need help, type -h, and you will get additional command line arguments.
```
> python createJobs.py -h
```

If you want to perform a dryrun(i.e. no permanent change is performed), run:
```
> python createJobs.py --dryrun load weinberg_config.py -bf ${base directory}
```
This will then print out every detail of the run, without submitting the job. Usefull for double check the settings without actually submitting.
