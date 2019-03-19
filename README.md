# GLAC - Gluonic Action

Lattice Quantum Chromo Dynamics program for generating pure gauge field configurations and applying gradient flow on them. Created by Hans Mathias Mamen Vege at UiO(University of Oslo).

<p align="center">
    <img src="docs/field_density_b62_b6.200000_N32_NT64_np512_config00600.gif" alt="Action density of the gauge field" width="600"/>
</p>

<center>
    *The action density of a gauge field. Created using [LatViz](https://github.com/hmvege/LatViz)*
</center>


#### Program folder structure
* `GLAC` contains main program. 
* `config_folder` contains parameter-files for different runs.
* `docs` contains the documentation of GLAC.
* `python_scripts` contains various Python scripts used during development.
* `JobRenamer.py` is a script for renaming Slurm `.out` output files.
* `multiJobSetup.py` is script for starting multiple jobs.
* `createJobs.py` is the main script for generating `.json` configuration files for GLAC, as well as submitting files to either Torque or Slurm.

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

`-s {system}` specifies the run system. Either Slurm, Torque or local.

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
This will then print out every detail of the run, without submitting the job. Useful for double check the settings without actually submitting.

#### Additional resources
* The analysis code used for this thesis program can be seen in [LatticeAnalyser](https://github.com/hmvege/LatticeAnalyser).
* To generate visualizations of the topological charge or energy, we used a program called [LatViz](https://github.com/hmvege/LatViz).
* A list of the Slurm commands can be found [here](https://slurm.schedmd.com/pdfs/summary.pdf). A list of the Torque commands can be found [here](https://gif.biotech.iastate.edu/torque-pbs-job-management-cheat-sheet). A Rosetta stone for going between Slurm and Torque can be seen [here](https://slurm.schedmd.com/rosetta.pdf).