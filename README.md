# GLAC - Gluonic Action

Lattice Quantum Chromo Dynamics program for generating pure gauge field configurations and applying gradient flow on them. Created by Hans Mathias Mamen Vege at UiO(University of Oslo).

<p align="center">
    <img src="docs/field_density_b62_b6.200000_N32_NT64_np512_config00600.gif" alt="Action density of the gauge field" width="600"/>
</p>
    
<p align="center">
    <i>The action density of a gauge field. Created using <a href="https://github.com/hmvege/LatViz">LatViz</a>, a program developed by Mathias M. Vege and Giovanni Pederiva.</i>
</p>

#### Program folder structure
* `GLAC` contains main program. 
* `config_folder` contains parameter-files for different runs.
* `docs` contains the documentation of GLAC.
* `scripts` contains various scripts, mostly written in Python, which has been used during development.
* `JobRenamer.py` is a script for renaming Slurm `.out` output files.
* `multiJobSetup.py` is script for starting multiple jobs.
* `createJobs.py` is the main script for generating `.json` configuration files for GLAC, as well as submitting files to either Torque or Slurm.

#### Installation
Compile with following libraries(these should be available at most HPC cluster, should run without any problems):
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

#### A (short) guide to GLAC
To run GLAC through its Python interface, you will for the moment being have to use python 2.7.

Example:
```
python createJobs.py load weinberg_config.py -s slurm -bf {base directory}
```
This will automatically create a job script and submit it to the cluster the script is running on.

`-s {system}` specifies the run system. Either Slurm, Torque or local.

`-bf {base_directory}` is the full path leading up to where to output directory is situated. This would be the scratch space or wherever you want to store the output data. `createJobs.py` will automatically setup all relevant folders in a following structure,
```
{base directory}/
├── output/
│   └── {run-name}/
│   │   ├── field_configurations/
│   │   │   └── {field configurations generated}.bin
│   │   ├── flow_observables/
│   │   │   └── {observables}/
│   │   │   │   └── {observable output files for each flow time}.dat
│   │   ├── observables/
│   │   │   └── {non-flowed observable data}.dat
│   │   └── scalar_fields/
│   │   │   └── {observables}/
│   │   │   │   └── {the observables for the entire lattice}.bin
└── input/
│   └── {run-name}.json
```

If you need help, type -h, and you will get additional command line arguments.
```
python createJobs.py -h
```

If you want to perform a dryrun(i.e. no permanent change is performed), run:
```
python createJobs.py --dryrun load weinberg_config.py -bf {base directory}
```
This will then print out every detail of the run, without submitting the job. Useful for double check the settings without actually submitting.

For a more extensive guide to GLAC, see the [examples folder](https://github.com/hmvege/GLAC/tree/master/examples). For more implementation specific details, see the [GLAC Documentation](http://hmvege.github.io/GLAC/html/index.html).


#### Additional resources
* The analysis code used for this thesis program can be seen in [LatticeAnalyser](https://github.com/hmvege/LatticeAnalyser).
* To generate visualizations of the topological charge or energy, we used a program called [LatViz](https://github.com/hmvege/LatViz).
* A list of the Slurm commands can be found [here](https://slurm.schedmd.com/pdfs/summary.pdf). A list of the Torque commands can be found [here](https://gif.biotech.iastate.edu/torque-pbs-job-management-cheat-sheet). A Rosetta stone for going between Slurm and Torque can be seen [here](https://slurm.schedmd.com/rosetta.pdf).
* To read `.json`-files, we make use of the highly recommended json reader provided [here](https://github.com/nlohmann/json).

## Todo list and future goals.
Aspirational goals and todo's that could/should be implemented.

#### Todo's
* [ ] Add check in `createJobs.py` such that when loading a single gauge configurations and generating new ones, we make sure that no previous gauge configurations of the names we will generates exist in the provided folder. Do this perhaps for observables as well.
* [ ] Add unit tests for the Action derivatives
* [ ] Add unit tests for the observables
* [ ] Add unit tests for remaining lattice methods and make sure they match their equivalent SU3 unit tests or similar.
* [ ] Add unit tests for the SU3 matrix exponential methods.
* [ ] Add performance timing for the flow method.
* [ ] Investigate ways of allocating the cubes for sharing at the beginning of a program. Will need to have the performance timing method implemented to check if this helps.

#### Future goals
* [ ] Implement the Luscher-Weisz action.
* [ ] Implement the option of calculating the field strength tensor with the plaquette.
* [ ] Implement O(a^4) error correction for the field strength tensor.
* [ ] Implement O(a^6) error correction for the field strength tensor.
* [ ] Implement a better method for selecting the observables to run for.
* [ ] Clean up the observables selection. Currently quite confusing names for the observables as well as a convoluted way of selecting observables.