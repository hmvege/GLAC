# Run examples for GLAC

This folder contains instructions for how to perform runs with GLAC. We will specify a run example for each of the different types of runs,
* `setup`
* `load`
* `utest`
* `perf_test`
* `field_density`

as well as a showing how to load an example configuration from Chroma to compare the output with. 

First, a quick note of the three types of systems we can run on,
* `local` generates run configuration for running on local computer. It does _not_ initiate a run. This will have to be done using e.g. ```mpirun -n {numprocs} build/GLAC input/{cfg-name}.json``` where `numprocs` is the number of processors, and `cfg-name` is the name of the `.json` configuration to pass to GLAC.
* `slurm` generates run configuration and initiates run for a system running the Slurm manager.
* `torque` generates run configuration and initiates run for a system running the Torque manager.


**NOTE:** the `sbatch` option in `createJobs.py` currently interacts with the Slurm system, but one can still see submitted jobs and their initiation parameters using the commands, as these is stored in a `.ids.json` file.

**NOTE:** if one specifies to load configurations with `-lcfg`, no additional configurations will be generated. In order to generate configurations, pass a single configuration file path to `-lcfgr` (load configuration and run) and specify how many configuration to generate with `-NCfgs`.




## Example: `setup`
The `setup` allows the user to setup a run while specifying all of the available parameters.

#### Setup example 1: generating and flowing
Example on how to start a local run generating configurations.
```
python2 createJobs.py setup local 8 -N 8 -NT 16 -fobs topct -NFlows 1000 -NCfgs 100 -NTh 100000 -fEps 0.01 -b 6.0 -rn GLACTestRun
mpirun -n 8 build/GLAC input/config_GLACTestRun.json
```

Parameters of run:
* `local` specifies that we will run on local computer. The next line is the run command.
* `8` specifies number of cores to utilize, which is needed in order to verify that the lattice size is valid.
* `-N 8` spatial lattice size of 16.
* `-NT 16` temporal lattice size of 32.
* `-fobs topct` specifies that we will use the topological charge sampled in euclidean time as our gauge observable in flow time.
* `-NFlows 1000` specifies GLAC to flow the configurations 1000 times.
* `-NCfgs 100` specifies GLAC to generate 1000 configurations.
* `-NTh 100000` specifies GLAC to use 100000 for thermalization.
* `-fEps 0.01` specifies the flow integration step to be 0.1.
* `-b 6.0` specifies the beta value.
* `-rn GLACTestRun` specifies the run name. Files will be stored with `GLACTestRun`in the name.

#### Setup example 2: loading configurations and flowing
An example of loading a folder of configurations and continuing flowing them.
```
python2 createJobs.py setup local 8 -N 8 -NT 16 -lcfg input/configs_to_run -fobs topct -NFlows 1000 -fEps 0.01 -b 6.0 -rn GLACLoadConfigsAndFlow
mpirun -n 8 build/GLAC input/config_GLACLoadConfigsAndFlow.json
```

Parameters of run:
* `local` specifies that we will run on local computer. The next line is the run command.
* `8` specifies number of cores to utilize, which is needed in order to verify that the lattice size is valid.
* `-N 8` spatial lattice size of 16.
* `-NT 16` temporal lattice size of 32.
* `-fobs topct` specifies that we will use the topological charge sampled in euclidean time as our gauge observable in flow time.
* `-NFlows 1000` specifies GLAC to flow the configurations 1000 times.
* `-fEps 0.01` specifies the flow integration step to be 0.1.
* `-b 6.0` specifies the beta value.
* `-rn GLACLoadConfigsAndFlow` specifies the run name. Files will be stored with `GLACTestRun`in the name.

#### Setup example 3: loading a configuration and continue generating gauge configurations
An example of loading a single configurations and continue to generate new gauge configurations from that configuration.




## Example: `load`
`load` allows the user to pass a simplified `.json` configuration(see [`config_folder`](https://github.com/hmvege/GLAC/tree/master/config_folder) for examples) for starting a run, while also being able to override run parameters.

#### Load example 1: Loading a simplified `.json` configuration for Slurm
An example of loading a `.json` configuration which initiates a job on Slurm.
```
python2 createJobs.py --dryrun load config_folder/weak_scaling/cfg_gen/weak_scaling_np2.py -lcfg output/testRunLocal/field_configurations/testRunLocal_beta6.000000_spatial8_temporal16_threads8_config00000.bin 
```

#### Load example 2: Loading a simplified `.json` configuration for Slurm and continue generating gauge configurations:
An example of loading a `.json` configuration which initiates a job on Slurm, and then continuing to generate gauge configurations.
```
python createJobs.py load config_folder/beta6_2.py -s abel -lcfg /work/users/hmvege/output/prodRunBeta6_2/field_configurations/ -lhr 25 -lmin 0 -bf /work/users/hmvege/ -NCfg 0 -cfgnum 395
```

python createJobs.py load config_folder/NUpdatesScaling/B60_NUP20_NCORR400.py -s abel -lcfgr output/B60_AC_SCALING/field_configurations/B60_AC_SCALING_b6.000000_N16_NT32_np512_config00000.bin -lhr 2 -lmin 0 -NCfgs 200




## Example: `utest`
The `setup` allows the user to setup a run while specifying all of the available parameters.

```
python2 createJobs.py  -v utest local 8 -vr -cgi input/activeProcTest_beta6.000000_spatial8_temporal16_threads4_config0.bin -N 8 -NT 16 -sq
mpirun -n 8 build/GLAC input/config_defaultTestRun.json
```




## Example: `perf_test`
The `setup` allows the user to setup a run while specifying all of the available parameters.

```
python2 createJobs.py -v perf_test local 8 -NExpTests 100 -TaylorPolDegree 16 -NDerivativeTests 10 -NRandTests 1000000 
```




## Example: `field_density`
The `setup` allows the user to setup a run while specifying all of the available parameters.

```
python2 createJobs.py --dryrun -v field_density -s local -nt 8 -lcfg testRunLocal/testRunLocal_beta6.000000_spatial8_temporal16_threads8_config00000.bin -rn lattice-density-test1 -N 8 -NT 16 -b 6.0 -sd 8 8 4 4 -NFlows 1000 -NCfgs 0
```




## Example: Chroma
Some LQCD programs saves the configurations with reversed byte order. An example of such a program is Chroma. To load and run configurations from Chroma, one can do the following,
```
python createJobs.py setup local 8 -rn chromaTest -N 16 -NT 32 -b 6.0 -lcfg input/chroma_config -fobs topct -obs topct -chroma -NFlows 10 -fEps 0.1
mpirun -n 8 build/GLAC input/config_chromaTest.json
```
Parameters of run:
* `local` specifies that we will run on local computer. The next line is the run command.
* `8` specifies number of cores to utilize, which is needed in order to verify that the lattice size is valid.
* `-rn chromaTest` specifies the run name. Files will be stored with `chromaTest`in the name.
* `-N 16` spatial lattice size of 16.
* `-NT 32` temporal lattice size of 32.
* `-b 6.0` specifies the beta value.
* `-lcfg input/chroma_config` is the folder containing the input configuration. Any `.bin` files in this folder will be passed to GLAC and loaded assuming it is a field configuration.
* `-fobs topct` specifies that we will use the topological charge sampled in euclidean time as our gauge observable in flow time.
* `-obs topct` specifies that we will use the topological charge sampled in euclidean time as our gauge observable for the zero flow time observable.
* The `-chroma` flag is passed on to GLAC in the `config_chromaTest.json` configuration and tells GLAC to load the configuration in reversed byte order. 
* `-NFlows 10` specifies GLAC to flow the configurations 10 times.
* `-fEps 0.1` specifies the flow integration step to be 0.1.