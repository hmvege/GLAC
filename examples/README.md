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




## Example: `setup`
The `setup` allows the user to setup a run while specifying all of the available parameters.

#### Setup example 1: generating gauge configurations
Example on how to start a local run generating configurations.
```
python2 createJobs.py setup local 8 -N 8 -NT 16 -NFlows 0 -NCfgs 1000 -NTh 100000 -b 6.0 -rn GLACConfigRun
mpirun -n 8 build/GLAC input/config_GLACConfigRun.json
```

Run parameters:
* `local` specifies that we will run on local computer. The next line is the run command.
* `8` specifies number of cores to utilize, which is needed in order to verify that the lattice size is valid.
* `-N 8` spatial lattice size of 16.
* `-NT 16` temporal lattice size of 32.
* `-NFlows 0` specifies GLAC to flow the configurations 0 times, thus only generating configurations.
* `-NCfgs 100` specifies GLAC to generate 1000 configurations.
* `-NTh 100000` specifies GLAC to use 100000 for thermalization.
* `-b 6.0` specifies the beta value.
* `-rn GLACConfigRun` specifies the run name. Files will be stored with `GLACConfigRun`in the name.


#### Setup example 2: generating and flowing
Example on how to start a local run generating and flowing configurations.
```
python2 createJobs.py setup local 8 -N 8 -NT 16 -fobs topct -NFlows 1000 -NCfgs 100 -NTh 100000 -fEps 0.01 -b 6.0 -rn GLACConfigGenAndFlow
mpirun -n 8 build/GLAC input/config_GLACConfigGenAndFlow.json
```

Run parameters:
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
* `-rn GLACConfigGenAndFlow` specifies the run name. Files will be stored with `GLACConfigGenAndFlow`in the name.


#### Setup example 3: loading configurations and flowing
In the case we want to flow the configurations we loaded in the first example, we can run the following,
```
python2 createJobs.py setup local 8 -N 8 -NT 16 -lcfg output/GLACConfigRun/field_configurations -fobs topct -NFlows 1000 -fEps 0.01 -b 6.0 -rn GLACLoadConfigsAndFlow
mpirun -n 8 build/GLAC input/config_GLACLoadConfigsAndFlow.json
```

Run parameters:
* `-lcfg input/configs_to_run` locates the folder `input/configs_to_run` and retrieves all of the `.bin` files in the folder, for then to pass them on to GLAC with the absolute file path in the configuration file `input/config_GLACLoadConfigsAndFlow.json`. No configurations will be generated since we are loading configurations.


#### Setup example 4: loading a configuration and continue generating gauge configurations
An example of loading a single configurations and continue to generate new gauge configurations from that configuration.

```
python2 createJobs.py setup local 8 -N 8 -NT 16 -lcfgr output/GLACConfigRun/field_configurations/GLACConfigRun_b6.000000_N8_N16_np8_config00099.bin -NFlows 0 -b 6.0 -rn GLACConfigRun -cfgnum 100
mpirun -n 8 build/GLAC input/config_GLACConfigRun.json
```

Run parameters:
* `-lcfgr output/GLACConfigRun/field_configurations/GLACConfigRun_b6.000000_N8_N16_np8_config00099.bin` will load the `.bin` configuration and continue generating configurations by updating the given configurations.
* `-cfgnum 100` specifies that the new configuration will by given a configuration number of 100, instead of the default 0. This is needed in order not to override any previous configurations. 
* `-rn GLACConfigRun` is the same run name as the first example, but we will not override any files due to specifying `-cfgnum 100`.




## Example: `load`
`load` allows the user to pass a simplified `.json` configuration(see [`config_folder`](https://github.com/hmvege/GLAC/tree/master/config_folder) for examples) for starting a run, while also being able to override run parameters.


#### Load example 1: Loading a simplified `.json` configuration for Slurm/Torque/local
An example of loading a `.json` configuration which initiates a job on Slurm/Torque.
```
python2 createJobs.py load examples/configExampleLoad/configExampleLoad.json -s {slurm | torque | local}
```
If one wants to run this locally, run
```
mpirun -n 2 build-release/GLAC input/config_configExampleLoad.json
```
The terminal output of this can be viewed in [`examples/test_config/configExampleLoad.out`](https://github.com/hmvege/GLAC/tree/master/examples/test_config/configExampleLoad.out). The configurations generated can be viewed in the same folder. The output can be viewed in [`examples/test_config_output`](https://github.com/hmvege/GLAC/tree/master/examples/test_config_output).


#### Load example 2: Flowing configurations on Slurm
An example of loading a `.json` configuration which initiates a job on Slurm.
```
python2 createJobs.py load examples/configExampleLoad/configExampleLoad.json -s {slurm | torque} -lhr 10 -lmin 0 -vr
```

Run parameters:
* `-s {slurm | torque}` submits job to either Slurm or Torque.
* `-lhr 10` run time human hours.
* `-lmin 0` run time human minutes. Added to number of hours.
* `-vr` specifies a verbose run. Additional information will be printed. Might be useful for systems such as Slurm where it can provide insight to how fare the run has progressed.


#### Load example 3: Loading a simplified `.json` configuration for Slurm and continue generating gauge configurations:
An example of loading a `.json` configuration which initiates a job on Slurm, and then continuing to generate gauge configurations.
```
python createJobs.py load config_folder/configExampleLoad -s slurm -lcfg /work/users/hmvege/output/configExampleLoad/field_configurations/ -lhr 25 -lmin 0 -bf /work/users/hmvege/ -NCfg 0 -cfgnum 395
```

Run parameters:
* `-lcfg /work/users/hmvege/output/configExampleLoad/field_configurations/` will load configurations from given folder. Use absolute path.
* `-bf /work/users/hmvege/` specifies base folder. The output folders ect. will all use this as their base folder. Useful if you need to write to a scratch space. Will not override the location of the input configurations.
* `-cfgnum 395` will specify the program to only load and flow configurations of a number higher than 395.


## Example: `utest`
The `utest` tells GLAC to perform unit, integration and validation tests. An example of the test output can be viewed in [`examples/test_config_output`](https://github.com/hmvege/GLAC/tree/master/examples/test_config_output).

```
python2 createJobs.py utest local 2 -cgi output/configExampleLoad/field_configurations/configExampleLoad_b6.000000_N4_NT8_np2_config00000.bin -N 4 -NT 8 -sq -vr
```

Run parameters:
* `-cgi output/configExampleLoad/field_configurations/configExampleLoad_b6.000000_N4_NT8_np2_config00000.bin` will load the provided configuration and check gauge invariance.
* `-N 4` and `-NT 8` specifies the lattice dimensions to check gauge invariance for. The same value is used for other tests, such as IO and lattice shift tests even if no gauge configuration is provided.
* `-vr` tells the program to print all tests being passed as well as some additional diagnostics.




## Example: `perf_test`
The `perf_test` tells GLAC to perform performance tests. An example is
```
python2 createJobs.py -v perf_test local 8 -NExpTests 100 -TaylorPolDegree 16 -NDerivativeTests 10 -NRandTests 1000000 
```





## Example: `field_density`
The `field_density` loads a configurations and flows it, storing the observable at every lattice point at a given frequency.

```
python2 createJobs.py -v field_density -s local -nt 8 -lcfg testRunLocal/testRunLocal_beta6.000000_spatial8_temporal16_threads8_config00000.bin -rn lattice-density-test1 -N 8 -NT 16 -b 6.0 -sd 8 8 4 4 -NFlows 1000 -NCfgs 0
```




## Example: Chroma
Examples on how to load a Chroma configuration. The data produced with the `topct` flow observable should match those files in `chroma_output` that starts with `_Qt`. The files starting with `_Wt` should match that of `weinberg` flow observable. Both should match down to the \~15 decimal, thus serving as good validation test.

#### Chroma load example
Some LQCD programs saves the configurations with reversed byte order. An example of such a program is Chroma. To load and run configurations from Chroma, one can do the following,
```
python createJobs.py setup local 8 -rn chromaTest -N 16 -NT 32 -b 6.0 -lcfg input/chroma_config -fobs topct -obs topct -chroma -NFlows 10 -fEps 0.1
mpirun -n 8 build/GLAC input/config_chromaTest.json
```
Run parameters:
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




## Managing runs on Slurm and Torque
Following is a few HPC specific commands.

#### Setting the account name
On Slurm and Torque systems, additional arguments such as
```
--account_name {your-account-name}
```
will have to be passed.

#### Allocating CPU time
Further, one will also either need to specify the estimated computation time in the simplified `.json` parameter file, or pass arguments
```
-lhr {estimated-human-hours}
-lmin {estimated-human-minutes}
```
Note that these arguments are not in CPU hours.

#### Overriding default of 16 tasks per node.
Most systems have 16 cores per node, and thus `createJobs.py` will try to optimize for this when submitting the batch script. Passing
```
-igntsk
```
will override this requirement, which might be needed in some cases.

#### Specifying HPC partition
Passing the argument
```
-p {partition}
```
where `partition` is the name of the cluster parts can be used to specify which parts you want to run on. E.g. `-p a[1-8],b[1-8]`.

#### Ignoring certain nodes
To ignore certain nodes, pass
```
-ex {node-names}
```
with `node-names` as the name of the node to ignore. E.g. `-ex a[1,2,4]`.




## Useful tips and tricks

#### Checking the validity of a configuration
Passed 
```
--debug
```
will tell GLAC to check when loading a configuration if it contains `nan` values and quit if this is the case.

#### The `sbatch` option in `createJobs.py`
The `sbatch` option in `createJobs.py` currently interacts with the Slurm system, but one can still see submitted jobs and their initiation parameters using the commands, as these is stored in a `.ids.json` file.

#### The `-lcfg` option in `createJobs.py`
If one specifies to load configurations with `-lcfg`, no additional configurations will be generated. In order to generate configurations, pass a single configuration file path to `-lcfgr` (load configuration and run) and specify how many configuration to generate with `-NCfgs`.
