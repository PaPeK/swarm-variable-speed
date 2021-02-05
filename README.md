# PredatorPrey – simulation code for agent-based predator-prey interactions
Version 0.1 described and applied in <http://arxiv.org/abs/2009.02079>


(Python-wrapper package calling c++ running)

[Github](https://github.com/PaPeK/predatorPrey)

## General Notes

PredatorPrey is a C++ code which runs numerical simulations of the agent based stochastic differential equations.
The motivation of this code are group living animals, e.g. fish shoals in shallow water(2 dimensions).
In principle the code can be easily extended to three dimensions. However, the higher the dimension the smaller the volume if the number of agents stays constant.
For a full mathematical description of the stochastic differential equations and the social forces (repulsion, attraction and alignment) see <http://arxiv.org/abs/2009.02079>. 

## Required C++ libraries:

next to standard C libraries the following libraries need to be installed via ``apt-get install``'

- libhdf5-serial-dev
- libgsl-dev
- libcgal-dev

### Docker-alternative

The ```Dockerfile``` contains the build instruction for an image
which can be used to run the c++ code. In this case the above libraries do not need to be installed.
For building the docker:

- install docker following https://docs.docker.com/engine/install/.
- build the image by running 
```
docker build -t gcc_docker .
```
    in the directory of this git-repository
    - note that `gcc\_docker` is an arbitrary name for the image and can be chose freely
- check if the image exists by running `docker images` ("gcc_docker" should be listed)
- change in "RunSingle.py" the line `dockerName = None` to `dockerName = 'gcc_docker'`

## Required python packages

The code runs in the anaconda environment specified in `environment.yml` which can be created from the projects root directory via
```shell
conda env create -f environment.yml
```
which creates the environment `animateSwarm` (name specified in `environment.yml`) which you activate by
```shell
conda activate animateSwarm
```
In case you want to update your python packages to enable the code to run the most important packages are:

- numpy
- h5py
- matplotlib
- pathlib

## Compilation(MacOS/Linux):

On MacOS/Linux system:
The included Makefile is for compilation of the source code (see tutorial: http://mrbook.org/blog/tutorials/make/). Just run 

```
make
``` 

in the terminal in the directory of this repository.

#### Possible compilation problem

If you are using anaconda and h5py is installed the linking to the libhdf5-serial-dev library might not work.
A simple work around is to comment in your `~/.bashrc` (or `~/.shrc` or `~/.zshrc` depending on the shell you are using) the anaconda initialization.
Than you start a new terminal and retry the compilation (run `make`).
If the compilation was succesfull you can uncomment the anaconda initialization again and run the python-scripts.

### Docker-alternative

run in the directory of this repository (assuming you followed the instructions above and called the docker image "gcc_docker" )

```
docker run --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp gcc_docker make
```

## Running

After successful compilation the executable file swarmdyn will be generated.
The executable is meant to be called via running:

```
python RunSingle.py
```


- RunSingle.py - starting a single simulation of swarmdyn and defining parameters.

RunSingle.py will also call a visualization of the simulations by calling:

- AnimateRun.py - animating the results using pythons matplotlib

### Docker-alternative

Remember to change in "RunSingle.py" the line `dockerName = None` to `dockerName = 'gcc_docker'` (or the docker-name you have chosen)

## User Agreement

By downloading PredatorPrey you agree with the following points: PredatorPrey is provided without any warranty or conditions of any kind. We assume no responsibility for errors or omissions in the results and interpretations following from application of PredatorPrey.

## License

Copyright (C) 2016-2020 Pascal Klamser, Pawel Romanczuk

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Implementation details

The code consists of a `c++` part (costly numerical simulations) which is accessed via `python`-wrapper scripts.

## C++ simulation

The c++ code is in the `src`-folder and split in multiple logic (more or less) parts.
The splitting in different files and headers is usefull for developing, because it strongly reduces compilation time if only part of the code is changed.
The following parts exists (sorted according their "conceptual importance"):

* main (main loop, output):
    * `swarmdyn.cpp`
* agent- and system-parameter: definitions, initialization and parsing:
    * `agents.cpp/.h` (definition)
    * `settings.cpp` (initialization)
    * `input_output.cpp` (parsing)
* agents: individual dynamics, interactions(network), social forces, operations
    * `agents_dynamics.cpp`
    * `agents_interact.cpp`
    * `social_forces.cpp`
    * `agents_operation.cpp`
* tools: hdf5 writing, maths, raycasting (for visual field ... actually not needed):
    * `h5tools.cpp`
    * `mathtools.cpp`
    * raycasting (based on Colin Twomey raycasting-repository: fovea link??):
        * `fov.cpp`
        * `geometry.cpp`
        * `ray_casting.cpp`

### Hard coded peculiarities

* if `fileID = 'xx'` the hdf5-file _out_xx.h5_ is created by the c++ code in the current directory.
    * otherwise the c++ code expects that an hdf5-file exist with the name _out_fileID.h5_ 
* if a anything in the header `src/agents.h` is changed, the whole code needs a recompilation 
    * do `rm -f obj/*`          (removes all obj and forces a recompilation)
    * run `make`
    * Todo: maybe place the class definitions in `agents.cpp` to prevent this

## Python wrapper


### Single Run (RunSingle.py)

#### Animation (AnimateRun.py)

* The code can be animated if the dataset `part` (for partiticles) exists via `python AnimateRun.py` assuming that the file `out_xx.h5` exists (or specify another file).
    * the script `AnimateRun.py` is based on the github gitrepository [animateswarm](https://github.com/PaPeK/animateswarm)

### Parameter Scans (Scan2D.py)
This script can call mutliple 2D-scans sequentially.
* Each scan is parallelized by running each point in parameter-space on a seperate process.
    * each process saves the output to a seperate hdf5-file (with the parameter-values in its name)
    * for each parameter-pair multiple samples can be computed and the scan can get extended
        * to extend: run the Scan2D.py with the same parameters but a larger value for _runs_
        * Problem: if _para_values0_ or _para_values0_ is changed the code checks only if already runs exist for the first parameter-pair -> if they exist, the total number of runs is reduced by the number of existing runs (no parameter)
    * if the number of parameter-pairs is smaller than the number of available processors: the samples are split and later reunited

#### Evaluation of Scans (EvaluateScan.py)

The script evaluates each parameter-point and returns the output (mean or std or more specific processing) in an array which whose rows and columns correspond to the respective paremter-points.
These arrays and a dictionary which links keywords to the entry in the last dimension of the array is saved in `Mean.pkl`.


