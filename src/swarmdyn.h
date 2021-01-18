/*  SwarmDynamics 
    Stochastic differential equation model of agents interacting via attraction/repulsion/alignment.
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser, Pawel Romanczuk

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
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <time.h>
#include <getopt.h> // for required and optional argument
#include <unistd.h>
#include <string.h>
#include <assert.h>

#ifdef __cplusplus

#include <chrono>
#include <list>
#include <vector>
#include <set>
#include <algorithm>    // to use std::set_difference, std::find, st::iota, std::shuffle
#include <iterator>     // to use std::begin
#include <utility>                              // for nearest neighbor search needed
// #include <limits>       // for accessing numerical limits
#include <fstream>      // for writing vector to file

#endif

// DIFFERERNT MODULES NEEDED TO RUN CODE------------------------------
// modules for agent structures as particle, predator, params, ...
#include "common_defines.h"
#include "agents.h"
// interactions of agents:
#include "agents_interact.h"
// dynamics of agents:
#include "agents_dynamics.h"
// operations on agents: find specific agents, get COM, ...
#include "agents_operation.h"
// has some math functions and operations in vectors and arrays
#include "mathtools.h"
// to write hdf5-files
#include <H5Cpp.h>      // for hdf5 output
#include "h5tools.h"
// settings: parameters, initialization, reset
#include "settings.h"
// input and outputs:
#include "input_output.h"

// FUNCTION DEFINITION
void InitRNG();             // initializes the random number generation
void Step(int s, std::vector<particle> &a, params *);      // numerical step
// fctns. for Output:
void Output(int s, std::vector<particle> &a, params &SP);
// will be initialized for particle and predator
template<class part>
void WriteParticles(std::vector<part> &a, params &SP, 
                    std::string name, double outstep, int Nori);
std::vector<double> Out_swarm(std::vector<particle> &a, params &SP);
double corr_velFluc_interact(std::vector<particle> &a,
                             std::vector<double> &mean_v);
template<class T, class MembType>
double get_varOFVecFluctuation(std::vector<particle> &a,
                               std::vector<double> &mean,
                               MembType vec,
                               std::vector<T> &nodes);
template<class T>
std::vector<double> basic_particle_averages(std::vector<particle> &a, std::vector<T> &nodes);
void DataCreateSaveWrite(std::vector< std::vector<double> > &data,
                     std::vector<double> &out, params &SP,
                     std::string name, bool forceSave=false);
