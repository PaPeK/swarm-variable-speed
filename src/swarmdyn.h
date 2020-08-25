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
void Step(int s, std::vector<particle> &a, params *, std::vector<predator> &preds,
          bool dummy=false);      // numerical step
// fctns. for Output:
void Output(int s, std::vector<particle> &a, params &SP, std::vector<predator> &pred,
            std::vector<particle> &d, std::vector<predator> &predD, bool forceSave=false);
void ErrorOutput(std::vector<particle> &a, int err, params *ptrSP, predator *pred);
std::vector<double> Out_swarm(std::vector<particle> &a, params &SP);
std::vector<double> Out_swarm_pred(std::vector<particle> &a, predator &pred, params &SP);
void Out_end_Save(std::vector<particle> &a, std::string name, params &SP);
void Out_End_PosVel(std::vector<particle> &a, predator &P, std::string name, params &SP);
void Out_start_Save(std::vector<particle> &a, std::string name, params &SP);
void Out_startPred_Save(std::vector<predator> &p, std::string name, params &SP);
void DataCreateSaveWrite(std::vector< std::vector<double> > &data,
                     std::vector<double> &out, params &SP,
                     std::string name, bool forceSave=false);
// fctns. for collecting means
template<class T>
std::vector<double> basic_particle_averages(std::vector<particle> &a, std::vector<T> &nodes);
template<class T>
double get_avg_pred_dist(std::vector<particle> &a,
                         std::vector<T> &nodes, predator &pred,
                         params &SP);
template<class T>
double get_avg_deg(std::vector<particle> &a, std::vector<T> &nodes);
void Write_out(std::vector<double> &out, params &SP,
               std::vector<hsize_t> & vec_dim,
               H5::DataSet *h5dset, std::string name);
std::vector<double> corr_velFlucTS_interact(std::vector<particle> &a,
                                            std::vector<double> &mean_v,
                                            std::vector<unsigned int> Fnn,
                                            std::vector<unsigned int> Fnn2);
std::vector<double> corr_velFluc_interact(std::vector<particle> &a,
                                          std::vector<double> &mean_v,
                                          std::vector<unsigned int> Fnn,
                                          std::vector<unsigned int> Fnn2);
double corr_velFlucDirection_interact_all(std::vector<particle> &a,
                                          std::vector<double> &mean_v);
double corr_velFluc_interact_all(std::vector<particle> &a,
                                 std::vector<double> &mean_v);
// will be initialized for particle and predator
template<class part>
void WriteParticles(std::vector<part> &a, params &SP, 
                    std::string name, double outstep, int Nori);
