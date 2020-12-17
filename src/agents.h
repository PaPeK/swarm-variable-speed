/*  agents
    defines structures/classes of agents used in SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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

#ifndef agents_H
#define agents_H
#include "common_defines.h"
#include <H5Cpp.h>      // for hdf5 output
#include <set>
#include <vector>

struct particle{
    std::vector<double> x;            // agent position vector
    std::vector<double> v;            // agent velocity vector
    std::vector<double> u;            // agent direction unit vector
    double phi;             // agent direction polar angle [0,2Pi]
    double vproj;           // vel along heading direction
    std::vector<double> force;        // total social force vector
    std::vector<double> force_alg;    // alignment vector
    std::vector<double> force_rep;    // repulsion vector
    unsigned int id;        // ID of each agent (needed to split and merge agents correctly)
    std::vector<unsigned int> NN;    // vector containing all NN

    // counters of interaction partners - important only for local metric coupling (not global)
    int counter_rep;        // counter repulsion partners 
    // parameter-variations in between the agents
    double alg_strength;
    double rep_steepness;
    std::vector<double> out(void);  // for output
};
typedef struct particle particle;


// data structure for system parameters, and auxiliary variables
struct params{

    std::string location;   // path of output-files
    std::string fileID;     // name of output-files
    unsigned int N;                  // number of agents
    double sizeL;           // system size
    double Dphi;            // angular noise intensity  (= standard deviation)
    double noisep;          // auxiliary var. - std. deviation of the angular noise

    double alg_strength;    // alignment strength
    double rep_strength;    // repulsion strength
    double rep_range;       // repulsion range
    double rep_steepness;   // repulsion - steepness of the sigmoid function


    double sim_time;        // simulation time in natural time units
    double speed0;          // stationary speed of single individuals
    double dt;              // time step of integration
    double output;          // output time in natural units

    int BC;                 // switch for boundary conditions
    int IC;                 // switch for initial conditions

    double trans_time;      // transient time before output starts

    int sim_steps;          // auxiliary var. - total number of integration steps int(sim_time/dt)
    int step_output;        // auxiliary var. - integrations steps between outputs int(output/dt)

    int out_h5;             // switch for ouput data format (txt, HDF5)
    unsigned int outstep;       // current output step (needed for hdf5)
    unsigned int total_outstep; // total Nr of output steps


    double cludist;         // criterion for max. particle dist. in cluster
    unsigned int Sclu; // # of particles in largest cluster
    unsigned int MinCluster;    // # of particles necessary, otherwise repeat 
};
typedef struct params params;

#endif
