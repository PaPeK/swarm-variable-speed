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
    std::vector<double> force_flee;   // flee vector
    double fitness;         // current fitness of particle
    double fit_decrease;    // how much fitness got decreased in current time_step
    bool dead;              // if killed by predator
    unsigned int id;        // ID of each agent (needed to split and merge agents correctly)
    std::vector<unsigned int> NN;    // vector containing all NN

    // counters of interaction partners - important only for local metric coupling (not global)
    int counter_rep;        // counter repulsion partners 
    int counter_flee;       // counter flee
    // parameter-variations in between the agents
    double alg_strength;
    double rep_steepness;
    std::vector<double> out(void);  // for output
};
typedef struct particle particle;


struct predator{
    unsigned int id;        // ID of each agent (needed to split and merge agents correctly)
    std::vector<double> x;            // agent position vector
    std::vector<double> v;            // agent velocity vector
    std::vector<double> u;            // agent direction unit vector
    double phi;             // agent direction polar angle
    double phi_start;       // phi at creation of predator
    double vproj;           // vel along heading direction
    std::vector<int> cluster;// contains indices of largest cluster hunted by pred.
    unsigned int kills; // counts how many prey got killed (only if pred_kill == 1)
    bool beforeCOM;         // needed for pred_attack == 2, true: follows COM, false: shoots straight
    unsigned int state;     // state of predator:0=approach, 1=hunt
    std::vector<unsigned int> NN;    // vector containing all F seen by P
    std::set<unsigned int> NNset;    // set of all F seeing P
    std::set<unsigned int> NN2set;   // set of all F seeing NNset which are not in NNset
    std::set<unsigned int> NN3set;   // set of all F seeing NN2set  which are not in NNset or NN2set
    std::vector<double> out(void);
};
typedef struct predator predator;

// data structure for system parameters, and auxiliary variables
struct params{

    std::string location;   // path of output-files
    std::string fileID;     // name of output-files
    unsigned int N;                  // number of agents
    int Npred;              // 1 if predator, 0 if not
    int Ndead;              // # of dead/killed agents
    double sizeL;           // system size
    double Dphi;            // angular noise intensity  (= standard deviation)
    double noisep;          // auxiliary var. - std. deviation of the angular noise

    double alg_strength;    // alignment strength
    double rep_strength;    // repulsion strength
    double rep_range;       // repulsion range
    double rep_steepness;   // repulsion - steepness of the sigmoid function

    double flee_strength;   // flee strength
    double pred_strength;   // force of predator

    double pred_time;       // time at which the predator comes in
    double pred_angle;      // angle relative to flock that the predator appears
    double pred_radius;     // how far away the predator starts from flock
    double pred_speed0;     // how fast the predator moves
    double pred_memo_t;     // memory time before pred_time from which on predator observes prey cluster
    double pred_angle_noise;// angular noise at creation of predator (directly behind or roughly directly behind)
    std::vector<double> pred_memo_av; // memory average velocity of estimated cluster, averaged in t \in [pred_time-pred_memo_t, pred_time] 
    // cluster which is first attacked by predator
    // -difference to pred.cluster: is not changing during simulation and computed only at creation of predator
    std::vector<int> cluster_attacked;


    double sim_time;        // simulation time in natural time units
    double speed0;          // stationary speed of single individuals
    double dt;              // time step of integration
    double output;          // output time in natural units

    int BC;                 // switch for boundary conditions
    int IC;                 // switch for initial conditions
    int pred_attack;// determines if the predator 1: strictly moves relative to COM, 2: follows the COM, 3: adjust start position to hit com with straight move
    int pred_kill;         // if predator kills prey or flies through
    bool stopAtKill;       // if simulations is stopped after 1st kill

    double trans_time;      // transient time before output starts

    int sim_steps;          // auxiliary var. - total number of integration steps int(sim_time/dt)
    int step_output;        // auxiliary var. - integrations steps between outputs int(output/dt)

    bool out_dummy;         // derived from Npred
    int out_h5;             // switch for ouput data format (txt, HDF5)
    unsigned int outstep;       // current output step (needed for hdf5)
    unsigned int outstep_pred;  // current output step (needed for hdf5)
    unsigned int total_outstep; // total Nr of output steps
    unsigned int total_outstep_pred; // total Nr of output steps


    double cludist;         // criterion for max. particle dist. in cluster
    unsigned int Sclu; // # of particles in largest cluster
    unsigned int MinCluster;    // # of particles necessary, otherwise repeat 
    double kill_rate;        // kills per second if prob_selected = 1 and prob_catch = 1
    double kill_range;       // distance between pred and prey at which prob_kill > 0
    double predcirclerad;       // radius of circle arond predator = 3*kill_range(for fct. GetPredPreyCircle )
};
typedef struct params params;

#endif
