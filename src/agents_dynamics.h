/*  AgentsDynamics
    defines movement rules/differential equations of agents in
    SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser

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
#ifndef agents_dynamics_H
#define agents_dynamics_H
// OWN MODULES:
#include "common_defines.h"
#include "agents.h"
#include "agents_operation.h"
// #include "mathtools.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <math.h>
#include <random>       // std::default_random_engine

void MoveParticle(particle &a, params * ptrSP, gsl_rng *r, double rnp);
template<class agent>
void Boundary(agent &a, double sizeL,  int BC);             // calculate boundary conditions
#endif
