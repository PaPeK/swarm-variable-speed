/*  Settings
    defines functions to set parameters for
    SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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
#ifndef settings_H
#define settings_H
#include "common_defines.h"
#include "agents.h"
#include "agents_dynamics.h"
#include "mathtools.h"
#include "input_output.h"

#include <H5Cpp.h>      // for hdf5 output
#include <vector>
#include <limits>       // for accessing numerical limits
#include <algorithm>    // to use std::set_difference, std::find, st::iota
#include <chrono>       // std::chrono::system_clock
#include <random>       // std::default_random_engine

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

void SetCoreParameters(params* SysParams);
void InitSystemParameters(params* SysParams);
void InitSystem(std::vector<particle> &a, params SysParams);
void InitPredator(std::vector<predator> &preds);
void ResetSystem(std::vector<particle> &a, params *ptrSP, bool out, gsl_rng *r);
#endif
