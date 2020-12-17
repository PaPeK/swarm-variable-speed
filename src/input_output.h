/*  InputOutput
    defines basic in- and output operations for 
    SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser, Pawel Rormanczuk

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
#ifndef input_output_H
#define input_output_H
#include <fstream>      // for writing vector to file
#include <iostream>     // to use std::cout
#include <getopt.h>     // for required and optional argument
#include <math.h>       // probably not needed
#include <iterator>     // to use std::begin ....
#include <string.h>
#include <string>
#include <algorithm>

#include "agents.h" 

void ParseParameters(int argc, char **argv, params *SysParams);
void OutputParameters(params SysParams);
void LoadCoordinatesCPP(params * ptrSP, std::string name, std::vector<particle> &a);
void LoadCoordinates(params * ptrSP, const char *fn, std::vector<particle> &a, int N, double sizeL);
void LoadVector(params * ptrSP, std::string name, std::vector<double> &vec);
void WritePosVel(std::vector<particle> &a, params* ptrSP, std::string name, bool append=true);
template<class T>
void WriteVector(std::string file, std::vector<T> &vec, bool append=true);
template<class T>
void WriteVector2d(std::string file, std::vector< std::vector<T> > &vec, bool append=true);
void pava_load_set_paras(std::vector<particle> &a, std::string in_name);

#endif
