/*  H5Tools
    defines standard procedures to save data in hdf5-files
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
#ifndef h5tools_H
#define h5tools_H
#include <H5Cpp.h>      // for hdf5 output
#include <vector>
#include <sys/stat.h>
#include <string.h>
#include <iterator>     // to use std::begin ....
#include <stdio.h>
#include <stdlib.h>

H5::DataSet h5CreateDSet(H5::H5File *file, std::vector<hsize_t> dim,
                   std::string name, std::string type);
H5::DataSet h5CreateDSetExt(H5::H5File *file, std::vector<hsize_t> dim,
                     std::string name, std::string type);
void h5WriteInt(H5::DataSet *dset, std::vector<int> data, std::vector<hsize_t> offset);
void h5WriteDouble(H5::DataSet *dset, std::vector<double> data, std::vector<hsize_t> offset);
void h5readDimension(H5::DataSet *dataset, std::vector<hsize_t> &dim);
void h5read_extend_dim(H5::DataSet *dataset, std::vector<hsize_t> &dim);
void h5CreateWriteDset(std::string f_h5out, std::string n_group,
                       std::string n_dset,
                       std::vector< std::vector<double> > data,
                       bool extend);
inline bool exists (const std::string& name) {
    struct stat buffer;   
    return (stat (name.c_str(), &buffer) == 0); 
}


#endif
