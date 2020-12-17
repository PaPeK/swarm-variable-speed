'''
    SwarmDynByPy
    functions which define the parameters and provide the appropriate c-call
    to start the simulation in C++ via swarmdyn
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
'''
import numpy as np
import os
from pathlib import Path


def get_base_params(record_time, trans_time=None):
    if trans_time is None:
        trans_time = 0
    params = dict()

    # params["path"] = "/mnt/data2/PredPrey/DEFAULT/"
    params['path'] = "./"
    params["fileID"] = 'xx'  # not passed to the code, only there
    params["out_h5"] = 1 # output-format 0:txt 1:hdf5
    params["trans_time"] = trans_time    # time till output starts (transient)
    params["time"] = record_time + trans_time    # total time of simulation
    params["dt"] = 0.02
    params["output"] = .2
    # Initial Conditions
    # IC    0:v-dis,x-dis-box, 1:v-order, x-dis, 2:v-order, x-dis-box
    #       3:v-order, x-dis-circ  4: v-milling, x-cricle  5:v-dis, x-circle
    #       6:v-dis, x-dis 99:from file
    params["IC"] = 3     # recommended:2 global, 3 voronoi
    # Boundary Conditions
    # BC    -1:no, 0: periodic, 1:inelastic, 2:elastic, 3:x peri, y ela
    #       4:x peri, y inela 5:elastic circle
    params["BC"] = -1
    params["N"] = 400
    params["size"] = 150

    # #################### F behavior
    params["Dphi"] = 0.5    # angular noise (= variance/2) --> Pi^2/4=2.5 is very large
    params["speed0"] = 1.0
    params["beta"] = 4
    params["rep_range"] = 1
    params["rep_strength"] = 2.0 # 20. (voro)
    params["rep_steepness"] = -2    # for rep_range = 1 steepness should be between [-4, -2]
    params["alg_strength"] = 2
    params["flee_strength"] = 4


    # #################### Params for OUTPUT computation
    # cludist< alpha-value(=12, compute distance of P to alpha-shape-swarm)
    #   if |r_ij| > cludist -> not same cluster (unless close to others)
    params["cludist"] = 5 * params['rep_range']
    params["MinCluster"] = 0.9 # percentage of cluster necessary to pass the simu

    return params


def solve_parameter_dependencies(dic):
    '''
    some parameter values need to be adjusted if specific parameters
    are changed to be roughly in the "same" parameter space
    ("same" defined by the properties of the collective, e.g. NND)

    Other parameters as dic['MinCluster'] do depend on parameters
    which might be varied in between different parameter-runs
    -> before simulation is started dependencies are re
    INPUT:
        dic, dictionary
            keys = parameter-names
            values = parameter-values
    '''
    if dic['path'][-1] != os.sep:
        dic['path'] += os.sep

    # ensure consistent parameters:
    dic['MinCluster'] *= dic['N']
    if dic['MinCluster'] == 0:
        dic['MinCluster'] = 1


def dic2swarmdyn_command(dic):
    '''
    transforms dictionary of parameters to calls to swarmdyn.cpp
    INPUT:
        dic dictionary
            keys = parameter names
            values = parameter values
    OUTPUT:
        command string
    '''
    #consistency check:
    for key in dic.keys():
        if 'time' in key:
            assert np.any(dic[key] >= 0), '{} has negative value {}'.format(key, dic[key])
    command = './swarmdyn'
    command += ' -l %s' % dic['path']
    command += ' -E %s' % dic['fileID']
    command += ' -J %d' % dic['out_h5']
    command += ' -d %g' % dic['dt']
    command += ' -o %g' % dic['output']
    # simulation
    command += ' -L %g' % dic['size']
    command += ' -t %g' % dic['time']
    command += ' -T %g' % dic['trans_time']
    command += ' -B %d' % dic['BC']
    command += ' -I %d' % dic['IC']
    command += ' -M %d' % dic['MinCluster']
    command += ' -C %g' % dic['cludist']
    # particles
    command += ' -N %d' % dic['N']
    command += ' -h %g' % dic['rep_range']
    command += ' -Y %g' % dic['rep_steepness']
    command += ' -H %g' % dic['rep_strength']
    command += ' -A %g' % dic['alg_strength']
    command += ' -D %g' % dic['Dphi']
    command += ' -b %g' % dic['beta']
    command += ' -s %g' % dic['speed0']
    return command
