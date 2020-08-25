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


def get_base_params(pred_time, record_time, trans_time=None):
    if trans_time is None:
        trans_time = 0
    params = dict()

    # params["path"] = "/mnt/data2/PredPrey/DEFAULT/"
    params['path'] = "./"
    params["fileID"] = 'xx'  # not passed to the code, only there
    params["out_h5"] = 1 # output-format 0:txt 1:hdf5
    params["trans_time"] = pred_time + trans_time    # time till output starts (transient)
    params["time"] = pred_time + record_time + trans_time    # total time of simulation
    params["dt"] = 0.02
    params["output"] = .2
    params["output_mode"] = 1 # 0:mean, 1:full, 2:no-output AND no-dummies (only pava_out)
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
    params["Npred"] = 1
    params["size"] = 150

    # #################### F behavior
    params["Dphi"] = 0.5    # angular noise (= variance/2) --> Pi^2/4=2.5 is very large
    params["speed0"] = 1.0
    params["rep_range"] = 1
    params["rep_strength"] = 2.0 # 20. (voro)
    params["rep_steepness"] = -2    # for rep_range = 1 steepness should be between [-4, -2]
    params["alg_strength"] = 2 
    params["flee_strength"] = 4
    params["pred_strength"] = 2
    #################### Evolution parameter
    # pava_sig: standard deviation of parameters AND mutation strength
    params["pava_sig"] = 0.075
    params["discretize"] = 1  # int: np.round(pavas, discretize)



    # #################### P behavior
    # pred_com = hunting behavior of P
    #   P approachs to: <10 largest clu, <20 closest clu, <30 closest prey
    #                   >30 just go straight
    #   P hunts:  0: straigth, 1: foll. COM,
    #             2: foll. csg p_catch-weighted mean
    params["pred_com"] = 2
    # pred_kill: 0: no kill, 1: probabilistic kill, 2:select closest (still probabilistic)
    # if pred_kill>=10 -> stopAtKill = True and pred_kill -= 10
    # than if pred_kill>5 : no kill and pred_kill %= 5
    params["pred_kill"] = 1
    params["pred_time"] = pred_time
    params["pred_angle"] = 1 * np.pi
    params["pred_angle_noise"] = np.pi    # angular noise on pre_angle
    params["pred_speed0"] = 2 * params["speed0"]
    params["kill_range"] = 3 * params['rep_range']
    params["kill_rate"] = 1
    params["pred_radius"] = 1.5 * params["kill_range"]

    # #################### Params for OUTPUT computation
    # cludist< alpha-value(=12, compute distance of P to alpha-shape-swarm)
    #   if |r_ij| > cludist -> not same cluster (unless close to others)
    params["cludist"] = 2 * 3 * params['rep_range']  # 2 * params['kill_range']
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
    command += ' -m %d' % dic['output_mode']
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
    command += ' -F %g' % dic['flee_strength']
    command += ' -f %g' % dic['pred_strength']
    command += ' -D %g' % dic['Dphi']
    command += ' -s %g' % dic['speed0']
    # predator
    command += ' -n %d' % dic['Npred']
    command += ' -e %g' % dic['pred_time']
    command += ' -p %g' % dic['pred_angle']
    command += ' -P %g' % dic['pred_radius']
    command += ' -S %g' % dic['pred_speed0']
    command += ' -c %d' % dic['pred_com']
    command += ' -z %d' % dic['pred_angle_noise']
    command += ' -x %d' % dic['pred_kill']
    command += ' -G %g' % dic['kill_range']
    command += ' -O %g' % dic['kill_rate']
    return command


def random_positive_normal(mean, std, size, mean_tolerance=None):
    '''
    creates a positive normal distribution with numpy.random.normal
    and sets all negatives to 0
    if resulting mean deviates stronger than mean_tolerance:
        wald-distribution is used to create distribution
    '''
    if mean_tolerance is None:
        mean_tolerance = 0.1
    if mean == 0:
        mean = 0.05 * std
    dist = np.random.normal(mean, std, size)
    floored = len( dist[dist<0] )
    if floored > 0:
        dist[dist<0] = 0
        dist_mean = dist.mean()
        dist_std = dist.std()
        if ( dist_mean - mean ) / std > mean_tolerance: # if rel. difference larger 10%
            dist = np.random.wald(mean, mean**3 / std**2, size)
    return dist


def heteroPopulation(dic):
    '''
    makes the population heterogeneous
    the parameters have their mean around the corresponding values of "dic"
    the parameters are randomly disturbed by stds defined by "sigs"
    saves the array in "path" with file-name containing "fileID"
    INPUT:
        dic dictionary
            dictionary with parameters of simulation
    OUPUT:
        pavas.shape(N, 4)
            each row represnting paras of 1 agent
    '''
    N = int(dic['N'])
    pavas = np.zeros(N, dtype=float)
    sig = dic['pava_sig']
    phenotype = 'alg_strength'
    mean = dic[phenotype]
    if sig > 0:
        signMean = np.sign(mean)
        if signMean == 0:
            signMean = 1
        pavas[:] = (signMean *
                    random_positive_normal(np.abs(mean), sig, N))
        ############ DISCRETIZE ############
        if dic['discretize'] is not None:
            pavaTest = np.round(pavas, dic['discretize'])
            increase = 2
            # to ensure heterogenous population
            while pavaTest.sum() == 0:
                pavaTest = np.round(pavas * increase, dic['discretize'])
                increase += 1
            pavas = pavaTest
    else:
        pavas = mean
    pavas = pavas.reshape(N, 1)
    f_pavaInput = Path(dic['path']) / ('pava_in_' + dic['fileID'] + '.in')
    np.savetxt(str(f_pavaInput), pavas)
    return pavas
