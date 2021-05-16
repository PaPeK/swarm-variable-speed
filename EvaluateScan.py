import os
import sys
import configparser
import multiprocessing as mp
import numpy as np
import glob     # to count number of folders
import pdb
import time as pytime
from functools import partial
# from numba import njit
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable  # for colorbar
import h5py
import importlib.util
import pickle
import Scan2D as s2d
import general as gen
from TsTools import TSPosition as tsp    # needed for spatial_properties
from pathlib import Path
from numba import jit

def GetParticleMeasures(datas, R0=None):
    '''
    INPUT:
        datas.shape(samples, time, N, value)
            with values
    '''
    if R0 is None:
        R0 = 2
    out_names = ['Corr_speed_op_local']
    out_dic = s2d.names2dict(out_names)
    in_dic = s2d.outDics().get_out_swarm_dict()
    varis = len(out_names)

    samples, time, N, varis = datas.shape
    means = np.empty(varis) * np.nan
    polLocs, sLocs = np.empty(0), np.empty(0)
    for i, data in enumerate(datas):
        t_end = getTEnd(data)
        if t_end > 0:
            dat = data[:t_end]
            pos, vel = dat.T[:2].T, dat.T[2:4].T
            phi = np.arctan2(vel.T[1].T, vel.T[0].T)
            s = np.sqrt(np.sum(vel**2, axis=-1))
            polLoc, sLoc = localValues(pos, phi, s, R0)
            polLocs = np.append(polLocs, polLoc.flatten())
            sLocs = np.append(sLocs, sLoc.flatten())
    means[out_dic['Corr_speed_op_local']] = np.corrcoef(polLocs, sLocs)[0, 1]
    return means, out_dic

def localValues(pos, phi, s, R0):
    '''
    wrapper for "calculateLocalValues"
    INPUT:
        pos shape=[time, N, 2]
        phi, s shape=[time, N]
        R0 float
    '''
    ll_t, Npart = phi.shape
    x, y = pos[:, :, 0], pos[:, :, 1]
    return calculateLocalValues(ll_t, Npart, x, y, phi, s, R0)


@jit
def calculateLocalValues(ll_t, Npart, x, y, phi, s, R0):
    '''
    "Luis Alberto Gomez Nava"-code to compute local-order and -velocity
    -modifications: before "local group-velocity" was computed instead of local speed
    INPUT:
        ll_t, Npart: ints
            length of timeseries, particle number
        x, y, phi, s: floats
            position, heading angle, speed
        R0: float
            radius
    '''
    polLoc = np.zeros((ll_t, Npart))
    vLoc   = np.zeros((ll_t, Npart))
    for t in range(ll_t):
        for i in range(Npart):
            aux1, aux2  = 0.0, 0.0
            aux3  = 0.0
            N_nei = 0
            for j in range(Npart):
                dist = np.sqrt((x[t,j] - x[t,i])**2 + (y[t,j] - y[t,i])**2)
                if (dist <= R0):
                    # Local polarization
                    aux1  = aux1 + np.cos(phi[t,j])
                    aux2  = aux2 + np.sin(phi[t,j])
                    # Local speed
                    aux3  = aux3 + s[t,j]
                    N_nei += 1
            aux1, aux2 = aux1/N_nei, aux2/N_nei
            aux3 = aux3/N_nei
            polLoc[t,i] = np.sqrt(aux1**2 + aux2**2)
            vLoc[t,i]   = aux3
    return polLoc, vLoc


def getTEnd(data):
    '''
    returns end-time: when all values are 0
    this is caused if the simulation is stopped before
    the simulation-time is over (due to break condition
                                 ...Pred leaves swarm)
    INPUT:
        data.shape(time, vars)
        OR
        data.shape(time)
    '''
    var_means = data
    while len(var_means.shape) >= 2:
        var_means = var_means.sum(axis=1)
    there = np.where(var_means == 0)[0]
    t_end = data.shape[0]
    if len(there) > 0:
        t_end = there[0]
    if t_end == 0:
        print('var_means[:3]: ', var_means[:3])
    return t_end


def GetSimpleMeanStd(datas):
    '''
    just computes mean and std
    INPUT:
        datas.shape(samples, time, value)
    '''
    samples, time, varis = datas.shape
    means = np.empty((samples, varis)) * np.nan
    weights = np.empty(samples)
    for i, data in enumerate(datas):
        t_end = getTEnd(data)
        weights[i] = t_end
        if t_end > 0:
            dat = data[:t_end]
            means[i] = np.nanmean(dat, axis=0)
    std = np.nanstd(means, axis=0)
    means = gen.NanAverage(means, weights)
    return means, std


def GetSwarmMeanStd(datas):
    '''
    returns averages of combinations of values (Correlations, ...)
    and a dictionary with appropriate name for value
    INPUT:
        datas.shape(samples, time, value)
            !!ASSUMES!!: h5dset=pred_swarm
    '''
    out_names = ['pol_order_Tgradient',
                 'pol_order_Moment2',
                 'pol_order_Moment4',
                 'suscept',
                 'grp_speed',
                 'std_grp_speed',
                 'std_pol',
                 'Corr_speed_op'
                ]
    out_dic = s2d.names2dict(out_names)
    in_dic = s2d.outDics().get_out_swarm_dict()
    varis = len(out_names)
    # preparing averaging
    samples = len(datas)
    means = np.empty((samples, varis)) * np.nan
    weights = np.empty(samples)
    for i, data in enumerate(datas):
        t_end = getTEnd(data)
        weights[i] = t_end
        if t_end > 0:
            dat = data[:t_end].T
            order = dat[in_dic['pol_order']]
            means[i, out_dic['pol_order_Tgradient']] = np.nanmean(np.diff((order)))
            means[i, out_dic['pol_order_Moment2']] = np.nanmean(order**2)
            means[i, out_dic['pol_order_Moment4']] = np.nanmean(order**4)
            means[i, out_dic['suscept']] = dat[in_dic['N'], -1] * np.var(order)
            means[i, out_dic['std_pol']] = np.nanstd(order)
            grp_speed = np.sqrt( dat[in_dic['avg_v[0]']]**2 +
                                 dat[in_dic['avg_v[1]']]**2)
            means[i, out_dic['grp_speed']] = np.nanmean(grp_speed)
            means[i, out_dic['std_grp_speed']] = np.nanstd(grp_speed)
            speed = dat[in_dic['avg_speed']]
            means[i, out_dic['Corr_speed_op']] = np.corrcoef(speed, order)[0, 1]
    means = gen.NanAverage(means, weights)
    stds = np.nanstd(means)
    out_dic = s2d.names2dict(out_names)
    return means, stds, out_dic


def evaluate_scan2d(paras, verb=None, deleteIfDone=None, paraRun=None):
    '''
    Goes sequentially or parallel (if paraRun is True) through all output files and evaluate the data
    '''
    if verb is None:
        verb = False
    if deleteIfDone is None:
        deleteIfDone = False
    if paraRun is None:
        paraRun = True
    if paraRun:
        para_Id_pairs = []
    para_values0 = paras['para_values0']
    para_values1 = paras['para_values1']
    xdim = len(para_values0)
    ydim = len(para_values1)
    for i in range(xdim * ydim):
        if verb:
            print('processing folder ', i, 'finished from total: ',
                  float(i) / (xdim * ydim) * 100, '%', end='\r')
        i1 = int(i % xdim)
        i2 = int(i / xdim)
        para_vals = [para_values0[i1], para_values1[i2]]
        # initialize
        if i == 0:
            means0, means_dic, fileInfo = AnalyseDataH5(paras, para_vals, verb=verb)
            means = np.empty((xdim, ydim, len(means0)))
            fileInfos = np.empty((xdim, ydim, 2), dtype=object) # file-name and group-name
            means[i1, i2] = means0
            fileInfos[i1, i2] = fileInfo
            if verb:
                print(means_dic.keys())
        else:
            if paraRun:
                para_Id_pairs += [[i1, i2]]
                continue
            means[i1, i2], means_dic, fileInfos[i1, i2] = AnalyseDataH5(paras, para_vals, verb=verb)
        if deleteIfDone:
            f_name = s2d.get_out_name(paras, para_vals)
            gen.silentRemove(f_name)
    if paraRun:
        partAnaDatH5 = partial(AnalyseDataH5, paras, verb=verb, paraRun=paraRun)
        comp_pool = mp.Pool(mp.cpu_count()-3)
        out = comp_pool.map(partAnaDatH5, para_Id_pairs, chunksize=1)   # para_values passed as tuple
        for i in range(len(out)):
            index = out[i][0]
            means[index[0], index[1]] = out[i][1]
            fileInfos[index[0], index[1]] = out[i][2]
    return means, means_dic, fileInfos


def AnalyseDataH5(para, para_vals, verb=None, paraRun=None):
    '''
    INPUT:
        para dict
            dictionary of parameters used to run the simulation
        para_vals [para0, para1]
            values of parameters of the 2D-scan
    '''
    if verb is None:
        verb = True
    if paraRun is None:
        paraRun = False
    if paraRun:
        para_ids = np.array(para_vals)
        para_vals = [para['para_values0'][para_vals[0]],
                     para['para_values1'][para_vals[1]]]
    f_name = s2d.get_out_name(para, para_vals)
    g_name = s2d.get_groupname(para, para_vals)

    def h5LoadIfExists(f, d_name):
        items = []
        f.visit(items.append)
        dset = None
        if d_name in items:
            dset = np.array(f[d_name])
        return dset

    with h5py.File(f_name, 'r+') as f:   # open existing file, must exist
        dset_swarm = h5LoadIfExists(f, g_name + '/swarm')
        dset_part = h5LoadIfExists(f, g_name + '/part')

    means = np.empty(0, dtype=float)
    out_dic = dict()

    if dset_swarm is not None:
        means0, stds = GetSimpleMeanStd(dset_swarm)
        out_dic0 = s2d.outDics().get_out_swarm_dict()
        means = np.append(means, means0)
        out_dic = s2d.join_enumerated_dicts(out_dic, out_dic0)

        means0, stds, out_dic0 = GetSwarmMeanStd(dset_swarm)
        means = np.append(means, means0)
        out_dic = s2d.join_enumerated_dicts(out_dic, out_dic0)

    if dset_part is not None:
        means0, out_dic0 = GetParticleMeasures(dset_part)
        means = np.append(means, means0)
        out_dic = s2d.join_enumerated_dicts(out_dic, out_dic0)

    if paraRun:
        return para_ids, means, [f_name, g_name]
    return means, out_dic, [f_name, g_name]


# Main function:
def Evaluate(sourcepath):
    f_paras = Path(sourcepath) /  'paras_independent.pkl'
    paras = pickle.load(f_paras.open("rb"))
    para_name0 = paras['para_name0']
    para_name1 = paras['para_name1']
    para_values0 = paras['para_values0']
    para_values1 = paras['para_values1']
    # parameters used in simulation:
    out_h5 = paras['out_h5']
    assert out_h5 == 1, "data-analysis only supported for hdf5"

    f_result = Path(sourcepath) / "MeanData.pkl"
    if f_result.exists():
        f_result.unlink()
    print('start evaluation of ', len(para_values0) * len(para_values1),
          ' folders: ')
    means, means_dic, fileInfos = evaluate_scan2d(paras, verb=True)
    pickle.dump([means, means_dic, fileInfos], f_result.open('wb'))
    f_current = os.path.realpath(__file__)
    gen.copyIfNeeded( f_current, sourcepath )

if __name__ == '__main__':
    Evaluate("./")
