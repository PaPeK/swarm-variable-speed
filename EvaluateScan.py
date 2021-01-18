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
from TsTools import TSPositionPlus as tspp    # needed for spatial_properties


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


def get_UpdateMeans(datas, pre=None):
    '''
    INPUT:
        datas.shape(samples, time, value)
            assumme
    OUTPUT:
        up_rate in units [1/frame]
        up_len in units [frame]
    '''
    samples, time, N, varis = datas.shape
    out_names = ['burst_rate',
                 'burst_duration',
                ]
    out_dic = s2d.names2dict(out_names)
    if datas is None:  # only need dictionary
        out_dic = s2d.names2dict(out_names, pre)
        return out_dic
    means = np.empty((samples, len(out_names))) * np.nan
    for i, data in enumerate(datas):
        t_end = getTEnd(data)
        if t_end > 0:
            simuTSPP = tspp.TSplus(data[:t_end, :, :2], None, False, None, burstWithTurn=False)
            means[i, out_dic['burst_rate']] = np.mean(simuTSPP.up_rate)
            uplens = []
            for j in range(8):
                up_len = np.mean(simuTSPP.up_len[j])
                uplens.append(up_len)
            means[i, out_dic['burst_duration']] = np.mean(uplens)
    means = np.nanmean(means, axis=0)
    return means, out_dic 


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


def get_fishNet_meansStd(datas, outRawData=None):
    '''
    just computes mean and std
    INPUT:
        datas.shape(samples, time, value)
    '''
    gen.setDefault(outRawData, False)
    samples, time, varis = datas.shape
    out_names = ['kill_rate',
                 'kill_rate2',
                ]
    if outRawData:
        out_names += ['N_dead', 'simu_time']
    out_dic = s2d.names2dict(out_names)
    if datas is None:  # only need dictionary
        out_dic = s2d.names2dict(out_names, pre)
        return out_dic
    in_dic = s2d.get_out_swarm_fishNet_dict()
    means = np.empty((samples, len(out_names))) * np.nan
    for i, data in enumerate(datas):
        t_end = getTEnd(data)
        if t_end > 0:
            means[i, out_dic['kill_rate']] = (data[t_end-1, in_dic['Ndead']] /
                                              t_end)
            if outRawData:
                means[i, out_dic['N_dead']] = data[t_end-1, in_dic['Ndead']]
                means[i, out_dic['simu_time']] = t_end
        there = np.where(data[:, in_dic['Ndead']] == 15)[0]
        if len(there) != 0:
            t_end = there[0] + 1
        if t_end > 0:
            means[i, out_dic['kill_rate2']] = (data[t_end-1, in_dic['Ndead']] /
                                               t_end)
    if outRawData:
        return means, out_dic 
    std = np.nanstd(means, axis=0)
    means = np.nanmean(means, axis=0)
    return means, std, out_dic 


def GetSwarmMeanStd(datas, pre=None):
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
                ]
    out_dic = s2d.names2dict(out_names)
    if datas is None:  # only need dictionary
        out_dic = s2d.names2dict(out_names, pre)
        return out_dic
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
            dat = data[:t_end]
            means[i, out_dic['pol_order_Tgradient']] = np.nanmean(np.diff((dat[1:, in_dic['pol_order']])))
            means[i, out_dic['pol_order_Moment2']] = np.nanmean(dat[:, in_dic['pol_order']]**2)
            means[i, out_dic['pol_order_Moment4']] = np.nanmean(dat[:, in_dic['pol_order']]**4)
            means[i, out_dic['suscept']] = dat[-1, in_dic['N']] * np.var(dat[:, in_dic['pol_order']])
            means[i, out_dic['grp_speed']] = np.nanmean(np.sqrt( dat[:, in_dic['avg_v[0]']]**2 +
                                                                 dat[:, in_dic['avg_v[1]']]**2) )
    means = gen.NanAverage(means, weights)
    stds = np.nanstd(means)
    out_dic = s2d.names2dict(out_names, pre)
    return means, stds, out_dic


def GetSwarmPredMeanStd(datas, paras, pre=None):
    '''
    returns averages of combinations of values (Correlations, ...)
    and a dictionary with appropriate name for value
    INPUT:
        datas.shape(samples, time, value)
            shape = (samples, time, value)
            !!ASSUMES!!: h5dset=pred_swarm
    '''
    out_names = ['END<fit>_f<0',
                 'END<fit>',
                 'ENDN_f<0',
                 'ENDpred.kills',
                 'ENDN_f<0/time',
                 'END<fit>/time',
                 't_1stdead',
                ]
    out_dic = s2d.names2dict(out_names)
    if datas is None:
        out_dic = s2d.names2dict(out_names, pre)
        return out_dic
    samples, time, varis = datas.shape
    means = np.empty((samples, len(out_dic.keys()))) * np.nan
    weights = np.empty(samples)
    dic = s2d.get_out_swarm_pred_dict()
    for i, data in enumerate(datas):
        t_end = getTEnd(data)
        weights[i] = t_end
        if t_end > 0:
            dat = data[:t_end]
            # final measures
            means[i, out_dic['END<fit>_f<0']] = (dat[-1, dic['fit_sum']] /
                                                  dat[-1, dic['fit_count']])
            means[i, out_dic['END<fit>']] = (dat[-1, dic['fit_sum']] /
                                              dat[-1, dic['N_clu']])
            means[i, out_dic['ENDN_f<0']] = dat[-1, dic['fit_count']]
            means[i, out_dic['ENDpred.kills']] = dat[-1, dic['pred.kills']]
            means[i, out_dic['ENDN_f<0/time']] = (dat[-1, dic['fit_count']] /
                                                   (paras['output'] * t_end))
            means[i, out_dic['END<fit>/time']] = (dat[-1, dic['fit_sum']] /
                                              (dat[-1, dic['N_clu']] * t_end))
            t_1stdead = np.where(dat[:, dic['pred.kills']] > 0)[0]
            if len(t_1stdead) > 0:
                t_1stdead = t_1stdead[0]
            else:
                t_1stdead = t_end
            means[i, out_dic['t_1stdead']] = t_1stdead
    means = gen.NanAverage(means, weights)
    stds = np.nanstd(means)
    out_dic = s2d.names2dict(out_names, pre)
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
            means0, means_dic = AnalyseDataH5(paras, para_vals, verb=verb)
            means = np.empty((xdim, ydim, len(means0)))
            means[i1, i2] = means0
            if verb:
                print(means_dic.keys())
        else:
            if paraRun:
                para_Id_pairs += [[i1, i2]]
                continue
            means[i1, i2], means_dic = AnalyseDataH5(paras, para_vals, verb=verb)
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
    return means, means_dic


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

    means, stds = GetSimpleMeanStd(dset_swarm)
    out_dic = s2d.outDics().get_out_swarm_dict()
    means0, stds, out_dic0 = GetSwarmMeanStd(dset_swarm)
    means = np.append(means, means0)
    out_dic = s2d.join_enumerated_dicts(out_dic, out_dic0)

    if dset_part is not None:
        means0, out_dic0 = get_UpdateMeans(dset_part)
        means = np.append(means, means0)
        out_dic = s2d.join_enumerated_dicts(out_dic, out_dic0)
    if paraRun:
        return para_ids, means
    return means, out_dic


# Main function:
def Evaluate(sourcepath):
    f_paras = os.path.join(sourcepath, "paras_independent.pkl")
    paras = pickle.load(open(f_paras, "rb"))
    para_name0 = paras['para_name0']
    para_name1 = paras['para_name1']
    para_values0 = paras['para_values0']
    para_values1 = paras['para_values1']
    # parameters used in simulation:
    out_h5 = paras['out_h5']
    assert out_h5 == 1, "data-analysis only supported for hdf5"

    f_result = os.path.join(sourcepath, "MeanData.pkl")
    os.system('rm ' + f_result)
    if os.path.isfile(f_result):
        print('loading:', f_result)
        means, means_dic = pickle.load(open(f_result, 'rb'))
    else:
        print('start evaluation of ', len(para_values0) * len(para_values1),
              ' folders: ')
        means, means_dic = evaluate_scan2d(paras, verb=True)
        pickle.dump([means, means_dic], open(f_result, 'wb'))
    f_current = os.path.realpath(__file__)
    gen.copyIfNeeded( f_current, sourcepath )

if __name__ == '__main__':
    Evaluate("./")
