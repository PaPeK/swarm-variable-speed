import os
import multiprocessing as mp
import numpy as np
from functools import partial
import glob
import sys
import time as pytime
import h5py
import shutil
import matplotlib
from pathlib import Path
import SwarmDynByPy as sd
import general as gen
import pdb
import pickle


if __name__ == '__main__':
    matplotlib.use('Agg')

def increaseDictValues(dic, inc):
    for key in dic.keys():
        dic[key] += inc

def get_out_name(dic, para_vals):
    groupname = get_groupname(dic, para_vals)
    return str(Path(dic['path']) / ('out_' + groupname + '.h5'))

def get_groupname(dic, para_vals):
    '''
    checks in path if the value is close to an already existing value
    if yes: the string of this close value is used in groupname
        no: the string representation is used
    '''
    files = glob.glob( os.path.join(dic['path'], 'out*h5') )
    para_valsStr = []
    for i, para in enumerate(para_vals):
        para_nam = dic['para_name' + str(i)]
        files_filtered = [f for f in files if para_nam in f] # necessary
        there = []
        if len(files_filtered) > 0:
            paravals = np.array([gen.string_NrAfterWord(fil, para_nam) for fil in files_filtered])
            there = np.where(paravals == None)
            paravals = np.delete(paravals, there)
            paravals = np.unique(paravals)
            there = np.where( np.isclose(para, np.array(paravals, dtype=float)) )[0]
        if len(there) > 0:
            para_valsStr.append(paravals[there[0]])
        else:
            para_valsStr.append(str(para))

    name = ''
    for i, para in enumerate(para_vals):
        key_name = 'para_name' + str(i)
        if key_name in dic.keys():
            name += dic[key_name] +  '{}_'.format(para_valsStr[i])
        else:
            raise ValueError('para_name{} not in params.keys()'.format(i))
    if name == '':
        name = 'run_base_paras'
    else:
        name = name[:-1]  # exclude last "_" from file_name
    return name


def names2dict(names, pre = None):
    if pre is None:
        pre = ''
    keys = [pre + key for key in names]
    vals = list(range(len(keys)))
    dic = dict(zip(keys, vals))
    assert len(dic.keys()) == len(keys), 'len(dic.keys()) != len(keys) => probably redundant key'
    return dic


# Dictionaries
class outDics:
    '''
    contains all the variables names for the data saved in the hdf5-datasets
    numbers correspond to entries of last dimension
    '''
    def __init__(self, pre=None):
        self.pre = pre
        self.swarm = self.get_out_swarm_dict()
        self.part = self.get_out_part_dict()


    def get_out_swarm_dict(self):
        keys0 = ['N',
                 'pol_order',
                 'L_norm',
                 'avg_x[0]',
                 'avg_x[1]',
                 'avg_v[0]',
                 'avg_v[1]',
                 'avg_speed',
                 'avg_vsquare',
                 'elongation',
                 'aspect_ratio',
                 'a_com_maxIID',
                 'Area_ConHull',
                 'NND',
                 'IID',
                 'ND',
                 'varVelFluc',
                 'C_velFluc',
                 'av_eddi',
                 'inv_NND2',
                 'var_inv_NND2',
                ]
        dic = gen.list2dict(keys0, pre=self.pre)
        return dic


    def get_out_part_dict(self):
        keys0 = ['x0',
                 'x1',
                 'v0',
                 'v1',
                 'force0',
                 'force1',
                 'eddi']
        dic = gen.list2dict(keys0, pre=self.pre)
        return dic


def join_enumerated_dicts(dic0, dic1):
    Nkeys0 = len(dic0.keys())
    if Nkeys0 > 0:
        assert Nkeys0 - 1 == max(dic0.values()), 'len(dic0.keys()) - 1 != dic0.values() => Wrong values'
    dic_new = dic0.copy()
    dic11 = dic1.copy()
    for key in dic11.keys():
        dic11[key] += Nkeys0
    dic_new.update(dic11)
    return dic_new


def h5CompareGrAttr2Dict(h5file, gname, dicce, verb=False):
    '''
    compares group values of group attributes of h5-file with
    values of keys of given dictionary 'dicce'
    '''
    if gname in list(h5file.keys()):
        grp = h5file.require_group(gname)
        # compare attributes:
        count = 0
        for i in list(grp.attrs.keys()):
            if i in dicce.keys():
                if type(dicce[i]) in [float, int] and i not in ["output_mode"]:
                    if dicce[i] != grp.attrs[i]:
                        print('MISMATCH OF ', i, ': ', dicce[i], '!=', grp.attrs.get(i))
                        count += 1
        if count == 0:
            if verb:
                print('all relevant attribtes MATCH')
        else:
            sys.exit()
    else:
        print('GROUP not in FILE')
        sys.exit()


def create_or_check_swarmdyn_hdf5(path, dirname, dic, verb=False):
    path = Path(path)
    f_h5 = str(path / ('out_' + dirname + '.h5'))
    outDic = outDics()  # swadyDics = SwarmdynDataDics()
    stepsBetween, dt = int(dic["output"]/dic["dt"]), dic["dt"]
    Steps = (int((int(dic["time"]/dt)-1)/stepsBetween) -
                 int(int(dic["trans_time"]/dt)/stepsBetween) +
                 int(0 == dic["trans_time"]/dt/stepsBetween % 1))
    try:    # try creation of file and datasets
        with h5py.File(f_h5, 'w-') as f:
            grp = f.create_group(dirname)
            for i in list(dic.keys()):      # create attributes of group
                if type(dic[i]) != str:     # using str as attribute value is complicated
                    grp.attrs.create(i, dic[i])
            dims = np.array([0, Steps, len(outDic.swarm)])
            part_h5createDataset = partial(gen.h5createDataset, f, dirname)
            part_h5createDataset('swarm', dims)
            if dic["output_mode"] > 0: # CREATE PARTICLE DSET:
                dims = np.array([0, Steps, dic["N"], len(outDic.part)])
                gen.h5createDataset(f, dirname, 'part', dims)
        existing_samples = 0
    except:  # load file and check if same attributes where used
        if verb:
            print('PPK: try opening')
        with h5py.File(f_h5, 'r+') as f:
            h5CompareGrAttr2Dict(f, dirname, dic)
            n_dsets = gen.h5AllDatasetNames(f, verbose=False)
            sizes = [f[dset].shape[0] for dset in n_dsets]
            existing_samples = np.nanmax(sizes)
    return existing_samples


def MultipleRunSwarmdyn(params, runs, paratuple, verb=False):
    '''
    sequantially runs swarmdyn with same parameters
    ?params[out_h5] == 1?
        -yes: ?existing hdf5-datasets with same name?
            ->yes: ?same paras used?
               ->yes: use same dataset
               ->no: stop simulation
            -> no: create new hdf5-dataset
        -> no: save in new txt-file (with index corresponding to run)
    INPUT:
        params dict()
            contains all parameters needed for simulation
        runs int
           # of repeats of same parameters 
        paratuple list
            para_values which shall be changed
            corresponding para_names are defined in params['para_name0']
            params['para_name1']
    OUTPUT:
        [time, paratuple]
            only needed to estimate computation time for specific
            parameters
    '''
    dic = params.copy()  # to not modify original part
    dirname = get_groupname(dic, paratuple)
    for i, val in enumerate(paratuple):
        dic[dic['para_name{}'.format(i)]] = paratuple[i]
    sd.solve_parameter_dependencies(dic)
    path = dic['path']
    dic['fileID'] = dirname
    existingRuns = 0
    if verb:
        print('Running {} runs for '.format(runs) + dirname)
        print(dic['path'])

    if dic['out_h5'] == 0:     # txt output, runs sorted in directories
        path = path + dirname + os.sep
        add = len(glob.glob(path + "*"))    # glob.glob(pattern) returns list of paths matching pattern
        if not os.path.exists(path):
            os.mkdir(path)
            add = 0
    else:     # hdf5 output, file extended for new run
        existingRuns = create_or_check_swarmdyn_hdf5(path, dirname, dic)

    runs = max(0, runs - existingRuns)
    t0 = pytime.time()
    for i in range(runs):
        if dic['out_h5'] == 0 and dic['output_mode'] != 2: # txt output, runs sorted in directories
            dic['path'] = path + str(i+add) + os.sep
            if verb:
                print("run: ", i, dic['path'])
            if not os.path.exists(dic['path']):
                os.mkdir(dic['path'])
        command = sd.dic2swarmdyn_command(dic)
        os.system(command)
    t1 = pytime.time()
    return [t1 - t0, paratuple]


def DoubleSimulationScan(params, runs=20,
                         scantype='seq', no_processors=24, verb=None):
    outpath = params['path']
    if verb is None:
        verb = False
    dic = params.copy()  # to not modify original part
    # ensure no redundant computation: also causing Bug in multiprocessing
    dic['para_values0'] = np.unique(dic['para_values0'])
    dic['para_values1'] = np.unique(dic['para_values1'])
    #  Runs Simulation and Saves it in Folder
    checkRuns = True
    if not os.path.exists(outpath):
        os.mkdir(outpath)
        checkRuns = False
    # check how many runs have already been done
    if (scantype == 'seq'):
        for p in dic['para_values0']:
            for q in dic['para_values1']:
                MultipleRunSwarmdyn(dic, runs, [p, q])
    else:
        parallel_run_func = partial(MultipleRunSwarmdyn, dic,
                                    runs)
        comp_pool = mp.Pool(no_processors)
        para_values = get_scanTuples(dic)
        print(para_values)
        out = comp_pool.map(parallel_run_func, para_values, chunksize=1)   # para_values passed as tuple
        comp_pool.terminate()
    return


def get_scanTuples(dic):
    '''
    Create parameter-tuples which are passed to parallel_run_func
    for 1D scan: returns a list
        2D scan: returns a list of pair-tuples
        3D scan: returns a list of triplet-tuples
    '''
    para_values = []
    for i in dic['para_values0']:
        pa_tuples = list(zip([i]*len(dic['para_values1']), dic['para_values1']))
        if 'para_name2' in dic.keys():
            for j in dic['para_values2']:
                para_values += [list(tup) + [j] for tup in pa_tuples]
        else:
            para_values += pa_tuples
    return para_values


def Scanner(params, scantype, no_processors, runs):
    outpath = params['path']
    # ensure reproducability
    sd.reproducible(outpath)

    t_start = pytime.time()
    # params used:
    for key in params:
        print(key + " = " + str(params[key]))
    # Run Simulation and Move Files
    ############################################################
    DoubleSimulationScan(params, runs=runs, scantype=scantype,
                         no_processors=no_processors)

    t_end = pytime.time()
    print(outpath + ' finished in ', (t_end - t_start)/60, ' minutes')


def run_func4list(func, outpaths, base_para, para_changes,
                  scantype, no_processors, runs):
    '''
    Helper functions which just repeatedly calls the
    same functions with modified parameters
    1. "base_para" contains the default parameters
    1.1. "base_para" is copied so the original is not modified
    2. "para_changes" is a dictionary with
        keys = parameters of "base_para" to be changed
        and
        values = list of the parameter values to substitute
    3. "outpaths" list of dictionary names where the repspective results will be
        to be saved
    '''
    for i in range(len(outpaths)):
        if not os.path.exists(outpaths[i]):
            os.mkdir(outpaths[i])
        newpara = base_para.copy()
        for k in para_changes.keys():
            assert len(para_changes[k]) == len(outpaths), (
                    'len(para_changes[{}]) != len(outpaths), {}'.format(k, para_changes[k]))
            if para_changes[k][i] is not None:
                newpara[k] = para_changes[k][i]
        f_paras = Path(outpaths[i]) / 'paras_independent.pkl'
        time_stamp = pytime.strftime('%y%m%d%H')
        f_para_stamped = f_paras.parent / ('paras_independent' + time_stamp + '.pkl')
        pickle.dump(newpara, f_paras.open('wb'))    # paras which will be used
        pickle.dump(newpara, f_para_stamped.open('wb'))
        print((("Outpath = {}, Processors = {}, Runs = {}, "
               ).format(outpaths[i], no_processors, runs)))
        func(newpara, scantype,
             no_processors, runs)


possible_modes = ['burst_coast', 'sinFisher', 'mulFisher', 'natPred', 'natPredNoConfu', 'exploration']
# This function is executed if the script is called from terminal e.g. using > python RunSimulationScan.py
if __name__ == '__main__':
    ###########################################################
    #
    #How this program works:
    #
    #The ./swarmdyn is run while scanning two parameters of your choice. The output is saved to a local
    #folder and then copied to another (network) folder after the scan is completed. You can choose how
    #many runs you would like to do, and at what number you would like to start at (default is 0).
    #
    #Additionally you may choose to run the program in parallel, using as many cores as you would like.
    ###########################################################

    # to make sure new code is used!!!
    os.system("make")

    t0 = pytime.time()
    #Scan Values
    scantype = 'para'# 'para': parallel; 'seq':sequential (Better for debugging)
    no_processors = 35  # 66: for itb cluster, 40: for CHIPS
    runs = 20   # 40 

    # Simulation BASE-Parameter Values
    equi_time = 100
    record_time = 300
    # AllOptiModes = ['WithAlphaSimplestI', 'WithAlphaSimplestII', 'WithAlphaAsPNAS']
    errOrder = False

    para_changes = dict()
    base_para = sd.get_base_params(record_time, trans_time=equi_time)
    rep = 1 # repeat: create the correct length of lists
    base_para['output_mode'] = 1 # output: 0=only mean, 1=mean+particle, 2=only particle
    # non-explorativ: parameters based on fitting and optimization
    para_changes['para_name0']   = rep * ['N']
    ra0 = np.arange(1, 9, step=1)
    ra1 = 2**np.arange(0, 8, step=1)
    para_changes['para_values0'] = rep * [np.unique(np.append(ra0, ra1))]
    para_changes['para_name1'] = rep * ['beta']
    para_changes['para_values1'] = rep * [0.025 * 2**np.arange(9, step=0.25)]
    para_changes['speed0'] = rep * [2]

    # OUTPATH-name
    # d_save = Path(os.path.realpath(__file__)).parent
    d_extra = 'Out_S__2D_Beta_N'
    d_save = Path('/home/klamser/Seafile/PoSI_EmergentStructure/Data/SimuPPK')
    d_save = Path('/mnt/data3/PoSI_VariableSpeed/') / d_extra
    d_save.mkdir(exist_ok=True)
    keys = [k for k in para_changes.keys() if 'values' not in k]
    keys.sort()
    para_changes['path'] = []
    for i in range(rep):
        f_name = 'Scan'
        for k in keys:
            f_name += '__{}_{}'.format(k, para_changes[k][i])
        para_changes['path'].append(str(d_save / f_name))
    outpaths = para_changes['path']
    # BELOW para_cahnges are not part of the name
    para_changes['motivation'] = ['Check if the transition at N=4 also exists for alg_strength=0!']

    # create Data:
    run_func4list(Scanner, outpaths, base_para, para_changes,
                  scantype, no_processors, runs)

    t1 = pytime.time()
    print('total time:', t1-t0)

    # Analyze Data, create Plots
    import EvaluateScan as es
    # TODO: pass number of processors
    for path in outpaths:
        es.Evaluate(path)
