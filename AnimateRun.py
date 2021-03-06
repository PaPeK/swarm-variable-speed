#!/usr/bin/env python
import numpy as np
import matplotlib
if __name__ == '__main__':
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import patches
import h5py
import pandas as pd
from pathlib import Path
from functools import partial
from optparse import OptionParser
import os


def smooth2D(dat, std, min_periods=None, smotype=None):
    '''
    only wrapper for smooth1D for 2 dimensions
    -> see smooth1D for details
    '''
    assert len(dat.shape) == 2, 'data has wrong dimension'
    smodat = np.empty(dat.shape, dtype=float)
    for i, da in enumerate(dat.T):
        smodat[:, i] = smooth1D(da, std, min_periods=min_periods,
                                smotype=smotype)
    return smodat


def smooth1D(dat, std=1, window=None, min_periods=1, smotype='gaussian'):
    '''
    INPUT:
        dat.shape(T)
        std double
            - standard deviation of gaussian kernel used
            - if window=None: int(6*std) = window size
        window int
            number of datapoints in a moving window
        min_periods int
            - minimum # of datapoints in window
            - if less data than min_periods -> None
        smotype string
            up to now only 'gaussian' is implemented
    '''
    if window is None:
        window = int(np.round(6*std))
    # use pandas to smooth
    smodat = pd.Series(dat)
    smodat = smodat.rolling(window=window, win_type=smotype,
                            center=True, min_periods=int(np.round(min_periods))
                            ).mean(std=std)
    return smodat


def pavas2colors(pavas):
    colors = np.squeeze(pavas)
    colors -= colors.min()
    colors /= colors.max()
    return colors


def setDefault(x, val):
    if x is None:
        x = val
    return x


def fixedLimits(dat, ax, extra=0):
    pos = dat.dat.T
    ax.set_xlim(np.nanmin(pos[0])-extra, np.nanmax(pos[0])+extra)
    ax.set_ylim(np.nanmin(pos[1])-extra, np.nanmax(pos[1])+extra)


def updateDelay(func):
    '''
    To use this decorator, the class must have
        create, delay, keyVar
        HOWEVER not working in ipython notebook -> do everything without header
    '''
    def dum(selfi, s):
        pass
    def inner(selfi, s):
        s *= 1
        s -= selfi.delay
        if s >= 0:
            if not hasattr(selfi, selfi.keyVar):
                selfi.create()
            return func(selfi, s)
        else:
            return dum(selfi, s)
    return inner


class datCollector:
    '''
    cleans the data by setting data to nan if only contains 0's
        nans are not plotted
    - collects also properties relevant for plotting/animation
        (delay, tail_length, colors)
    INPUT:
        dat.shape(T, N, M)
        OPTIONAL:
            colors
                either color-string
                OR
                2- OR 3-D array of matplotlib.pyplot.cmap colors
                if 2D: agents keep color over whole simulation
                if 3D: agents change color at every time-step
    '''
    def __init__(self, dat, delay=None, colors=None, tail_length=None, radius=None):
        self.delay = setDefault(delay, 0)
        self.dat = setDefault(dat, None)
        if dat is not None:
            dat = dat.T
            self.dat = dat[:2].T
            dead = np.where( np.sum(self.dat, axis=-1) == 0 )
            self.dat[dead] = np.nan
        # some default parameter
        self.colors = setDefault(colors, 'k')
        self.tail_length = setDefault(tail_length, 2)
        self.radius = setDefault(radius, 1/4) # radius of marker


class headAndTail:
    def __init__(self, dat, ax, marker=None):
        '''
        INPUT:
            data datCollector-object
            ax matplotlib.axes.AxesSubplot object
                Subplot where the data is plotted
        OUTPUT:
            heads
            tails
        '''
        self.scatter = False # scatter is slow BUT allows varying color
        self.delay = dat.delay
        colors = dat.colors
        if type(dat.colors[0]) == np.ndarray:
            self.scatter = True
            if len(dat.colors.shape) == 3:
                self.colors = dat.colors # of self.colors exists: update color every time-step
                colors = colors[0]
        self.pos = dat.dat
        self.tail_length = dat.tail_length
        self.ax = ax
        self.sizeRaw = dat.radius
        self.keyVar = 'heads'
        marker = setDefault(marker, 'o')
        if self.scatter:
            self.createFunc = partial(ax.scatter, marker=marker, c=colors)
        else:
            self.createFunc = partial(ax.plot, c=colors, marker=marker, linestyle='')


    def sizeAdjust(self):
        trans = self.ax.transData
        pixels = trans.transform([(0, 1), (1, 0)]) - trans.transform((0, 0))
        xpix = pixels[1].max()
        self.size = self.sizeRaw * xpix / 77. * 100
        if self.scatter:
            self.size = self.size ** 2


    def create(self):
        self.sizeAdjust()
        posInit = self.pos[0]    # use Time=1 for initialization
        if self.scatter:
            self.heads = self.createFunc(posInit.T[0], posInit.T[1], s=self.size)
            self.tails = self.createFunc(posInit.T[0], posInit.T[1],
                                         s=self.size/4, alpha=0.5)
        else:
            self.heads = self.createFunc(posInit.T[0], posInit.T[1], ms=self.size)[0]
            self.tails = self.createFunc(posInit.T[0], posInit.T[1],
                                         ms=self.size/2, alpha=0.5)[0]


    @updateDelay
    def update(self, s):
        self.sizeAdjust()
        endtail = s - self.tail_length
        if endtail < 0:
            endtail = 0
        if self.scatter:
            tail = self.pos[endtail:s].reshape(-1, 2)
            self.heads.set_offsets(self.pos[s])
            self.tails.set_offsets(tail)
            self.heads.set_sizes(len(self.pos[s]) * [self.size])
            self.tails.set_sizes(len(tail) * [self.size / 4])
            if hasattr(self, 'colors'):
                self.heads.set_color(self.colors[s])
                self.tails.set_color(self.colors[s])
        else:
            self.heads.set_data(self.pos[s].T[0], self.pos[s].T[1])
            self.tails.set_data(self.pos[endtail:s].T[0],
                                self.pos[endtail:s].T[1])
            self.heads.set_markersize(self.size)
            self.tails.set_markersize(self.size / 2)
        return self.heads, self.tails


class Limits4Pos:
    '''
    adjust the limits of the axis to the positions of the agents
    -> always see all agents and not more
    '''
    def __init__(self, positions, ax):
        '''
        Assumes that the start-positions can be different
        but all end synchrounously
        INPUT:
            positions [posA, posB, ....]
                list of datCollector objects
                each datCollector has
                    datCollector.dat.shape = [Time, Nagents, 2] OR [Time, 2]
        '''
        self.delay = np.min([pos.delay for pos in positions])
        self.ax = ax
        times = [len(pos.dat) for pos in positions]
        mins = np.zeros((max(times), 2)) * np.nan
        maxs = np.zeros((max(times), 2)) * np.nan
        lims = [mins, maxs]
        for i, limitFunc in enumerate([np.nanmin, np.nanmax]):
            lim = lims[i] # either mins or maxs
            for k, pos in enumerate(positions):
                dat, delay = pos.dat * 1, pos.delay - self.delay
                if len(dat.shape) < 3:
                    dat = np.expand_dims(dat, axis=1)
                compare = np.concatenate((dat, np.expand_dims(lim[delay:],
                                                              axis=1)), axis=1)
                lim[delay:] = limitFunc(compare, axis=1)
        extra = 1
        self.xMin = mins[:, 0] - extra
        self.yMin = mins[:, 1] - extra
        self.xMax = maxs[:, 0] + extra
        self.yMax = maxs[:, 1] + extra
        self.keyVar = 'xMin'


    @updateDelay
    def update(self, s):
        self.ax.set_xlim(self.xMin[s], self.xMax[s])
        self.ax.set_ylim(self.yMin[s], self.yMax[s])
        return self.ax


class taskCollector:
    '''
    collects tasks (=objects with update(s) method returning artist objects)
    '''
    def __init__(self):
        self.taskList = []

    def append(self, task):
        self.taskList.append(task)

    def update(self, s):
        '''
        returns all output summarized
        '''
        artistObjs = ()
        for task in self.taskList:
            artist = task.update(s)
            if type(artist) != tuple:
                artist = tuple([artist])
            # print(type(artistObjs), type(artist), task.__class__)
            artistObjs += artist
        return artistObjs


def UpdateViaAnimation(fig, tasks, tmin, tmax, fps=None, dpi=None,
                       mode=None, name=None, repeat=None):
    fps = setDefault(fps, 15)
    dpi = setDefault(dpi, 300)
    mode = setDefault(mode, 'normal')
    name = setDefault(name, 'Animation')
    repeat = setDefault(repeat, True)
    interval = 1000*(1/fps)
    anim = animation.FuncAnimation(fig, tasks.update, interval=interval,
                                   frames=range(tmin-1, tmax), repeat=repeat)
    if mode == 'movie':
        anim.save(name + '.mp4', writer='ffmpeg', dpi=dpi)
    elif mode == 'gif':
        anim.save(name + '.gif', writer='imagemagick', dpi=dpi)
    else:
        plt.show()


def UpdateViaDraw(fig, tasks, tmin, tmax, fps=None, dpi=None,
                  mode=None, name=None, repeat=None):
    '''
    TODO: repeat is only dummy (ensure same parameters as UpdateViaAnimation)
    '''
    fps = setDefault(fps, 15)
    dpi = setDefault(dpi, 300)
    mode = setDefault(mode, 'normal')
    name = setDefault(name, 'Animation')
    plt.show(block=False) # necessary for normal
    maxFps = 20 # this is system specific
    interval = 1/fps - 1/maxFps
    for s in range(tmin, tmax):
        _ = tasks.update(s)
        if(mode != 'movie'):
            fig.canvas.draw()
            if interval > 0:
                sleep(interval)
        if(mode != 'normal'):
            fig.savefig('f%06d.jpg' % s, dpi=dpi)
    if (mode == 'movie'):
        com0 = 'mencoder mf://f*.jpg -mf fps={}:type=jpg'.format(fps)
        com1 = ' -vf scale=-10:-1 -ovc x264 -x264encopts'
        com2 = ' bitrate=4000 -o {}.mp4'.format(name)
        os.system(com0+com1+com2)
        os.system("rm f*.jpg")
    elif(mode == 'gif'):
        print('gif-representation not implemented')


def main(mode=None, size=None):
    if size is not None:
        size = int(size)
    if mode is None:
        mode = 'normal' # 'normal', 'pictures', 'movie', 'gif'
    # # Definitions
    tail_length = 2 # 3
    fps = 15
    dpi = 200
    sizePrey = 1/8
    name = 'Animation'
    cmap = plt.get_cmap('bwr') # alternatives 'bwr', 'Reds'
    folder = Path.cwd()

    # # Load data
    f_h5 = folder / 'out_xx.h5'
    if f_h5.exists:
        with h5py.File(f_h5) as fh5:
            preys = datCollector( np.array(fh5['/part']), radius=sizePrey)
            s = np.array(fh5['/part'][:, :, 2:4])
    else:
        print('{} does not EXIST'.format(f_h5))
    # pdb.set_trace()
    pava_in_name = folder / 'pava_in_xx.in'
    colors = 'k'
    if pava_in_name.exists():
        pavas = np.loadtxt(str(pava_in_name))
        colors = pavas2colors(pavas)
    # color from velocity
    if False:
        s = np.diff(preys.dat, axis=0)
        s = np.sqrt(np.sum(s**2, axis=-1))
        s = np.vstack((s, s[-1]))
        s = smooth2D(s, 1)
        s = pavas2colors(s)
        print(s.shape, s.min(), s.max())
        colors = cmap(s)
        print('colors.shape', colors.shape)
        preys.colors = colors
    # get info from files
    time, N, _ = preys.dat.shape

    # # Animation
    f, ax = plt.subplots(1)
    # ax.axis('off')
    ax.set_aspect('equal')
    # Collect update-tasks
    tasks = taskCollector()
    if N == 1 or size is not None:
        fixedLimits(preys, ax, extra=2)
        if size is not None:
            rec = patches.Rectangle((0, 0), size, size, color='k', fill=False)
            ax.add_artist(rec)
    # elif size is not None:
    #     mine, maxe = 0, size
    #     if np.min(preys.dat) < 0:
    #         mine, maxe = -size, size
    #     ax.set_xlim(mine, maxe)
    #     ax.set_ylim(mine, maxe)
    else:
        positions = [preys]
        tasks.append( Limits4Pos(positions, ax) )
    tasks.append( headAndTail(preys, ax) )
    # animation
    if mode == 'movie':
        UpdateViaDraw(f, tasks, 0, time, fps=fps, repeat=False, mode=mode, name=name)
    else:
        UpdateViaAnimation(f, tasks, 0, time, fps=fps, repeat=False, mode=mode, name=name)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-s", "--size", dest="size", default=None,
                      help="size of area", metavar=float)
    options, args = parser.parse_args()
    main(None, size=options.size)
