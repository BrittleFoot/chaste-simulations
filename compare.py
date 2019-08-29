import matplotlib.pylab as plt
import numpy as np
import fire
import sys
import os
import re

from matplotlib.widgets import Slider
from os import path
from math import log, log10
from glob import glob


ischemia_procs = [
    0.0,  # unused. only for normalization
    0.2968,
    0.6342,
    0.5398,
    1.0092,
    0.7022,
    1.2489,
    0.0016,
    0.0429,
    0.000000001,
    0.0002,
    0.9852,
    1.7039,
    0.8434,
    1.4418,
    1.2413,
    2.2188,
    1.2224,
    2.3013,
    0.3693,
    0.8092,
    0.3284,
    1.0041,
    0.5488,
    1.0380,
    0.3015,
    0.6906,
    1.1789,
    2.0570,
    1.3322,
    2.4703,
    0.0387,
    0.1802,
    0.6024,
    1.1471
]

max_proc = max(ischemia_procs) + 0.01

normed_opacity = [(x+0.01)/max_proc for x in ischemia_procs]


def is_dat(file):
    return file.endswith('.dat')


def load(file_mask):
    for file in filter(is_dat, glob(file_mask)):
        yield np.loadtxt(file).transpose()


def filename(file):
    return path.split(file)[-1]

def power_up(alpha):
    ls = '-'
    if alpha < 0.01:
        alpha *= 10
        ls = '--'
    if alpha < 0.01:
        alpha *= 10
        ls = ':'
    return alpha, ls


def plot_by_electrode(file_mask):
    charts = {}
    for file in filter(is_dat, glob(file_mask)):
        charts.setdefault(filename(file), []).append(
            (file, np.loadtxt(file).transpose()))

    keys = list(charts.keys())[0:]

    # the most representative keys
    # keys = [
    #     keys[0],
    #     keys[1],
    #     keys[2],
    #     keys[6]
    # ]
    
    ptrn = re.compile("experiment\s(\d+)")

    mfig, maxes = plt.subplots(1, 1, figsize=(15, 8))
    maxes.set_visible(False)

    axes = mfig.subplots(3, 5, sharex='all')

    [a.set_visible(False) for a in axes[2, :]]

    ax = mfig.add_subplot(313)


    def update(experiment_filter):
        for key, axe in zip(keys, list(np.nditer(axes, flags=['refs_ok']))):
            axe = axe.any()
            axe.set_title(key[13:-4])

            for file, (x, y) in charts[key]:
                m = ptrn.search(file)
                if not m:
                    continue

                experiment_no = int(m.group(1))

                if experiment_filter and experiment_no != experiment_filter:
                    continue


                alpha = normed_opacity[experiment_no]
                
                alpha, ls = power_up(alpha)
                
                axe.plot(x, y, color=(0.4, 0.6, 0.4, alpha), ls=ls)



    def slider_update(value):
        list(map(plt.Axes.cla, axes[0, :]))
        list(map(plt.Axes.cla, axes[1, :]))
        update(int(value))

    global slider
    slider = Slider(ax, 'experiment #', 0, 34, 0, '%i', valstep=1)
    slider.on_changed(slider_update)
    # ax.set_visible(False)

    update(0)
    plt.show()



plot_by_electrode(r"C:\Users\ostanin.igor\d\3d\experiment * electrode_at_*\output\PseudoEcgFromElectrodeAt_*.dat")
# if __name__ == '__main__':
#     fire.Fire()
