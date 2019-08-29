import matplotlib.pylab as plt
import numpy as np
import fire
import sys
import os
import re

from os import path
from math import log, log10
from glob import glob


ischemia_procs = [
    3.0,  # unused. only for normalization
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

max_proc = max(ischemia_procs)
print(max_proc)

normed_opacity = [x/max_proc for x in ischemia_procs]


def set_ecg(axe):
    axe.set_xlabel('Время (мс)')
    axe.set_ylabel('Напряжение (мВ)')


def is_dat(file):
    return file.endswith('.dat')


def load(file_mask):
    for file in filter(is_dat, glob(file_mask)):
        yield np.loadtxt(file).transpose()


def filename(file):
    return path.split(file)[-1]


def plot_by_electrode(file_mask):
    charts = {}
    for file in filter(is_dat, glob(file_mask)):
        charts.setdefault(filename(file), []).append(
            (file, np.loadtxt(file).transpose()))

    keys = list(charts.keys())[0:]

    # the most representative keys
    keys = [
        keys[0],
        keys[1],
        keys[2],
        keys[6]
    ]
    
    ptrn = re.compile(r"experiment\s(\d+)")

    fig, axes = plt.subplots(2, 2, sharex='all')

    for key, axe in zip(keys, list(np.nditer(axes, flags=['refs_ok']))):
        axe = axe.any()
        set_ecg(axe)

        en = key[13:-4].replace('_', ', ').replace("ElectrodeAt, ", "")
        axe.set_title("Позиция электрода (" + en + ")")

        for file, (x, y) in charts[key]:
            m = ptrn.search(file)
            if not m:
                continue

            alpha = normed_opacity[int(m.group(1))]
            axe.plot(x, y, color=(0.4, 0.6, 0.4, alpha))

    plt.show()


plot_by_electrode(
    r"C:\Users\ostanin.igor\d\3d\experiment * electrode_at_*\output\PseudoEcgFromElectrodeAt_*.dat")
# if __name__ == '__main__':
#     fire.Fire()
