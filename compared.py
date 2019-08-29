import matplotlib.pylab as plt
import numpy as np
import math
import fire
import sys
import os
import re

from os import path
from math import log, log10
from glob import glob
from itertools import chain
from collections import namedtuple
from matplotlib.widgets import Slider


oned_glob = r"C:\Users\ostanin.igor\d\experiments\*i(*\output\*.dat"
tred_glob = r"C:\Users\ostanin.igor\d\3d\*\output\*.dat"


ischemia_procs_3d = [
    3.0,  # unused. only for normalization
    0.2968, 0.6342, 0.5398, 1.0092, 0.7022,
    1.2489, 0.0016, 0.0429, 0.0001, 0.0002,
    0.9852, 1.7039, 0.8434, 1.4418, 1.2413,
    2.2188, 1.2224, 2.3013, 0.3693, 0.8092,
    0.3284, 1.0041, 0.5488, 1.0380, 0.3015,
    0.6906, 1.1789, 2.0570, 1.3322, 2.4703,
    0.0387, 0.1802, 0.6024, 1.1471
]
opacity = sorted(enumerate(ischemia_procs_3d), key=lambda e: e[-1])

point_opacity = ischemia_procs_3d[:]

j = 0
for i, e in opacity:
    point_opacity[i] = j / len(ischemia_procs_3d)
    j += 1

print(point_opacity)

def divisors(n):
    divs = [1]
    for i in range(2, int(math.sqrt(n))+1):
        if n % i == 0:
            divs.extend([i, n/i])
    divs.extend([n])
    return list(set(divs))


def rect(size: int):
    while True:
        sq = int(math.sqrt(size))
        divs = sorted(divisors(size))
        d = {min(e, size // e) for e in divs}

        md = max(d)
        if md >= sq:
            return md, size // md

        size += 1


def set_ecg(axe):
    axe.set_xlabel('Время (мс)')
    axe.set_ylabel('Напряжение (мВ)')


def is_dat(file):
    return file.endswith('.dat')


def load(file_mask):
    for file in filter(is_dat, glob(file_mask)):
        yield np.loadtxt(file).transpose()[1]  # [0] always 1..2..3...


def filename(file):
    return path.split(file)[-1]


def name2coords(ecg_name):
    pecg = "PseudoEcgFromElectrodeAt_"
    return tuple(map(float, ecg_name[:-4].replace(pecg, "").split('_')))


def iterate(array):
    yield from (x.any() for x in np.nditer(array, flags=['refs_ok']))


def power_up_1d(ischemia_proc):
    return {
        0.01: 1,
        0.05: 0.6,
        0.1: 0.4,
        0.2: 0.2,
        0.3: 0.1
    }.get(ischemia_proc, 0)


"""
    d - dimension
    p - path,
    e - electrode,
    ip - ischemia position
    i - ischemia proc
    a - activation proc
    n - test number
"""
Params = namedtuple('Params', ['d', 'p', 'e', 'ip', 'i', 'a', 'n'])


def loadoned(dat_glob):

    dir_ptrn = re.compile(
        r'(?P<n>\d+)\.\sexperiment_i\((?P<x>\d\.\d{1,2}),\s(?P<y>\d)\)_a\((?P<z>\d\.\d)\)')

    for dat_file in filter(is_dat, glob(dat_glob)):
        _, y = np.loadtxt(dat_file).transpose()
        props = {}

        m = dir_ptrn.search(dat_file)
        if not m:
            print("warn: cannot parse 1d ecg from: " + dat_file)
            continue

        d = 1
        p = dat_file
        e = name2coords(path.split(dat_file)[-1])
        i, ip, a = map(float, (m.group('x'), m.group('y'), m.group('z')))
        n = int(m.group('n'))

        yield y, Params(d, p, e, ip, i, a, n)


def loadtred(dat_glob):
    dir_ptrn = re.compile(r"experiment\s(\d+)")

    for dat_file in filter(is_dat, glob(dat_glob)):
        _, y = np.loadtxt(dat_file).transpose()
        props = {}

        m = dir_ptrn.search(dat_file)
        if not m:
            print("warn: cannot parse 3d ecg from: " + dat_file)
            continue

        d = 3
        p = dat_file
        e = name2coords(path.split(dat_file)[-1])
        n = int(m.group(1))  # experiment №
        i = ischemia_procs_3d[n]

        yield y, Params(d, p, e, 0, i, 0, n)


def normalize(v):
    norm = np.linalg.norm(v, ord=1)
    if norm == 0:
        norm = np.finfo(v.dtype).eps
    return v / norm


def sign(y):
    return np.sum(y[20:50])


def fit(y):
    y = normalize(y)
    if sign(y) < 0:
        y = -y
    return y


def create_plot():
    fig = plt.figure(figsize=(16, 8))

    all_axes = fig.subplots(2, 4, sharex=True, sharey=False)

    axes = all_axes[0, :]

    for ax in iterate(all_axes[1:, :]):
        ax.remove()

    config_axes = fig.subplots(13, 2)

    for ax in iterate(config_axes[0:7, :]):
        ax.remove()

    config_axes = config_axes[7:, :]

    o_axes = config_axes[:, 0]
    t_axes = config_axes[:, 1]

    config_axes[0, 0].set_title('1D config')
    config_axes[0, 1].set_title('3D config')

    return axes, o_axes, t_axes


class SliderVariable:

    def __init__(self, ax, name, iterable):
        self._items = tuple(iterable)
        l = len(self._items)

        self.on_changed = None
        self.slider = Slider(
            ax, name,
            valmin=0,
            valmax=l - 1,
            valinit=0,
            valstep=1,
            valfmt='%i'
        )
        self.value = self._items[0]
        self.slider.on_changed(lambda x: self._handle_changed(int(x)))

    def _handle_changed(self, idx):
        self.value = self._items[idx]

        if self.on_changed:
            self.on_changed(self.value)

    def install(self, callback):
        self.on_changed = callback


def main():

    global od, td
    od, td = loadoned(oned_glob), loadtred(tred_glob)
    od = [(fit(x), p) for x, p in od]
    td = [(fit(x), p) for x, p in td]

    global axes, o_axes, t_axes
    axes, o_axes, t_axes = create_plot()

    o_ax1, o_ax2, o_ax3, o_ax4, o_ax5, o_ax6 = o_axes

    s1 = SliderVariable(o_ax1, "ischemia position", [2, 3, 1])
    s2 = SliderVariable(o_ax2, "ischemia proc", [0.01, 0.05, 0.1, 0.2, 0.3])
    s3 = SliderVariable(o_ax3, "activation proc", [0.1, 0.2, 0.3])

    o_ax4.remove()
    o_ax5.remove()
    o_ax6.remove()

    t_ax1, t_ax2, t_ax3, t_ax4, t_ax5, t_ax6 = t_axes

    t1 = SliderVariable(t_ax1, "experiment #", list(
        sorted(range(1, 35), key=lambda i: ischemia_procs_3d[i])))
    t2 = SliderVariable(t_ax2, "electrode", sorted({p.e for _, p in td}))

    t_ax3.remove()
    t_ax4.remove()
    t_ax5.remove()
    t_ax6.remove()

    d1_i_pos = s1
    d1_i_proc = s2
    d1_a_proc = s3

    d3_testno = t1
    d3_electrode = t2

    need_backup = False
    route_ax = {e: ax for e, ax in zip(sorted({p.e for _, p in od}), axes)}

    def draw(*args):
        nonlocal need_backup
        if need_backup:
            backup_camera = [ax.axis() for ax in axes]

        for ax in axes:
            ax.clear()

        for oy, op in od:
            ax = route_ax[op.e]
            if op.e != (2.1, 0, 0):
                continue


            if op.i <= d1_i_proc.value and op.ip == d1_i_pos.value and op.a == d1_a_proc.value:
                axes[1].plot(oy, color=(0, 0, 1, power_up_1d(op.i)))
                axes[2].plot(np.diff(oy, 1), color=(0, 0, 1, power_up_1d(op.i)))
                axes[3].plot(np.diff(oy, 2), color=(0, 0, 1, power_up_1d(op.i)))

            ass = []
            for ty, tp in td:
                if ischemia_procs_3d[tp.n] <= ischemia_procs_3d[d3_testno.value] and tp.e == d3_electrode.value:
                    alpha = 1 - point_opacity[tp.n]
                    axes[1].plot(ty,             color=(0, alpha, 0), zorder=1-alpha)
                    axes[2].plot(np.diff(ty, 1), color=(0, alpha, 0), zorder=1-alpha)
                    axes[3].plot(np.diff(ty, 2), color=(0, alpha, 0), zorder=1-alpha)

        if need_backup:
            for ax, restore in zip(axes, backup_camera):
                ax.axis(restore)
        need_backup = True

    s1.install(draw)
    s2.install(draw)
    s3.install(draw)
    t1.install(draw)
    t2.install(draw)

    draw()
    plt.show()


if __name__ == '__main__':
    main()
