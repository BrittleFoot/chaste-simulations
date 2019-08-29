import matplotlib.pylab as plt
import numpy as np
import fire
import sys
import os
import re

from os import path
from math import log
from glob import glob


def is_dat(file):
    return file.endswith('.dat')

def set_ecg(axe):
    axe.set_xlabel('Время (мс)')
    axe.set_ylabel('Напряжение (мВ)')


def load(file_mask):
    for file in filter(is_dat, glob(file_mask)):
        yield np.loadtxt(file).transpose()


def filename(file):
    return path.split(file)[-1]



def plot_by_electrode(file_mask):
    charts = {}
    for file in filter(is_dat, glob(file_mask)):
        charts.setdefault(filename(file), []).append((file, np.loadtxt(file).transpose()))

    keys = list(charts.keys())

    assert len(keys) == 4

    fig, axes = plt.subplots(2, 2, sharex='all')

    list(map(lambda a: set_ecg(a.any()), np.nditer(axes, flags=['refs_ok'])))

    ptrn = re.compile(r'.*i\((?P<x>\d\.\d{1,2}),\s(?P<y>\d)\)_a\((?P<z>\d\.\d)\).*')

    for key, axe in zip(keys, list(np.nditer(axes, flags=['refs_ok']))):
        axe = axe.any()
        en = key[13:-4].replace('_', ', ').replace("ElectrodeAt, ", "")
        axe.set_title("Позиция электрода (" + en + ")")

        for file, (x, y) in charts[key]:
            m = ptrn.match(file)
            if not m:
                continue
            i_proc, i_pos, a_proc = map(float, (m.group('x'), m.group('y'), m.group('z')))

            r = int(i_pos == 1)
            g = int(i_pos == 2)
            b = int(i_pos == 3)

            a = {
                0.01: 1,
                0.05: 0.6,
                0.1: 0.4,
                0.2: 0.2,
                0.3: 0.1
            }[i_proc]

            ls = {
                0.3: '-',
                0.2: '--',
                0.1: ':',
            }[a_proc]

            if ls != '-':
                continue

            if r != 0:
                continue


            axe.plot(np.diff(y), color=(r, g, b, a), ls=ls)
    plt.show()



plot_by_electrode(r"C:\Users\ostanin.igor\d\experiments\*i(*\output\*.dat")
# if __name__ == '__main__':
#     fire.Fire()
