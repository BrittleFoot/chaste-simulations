#! /usr/bin/python3
# -*- utf-8 -*-

import matplotlib.pylab as plt
import numpy as np
import fire
import sys
import os

from os import path


CHASTE_OUTPUT = os.environ['CHASTE_OUTPUT_DIR']
if not CHASTE_OUTPUT:
    print('Cannot find CHASTE_OUTPUT_DIR')
    sys.exit(1)

ECG_DEFAULT_OUTPUT = path.join(CHASTE_OUTPUT, "ChasteResults", "output")


def resolve_location(dat_path):
    if path.isabs(dat_path):
        print('Absolute path: ' + dat_path)
        return dat_path

    print('Path relatvie to %CHASTE_OUTPUT%: ' + dat_path)
    return path.join(ECG_DEFAULT_OUTPUT, dat_path)


def resolve(dat_path):
    full_path = resolve_location(dat_path)
    if not path.isfile(full_path):
        print('Expected file to exist: ' + full_path)
        sys.exit(2)

    if path.splitext(full_path)[1] != '.dat':
        print('File expected to be `.dat`')
        sys.exit(3)

    return full_path


def load(dat_file):
    ms_s = []
    ecg_s = []

    with open(dat_file, mode='r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue

            ms, ecg = tuple(map(float, line.split('\t')))

            ms_s.append(ms)
            ecg_s.append(ecg)

    return ms_s, ecg_s


def plot(dat_path):
    """
    :param dat_path: path to .dat file 
        either relatvie to %CHASTE_OUTPUT_DIR% directory or absolute

    :return: shows ecg plot with some comments
    """
    x, y = load(resolve(dat_path))
    plt.plot(x, y)

    plt.show()


def plot_all():
    """
    :param dat_path: plot all ecg-s from %CHASTE_OUTPUT_DIR%

    :return: shows ecg plot with some comments
    """

    legend = []

    with os.scandir(ECG_DEFAULT_OUTPUT) as it:
        for entry in it:
            if not (entry.is_file() and entry.name.endswith(".dat")):
                continue

            x, y = load(resolve(entry.name))

            title = path.splitext(entry.name)[0]
            if title.startswith("PseudoEcgFromElectrodeAt_"):
                title = title[25:]

            legend.append(title)

            plt.plot(x, y)

    plt.legend(legend)
    plt.show()


if __name__ == '__main__':
    fire.Fire({"single": plot, "all": plot_all})
