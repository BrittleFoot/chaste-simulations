import matplotlib.pylab as plt
import numpy as np

X = np.array([
0.6491439041430234,
0.0063419999999999995,
0.540217442697094,
0.008092,
0.4806414212383382,
0.005397999999999999,
0.5758145059843119,
0.010092,
0.41999999999999993,
0.0070220000000000005,
0.5161956318234955,
0.012488999999999998,
0.3547868569224586,
1.6e-05,
0.976905989232415,
0.000429,
0.880417392568986,
1e-06,
0.9942264973081038,
0.002968,
0.6854633036777213,
2e-06,
0.9918350341907227,
0.009852,
0.4269380487242238,
0.017039,
0.2463643674383047,
0.008434,
0.4697799199074658,
0.014418,
0.3067467995025194,
0.012413,
0.35675302824912836,
0.022188,
0.14,
0.012223999999999999,
0.36166884247542275,
0.023013,
0.12415754841409976,
0.003693,
0.6491439041430234,
0.0063419999999999995,
0.540217442697094,
0.008092,
0.4806414212383382,
0.005397999999999999,
0.5758145059843119,
0.010092,
0.9942264973081038,
0.41999999999999993,
0.0070220000000000005,
0.5161956318234955,
0.012488999999999998,
0.3547868569224586,
1.6e-05,
0.976905989232415,
0.000429,
0.880417392568986,
1e-06
])



plt.hist(X ** (1/4))

plt.show()