from paraview.simple import *


for ((name, port), proxy) in GetSources().items():
    if name != 'heart.vtu':
        Delete(proxy)