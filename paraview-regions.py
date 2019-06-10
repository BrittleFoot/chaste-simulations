from paraview.simple import *
from os import path

ischemia = [
    (-3, -248, -268),
    (-3, -248, -228),
    (-3, -208, -268),
    (-3, -168, -268),
    (-3, -168, -228),
    (37, -248, -268),
    (37, -248, -228),
    (37, -208, -268),
    (37, -208, -228),
    (37, -168, -268),
    (37, -168, -228),
    (77, -248, -268),
    (77, -248, -228),
    (77, -208, -268),
    (77, -208, -228),
    (77, -168, -268),
    (77, -168, -228)
]

rads = (13, 16)


heart = FindSource('heart.vtu')
ImplicitSphere = servermanager.implicit_functions.Sphere


def extract(center, radius):
    return ExtractCellsByRegion(heart, IntersectWith=ImplicitSphere(Radius=radius, Center=center))


folder = "C:\Users\ostanin.igor\d\paraview"

i = 0
for x, y, z in ischemia:
    for r in rads:
        i += 1

        print("Extracting", (x, y, z), r)

        ex = extract((x, y, z), r)
        # Show(ex)
        name = "%s.vtu" % i

        SaveData(path.join(folder, name), ex)

Render()
