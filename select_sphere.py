from paraview.simple import *
from collections import namedtuple
from itertools import permutations
from os import path

import operator as op

Point = namedtuple("Point", ["x", "y", "z"])


heart = FindSource('heart.vtu')
heart_center = Point(37.503335543345344, -205.7330789685239, -243.42471161095895)
full_radius = 65


ImplicitSphere = servermanager.implicit_functions.Sphere


def extract(center, radius):
    return ExtractCellsByRegion(heart, IntersectWith=ImplicitSphere(Radius=radius, Center=center))


def do(operator, v1, v2):
    return Point(operator(v1.x, v2.x), operator(v1.y, v2.y), operator(v1.z, v2.z))

def scale(v, scalar):
    return Point(v.x * scalar, v.y * scalar, v.z * scalar)


sizes = (15, 25)
step = 40

bounds = [-21.474571228027, 99.588241577148, -260.83679199219, -140.49102783203, -283.66015625, -183.43103027344]
box_scale = Point(0.7, 0.8, 0.7)


top = Point(*bounds[::2])
bot = Point(*bounds[1::2])

c = scale(do(op.add, top, bot), 0.5)
h = do(op.mul, do(op.sub, c, top), box_scale)

new_top = do(op.add, c, h)
new_bot = do(op.sub, c, h)


ignore = {7, 8, 9, 11, 33, 35}


folder = "C:\\Users\\ostanin.igor\\d\\paraview"

i = 0
for x in range(int(new_bot.x), int(new_top.x), step):
    for y in range(int(new_bot.y), int(new_top.y), step):
        for z in range(int(new_bot.z), int(new_top.z), step):

            curr_center = (x, y, z)

            for radius in sizes:
                i += 1

                if i in ignore:
                    continue

                print("Extracting", curr_center, radius)

                ex = extract(curr_center, radius)
                # Show(ex)
                filename = path.join(folder, "%0.2d. (%s, %s, %s) r=%s.vtu" % (i, x, y, z, radius))
                print("Saving: %s" % filename)
                SaveData(filename, ex)


Render()
