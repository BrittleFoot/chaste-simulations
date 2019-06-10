import json

define = "00000000"
ptrn_cx = "11111111"
ptrn_cy = "22222222"
ptrn_cz = "33333333"
ptrn_rad = "44444444"
ptrn_ichemia_applied = "true && true && true && true"



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


configuration = {}


with open("IschemiaTemplate.hpp") as input:
    text = input.read()

    i = 0
    for x, y, z in ischemia:
        for rad in rads:
            for apply in ["false"]:
                i += 1


                s_x = str(x / 10.0)
                s_y = str(y / 10.0)
                s_z = str(z / 10.0)
                s_rad = str(rad / 10.0)


                interpretation = text                         \
                    .replace(define, str(i))                  \
                    .replace(ptrn_cx, s_x)                    \
                    .replace(ptrn_cy, s_y)                    \
                    .replace(ptrn_cz, s_z)                    \
                    .replace(ptrn_rad, s_rad)                 \
                    .replace(ptrn_ichemia_applied, apply)


                filename = "%s_test.hpp" % (i,)

                configuration[i] = {
                    'id': i,
                    'x': x,
                    'y': z,
                    'z': y,
                    'r': rad,
                    'has_ischemia': apply == 'false'
                }
                
                with open("prepared/" + filename, mode='w') as output:
                    output.write(interpretation); 

with open('generated-test-config.json', mode='w') as config_fp:
    json.dump(configuration, config_fp, indent=2)


                


