import modfinder
from pprint import pprint

moddesc = [
    ('Lys', 128.1726, 2),
    ('Man5', 1217.0900, 2),
    ('G0', 1299.1940, 2),
    ('G0F', 1445.3355, 2),
    ('G1', 1461.3349, 2),
    ('G1F', 1607.4763, 2),
    ('G2', 1623.4757, 2),
    ('G2F', 1769.6172, 2),
    ('G2FSA1', 2060.8722, 2),
    ('G2FSA2', 2352.1273, 2)
]

masses = [
    2677.90428,
    2747.24798,
    2890.91988,
    3050.54488,
    3213.87298,
    3375.70108,
    3532.31048,
    3670.10738
]

mods = [(mod[1], mod[2]) for mod in moddesc]
for mass in masses:
    print(mass),
    result = modfinder.examine_modifications(mods, mass, 5.0)
    for r in result:
        pprint(list(zip([i[0] for i in moddesc],r)))
