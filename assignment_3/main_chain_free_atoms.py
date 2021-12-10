# Author: Stefano Ribes

import argparse
import math
from itertools import combinations, tee, groupby

def dot(a, b):
    x = a.x - b.x
    y = a.y - b.y
    z = a.z - b.z
    return x * x + y * y + z * z

def angle3d(a, b, c):
    v1 = (a.x - b.x, a.y - b.y, a.z - b.z)
    v2 = (c.x - b.x, c.y - b.y, c.z - b.z)

    v1_mag = math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])
    v1_norm = (v1[0] / v1_mag, v1[1] / v1_mag, v1[2] / v1_mag)

    v2_mag = math.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])
    v2_norm = (v2[0] / v2_mag, v2[1] / v2_mag, v2[2] / v2_mag)

    res = v1_norm[0] * v2_norm[0] + v1_norm[1] * v2_norm[1] + v1_norm[2] * v2_norm[2]
    angle_rad = math.acos(res)

    return math.degrees(angle_rad)

class Atom(object):
    """Atom class"""
    def __init__(self, id, x, y, z, label=' CA '):
        super(Atom, self).__init__()
        self.id = int(id)
        self.x = x
        self.y = y
        self.z = z
        self.label = label

    def distance_from(self, atom):
        assert type(atom) == Atom, f'ERROR. Input argument must be an Atom class (got {type(atom)} instead). Exiting'
        return math.sqrt(dot(self, atom))

    def __eq__(self, a):
        return self.id == a.id

    def __repr__(self):
        # return self.label.strip(' ') + f' atom n.{self.id} at ({self.x}, {self.y}, {self.z})'
        # return self.label.strip(' ') + f' atom n.{self.id}'
        return f'{self.id}'


def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def main():
    test_file = 'data_q2.txt'
    atoms = []
    with open(test_file) as f:
        for i, line in enumerate(f):
            atom_descr = tuple([float(x) for x in line.rstrip().split()])
            atoms.append(Atom(*atom_descr))

    dist_neighbours = 4.5
    dist_adjecent = 3.8

    within = {}
    within[dist_neighbours] = {}
    within[dist_adjecent] = {}

    for a, b in combinations(atoms, 2):
        dist = a.distance_from(b)
        for distance_threshold in within.keys():
            if dist < distance_threshold and dist > 0:
                within[distance_threshold][dist] = (a, b)

    adjacent_atoms = []
    for dist, (a, b) in within[dist_adjecent].items():
        # print(dist)
        adjacent_atoms.append(a.id)
        adjacent_atoms.append(b.id)

    # adjacent_atoms = list(adjacent_atoms)
    adjacent_atoms = list(dict.fromkeys(adjacent_atoms))

    print(adjacent_atoms)
    print(len(adjacent_atoms))

    within['angle'] = []
    for a, b, c in combinations(atoms, 3):
        if a.distance_from(b) < dist_adjecent and b.distance_from(c) < dist_adjecent:
            angle = angle3d(a, b, c)
            if 70 <= angle <= 180:
                # Give a score of 1 for a 30 degree from 110 degree.
                score = abs(110 - angle) / 30
                within['angle'].append((score, a, b, c))
    print(within['angle'])

if __name__ == '__main__':
    main()