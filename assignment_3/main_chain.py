import argparse
import math
from itertools import combinations, tee, groupby

import matplotlib.pyplot as plt

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
    test_file = 'test_q1.txt'
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
        for distance in within.keys():
            if dist <= distance and dist > 0:
                within[distance][dist] = (a, b)

    within['angle'] = []
    for a, b, c in combinations(atoms, 3):
        angle = angle3d(a, b, c)
        if 70 <= angle <= 180:
            # Give a score of 1 for a 30 degree from 110 degree.
            score = abs(110 - angle) / 30
            within['angle'].append((score, a, b, c))

    neighboring_atoms = list(within[dist_neighbours].values())
    # Init the chain with the first pair of atoms
    chain = [*neighboring_atoms.pop(0)]
    while len(neighboring_atoms):
        # While not all pairs are attached to the chain...
        for i, (x, y) in enumerate(neighboring_atoms):
            # Check extremes of the chain: if an atom in the pair matches the
            # extremes and the other one in the pair is not already in the
            # chain, then add the non-matching atom to the chain in the proper
            # position. Example:
            # 
            # chain: [2, 4, 3] -\
            #                    > updated chain: [2, 4, ( 3] + [5 )] 
            # pair: (3, 5) -----/
            # 
            # In this way, atoms that are close together stay close together.
            if x == chain[0] and y not in chain:
                chain = [y] + chain
                neighboring_atoms.pop(i)
            if x == chain[-1] and y not in chain:
                chain.append(y)
                neighboring_atoms.pop(i)
            if y == chain[0] and x not in chain:
                chain = [x] + chain
                neighboring_atoms.pop(i)
            if y == chain[-1] and x not in chain:
                chain.append(x)
                neighboring_atoms.pop(i)
            print(f'[DEBUG] Appending {(x, y)}, chain: {chain}')

    gold_order = [8, 10, 9, 7, 6, 5, 3, 1, 4, 2]
    print('chain:', *reversed(chain))
    print('gold: ', *gold_order)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for a, b in pairwise(chain):
        # Plot lines connecting two atoms in the chain
        ax.plot([a.x, b.x], [a.y, b.y], [a.z, b.z], c='r')
    for a in chain:
        ax.scatter(a.x, a.y, a.z, c='blue') # Plot actual atoms
    plt.show()

if __name__ == '__main__':
    main()