#!/usr/bin/env python

# this script needs 'plot.3d.vector' or 'plot.3d.scalar' file from **visual calculations
# and creates 'cube_grid' file with 3d grid for visualization in DIRAC (with '.3D_IMP' keyword)

import sys

with open(sys.argv[1], 'r') as f:
    with open('cube_grid', 'w') as g:
        lines = f.readlines()
        for line in lines:
            l = line.strip().split()
            if (len(l) > 3):
                g.write(" {0:>12} {1:>12} {2:>12}\n".format(l[0], l[1], l[2]))


