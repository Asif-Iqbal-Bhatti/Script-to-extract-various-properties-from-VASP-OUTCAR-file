#!/usr/bin/env python3

'''
#=====================================================
# AUTHOR:  Asif Iqbal
# GitHUB:  Asif_em
# PURPOSE: Read CHGCAR file and compute the integrated
#          charge density to sum up to # of electrons.
#          Dependency on pymatgen module
#          The following code is uploaded to my GitHub
#          
#=====================================================
'''

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import VolumetricData

poscar, data, data_aug = VolumetricData.parse_file('CHGCAR')
latt = np.array(poscar.structure.lattice.matrix)

# As per VASP the CHGCAR needs to be divided by Volume
# COuld be total or diff if ISPIN = 2
chgd = data["total"] / np.linalg.det(latt) 
pos = poscar.structure.cart_coords
# get volume per grid point
dv = np.linalg.det(latt)/(np.array(chgd).shape[0]*np.array(chgd).shape[1]*np.array(chgd).shape[2])
# get the step size in each direction (x, y, z)
dr = latt/np.array(chgd).shape
total_elect = np.sum(np.array(chgd))*dv
num_grid_points = np.prod(np.array(chgd).shape)
integrated_charge = np.cumsum(np.array(chgd)) * dv
# Measure the charge density sphere
distance = np.linspace(0, np.linalg.norm(dr), num_grid_points)

print(dv, np.linalg.norm(dr))
print(f"Shape of data in CHGCAR file:: {np.array(chgd).shape}")
print(f"Total integrated electron number = {total_elect:8.5f}")
# ploting 
plt.plot(distance, integrated_charge, '-o', color='m', linewidth=0.5, markersize=0.5, markerfacecolor='none')
plt.xlabel("Distance from the freestanding O atom (Å)")
plt.ylabel("Integrated charge")
plt.savefig('dvsQ.png', dpi=400, bbox_inches='tight')
plt.show()

