import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path

ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
ang2bohr   = 6.7483330371   # 1 A^3 = 6.7483330371 [a.u]^3

def poscar():
	while True:
		a = input ("Enter POSCAR filename: \t")
		if (a == 'POSCAR') or (a == 'CONTCAR'):	
			print('Reading POSCAR: \n')
			pos = []; kk = []; lattice = []; sum = 0
			file = open(a,'r')
			
			firstline   = file.readline()
			secondfline = file.readline()
			Latvec1 = file.readline()
			#print ("Lattice vector 1:", (Latvec1), end = '')
			Latvec2 = file.readline()
			#print ("Lattice vector 2:", (Latvec2), end = '')
			Latvec3 = file.readline()
			#print ("Lattice vector 3:", (Latvec3), end = '')
			elementtype=file.readline()
			elementtype = elementtype.split()			
			print ("Types of elements:", str(elementtype), end = '\n')
			numberofatoms=file.readline()
			Coordtype=file.readline()
			print ("Coordtype:", (Coordtype), end = '\n')	
##########################---------------------------------------------------------

			print ("*************-------------------# of Atoms--------------------")
			
			nat = numberofatoms.split()
			nat = [int(i) for i in nat]
			print (nat)
			for i in nat:
				sum = sum + i
			numberofatoms = sum
			print ("Number of atoms:", (numberofatoms), end = '\n')
##########################---------------------------------------------------------

			print ("*************---------------Atomic positions------------------")				
			for x in range(int(numberofatoms)):
				coord = file.readline().split()
				coord = [float(i) for i in coord]
				pos = pos + [coord]
			pos = np.array(pos)
			print (pos)
				
			file.close()	

##########################---------------------------------------------------------
			a=[]; b=[]; c=[];
			Latvec1=Latvec1.split()
			Latvec2=Latvec2.split()
			Latvec3=Latvec3.split()	
			for ai in Latvec1:
				a.append(float(ai))
			for bi in Latvec2:
				b.append(float(bi))
			for ci in Latvec3:
				c.append(float(ci))
				
##########################---------------------------------------------------------

			print ("//////---------------Lattice vectors-----------------")				
			lattice = np.array([a] + [b] + [c])
			print (lattice)
			print (" ")
			print ("//////---------------Space group-----------------")		
			print (" ")
			numbers = [1,1]			
			cell = (lattice, pos, numbers)
			print(spglib.get_spacegroup(cell, symprec=1e-5))
			#print(spglib.get_symmetry(cell, symprec=1e-5))
			#print(spglib.niggli_reduce(lattice, eps=1e-5))
			
			#mesh = [8, 8, 8]
			#mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])
			## All k-points and mapping to ir-grid points
			#for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
			#	print("%3d ->%3d %s" % (i, ir_gp_id, gp.astype(float) / mesh))
			#
			## Irreducible k-points
			#print("Number of ir-kpoints: %d" % len(np.unique(mapping)))
			#print(grid[np.unique(mapping)] / np.array(mesh, dtype=float))
			#
			##
			## With shift
			##
			#mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[1, 1, 1])
			#
			## All k-points and mapping to ir-grid points
			#for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
			#	print("%3d ->%3d %s" % (i, ir_gp_id, (gp + [0.5, 0.5, 0.5]) / mesh))
			#
			## Irreducible k-points
			#print("Number of ir-kpoints: %d" % len(np.unique(mapping)))
			#print((grid[np.unique(mapping)] + [0.5, 0.5, 0.5]) / mesh)

			print (" ")
			print ("/////------------------------------------------------")
			
			print ('a=', a)
			print ('b=', b)
			print ('c=', c)
			gamma = math.degrees(math.acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))))
			alpha = math.degrees(math.acos(np.dot(b,c) / (np.linalg.norm(b) * np.linalg.norm(c))))
			beta  = math.degrees(math.acos(np.dot(a,c) / (np.linalg.norm(a) * np.linalg.norm(c))))
			print ("ratio c/a = %2f" %(np.linalg.norm(c) / np.linalg.norm(a) ))
			print ("#####------------------------------------------------")
			print ('||a||=%2f, \u03B1= %2f' %(np.linalg.norm(a), alpha))
			print ('||b||=%2f  \u03B2= %2f' %(np.linalg.norm(b), beta))
			print ('||c||=%2f  \u03B3= %2f' %(np.linalg.norm(c), gamma))
			print ('Vol= %4.8f A^3; %4.8f [a.u]^3' %(volume(a,b,c,math.radians(alpha),math.radians(beta),math.radians(gamma) )))			
			break
		else:
			print ('NO file entered or wrong filename') 
			break	

#### math.sin function takes argument in radians ONLY
def volume(a,b,c,alpha,beta,gamma):
	length = np.linalg.norm(a) * np.linalg.norm(b) * np.linalg.norm(c) 
	volume = length * ( np.sqrt(1 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma) - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2) )
	vol_au = volume * ang2bohr
	return volume, vol_au
	