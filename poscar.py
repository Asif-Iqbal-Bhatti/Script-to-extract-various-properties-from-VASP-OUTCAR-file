import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path

ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
Ang32Bohr3 = 6.74833304162   # 1 A^3 = 6.7483330371 [a.u]^3

def poscar():
	if not os.path.exists('POSCAR'):
		print (' ERROR: POSCAR does not exist here.')
		sys.exit(0)
	print('Reading POSCAR/CONTCAR: \n')
	pos = []; kk = []; lattice = []; sum = 0
	file = open('POSCAR','r') or open('CONTCAR','r')
	
	firstline   = file.readline() # IGNORE first line comment
	secondfline = file.readline() # scale
	Latvec1 = file.readline()
	#print ("Lattice vector 1:", (Latvec1), end = '')
	Latvec2 = file.readline()
	#print ("Lattice vector 2:", (Latvec2), end = '')
	Latvec3 = file.readline()
	#print ("Lattice vector 3:", (Latvec3), end = '')
	elementtype=file.readline().split()
	if (str.isdigit(elementtype[0])):
		sys.exit("VASP 4.X POSCAR detected. Please add the atom types")
	print ("Types of elements:", str(elementtype), end = '\n')
	numberofatoms=file.readline()
	Coordtype=file.readline()
	print ("Coordtype:", (Coordtype), end = '\n')	
	
########################---------------------------------------------------------
	print (">>>>>>>>>-------------------# of Atoms--------------------")
	nat = numberofatoms.split()
	nat = [int(i) for i in nat]
	print (nat)
	for i in nat:
		sum = sum + i
	numberofatoms = sum
	print ("Number of atoms:", (numberofatoms), end = '\n')
########################---------------------------------------------------------

	#print (">>>>>>>>>---------------Atomic positions------------------")				
	for x in range(int(numberofatoms)):
		coord = file.readline().split()
		coord = [float(i) for i in coord]
		pos = pos + [coord]
	pos = np.array(pos)
	#print (pos)
		
	file.close()	

########################---------------------------------------------------------
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
		
########################---------------------------------------------------------

	print (">>>>>>>>>---------------Lattice vectors distortions-----------------")				
	lattice = np.array([a] + [b] + [c])
	#determinant = np.linalg.det(lattice)
	lld = local_lattice_distortion(a,b,c)
	print ("lattice distortion parameter g: {}".format(lld) )
	
	print (">>>>>>>>>---------------Space group-----------------")		
	print (" ")
	sp, symm = space_group_analyse(lattice, pos)
	print ( sp, symm )
	print (" ")
	print (">>>>>>>>>--------------------------------------------")
	print ('a=', a)
	print ('b=', b)
	print ('c=', c)
	gamma = math.degrees(math.acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))))
	alpha = math.degrees(math.acos(np.dot(b,c) / (np.linalg.norm(b) * np.linalg.norm(c))))
	beta  = math.degrees(math.acos(np.dot(a,c) / (np.linalg.norm(a) * np.linalg.norm(c))))
	print ("ratio c/a = %2f" %(np.linalg.norm(c) / np.linalg.norm(a) ))
	print ("-"*100)
	print ('||a||=%2f, \u03B1= %2f' %(np.linalg.norm(a), alpha))
	print ('||b||=%2f  \u03B2= %2f' %(np.linalg.norm(b), beta))
	print ('||c||=%2f  \u03B3= %2f' %(np.linalg.norm(c), gamma))
	print ('Vol= %4.8f A^3; %4.8f [a.u]^3' %(volume(a,b,c,math.radians(alpha),math.radians(beta),math.radians(gamma) )))			
###

def local_lattice_distortion(a1,b1,c1):
	#print ("The lattice distortion in paracrystals is measured by the lattice distortion parameter g")
	#print (Back.YELLOW + "Wang, S. Atomic structure modeling of multi-principal-element alloys by the principle")
	#print (Back.YELLOW + "of maximum entropy. Entropy 15, 5536–5548 (2013).")
	#print ("")
	a=np.linalg.norm(a1); b=np.linalg.norm(b1); c=np.linalg.norm(c1)
	d = np.array([a,b,c])
	d_mean = np.mean(d); d_std = np.std(d)
	d_square_mean = (a**2 + b**2 + c**2)/3
	g = np.sqrt( d_square_mean/(d_mean)**2 - 1 )
	return g
###

def space_group_analyse(lattice, pos):
	numbers = [1,2]			
	cell = (lattice, pos, numbers)
	sp=spglib.get_spacegroup(cell, symprec=1e-5)
	symm=spglib.get_symmetry(cell, symprec=1e-5)
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
	return sp, symm
###

def poscar_VASP42VASP5():
	if not os.path.exists('POSCAR' and 'POTCAR'):
		print (' ERROR: POSCAR does not exist here.')
		sys.exit(0)	
	file1 = open("POSCAR",'r')
	line1 = file1.readlines()		
	file1.close()
	
	file2 = open("POTCAR",'r')
	line2 = file2.readlines()		
	file2.close()
	
	atom_number=[]
	for i in line1:
		if ("Direct" or "direct" or "d" or "D") in i:
			PP=line1.index(i)
	atom_number = line1[5].split()
	print(atom_number)
	
	elementtype=[]; count=0
	for i in line2:
		if ("VRHFIN") in i:
			count+=1
			#print (i.split('=')[1].split(':')[0])
			elementtype.append(i.split('=')[1].split(':')[0])
	
	test = open("POSCAR_W",'w')
	
	for i in range(5):
		test.write( line1[i] )
		
	for j in elementtype:
		test.write("\t" +  j)
	test.write("\n" )
	
	for j in atom_number:
		test.write("\t" +  j )
	test.write("\n" )
	
	test.write("Selective dynamics")
	test.write("\n" )
	
	for i in range(len(line1)-PP):
		test.write(line1[PP+i] )
		
	test.close()
	
	print ("                        File is converted: POSCAR_W")
###


#### math.sin function takes argument in radians ONLY
def volume(a,b,c,alpha,beta,gamma):
	ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
	Ang32Bohr3 = 6.74833304162   # 1 A^3 = 6.7483330371 [a.u]^3
	
	length = np.linalg.norm(a) * np.linalg.norm(b) * np.linalg.norm(c) 
	volume = length * ( np.sqrt(1 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma) - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2) )
	vol_au = volume * Ang32Bohr3
	return volume, vol_au
	
