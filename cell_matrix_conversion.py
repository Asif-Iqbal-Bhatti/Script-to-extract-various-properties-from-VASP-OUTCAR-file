#!/usr/bin/python3
#####---------------------------------------------------------
#####---------------------------------------------------------
#    Credit: Asif Iqbal BHATTI
#    CODE to convert Cell Matrix to Cell Parameters
#    VERSION: This script runs with python3 or later
#    FORMAT of input file POSCAR VASP5 format
#    DATE: 23/12/2019
#    USAGE: python3 sys.argv[0]
#####---------------------------------------------------------
#####---------------------------------------------------------

import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path
from termcolor import colored

#####---------------------------------------------------------
# Reading just a POSCAR file from the command prompt
#####--------------------------------------------------------- 

def poscar():
	while True:
		if (sys.argv[1] == "-h"):
			print('HELP: execute by typing python3', sys.argv[0])
			break
		elif (sys.argv[1] == 'POSCAR'):
			print('Reading a POSCAR file:', sys.argv[1])
			pos = []; kk = [];
			lattice = [];
			file = open(sys.argv[1],'r')
			
			firstline   = file.readline()
			secondfline = file.readline()
			Latvec1 = file.readline()
			#print ("Lattice vector 1:", (Latvec1), end = '')
			Latvec2 = file.readline()
			#print ("Lattice vector 2:", (Latvec2), end = '')
			Latvec3 = file.readline()
			#print ("Lattice vector 3:", (Latvec3), end = '')
			elementtype=file.readline()
			print ("Types of elements:", str(elementtype), end = '')
			numberofatoms=file.readline()
			print ("Number of atoms:", (numberofatoms), end = '')
			Coordtype=file.readline()
			print ("Coordtype:", (Coordtype), end = '')	
			
			print ("//////---------------Atomic positions-----------------")				
			for x in range(int(numberofatoms)):
				coord = file.readline().split()
				coord = [float(i) for i in coord]
				pos = pos + [coord]
			pos = np.array(pos)
			print (pos)
				
			file.close()	

	#####---------------------------------------------------------			
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
	#####---------------------------------------------------------				
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
			print(spglib.niggli_reduce(lattice, eps=1e-5))
			
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
			print ('Vol= %3.5f' %(volume(a,b,c,math.radians(alpha),math.radians(beta),math.radians(gamma) )))			
			break
		else:
			print ('NO file entered or wrong filename') 
			break			
#####---------------------------------------------------------
# Looping over all directories in a current folder containing POSCARS files
#####--------------------------------------------------------- 

def main_poscar():
	count = 0
	os.system("rm out.dat")
	VOL_P = [];
	mypath = os.getcwd()
	#print (mypath)
	print ("               >>>>> Converting Cell Matrix to Cell Parameters <<<<<<")
	for entry in os.listdir(mypath):
		if os.path.isdir(os.path.join(mypath, entry)):
			#print (entry)
			for file in os.listdir(entry):
				#print (file)
				if file == "POSCAR":
					count+=1
					filepath = os.path.join(entry, file)
					#f = open(filepath, 'r')
					#print (f.read())
					#f.close()	
					fo = open(filepath, 'r')
					ofile=open('out.dat','a+')
					print (colored('>>>>>>>>  Name of the file: ','red'), fo.name, end = '\n', flush=True)
					ofile.write (fo.name + '\n')
					ofile.write ("")
					firstline   = fo.readline()
					secondfline = fo.readline()
					Latvec1 = fo.readline()
					#print ("Lattice vector 1:", (Latvec1), end = '')
					#ofile.write (Latvec1)
					Latvec2 = fo.readline()
					#print ("Lattice vector 2:", (Latvec2), end = '')
					#ofile.write (Latvec2)
					Latvec3 = fo.readline()
					#print ("Lattice vector 3:", (Latvec3), end = '')
					#ofile.write (Latvec3)
					elementtype=fo.readline()
					#print ("Types of elements:", str(elementtype), end = '')
					#ofile.write (str(elementtype))
					numberofatoms=fo.readline()
					#print ("Number of atoms:", (numberofatoms), end = '')
					#ofile.write ((numberofatoms))
					ofile.write ("\n")
					
					fo.close()
		#####---------------------------------------------------------
					a=[]; b=[]; c=[];
					Latvec1=Latvec1.split()
					Latvec2=Latvec2.split()
					Latvec3=Latvec3.split()
		#####---------------------------------------------------------
					for ai in Latvec1:
						a.append(float(ai))
					for bi in Latvec2:
						b.append(float(bi))
					for ci in Latvec3:
						c.append(float(ci))	
					print ("////------------------------------------------------")
					print ('a=', a)
					ofile.write ("'a=' {}\n".format(a))
					print ('b=', b)
					ofile.write ("'b=' {}\n".format(b))
					print ('c=', c)
					ofile.write ("'c=' {}\n".format(c))			
		#####---------------------------------------------------------
		
					alpha, beta, gamma = lattice_angles(a,b,c)
					VOL_POS = np.dot(a, np.cross(b,c))	
					VOL_P.append(VOL_POS)	
					
					print ('\u03B1=', alpha, '\u03B2=', beta, '\u03B3=', gamma)
					ofile.write ("'\u03B1=' {} '\u03B2=' {} '\u03B3=' {}\n".format(alpha,beta,gamma))
					print ("#####------------------------------------------------")
					print ('||a||=', np.linalg.norm(a))
					ofile.write ("'||a||=' {}\n".format(np.linalg.norm(a)))
					print ('||b||=', np.linalg.norm(b))
					ofile.write ("'||b||=' {}\n".format(np.linalg.norm(b)))			
					print ('||c||=', np.linalg.norm(c)) 
					ofile.write ("'||c||=' {}\n".format(np.linalg.norm(c)))
					print ('Vol= %3.5f' %(VOL_POS))						
					ofile.write ("***************************************************\n")
					ofile.close()
	print ("Number of folders detected: ", count)
	return VOL_P
	
def main_contcar():
	count = 0
	os.system("rm out_contcar.dat")
	VOL_C = [];	
	mypath = os.getcwd()
	print ("               >>>>> Converting Cell Matrix to Cell Parameters <<<<<<")
	for entry in os.listdir(mypath):
		if os.path.isdir(os.path.join(mypath, entry)):	
			for file in os.listdir(entry):
				if file == "CONTCAR":
					count+=1
					filepath = os.path.join(entry, file)
					fo = open(filepath, 'r')
					
					ofile=open('out_contcar.dat','a+')
					
					print (colored('>>>>>>>>  Name of the file: ','yellow'), fo.name, end = '\n', flush=True)
					ofile.write (fo.name + '\n')
					ofile.write ("")
					firstline   = fo.readline()
					secondfline = fo.readline()
					Latvec1 = fo.readline()
					Latvec2 = fo.readline()
					Latvec3 = fo.readline()
					elementtype=fo.readline()
					numberofatoms=fo.readline()
					
					ofile.write ("\n")
					
					fo.close()
		#####---------------------------------------------------------
					a=[]; b=[]; c=[];
					Latvec1=Latvec1.split()
					Latvec2=Latvec2.split()
					Latvec3=Latvec3.split()
		#####---------------------------------------------------------
					for ai in Latvec1:
						a.append(float(ai))
					for bi in Latvec2:
						b.append(float(bi))
					for ci in Latvec3:
						c.append(float(ci))	
					print ("#####------------------------------------------------")
					print ('a=', a)
					ofile.write ("'a=' {}\n".format(a))
					print ('b=', b)
					ofile.write ("'b=' {}\n".format(b))
					print ('c=', c)
					ofile.write ("'c=' {}\n".format(c))			
		#####---------------------------------------------------------

					alpha, beta, gamma = lattice_angles(a,b,c)
					VOL_CON = np.dot(a, np.cross(b,c))
					VOL_C.append(VOL_CON)
						
					print ('\u03B1=', alpha, '\u03B2=', beta, '\u03B3=', gamma)
					ofile.write ("'\u03B1=' {} '\u03B2=' {} '\u03B3=' {}\n".format(alpha,beta,gamma))
					print ("#####------------------------------------------------")
					print ('||a||=', np.linalg.norm(a))
					ofile.write ("'||a||=' {}\n".format(np.linalg.norm(a)))
					print ('||b||=', np.linalg.norm(b))
					ofile.write ("'||b||=' {}\n".format(np.linalg.norm(b)))			
					print ('||c||=', np.linalg.norm(c)) 
					ofile.write ("'||c||=' {}\n".format(np.linalg.norm(c)))
					print ('Vol= %4.6f ' %(VOL_CON)		)				
					ofile.write ("***************************************************\n")
					ofile.close()
			#print (VOL_C)		
	print ("Number of folders detected: ", count)
	return VOL_C
	
	
	
#### math.sin fucntion takes arguments in radians ONLY
def volume(a,b,c,alpha,beta,gamma):
	length = np.linalg.norm(a) * np.linalg.norm(b) * np.linalg.norm(c) 
	volume = length * ( np.sqrt(1 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma) - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2) )
	return volume

# Ordering of angles does matter
def lattice_angles(a,b,c):
	### gamma = Cos-1( (a.b)/||a||.||b|| )
	### alpha = Cos-1( (b.c)/||b||.||c|| )
	### beta  = Cos-1( (a.c)/||a||.||c|| )
	gamma = math.degrees(math.acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))))
	alpha = math.degrees(math.acos(np.dot(b,c) / (np.linalg.norm(b) * np.linalg.norm(c))))
	beta  = math.degrees(math.acos(np.dot(a,c) / (np.linalg.norm(a) * np.linalg.norm(c))))
	return alpha, beta, gamma

def volume_diff(VOL_P, VOL_C):
	n=os.popen("find . -mindepth 1 -maxdepth 1 -type d | wc -l").read()
	print ("VOLUME Diff %18s %12s %15.15s" %("CONTCAR",  "POSCAR",  "contcar-poscar"))
	for i in range(int(n)):
		print ("The difference is: %12.6f %12.6f %15.8f " %(VOL_C[i], VOL_P[i], VOL_C[i] - VOL_P[i]) )
	
	
if __name__ == "__main__":

	#poscar()

	VOL_P = main_poscar()
	VOL_C = main_contcar()
	print (" ----------------------------------------------------       ")
	print (" ----------------------------------------------------       ")
	print (" ----------------------------------------------------       ")
	volume_diff(VOL_P, VOL_C)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
