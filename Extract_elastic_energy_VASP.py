#!/usr/bin/env python3

'''
#####---------------------------------------------------------
#####---------------------------------------------------------
#    Credit	: 	Asif Iqbal BHATTI
#    CODE to: 	OBTAIN Elastic properties from OUTCAR file,
#               compare POSCAR and CONTCAR volume deformation
#               upon minimization, and extract energy from a number
#				of directories. extract lattice distortion  
#    VERSION: 	This script runs with python3 or later
#    FORMAT	:	POSCAR VASP5 format
#    DATE	: 	28/12/2019
#    USAGE	: 	python3 sys.argv[0]
#	 VERSION:   4.0	
#####---------------------------------------------------------
#####---------------------------------------------------------
'''

import os, sys, spglib, scipy
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path
from termcolor import colored
import multiprocessing as mp
from colorama import Fore, Back, Style, init
from   pylab import *
init(autoreset=True)


'''
##########################---------------------------------------------------------
#		1-	Read only POSCAR file in a given directory & obtain local lattice distortion
##########################---------------------------------------------------------
'''

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
	lld = lattic_distortion.local_lattice_distortion(a,b,c)
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
class lattic_distortion():
	@classmethod
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
	
	def local_lattice_distortion_DEF1():
		#print ("The lattice distortion in paracrystals is measured by the lattice distortion parameter g")
		#print (Back.YELLOW + "Wang, S. Atomic structure modeling of multi-principal-element alloys by the principle")
		#print (Back.YELLOW + "of maximum entropy. Entropy 15, 5536–5548 (2013).")
		print ("+"*40,"HUME ROTHERY RULE","+"*40)
		C_i=C=0.2 ; r_avg = 0.0; del_sum=0.0
		elements = ["Nb", "Hf", "Ta", "Ti", "Zr"]
		eta = {
		"Nb" : 1.98,
		"Hf" : 2.08,
		"Ta" : 2.00,
		"Ti" : 1.76,
		"Zr" : 2.06, }
		
		print ("                      {element: atomic radius}")
		print (eta)
		
		for i in elements: 
			r_avg = r_avg + C * eta[i] 
		
		for j in elements:
			del_sum = del_sum + C * ( 1 - float(eta[j]) / r_avg )**2
		del_sum = 100 * np.sqrt(del_sum) 	
		print("HEA_atomic_size_mismatch: \u03B4={}".format(del_sum))
	###
		
def local_lattice_distortion_DEF2():
	print (">"*10,"local_lattice_distortion_DEF2")
	print ("	Song, H. et al. Local lattice distortion in high-entropy alloys.")
	print ("	Phys. Rev. Mater. 1, 23404 (2017).")
	print ("	(***) Different definition of the atomic radius for the description ")
	print ("	of the local lattice distortion in HEAs")
	
	if not os.path.exists('POSCAR' and 'CONTCAR'):
		print ('>>> ERROR: POSCAR & CONTCAR does not exist (Both should be in the same directory)')
		sys.exit(0)
	print('Reading POSCAR and CONTCAR ... \n')
	
	x = []; y =[]; z=[]
	xp =[]; yp = []; zp = []; temp=0
	
	f = open('POSCAR','r')
	lines_poscar = f.readlines()
	f.close()
	
	f = open('CONTCAR','r')
	lines_contcar = f.readlines()
	f.close()
	
	file_P = ase.io.read('POSCAR')
	pos = file_P.get_cell_lengths_and_angles()
	print (CRED + "POSCAR=>Length&Angles->{}".format(pos) + CEND)
	file_C = ase.io.read('CONTCAR')
	con = file_C.get_cell_lengths_and_angles() 
	print (CRED + "CONTCAR=>Length&Angles->{}".format(con) + CEND)
	print ("Cell vectors difference:: ",con-pos)
	
	sum_atoms = lines_poscar[6].split()  ### reading 7th lines for reading # of atoms
	sum_atoms = [int(i) for i in sum_atoms]
	sum_atoms = sum(sum_atoms)
	
	for i in lines_poscar:
		if "Direct" in i:
			lp=lines_poscar.index(i)
	for j in lines_contcar:
		if "Direct" in j:
			lc=lines_contcar.index(j)
			
	for i in range(sum_atoms):
		x, y, z    = lines_poscar[lp+1+i].split()
		xc, yc, zc = lines_contcar[lp+1+i].split()
		x = float(x); y = float(y); z = float(z)
		xc = float(xc); yc = float(yc); zc = float(zc)
		temp = temp + np.sqrt( (x-xc)**2 + (y-yc)**2 + (z-zc)**2 )
	temp = temp/sum_atoms
	print("local lattice distortion (LLD): \u0394d={}".format(temp))	

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

'''
#####---------------------------------------------------------
#        2- Looping over all directories containing POSCAR & CONTCAR files
#####--------------------------------------------------------- 
'''

def main_poscar():
	count = 0
	os.system("rm out_POSCARS.dat")
	VOL_P = []; pos = []; kk = []; lattice = [];
	mypath = os.getcwd()
	print ("-"*100)	
	print (Back.YELLOW + "{:15s} {:15s} {:15.6s} {:15.6s} {:15.6s} {:15.15s}".format("Directory", "# of atoms", "||a||", \
	"||b||", "||c||", "VOL_POS[A^3]"),end="\n" )	
	print ("-"*100)	

	for entry in os.listdir(mypath):
		if os.path.isdir(os.path.join(mypath, entry)):
			#print (entry)
			for file in os.listdir(entry):
				if file == "POSCAR":
					count+=1; sum = 0
					filepath = os.path.join(entry, file)
					#f = open(filepath, 'r')
					#print (f.read())
					#f.close()	
					fo = open(filepath, 'r')
					ofile=open('out_POSCARS.dat','a')
					
					#print (colored('>>>>>>>>  Name of the file: ','red'), fo.name, end = '\n', flush=True)
					
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
					elementtype = elementtype.split()						
					#print ("Types of elements:", str(elementtype), end = '')
					#ofile.write (str(elementtype))
					numberofatoms=fo.readline()
					#print ("Number of atoms:", (numberofatoms), end = '')
					#ofile.write ((numberofatoms))
					Coordtype=fo.readline()

##########################---------------------------------------------------------
					#print ("**********-------------------# of Atoms--------------------")
					
					nat = numberofatoms.split()
					nat = [int(i) for i in nat]
					for i in nat:
						sum = sum + i
					numberofatoms = sum
					#print ("{} :: Number of atoms: {}".format(nat, numberofatoms) )
##########################---------------------------------------------------------					
					#print ("-----------------------Atomic positions-----------------")
					#print ("Coordtype:", (Coordtype), end = '')						
					for x in range(int(numberofatoms)):
						coord = fo.readline().split()
						coord = [float(i) for i in coord]
						pos = pos + [coord]
					pos = np.array(pos)
					#print (pos)
					
					ofile.write ("\n")			
					fo.close()
##########################---------------------------------------------------------

					a=[]; b=[]; c=[];
					Latvec1=Latvec1.split()
					Latvec2=Latvec2.split()
					Latvec3=Latvec3.split()
					
##########################---------------------------------------------------------
					for ai in Latvec1: a.append(float(ai))
					for bi in Latvec2: b.append(float(bi))
					for ci in Latvec3: c.append(float(ci))	
					#print ('a=', a)
					ofile.write ("'a=' {}\n".format(a))
					#print ('b=', b)
					ofile.write ("'b=' {}\n".format(b))
					#print ('c=', c)
					ofile.write ("'c=' {}\n".format(c))		
					lld = lattic_distortion.local_lattice_distortion(a,b,c)
##########################---------------------------------------------------------
		
					alpha, beta, gamma = lattice_angles(a,b,c)
					VOL_POS = np.dot(a, np.cross(b,c))	
					VOL_P.append(VOL_POS)	

					ofile.write ("'\u03B1=' {} '\u03B2=' {} '\u03B3=' {}\n".format(alpha,beta,gamma))
					ofile.write ("'||a||=' {}\n".format(np.linalg.norm(a)))
					ofile.write ("'||b||=' {}\n".format(np.linalg.norm(b)))		
					ofile.write ("'||c||=' {}\n".format(np.linalg.norm(c)))					
					#print ("#####------------------------------------------------")
					
					#print ("a={} \t ||a||={:10.6f}".format(a, np.linalg.norm(a)) )
					#print ("b={} \t ||b||={:10.6f}".format(b, np.linalg.norm(b)) )
					#print ("c={} \t ||c||={:10.6f}".format(c, np.linalg.norm(c)) )
					print ("{:15s} {:6d} {:15.6f} {:15.6f} {:15.6f} {:15.6f}".format(fo.name, numberofatoms, np.linalg.norm(a), \
					np.linalg.norm(b), np.linalg.norm(c), VOL_POS) )
					
					print ("'\u03B1=' {:6.6f} '\u03B2=' {:6.6f} '\u03B3=' {:6.6f} g={:6.6f}".format(alpha,beta,gamma,lld))
					print ("."*5)
					#print ('Vol= {:6.6f} A^3'.format(VOL_POS))						
					ofile.write ("***************************************************\n")
					ofile.close()
	print ("_"*30)				
	print ("Number of folders detected: ", count)
	return VOL_P

####
def main_contcar():
	count = 0
	os.system("rm out_CONTCARS.dat")
	VOL_C = [];	pos = []; kk = []; lattice = []; sum = 0
	mypath = os.getcwd()	
	print ("-"*100)	
	print (Back.GREEN + "{:15s} {:15s} {:15.6s} {:15.6s} {:15.6s} {:15.15s}".format("Directory", "# of atoms", "||a||", \
	"||b||", "||c||", "VOL_CON[A^3]"),end="\n" )	
	print ("-"*100)	
	print(Style.RESET_ALL)	
	for entry in os.listdir(mypath):
		if os.path.isdir(os.path.join(mypath, entry)):	
			for file in os.listdir(entry):
				
				if file == "CONTCAR":
					count+=1; sum = 0
					filepath = os.path.join(entry, file)
					fo = open(filepath, 'r')
					
					ofile=open('out_CONTCARS.dat','a')
					
					#print (colored('>>>>>>>>  Name of the file: ','yellow'), fo.name, end = '\n', flush=True)
					ofile.write (fo.name + '\n')
					ofile.write ("")
					firstline   = fo.readline()
					secondfline = fo.readline()
					Latvec1 = fo.readline()
					Latvec2 = fo.readline()
					Latvec3 = fo.readline()
					elementtype=fo.readline()
					elementtype = elementtype.split()					
					numberofatoms=fo.readline()
					Coordtype=fo.readline()
					#print ("Coordtype:", (Coordtype), end = '')
##########################---------------------------------------------------------
					#print ("**********-------------------# of Atoms--------------------")
					
					nat = numberofatoms.split()
					nat = [int(i) for i in nat]
					for i in nat:
						sum = sum + i
					numberofatoms = sum
					#print ("{} :: Number of atoms: {}".format(nat, numberofatoms) )
##########################---------------------------------------------------------						
					#print ("//////---------------Atomic positions-----------------")				
					for x in range(int(numberofatoms)):
						coord = fo.readline().split()
						coord = [float(i) for i in coord]
						pos = pos + [coord]
					pos = np.array(pos)
					#print (pos)
					
					ofile.write ("\n")				
					fo.close()
##########################---------------------------------------------------------
					a=[]; b=[]; c=[];
					Latvec1=Latvec1.split()
					Latvec2=Latvec2.split()
					Latvec3=Latvec3.split()					
##########################---------------------------------------------------------
					for ai in Latvec1: a.append(float(ai))
					for bi in Latvec2: b.append(float(bi))
					for ci in Latvec3: c.append(float(ci))	
					#print ('a=', a)
					ofile.write ("'a=' {}\n".format(a))
					#print ('b=', b)
					ofile.write ("'b=' {}\n".format(b))
					#print ('c=', c)
					ofile.write ("'c=' {}\n".format(c))		
					lld = lattic_distortion.local_lattice_distortion(a,b,c)					
##########################---------------------------------------------------------

					alpha, beta, gamma = lattice_angles(a,b,c)
					VOL_CON = np.dot(a, np.cross(b,c))
					VOL_C.append(VOL_CON)
						
					ofile.write ("'\u03B1=' {} '\u03B2=' {} '\u03B3=' {}\n".format(alpha,beta,gamma))
					ofile.write ("'||a||=' {}\n".format(np.linalg.norm(a)))
					ofile.write ("'||b||=' {}\n".format(np.linalg.norm(b)))		
					ofile.write ("'||c||=' {}\n".format(np.linalg.norm(c)))					
					#print ("-"*80)
					
					#print ("a={} \t ||a||={:10.6f}".format(a, np.linalg.norm(a)) )
					#print ("b={} \t ||b||={:10.6f}".format(b, np.linalg.norm(b)) )
					#print ("c={} \t ||c||={:10.6f}".format(c, np.linalg.norm(c)) )
					print ("{:15s} {:6d} {:15.6f} {:15.6f} {:15.6f} {:15.6f}".format(fo.name, numberofatoms, np.linalg.norm(a), \
					np.linalg.norm(b), np.linalg.norm(c), VOL_CON) )

					print ("'\u03B1=' {:6.6f} '\u03B2=' {:6.6f} '\u03B3=' {:6.6f} g={:6.6f}".format(alpha,beta,gamma,lld))
					print ("."*5)
					#print ('Vol= {:6.6f} A^3'.format(VOL_CON), end="\n")						
					ofile.write ("***************************************************\n")
					ofile.close()
			#print (VOL_C)	
	print ("-"*80)			
	print ("Number of folders detected: ", count)
	return VOL_C
	
#### math.sin function takes argument in radians ONLY
def volume(a,b,c,alpha,beta,gamma):
	ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
	Ang32Bohr3 = 6.74833304162   # 1 A^3 = 6.7483330371 [a.u]^3
	
	length = np.linalg.norm(a) * np.linalg.norm(b) * np.linalg.norm(c) 
	volume = length * ( np.sqrt(1 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma) - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2) )
	vol_au = volume * Ang32Bohr3
	return volume, vol_au

#### Ordering of returning angles variables does matter
def lattice_angles(a,b,c):
	### gamma = Cos-1( (a.b)/||a||.||b|| )
	### alpha = Cos-1( (b.c)/||b||.||c|| )
	### beta  = Cos-1( (a.c)/||a||.||c|| )
	gamma = math.degrees(math.acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))))
	alpha = math.degrees(math.acos(np.dot(b,c) / (np.linalg.norm(b) * np.linalg.norm(c))))
	beta  = math.degrees(math.acos(np.dot(a,c) / (np.linalg.norm(a) * np.linalg.norm(c))))
	return alpha, beta, gamma

####### BASH way of finding the # of directories in a working directory ###
def volume_diff(VOL_P, VOL_C):
	n=os.popen("find . -mindepth 1 -maxdepth 1 -type d | wc -l").read()
	print ("-"*80)	
	print ("VOL Diff A^3 %18s %12s %15.15s" %("CONTCAR",  "POSCAR",  "contcar-poscar"))
	print ("-"*80)		
	for i in range(int(n)):
		print ("The difference is: %12.6f %12.6f %15.8f " %(VOL_C[i], VOL_P[i], VOL_C[i] - VOL_P[i]) )

####

'''
#####---------------------------------------------------------
#            3- ELASTIC PROPERTIES from VASP OUTCAR file "STRESS APPROACH"
#####--------------------------------------------------------- 
'''

class Elastic_Matrix:
	def print_Cij_Matrix(): ###EXERCISE
		Bij = []; C = "C"
		for i in range(0, 6, 1): 
			Bij.append([])
			for j in range(1,7, 1): 
				Bij[i].append((C + str(i+1) + str(j)))
		l = np.matrix(Bij)
		return l

	def langragian_strain():
		kk = ('x', 'y', 'z'); C = "E"; Eij=[]
		eta = {
		"xx" : 11,
		"yy" : 22,
		"zz" : 33,
		"yz" : 23,
		"xz" : 13,
		"xy" : 12  }
		print ( "The Langragian strain in Voigt notation: ")
		print ( eta.items())
		
		for n in kk: 
			for m in kk: 
				Eij.append( (C + str(n) + str(m)) )
		ls = np.matrix(Eij).reshape(3,3)
		return ls

###
def elastic_matrix_VASP_STRESS():
	while True:		
		if not os.path.exists('OUTCAR'):
			print (' ERROR: OUTCAR does not exist here.')
			sys.exit(0)
			
		s=np.zeros((6,6))
		c=np.zeros((6,6))	
		file = open("OUTCAR",'r')
		lines = file.readlines()		
		file.close()
		
		for i in lines:
			if "TOTAL ELASTIC MODULI (kBar)" in i:
				ll=lines.index(i)
			if "LATTYP" in i:
				crystaltype=str(i.split()[3])
		print ("DETECTED CRYSTAL FROM OUTCAR:", crystaltype)			
		print (" ")			
		for i in range(0,6):
			l=lines[ll+3+i] # indexing a line in huge file
			word = l.split()
			s[i][:] = word[1:7]
			for j in range(0,6):
				c[i][j] = float(s[i][j])/10.0
				
##########-------------------stress tensor------------------	
			
		Cij = np.matrix(c)
		np.set_printoptions(precision=4, suppress=True)
		print (Cij)
		print(Elastic_Matrix.print_Cij_Matrix()  )
		print(Elastic_Matrix.langragian_strain() )
		print ("\nEigen Values of the matrix Cij:")
		evals = np.linalg.eigvals(Cij)
		if evals.all() > 0:
			print(evals)
			print("All Eigen values are positive")
		
##########-------------------Compliance tensor------------------
##########-------------------  s_{ij} = C_{ij}^{-1}

		Sij = np.linalg.inv(Cij)
		
#---------------------------- ELASTIC PROPERTIES -----------------------------------

		stability_test(Cij, crystaltype)
		
#------------------------------- Voigt bulk modulus  K_v  $(GPa)$---------------
#------------------ 9K_v = (C_{11}+C_{22}+C_{33}) + 2(C_{12} + C_{23} + C_{31}) 
		Kv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[1,2] + Cij[2,0]))/9.0
#------------------------------- Reuss shear modulus  G_v  $(GPa)$------------------
#------------------ 15/G_R = 4(s_{11}+s_{22}+s_{33}) - 4(s_{12} + s_{23} + s_{31}) + 3(s_{44} + s_{55} + s_{66})$
		Gv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) - (Cij[0,1] + Cij[1,2] + Cij[2,0]) + 3 * (Cij[3,3] + Cij[4,4] + Cij[5,5]))/15.0
#------------------------------- Reuss bulk modulus  K_r  $(GPa)$----------------
#------------------  1/K_R = (s_{11}+s_{22}+s_{33}) + 2(s_{12} + s_{23} + s_{31})$
		Kr = 1/((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) )
#------------------------------- Reuss shear modulus  G_r  $(GPa)$------------------
		Gr = 15/(4 * (Sij[0,0] + Sij[1,1] + Sij[2,2]) - 4 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) + 3 * (Sij[3,3] + Sij[4,4] + Sij[5,5]))
		
#-----------------------------------------------------------------------------------

		## Young's Modulus "E": Voigt
		Ev = (9*Kv*Gv)/(3*Kv + Gv)				
		## Young's: Reuss
		Er = (9*Kr*Gr)/(3*Kr + Gr)		
		
		## Poisson's ratio: Voigt
		Nu_V = (3*Kv - Ev)/(6*Kv)				
	       ## Poisson's ratio: Reuss
		Nu_R = (3*Kr - Er)/(6*Kr)
		
		## P-wave modulus, M: Voigt
		MV = Kv + (4*Gv/3.0)	        
	       ## P-wave modulus, M: Reuss
		MR = Kr + (4*Gr/3.0)	
			
#-----------------------------------------------------------------------------------
#-------------- Voigt-Reuss-Hill Approximation: average of both methods
		Kvrh = (Kv + Kr)/2.0
		Gvrh = (Gv + Gr)/2.0
		Mvrh = (MV + MR)/2.0
		Evrh = (Ev + Er)/2.0
		Nu_vrh = (Nu_V + Nu_R)/2.0
		KG_ratio_V = Kv/Gv
		KG_ratio_R = Kr/Gr
		KG_ratio_vrh = Kvrh/Gvrh			
#-------------- Isotropic Poisson ratio $\mu 
#-------------- $\mu = (3K_{vrh} - 2G_{vrh})/(6K_{vrh} + 2G_{vrh})$
		mu = (3 * Kvrh - 2 * Gvrh) / (6 * Kvrh + 2 * Gvrh )
		
#----------------------------------------------------------------------
#-----------------------------------------------------------------------------------
		
		print ("\n \n                         Voigt     Reuss    Average")
		print ("-------------------------------------------------------")
		print ("Bulk modulus   (GPa)  %9.3f %9.3f %9.3f " % (Kv, Kr, Kvrh))
		print ("Shear modulus  (GPa)  %9.3f %9.3f %9.3f " % (Gv, Gr, Gvrh))
		print ("Young modulus  (GPa)  %9.3f %9.3f %9.3f " % (Ev, Er, Evrh))
		print ("Poisson ratio         %9.3f %9.3f %9.3f " % (Nu_V, Nu_R, Nu_vrh))
		print ("P-wave modulus  (GPa) %9.3f %9.3f %9.3f " % (MV, MR, Mvrh))
		print ("Bulk/Shear ratio      %9.3f %9.3f %9.3f (%s) " %(KG_ratio_V, KG_ratio_R, KG_ratio_vrh,  ductile_test(KG_ratio_vrh) ))
		print ("-------------------------------------------------------")
		print("Isotropic Poisson ratio: ", mu)			
		break
		
###		
def ductile_test(ratio):
	if(ratio > 1.75):
		return "ductile"
	else:
		return "brittle"
###
	
def stability_test(matrix, crystaltype):
	c = np.copy(matrix)

	if(crystaltype =="cubic"):
		print ("Cubic crystal system \n")
		print ("Born stability criteria for the stability of cubic system are : \ [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
		print ("(i) C11 - C12 > 0;    (ii) C11 + 2C12 > 0;   (iii) C44 > 0 \n ")
		
		## check (i)   keep in mind list starts with 0, so c11 is stored as c00
	if(c[0][0] - c[0][1] > 0.0):
		print ("Condition (i) satisfied.")
	else:
		print ("Condition (i) NOT satisfied.")
	
	if(c[0][0] + 2*c[0][1] > 0.0):
		print ("Condition (ii) satified.")
	else:
		print ("Condition (ii) NOT satisfied.")
	
	if(c[3][3] > 0.0):
		print ("Condition (iii) satified.")
	else:
		print ("Condition (iii) NOT satisfied.")

	if(crystaltype =="hexagonal"):
		print ("Hexagonal crystal system \n")
		print ("Born stability criteria for the stability of hexagonal system are \: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
		print ("(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0 \n ")
		
		## check (i)   keep in mind list starts with 0, so c11 is stored as c00
		if(c[0][0] - c[0][1] > 0.0):
			print ("Condition (i) satisfied.")
		else:
			print ("Condition (i) NOT satisfied.")
		
		if(2*(c[0][2]*c[0][2]) < c[2][2]*(c[0][0] + c[0][1])):
			print ("Condition (ii) satified.")
		else:
			print ("Condition (ii) NOT satisfied.")
		
		if(c[3][3] > 0.0):
			print ("Condition (iii) satified.")
		else:
			print ("Condition (iii) NOT satisfied.")
###

def born_stability_criterion():
	print ("Born stability criteria for the stability of following systems \n")
	print ("Cubic crystal system.... \n")	
	print ("(i) C11 - C12 > 0;    (ii) C11 + 2C12 > 0;   (iii) C44 > 0 \n ")
	print ("Hexagonal crystal system.... \n")	
	print ("(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0 \n ")
	print ("Tetragonal crystal system.... \n")	
	print ("(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0;   (iv) C66 > 0;    (v) 2C16^2 < C66*(C11-C12) \n ")
	print ("rhombohedral crystal system.... \n")	
	print ("(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0; \n ")
	print ("orthorhombic crystal system.... \n")	
	print ("(i) C11 > 0;   (ii) C11*C22 > C12^2;   (iii) C11*C22*C33 + 2C12*C13*C23 - C11*C23^2 - C22*C13^2 - C33*C12^2 > 0;   (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0 \n")
	print ("Monoclinic crystal system.... \n")
	print ("[Ref- Mouhat and Coudert, PRB 90, 224104 (2014), and Wu et al. PRB 76, 054115 (2007)]  \n")
	print ("(i) C11 > 0;  (ii)  C22 > 0; (iii)  C33 > 0; (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0  ")
	print ("(vii) [C11 + C22 + C33 + 2*(C12 + C13 + C23)] > 0;    (viii)  C33*C55 - C35^2 > 0;   (ix)  C44*C66 - C46^2 > 0;   (x) C22 + C33 - 2*C23  > 0 ")
	print ("(xi) C22*(C33*C55 - C35^2) + 2*C23*C25*C35 - (C23^2)*C55 - (C25^2)*C33   > 0  ")
	print ("(xii)  2*[C15*C25*(C33*C12 - C13*C23) + C15*C35*(C22*C13 - C12*C23) + C25*C35*(C11*C23 - C12*C13)] - [C15*C15*(C22*C33 - C23^2) + C25*C25*(C11*C33 - C13^2) + C35*C35*(C11*C22 - C12^2)] + C55*g > 0  ")
	print (" where, g = [C11*C22*C33 - C11*C23*C23 - C22*C13*C13 - C33*C12*C12 + 2*C12*C13*C23 ] "	)	
###

'''
#####---------------------------------------------------------
#            4- CREATE ENERGY vs VOLUME file from VASP OUTCAR file "STRAIN APPROACH"
#####--------------------------------------------------------- 
'''

def create_energy_vs_volume():
	import fnmatch
	mypath = os.getcwd()
	os.system("rm energy-vs-volume energy-vs-strain")
	eV2Hartree=0.036749309
	Ang32Bohr3=6.74833304162
	
	E=[]; dir_list=[]; count = 0; dir_E=[];
	vol_cell=[]; strain_file=[]; strain_value=[] # strain_value is deformation
	
	print ("               >>>>> Extracting Energies from directories  <<<<<<")
	for entry in os.listdir(mypath):
		if not os.path.exists('strain-01'):
			print (' ERROR: strain-* files does not exist. create a strain file for each deformation values.')
			sys.exit(0)	
		if fnmatch.fnmatchcase(entry,'strain-*'):		
			f = open(entry,'r')
			lines = f.readline()  #Read  first line only
			strain_value.append( float(lines) )
			f.close()
			if os.path.isfile(os.path.join(mypath, entry)):
				strain_file.append(entry)

		if os.path.isdir(os.path.join(mypath, entry)):
			dir_list.append(entry)
			
			for file in os.listdir(entry):
				if file == "OUTCAR":
					filepath = os.path.join(entry, file)
					if not os.path.exists(filepath):
						print (' ERROR: OUTCAR does not exist here.')
						sys.exit(0)
					f = open(filepath,'r')
					lines = f.readlines()
					f.close()
					
					for i in lines:
						if "  free  energy   TOTEN  =" in i:
							m=float(i.split()[4])
						if "  volume of cell :"	in i:
							v=float(i.split()[4])
					vol_cell.append(v)		
					E.append(m)
					count+=1	
	print("# of folders detected: ", count)
	print ("Directory :%10.6s %14s %18s %25.20s " % ("Folder", "Energy(eV)", "Vol_of_cell(A^3)", "strain_deformation" ))
		
	for i in range(math.floor(count/2)): # 0 to 4
		print ("Folder name: %10.10s %16.8f %16.8f %16.12s %14.4f" %(dir_list[i], E[i], vol_cell[i], strain_file[i], strain_value[i] ))
	if (bool(math.floor(count/2))):
		i = math.floor(count/2)
		print(Back.GREEN + 'Folder name: %10.10s %16.8f %16.8f %16.12s %14.4f <--Ref' % (dir_list[i], E[i], vol_cell[i], strain_file[i], strain_value[i] ))		
		print(Style.RESET_ALL, end="")
	for i in range(math.ceil(count/2), count, 1):
		print ("Folder name: %10.10s %16.8f %16.8f %16.12s %14.4f" %(dir_list[i], E[i], vol_cell[i], strain_file[i], strain_value[i] ))		
	#rc = subprocess.Popen(['bash', 'extract_energy.sh'])
		
	print (colored('ENERGIES & VOLUMES ARE WRITTEN IN ATOMIC UNITS TO A FILE <energy-vs-volume>','yellow'), end = '\n', flush=True)
	print (colored('IT WILL BE READ BY ELASTIC SCRIPTS FOR POSTPROCESSING eV-->Ha; A^3-->Bohr^3','yellow'), end = '\n', flush=True)	
	
	file = open("energy-vs-volume",'w')
	for i in range(count):
		file.write ("%14.6f %14.6f\n" %(vol_cell[i] * Ang32Bohr3, E[i] * eV2Hartree))	
	file.close()

	file = open("energy-vs-strain",'w')
	for i in range(count):
		file.write ("%12.6f %14.6f\n" %(strain_value[i], E[i] * eV2Hartree))	
	file.close()
###

'''
#####---------------------------------------------------------
#            5- FITTING ENERGY vs VOLUME CURVE from VASP OUTCAR file
#####--------------------------------------------------------- 
'''

###
def fitting_energy_vs_volume_curve_ELASTIC():
	from   sys   import stdin
	import matplotlib.pyplot as plt
	import matplotlib.ticker as ptk 
	import pylab             as pyl
	import matplotlib.style
	
	bohr_radius     = 0.529177
	bohr32ang3		= 0.14818474347
	joule2hartree	= 4.3597482
	joule2rydberg   = joule2hartree/2.
	unitconv        = joule2hartree/bohr_radius**3*10.**3
	
	if (str(os.path.exists('energy-vs-volume'))=='False'): 
		sys.exit("ERROR: file energy-vs-volume not found!\n")
	energy = []; volume = []
	read_energy = open('energy-vs-volume',"r")
	
	while True:
		line = read_energy.readline()
		line = line.strip()
		if len(line) == 0: break
		energy.append(float(line.split()[1]))
		volume.append(float(line.split()[0]))
	volume,energy=sortvolume(volume,energy)

	print ("===============================")
	print ("Lattice symmetry codes"         )
	print ("-------------------------------")
	print ("1 --> Simple cubic (sc)"        )
	print ("2 --> Body-centered cubic (bcc)")
	print ("3 --> Face-centered cubic (fcc)")
	print ("-------------------------------")
	print ("0 --> Others")
	print ("===============================\n")
	scheck = input("Enter lattice symmetry code [default 0] >>>> ").replace(" ", "") 
		
	'''
	These factors are for the conversion from the conventional cell to primitive cell
	BCC: (a^3)/2 primitive cell volume
	FCC: (a^3)/4 primitive cell volume
	
	'''
	
	isym   = 0; factor = 1
	if ( scheck == "1" ): isym = 1 ; factor=1 ; slabel = "(sc) "
	if ( scheck == "2" ): isym = 2 ; factor=2 ; slabel = "(bcc)"
	if ( scheck == "3" ): isym = 3 ; factor=4 ; slabel = "(fcc)"
	print ("Verification lattice symmetry code      >>>> %d " %(isym) )
#-------------------------------------------------------------------------------
	print ('%s' %('-'*105) )
	print ('%20.25s %29.30s %21.30s %11.12s %18.30s' %("Opt_vol Bohr^3 (Ang^3)", "Lattice_const Bohr (A)", "Bulk_modulus [GPa]", "Log(chi)", "Polynomial_order"))
	print ('%s' %('-'*105) )
	for order_of_fit in range(2, 11): #order of polynomial fitting
		if order_of_fit % 2 == 0: 
			order_of_fit = int(order_of_fit)
			fitr = np.polyfit(volume,energy,order_of_fit)
			curv = np.poly1d(fitr)
			oned = np.poly1d(np.polyder(fitr,1)) #
			bulk = np.poly1d(np.polyder(fitr,2))
			bpri = np.poly1d(np.polyder(fitr,3)) #
			vmin = np.roots(np.polyder(fitr))
			dmin=[]
			for i in range(len(vmin)):
				if (abs(vmin[i].imag) < 1.e-10): 
					if (volume[0] <= vmin[i] and vmin[i] <= volume[-1]): 
						if(bulk(vmin[i]) > 0): dmin.append(vmin[i].real)
			
			xvol = np.linspace(volume[0],volume[-1],100)
			if (len(dmin) > 1): print ("WARNING: Multiple minima are found!\n")
			if (len(dmin) == 0): print ("WARNING: No minimum in the given xrange!\n")
			
			chi = 0
			for i in range(len(energy)): 
				chi=chi+(energy[i]-curv(volume[i]))**2
			chi=math.sqrt(chi)/len(energy)
#-------------------------------------------------------------------------------			
			for i in range(len(dmin)):
				x0=dmin[len(dmin)-1-i]
				v0=dmin[len(dmin)-1-i]
				a0=(factor*v0)**(0.33333333333)
				b0=bulk(v0)*v0*unitconv
				#if (isym > 0): print ('Lattice constant = %5s %12.6f %12.6f' %(slabel,a0, a0*bohr_radius), '[Bohr, Angstrom]')	
				if (isym > 0): 
					print("%12.6f (%11.6f) %12.6f (%9.6f) %17.6f %13.2f %10d\n" %(v0, v0*bohr32ang3, a0, a0*bohr_radius, b0, math.log10(chi), order_of_fit), end="") 				
				else: 
					print("%12.6f(%12.6f) %12.6f(%12.6f) %17.6f %13.2f %10d\n" %(v0, v0*bohr32ang3, a0, a0*bohr_radius, b0, math.log10(chi), order_of_fit), end="")
	print ('%s' %('-'*105) )
#-------------------------------------------------------------------------------

	xlabel = u'Volume [Bohr\u00B3]'; ylabel = r'Energy [Ha]'

	fontlabel=20
	fonttick=14
	
	params = {'ytick.minor.size': 6,
			'xtick.major.pad': 8,
			'ytick.major.pad': 4,
			'patch.linewidth': 2.,
			'axes.linewidth': 2.,
			'lines.linewidth': 1.8,
			'lines.markersize': 8.0,
			'axes.formatter.limits': (-4, 6)}
	
	plt.rcParams.update(params)
	plt.subplots_adjust(left=0.21, right=0.93,
						bottom=0.18, top=0.88,
						wspace=None, hspace=None)
							
	yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)                           
							
	figure = plt.figure(1, figsize=(9,9))  
	ax     = figure.add_subplot(111)
	ax.text(0.5,-0.18,xlabel,size=fontlabel, transform=ax.transAxes,ha='center',va='center')
	ax.text(-0.25,0.5,ylabel,size=fontlabel, transform=ax.transAxes,ha='center',va='center',rotation=90)
	for line in ax.get_xticklines() + ax.get_yticklines():
		line.set_markersize(6)
		line.set_markeredgewidth(2)
	plt.xticks(size=fonttick)
	plt.yticks(size=fonttick)
	pyl.grid(True)
	############***************************************************************************

	file=open("tmp","w")
	for j in range( len(xvol) ):
		file.write("{:12.8f} {} {:12.8f}\n".format(xvol[j]," ",curv(xvol[j]) ) )
	#for i in range( len(strain) ):
	#	file.write( "{} {}\n".format(strain[i],energy[i]) )
	file.close()
	
	os.system("paste -d ' ' tmp energy-vs-volume > volume.dat ")
	#subprocess.call(("paste -d  tmp energy-vs-volume > volume.dat "), shell = True)
	os.system("rm tmp")
			
	############***************************************************************************
	plt.plot(xvol,curv(xvol),'b-',label='n='+str(order_of_fit)+' fit')
	plt.plot(volume,energy,'go',label='calculated')
	plt.plot(dmin,curv(dmin),'ro')
	plt.legend(loc=9,borderaxespad=.8,numpoints=1)
	
	ymax  = max(max(curv(xvol)),max(energy))
	ymin  = min(min(curv(xvol)),min(energy))
	dxx   = abs(max(xvol)-min(xvol))/18
	dyy   = abs(ymax-ymin)/18
	ax.yaxis.set_major_formatter(yfmt)
	ax.set_xlim(min(xvol)-dxx,max(xvol)+dxx)
	ax.set_ylim(ymin-dyy,ymax+dyy)
	
	ax.xaxis.set_major_locator(MaxNLocator(7))
	
	plt.savefig('PLOT.png',orientation='portrait',format='png')
###

def sortvolume(s,e):
    ss=[]; ee=[]; ww=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee
###

'''
#####---------------------------------------------------------
#            6- EVALUATE Mechanical properties from Cij Matrix
#####--------------------------------------------------------- 
'''


def mechanical_properties():
	import numpy as np
	import math, scipy, os, sys
	from numpy import linalg as LA
	import statistics as st
	import ase.io
	
	np.set_printoptions(precision=3)
	CRED = '\033[91m';CEND = '\033[0m'
	CYEL = '\033[33m'; CEND = '\033[0m'
	CPIN = '\033[46m';
	# Elastic constants have been obtained from SOEC approach mentioned in the above 
	# paper. For three distortions of the crystal three constants are extracted: 
	# C11, C44, C12. For cubic system the matrix is symmetric.

	print (CRED + "Values has to be supplied in the script::" + CEND)

	########################## INPUT PARAMETERS assuming cubic #######################
	
	c11=c22=c33=123.22 ; 
	c44=c55=c66=47
	B = 120
	#c12=c21=c13=c31=c23=c32=(1/6) * (9 * B - 3 * c11)
	c12=c21=c13=c31=c23=c32=94.22
	c14=c15=c16=c24=c25=c26=c34=c35=c36=0
	c45=c46=c56=0
	################################################################################

	print (CRED +"{:_^80s}".format("The STIFNESS MATRIX Cij is:")+ CEND)
	Cij=[   [c11,c12,c13,c14,c15,c16], 
			    [c21,c22,c23,c24,c25,c26], 
			    [c31,c32,c33,c34,c35,c36],
			    [0  ,0  ,0  ,c44,c45,c46],
			    [0  ,0  ,0  ,0  ,c55,c56],
			    [0  ,0  ,0  ,0  ,0  ,c66] ]
	Cij=np.matrix(Cij).reshape((6,6))
	print (Cij)
	
	########------------Calculate variance and Mean of the Cij -----------
	# ONLY for C11 , C22 , C33  
	Cij_stat =  [ Cij[0,0] , Cij[1,1] , Cij[2,2] ]
	print (">>>", (Cij_stat) )
	STD = np.std(Cij_stat); Var = np.var(Cij_stat) ; Mean = np.mean(Cij_stat)
	print ("{:6.3s} {:6.3s} {:6.3s}".format("C", "STD", "Var") )
	print ("{:6.3f} {:6.3f} {:6.3f}".format(Mean, STD/np.sqrt(3), Var) )

	######## ------------------------------------------------------------
	
	evals, eigenvec = LA.eig(Cij)
	print ("-"*40)
	print("Eigenvalues are: ", evals>0)
	print ("-"*40)
	# ### Compliance tensor  s_{ij}$ $(GPa^{-1})$
	#  s_{ij} = C_{ij}^{-1}$
	
	print (CRED +"{:_^80s}".format("The COMPLIANCE MATRIX Sij is:")+ CEND)
	Sij = np.linalg.inv(Cij)

	print ("{} ".format(Sij))
	print ("-"*80)
	
	######## -----------------------------VOIGT-------------------------------
  
	'''Voigt bulk modulus  (GPa)'''
	
	#9K_v = (C_{11}+C_{22}+C_{33}) + 2(C_{12} + C_{23} + C_{31}) 
	Bv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[1,2] + Cij[2,0])) / 9.0
	
	'''Voigt shear modulus  (GPa)'''
	
	#15*G_v = (C_{11}+C_{22}+C_{33}) - (C_{12} + C_{23} + C_{31}) + 3(C_{44} + C_{55} + C_{66})$
	Gv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) - (Cij[0,1] + Cij[1,2] + Cij[2,0]) 
	+ 3 * (Cij[3,3] + Cij[4,4] + Cij[5,5]))/15.0
	
	## Young's: Voigt
	Ev = (9*Bv*Gv)/(3*Bv + Gv)
		
	## Poisson's ratio: Voigt
	NuV = (3*Bv - Ev)/(6*Bv)
	
	######## -----------------------------REUSS-------------------------------
	
	# Reuss bulk modulus  K_R  $(GPa)$
	#  1/K_R = (s_{11}+s_{22}+s_{33}) + 2(s_{12} + s_{23} + s_{31})$
	Br = 1/((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2*(Sij[0,1] + Sij[1,2] + Sij[2,0])) 
	
	# Reuss shear modulus  G_v  $(GPa)$
	# 15/G_R = 4*(s_{11}+s_{22}+s_{33}) - 4*(s_{12} + s_{23} + s_{31}) + 3(s_{44} + s_{55} + s_{66})$
	Gr = (4 * (Sij[0,0] + Sij[1,1] + Sij[2,2]) - 4*(Sij[0,1] + Sij[1,2] + Sij[2,0]) 
	+ 3 * (Sij[3,3] + Sij[4,4] + Sij[5,5]))
	Gr = 15.0/Gr
	
	## Young's: Reuss
	Er = (9*Br*Gr)/(3*Br + Gr)	
	
	## Poisson's ratio: Reuss
	NuR = (3*Br - Er)/(6*Br)

	##########################################################################

	######## -----------------------------Averages-------------------------------
	
	# #Hill bulk modulus  K_{VRH}$ $(GPa)$
	#  K_{VRH} = (K_R + K_v)/2 
	B_H = (Bv + Br)/2
	#print ("VRH bulk modulus  (GPa): %20.8f " %(B_H) )
	
	# Hill shear modulus  G_{VRH}$ $(GPa)$
	#  G_{VRH} = (G_R + G_v)/2 
	G_H = (Gv + Gr)/2
	#print ("VRH shear modulus (GPa): %20.8f " %(G_H) )
	
	# Young modulus E = 9BG/(3B+G)
	#E_H = (9 * B_H * G_H) / (3 * B_H + G_H)
	E_H = (Ev + Er)/2
	#print ("Young modulus E : {:1.8s} {:20.8f}".format(" ",E_H) )
	
	# ### Isotropic Poisson ratio $\mu 
	# $\mu = (3K_{VRH} - 2G_{VRH})/(6K_{VRH} + 2G_{VRH})$
	#nu_H = (3 * B_H - 2 * G_H) / (6 * B_H + 2 * G_H )
	nu_H = (NuV + NuR) / 2
	#print ("Isotropic Poisson ratio: {:15.8f} ".format(nu_H) )
	
	## Elastic Anisotropy
	## Zener anisotropy for cubic crystals only
	A = 2*(c44)/(c11-c12)
	
	# Universal Elastic Anisotropy AU
	AU = (Bv/Br) + 5*(Gv/Gr) - 6.0
	
	# C' tetragonal shear modulus
	C = (c11-c12)/2
	
	ratio_V = Bv/Gv
	ratio_R = Br/Gr
	
	print ("{:30.8s} {:20.8s} {:20.8s} {:20.8s}".format(" ","Voigt", "Reuss ", "Hill") )
	print ("{:_^80}".format("GPa"))	
	print ("{:16.20s} {:20.3f} {:20.3f} {:20.3f}".format("Bulk Modulus",Bv, Br, B_H) )
	print ("{:16.20s} {:20.3f} {:20.3f} {:20.3f}".format("Shear Modulus",Gv, Gr, G_H) )
	print ("{:16.20s} {:20.3f} {:20.3f} {:20.3f}".format("Young Modulus",Ev, Er, E_H) )
	print ("{:16.20s} {:20.3f} {:20.3f} {:20.3f}".format("Poisson ratio ", NuV, NuR, nu_H) )
	print ("{:16.20s} {:20.3f} {:20.3f} {:20.3f}({:5.3f})".format("B/G ratio ",Bv/Gv,Br/Gr, B_H/G_H, G_H/B_H) )
	print ("{:16.20s} {:20.3s} {:20.3s} {:20.3f}".format("Avr ratio ",'','', (Gv-Gr)/(Gv+Gr)) )
	print ("{:16.20s} {:20.3s} {:20.3s} {:20.3f}".format("Zener ratio ",'','', A) )
	print ("{:16.20s} {:20.3s} {:20.3s} {:20.3f}".format("AU ",'','', AU) )
	print ("{:16.20s} {:20.3s} {:20.3s} {:20.3f}".format("Cauchy pressure ",'','', (c12-c44)) )
	print ("{:16.20s} {:20.3s} {:20.3s} {:20.3f}".format("C'tetra Shear ",'','',  C) )
	
	print ("-"*80)
	return Sij
###


'''
#####---------------------------------------------------------
#            				 INTRODUCTION
#####--------------------------------------------------------- 
'''

####
def Introduction():

	print(colored('@'*80,'yellow'), end = '\n', flush=True)
	print("Credit	:	Asif Iqbal BHATTI")
	print("USAGE	:	Extract essential properties from VASP outputfiles.") 
	print("VERSION	:	python3 or above ")
	print("FORMAT	:	POSCAR VASP5 format                    ")
	print("DATE	:	28/12/2019                             ")
	print('USAGE	:	execute by typing python3 sys.argv[0]')
	print("VERSION	:	4.0	                                   ")
	print(colored('@'*80,'yellow'), end = '\n', flush=True)
	print("              ____| Python script to process various properties |____")

'''
#####---------------------------------------------------------
#            				 MAIN ENGINE
#####--------------------------------------------------------- 
'''

####
if __name__ == "__main__":
	
	Introduction()

	print("Number of processors Detected: ", mp.cpu_count())
	print(Back.MAGENTA + ' NB: POSCAR should be in VASP 5 format & without selective dynamics', end = '\n', flush=True)
	print(Style.RESET_ALL)	
	print(colored('~'*80,'red'), end = '\n', flush=True)
	print("**** Following are the options: ")
	print(colored('-'*80,'red'), end = '\n', flush=True)

	print("(01) To execute only POSCAR file (local lattice distortion DEF 1-3)")
	print("(02) To execute POSCAR CELL VOLUME DIFFERENCE with final CONTCAR file")
	print("(03) To extract ENERGY from directories")
	print("(04) To extract ELASTIC CONSTANTS from OUTCAR file (IBRION=6,ISIF=3)")
	print("(05) To Fit energy vs volume curve to extract Elastic Moduli: B0 ... ")
	print("(06) CONVERT POSCAR file from VASP4 to VASP5 format")
	print("(07) Calculate manually Elastic properties by entering Cij values (Energy-vs-Strain)")
	print("(08) Calculate lattice distortion of the structure")
	print("(09) Create directories for PHONON calculations generated with PHONOPY code")
	
	print (colored('~'*80,'red'), end = '\n', flush=True)

	option = input("Enter the option as listed above: ")
	option = int(option)
	
	if (option == 1):
		poscar()
		
	elif (option == 2):
		VOL_P = main_poscar()
		VOL_C = main_contcar()
		volume_diff(VOL_P, VOL_C)
		
	elif (option == 3):
		create_energy_vs_volume()
		
	elif (option == 4):
		print("OUTCAR should be in the same directory from which this script is run ")		
		pool = mp.Pool(mp.cpu_count())
		elastic_matrix_VASP_STRESS()
		pool.close()
		
	elif (option == 5):
		fitting_energy_vs_volume_curve_ELASTIC()
		
	elif (option == 6):
		poscar_VASP42VASP5()	

	elif (option == 7):
		mechanical_properties()

	elif (option == 8):
		#lattic_distortion.local_lattice_distortion(a,b,c)
		lattic_distortion.local_lattice_distortion_DEF1()
		#lattic_distortion.local_lattice_distortion_DEF2()

	elif (option == 9):
		pass
		
	else:
		print ("INVALID OPTION")






