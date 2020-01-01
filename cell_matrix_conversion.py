#!/usr/bin/python3

'''
#####---------------------------------------------------------
#####---------------------------------------------------------
#    Credit	: 	Asif Iqbal BHATTI
#    CODE to: 	convert Cell Matrix to Cell Parameters
#    VERSION: 	This script runs with python3 or later
#    FORMAT	:	POSCAR VASP5 format
#    DATE	: 	28/12/2019
#    USAGE	: 	python3 sys.argv[0]
#####---------------------------------------------------------
#####---------------------------------------------------------
'''

import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path
from termcolor import colored
import multiprocessing as mp



'''
##########################---------------------------------------------------------
#		1-		 Reading only POSCAR file from the command prompt in a given directory
##########################---------------------------------------------------------
'''

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

'''
#####---------------------------------------------------------
#        2- Looping over all directories containing POSCAR & CONTCAR files
#####--------------------------------------------------------- 
'''

def main_poscar():
	count = 0
	os.system("rm out.dat")
	VOL_P = []; pos = []; kk = []; lattice = [];
	mypath = os.getcwd()
	#print (mypath)

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
					elementtype = elementtype.split()						
					#print ("Types of elements:", str(elementtype), end = '')
					#ofile.write (str(elementtype))
					numberofatoms=fo.readline()
					#print ("Number of atoms:", (numberofatoms), end = '')
					#ofile.write ((numberofatoms))
					Coordtype=fo.readline()

##########################---------------------------------------------------------
					print ("**********-------------------# of Atoms--------------------")
					
					nat = numberofatoms.split()
					nat = [int(i) for i in nat]
					print (nat)
					for i in nat:
						sum = sum + i
					numberofatoms = sum
					print ("Number of atoms:", (numberofatoms), end = '\n')
##########################---------------------------------------------------------					
					print ("//////---------------Atomic positions-----------------")
					print ("Coordtype:", (Coordtype), end = '')						
					for x in range(int(numberofatoms)):
						coord = fo.readline().split()
						coord = [float(i) for i in coord]
						pos = pos + [coord]
					pos = np.array(pos)
					print (pos)
					
					ofile.write ("\n")			
					fo.close()
##########################---------------------------------------------------------

					a=[]; b=[]; c=[];
					Latvec1=Latvec1.split()
					Latvec2=Latvec2.split()
					Latvec3=Latvec3.split()
					
##########################---------------------------------------------------------
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
					
##########################---------------------------------------------------------
		
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
	VOL_C = [];	pos = []; kk = []; lattice = []; sum = 0
	mypath = os.getcwd()
	
	for entry in os.listdir(mypath):
		if os.path.isdir(os.path.join(mypath, entry)):	
			for file in os.listdir(entry):
				
				if file == "CONTCAR":
					count+=1; sum = 0
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
					elementtype = elementtype.split()					
					numberofatoms=fo.readline()
					Coordtype=fo.readline()
					print ("Coordtype:", (Coordtype), end = '')
##########################---------------------------------------------------------
					print ("**********-------------------# of Atoms--------------------")
					
					nat = numberofatoms.split()
					nat = [int(i) for i in nat]
					print (nat)
					for i in nat:
						sum = sum + i
					numberofatoms = sum
					print ("Number of atoms:", (numberofatoms), end = '\n')
##########################---------------------------------------------------------						
					print ("//////---------------Atomic positions-----------------")				
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
##########################---------------------------------------------------------

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
	
	
#### math.sin function takes argument in radians ONLY
def volume(a,b,c,alpha,beta,gamma):
	length = np.linalg.norm(a) * np.linalg.norm(b) * np.linalg.norm(c) 
	volume = length * ( np.sqrt(1 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma) - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2) )
	vol_au = volume * ang2bohr
	return volume, vol_au

#### Ordering of angles does matter
def lattice_angles(a,b,c):
	### gamma = Cos-1( (a.b)/||a||.||b|| )
	### alpha = Cos-1( (b.c)/||b||.||c|| )
	### beta  = Cos-1( (a.c)/||a||.||c|| )
	gamma = math.degrees(math.acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))))
	alpha = math.degrees(math.acos(np.dot(b,c) / (np.linalg.norm(b) * np.linalg.norm(c))))
	beta  = math.degrees(math.acos(np.dot(a,c) / (np.linalg.norm(a) * np.linalg.norm(c))))
	return alpha, beta, gamma
####
def volume_diff(VOL_P, VOL_C):
	n=os.popen("find . -mindepth 1 -maxdepth 1 -type d | wc -l").read()
	print ("VOL Diff A^3 %18s %12s %15.15s" %("CONTCAR",  "POSCAR",  "contcar-poscar"))
	for i in range(int(n)):
		print ("The difference is: %12.6f %12.6f %15.8f " %(VOL_C[i], VOL_P[i], VOL_C[i] - VOL_P[i]) )
	
####

def energy():
	mypath = os.getcwd()
	E=[]; dir_list=[]; count = 0; dir_E=[]
	print ("               >>>>> Extracting Energy from directories  <<<<<<")
	for entry in os.listdir(mypath):
		if os.path.isdir(os.path.join(mypath, entry)):
			dir_list.append(entry); 
			
			for file in os.listdir(entry):
				if file == "OUTCAR":
					filepath = os.path.join(entry, file)
					f = open(filepath,'r')
					lines = f.readlines()
					f.close()
					for i in lines:
						if "  free  energy   TOTEN  =" in i:
							m=float(i.split()[4])
					E.append(m)
					count+=1
	#print (dir_list); print (E)
	
	for i in range(count):
		print ("_______|		", dir_list[i], "-->" , E[i] )

#########

def print_Cij_Matrix():
	Bij = []
	C = "C"
	for i in range(0,6):
		Bij.append([])
		for j in range(0,6):
			Bij[i].append((C + str(i) + str(j)))
	l = np.matrix(Bij)
	return l
	
def elastic_matrix():
	while True:
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
					
##########--------------------stress tensor------------------	
				
			Cij = np.matrix(c)
			np.set_printoptions(precision=4, suppress=True)
			print (Cij)
			print(print_Cij_Matrix())
			print ("\nEigen Values of the matrix Cij:")
			evals = np.linalg.eigvals(Cij)
			if evals.all() > 0:
				print(evals)
				print("All Eigen values are positive")
			
##########--------------------Compliance tensor------------------
##########--------------------  s_{ij} = C_{ij}^{-1}

			Sij = np.linalg.inv(Cij)
			
#----------------------------- ELASTIC PROPERTIES -----------------------------------

			stability_test(Cij, crystaltype)
			
#-------------------------------- Voigt bulk modulus  K_v  $(GPa)$---------------
#------------------- 9K_v = (C_{11}+C_{22}+C_{33}) + 2(C_{12} + C_{23} + C_{31}) 
			Kv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[1,2] + Cij[2,0]))/9.0
#-------------------------------- Reuss shear modulus  G_v  $(GPa)$------------------
#------------------- 15/G_R = 4(s_{11}+s_{22}+s_{33}) - 4(s_{12} + s_{23} + s_{31}) + 3(s_{44} + s_{55} + s_{66})$
			Gv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) - (Cij[0,1] + Cij[1,2] + Cij[2,0]) + 3 * (Cij[3,3] + Cij[4,4] + Cij[5,5]))/15.0
#-------------------------------- Reuss bulk modulus  K_r  $(GPa)$----------------
#-------------------  1/K_R = (s_{11}+s_{22}+s_{33}) + 2(s_{12} + s_{23} + s_{31})$
			Kr = 1/((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) )
#-------------------------------- Reuss shear modulus  G_r  $(GPa)$------------------
			Gr = 15/(4 * (Sij[0,0] + Sij[1,1] + Sij[2,2]) - 4 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) + 3 * (Sij[3,3] + Sij[4,4] + Sij[5,5]))
		
#------------------------------------------------------------------------------------

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
				
#------------------------------------------------------------------------------------
#--------------- Voigt-Reuss-Hill Approximation: average of both methods
			Kvrh = (Kv + Kr)/2.0
			Gvrh = (Gv + Gr)/2.0
			Mvrh = (MV + MR)/2.0
			Evrh = (Ev + Er)/2.0
			Nu_vrh = (Nu_V + Nu_R)/2.0
			KG_ratio_V = Kv/Gv
			KG_ratio_R = Kr/Gr
			KG_ratio_vrh = Kvrh/Gvrh			
#--------------- Isotropic Poisson ratio $\mu 
#--------------- $\mu = (3K_{vrh} - 2G_{vrh})/(6K_{vrh} + 2G_{vrh})$
			mu = (3 * Kvrh - 2 * Gvrh) / (6 * Kvrh + 2 * Gvrh )
			
#-----------------------------------------------------------------------
#------------------------------------------------------------------------------------
		
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
	
			
def ductile_test(ratio):
	if(ratio > 1.75):
		return "ductile"
	else:
		return "brittle"
		
	
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


def Introduction():
	global message
	message = "              ____| Python script to process various properties |____"
	print(message)
	
####
if __name__ == "__main__":

	Introduction()
	print("Number of processors Detected: ", mp.cpu_count())    
	print (colored(' ----------------------------------------------------       ','red'), end = '\n', flush=True)
	print (colored(' ----------------------------------------------------       ','red'), end = '\n', flush=True)
	print (colored(' ----------------------------------------------------       ','red'), end = '\n', flush=True)
	print ('>>> USAGE: execute by typing python3 sys.argv[0]')

	print ("***************************** Following are the options ... ")
	print ("(1) To process only POSCAR file (Convert Lattice Matrix to Lattice parameter)")
	print ("(2) To process only ENERGY from directories")
	print ("(3) To process only CELL VOLUME DIFFERENCE from directories by comparing with final CONTCAR file")
	print ("(4) To process ELASTIC CONSTANTS from OUTCAR file")
	
	print (colored(' ----------------------------------------------------       ','red'), end = '\n', flush=True)
	print (colored(' ----------------------------------------------------       ','red'), end = '\n', flush=True)
	print (colored(' ----------------------------------------------------       ','red'), end = '\n', flush=True)

	option = input("Enter the option as listed above: ")
	option = int(option)
	if (option == 1):
		poscar()
	elif (option == 2):
		energy()
	elif (option == 3):
		VOL_P = main_poscar()
		VOL_C = main_contcar()
		volume_diff(VOL_P, VOL_C)		
	elif (option == 4):
		print("Reading OUTCAR. OUTCAR should be in the same directory from which this script is run ")
		pool = mp.Pool(mp.cpu_count())
		elastic_matrix()
		pool.close()
	else:
		print ("INVALID OPTION")
	
	
	
	
	
	
