import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path
from colorama import Fore, Back, Style, init

ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
ang2bohr   = 6.7483330371   # 1 A^3 = 6.7483330371 [a.u]^3

def main_poscar():

	count = 0
	os.system("rm out_POSCARS.dat")
	VOL_P = []
	pos = []
	kk = []
	lattice = [];
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
					count+=1
					sum = 0
					filepath = os.path.join(entry, file)
					with open(filepath, 'r') as fo:
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
						for _ in range(int(numberofatoms)):
							coord = fo.readline().split()
							coord = [float(i) for i in coord]
							pos = pos + [coord]
						pos = np.array(pos)
						#print (pos)

						ofile.write ("\n")
					Latvec1=Latvec1.split()
					Latvec2=Latvec2.split()
					Latvec3=Latvec3.split()

					a = [float(ai) for ai in Latvec1]
					b = [float(bi) for bi in Latvec2]
					c = [float(ci) for ci in Latvec3]
					#print ('a=', a)
					ofile.write(f"'a=' {a}\n")
					#print ('b=', b)
					ofile.write(f"'b=' {b}\n")
					#print ('c=', c)
					ofile.write(f"'c=' {c}\n")
					lld = local_lattice_distortion(a,b,c)
##########################---------------------------------------------------------

					alpha, beta, gamma = lattice_angles(a,b,c)
					VOL_POS = np.dot(a, np.cross(b,c))
					VOL_P.append(VOL_POS)	

					ofile.write(f"'\u03B1=' {alpha} '\u03B2=' {beta} '\u03B3=' {gamma}\n")
					ofile.write(f"'||a||=' {np.linalg.norm(a)}\n")
					ofile.write(f"'||b||=' {np.linalg.norm(b)}\n")
					ofile.write(f"'||c||=' {np.linalg.norm(c)}\n")
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
	VOL_C = []
	pos = []
	kk = []
	lattice = []
	sum = 0
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
					count+=1
					sum = 0
					filepath = os.path.join(entry, file)
					with open(filepath, 'r') as fo:
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
						for _ in range(int(numberofatoms)):
							coord = fo.readline().split()
							coord = [float(i) for i in coord]
							pos = pos + [coord]
						pos = np.array(pos)
						#print (pos)

						ofile.write ("\n")
					Latvec1=Latvec1.split()
					Latvec2=Latvec2.split()
					Latvec3=Latvec3.split()
					a = [float(ai) for ai in Latvec1]
					b = [float(bi) for bi in Latvec2]
					c = [float(ci) for ci in Latvec3]
					#print ('a=', a)
					ofile.write(f"'a=' {a}\n")
					#print ('b=', b)
					ofile.write(f"'b=' {b}\n")
					#print ('c=', c)
					ofile.write(f"'c=' {c}\n")
					lld = local_lattice_distortion(a,b,c)
##########################---------------------------------------------------------

					alpha, beta, gamma = lattice_angles(a,b,c)
					VOL_CON = np.dot(a, np.cross(b,c))
					VOL_C.append(VOL_CON)

					ofile.write(f"'\u03B1=' {alpha} '\u03B2=' {beta} '\u03B3=' {gamma}\n")
					ofile.write(f"'||a||=' {np.linalg.norm(a)}\n")
					ofile.write(f"'||b||=' {np.linalg.norm(b)}\n")
					ofile.write(f"'||c||=' {np.linalg.norm(c)}\n")
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
	print ("-"*100)
	print ("VOL Diff A^3 %18s %12s %15.15s" %("CONTCAR",  "POSCAR",  "contcar-poscar"))
	print ("-"*100)
	for i in range(int(n)):
		print ("The difference is: %12.6f %12.6f %15.8f " %(VOL_C[i], VOL_P[i], VOL_C[i] - VOL_P[i]) )
		
		