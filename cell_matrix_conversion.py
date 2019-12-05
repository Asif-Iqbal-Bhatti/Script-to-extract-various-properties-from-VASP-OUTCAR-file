#!/usr/bin/python3
#####---------------------------------------------------------
#####---------------------------------------------------------
#    Credit: Asif Iqbal BHATTI
#    CODE to convert Cell Matrix to Cell Parameters
#    python3 this_file
#    POSCAR VASP5 format
#    DATE: 05/12/2019
#####---------------------------------------------------------
#####---------------------------------------------------------

import os, sys
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
			#print ("Types of elements:", str(elementtype), end = '')
			numberofatoms=file.readline()
			#print ("Number of atoms:", (numberofatoms), end = '')
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
			print ('b=', b)
			print ('c=', c)
			gamma = math.degrees(math.acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))))
			alpha = math.degrees(math.acos(np.dot(b,c) / (np.linalg.norm(b) * np.linalg.norm(c))))
			beta  = math.degrees(math.acos(np.dot(a,c) / (np.linalg.norm(a) * np.linalg.norm(c))))
			print ('\u03B1=', alpha, '\u03B2=', beta, '\u03B3=', gamma)
			print ("#####------------------------------------------------")
			print ('||a||=', np.linalg.norm(a))
			print ('||b||=', np.linalg.norm(b))
			print ('||c||=', np.linalg.norm(c)) 
			break
		else:
			print ('NO file eneterd or wrong file') 
			break			
#####---------------------------------------------------------
# Looping over all directories in a current folder containing POSCARS files
#####--------------------------------------------------------- 

def main():
	os.system("rm out.dat")
	mypath = os.getcwd()
	#print (mypath)
	print ("               >>>>> Converting Cell Matrix to Cell Parameters <<<<<<")
	for entry in os.listdir(mypath):
		if os.path.isdir(os.path.join(mypath, entry)):
			#print (entry)
			for file in os.listdir(entry):
				#print (file)
				filepath = os.path.join(entry, file)
				#f = open(filepath, 'r')
				#print (f.read())
				#f.close()	
				fo = open(filepath, 'r')
				ofile=open('out.dat','a+')
				print (colored('>>>>>>>>  Name of the file: ','yellow'), fo.name, end = '\n', flush=True)
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
				print ("#####------------------------------------------------")
				print ('a=', a)
				ofile.write ("'a=' {}\n".format(a))
				print ('b=', b)
				ofile.write ("'b=' {}\n".format(b))
				print ('c=', c)
				ofile.write ("'c=' {}\n".format(c))			
	#####---------------------------------------------------------
	### gamma = Cos-1( (a.b)/||a||.||b|| )
	### alpha = Cos-1( (b.c)/||b||.||c|| )
	### beta  = Cos-1( (a.c)/||a||.||c|| )
	
				gamma = math.degrees(math.acos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))))
				alpha = math.degrees(math.acos(np.dot(b,c) / (np.linalg.norm(b) * np.linalg.norm(c))))
				beta  = math.degrees(math.acos(np.dot(a,c) / (np.linalg.norm(a) * np.linalg.norm(c))))
				print ('\u03B1=', alpha, '\u03B2=', beta, '\u03B3=', gamma)
				ofile.write ("'\u03B1=' {} '\u03B2=' {} '\u03B3=' {}\n".format(alpha,beta,gamma))
				print ("#####------------------------------------------------")
				print ('||a||=', np.linalg.norm(a))
				ofile.write ("'||a||=' {}\n".format(np.linalg.norm(a)))
				print ('||b||=', np.linalg.norm(b))
				ofile.write ("'||b||=' {}\n".format(np.linalg.norm(b)))			
				print ('||c||=', np.linalg.norm(c)) 
				ofile.write ("'||c||=' {}\n".format(np.linalg.norm(c)))			
				ofile.write ("***************************************************\n")
				ofile.close()

if __name__ == "__main__":

	poscar()

	#main()

