import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path

ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
ang2bohr   = 6.7483330371   # 1 A^3 = 6.7483330371 [a.u]^3
eV2Hartree = 0.036749309
Ang32Bohr3 = 6.74833304162
def energy_vs_volume():
	import fnmatch
	mypath = os.getcwd()
	os.system("rm energy-vs-volume energy-vs-strain")
	E=[]; dir_list=[]; count = 0; dir_E=[]; vol_cell=[]; strain=[]; a=[]
	print ("               >>>>> Extracting Energy from directories  <<<<<<")
	for entry in os.listdir(mypath):
		if fnmatch.fnmatchcase(entry,'strain-*'):
			f = open(entry,'r')
			lines = f.readline()
			a.append( float(lines) )
			f.close()
			if os.path.isfile(os.path.join(mypath, entry)):
				strain.append(entry)		
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
						if "  volume of cell :"	in i:
							v=float(i.split()[4])
					vol_cell.append(v)		
					E.append(m)
					count+=1	
	#print (dir_list); print (E); print (vol_cell)
	print ("Directory :%10.6s %14s %16s %25.20s " % ("Folder", "Energy (eV)", "Vol_of_cell", "strain_deformation" ))
	
	for i in range(count):
		print ("Folder name: %10.10s %16.8f %16.8f %16.12s %14.4f" % (dir_list[i], E[i], vol_cell[i], strain[i], a[i] ))
	#rc = subprocess.Popen(['bash', 'extract_energy.sh'])
		
	print (colored('ENERGIES & VOLUMES ARE WRITTEN IN ATOMIC UNITS TO A FILE <energy-vs-volume>','yellow'), end = '\n', flush=True)
	print (colored('IT WILL BE READ BY ELASTIC SCRIPTS FOR POSTPROCESSING eV-->Ha; A^3-->Bohr^3','yellow'), end = '\n', flush=True)	
	
	file = open("energy-vs-volume",'w')
	for i in range(count):
		file.write ("%12.6f %14.6f\n" %(vol_cell[i] * Ang32Bohr3, E[i] * eV2Hartree))	
	file.close()

	file = open("energy-vs-strain",'w')
	for i in range(count):
		file.write ("%12.6f %14.6f\n" %(a[i], E[i] * eV2Hartree))	
	file.close()
