import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path

ang2atomic 	= 	1.889725988579  # 1 A = 1.889725988579 [a.u]
Ang32Bohr3	=	6.74833304162   # 1 A^3 = 6.7483330371 [a.u]^3
eV2Hartree	=	0.036749309
	
def energy_vs_volume():
	import fnmatch
	mypath = os.getcwd()
	os.system("rm energy-vs-volume energy-vs-strain")
	eV2Hartree=0.036749309
	Ang32Bohr3=6.74833304162
	
	E=[]; dir_list=[]; count = 0; dir_E=[];
	vol_cell=[]; strain_file=[]; strain_value=[] # strain_value is deformation
	
	print ("               >>>>> Extracting Energy from directories  <<<<<<")
	for entry in os.listdir(mypath):
		if not os.path.exists('strain-01'):
			print (' ERROR: strain-* does not exist here.')
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
