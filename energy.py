import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path

ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
ang2bohr   = 6.7483330371   # 1 A^3 = 6.7483330371 [a.u]^3

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