#!/usr/bin/env python3

'''
# *****************************************************************************
#   USAGE :: python3 grad2.py
#   AUTHOR:: ADAPTATION OF Peter Larsson script to new version @ASIFIQBAL
#   please see => https://www.nsc.liu.se/~pla/vasptools/
#   ABOUT THE PROGRAM:
#   It prints out the forces at each ionic steps during simulation
# *****************************************************************************
'''

import os, sys, re
import numpy as np
import subprocess

def get_number_of_atoms(out_car):
  return int(subprocess.Popen('grep "NIONS" ' + out_car, stdout=subprocess.PIPE, shell=True).communicate()[0].split()[11])
  #print( int(subprocess.Popen('grep "NIONS" ' + where, stdout=subprocess.PIPE, shell=True).communicate()[0].split()[11]) )

def get_ediff(out_car):
  return float(subprocess.Popen('grep "  EDIFF" ' + out_car, stdout=subprocess.PIPE, shell=True).communicate()[0].split()[2])
  #print( float(subprocess.Popen('grep "  EDIFF" ' + where, stdout=subprocess.PIPE, shell=True).communicate()[0].split()[2]) )

OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'

# *****************************************
#             Main program
# *****************************************
try:
	outcar = open(sys.argv[1],"r")
except IOError:
	sys.stderr.write(FAIL)
	sys.stderr.write("No OUTCAR file")
	sys.stderr.write(ENDC+"\n")
	sys.exit(1)

if outcar != None:
	outcarfile = sys.argv[1]
	outcarlines = outcar.readlines()
	
	#Find max iterations
	nelmax = int(
			subprocess.Popen(
					'grep "NELM " ' + outcarfile, stdout=subprocess.PIPE,
					shell=True).communicate()[0].split()[2][:-1])
	natoms = get_number_of_atoms(outcarfile)
	ediff = np.log10(float(get_ediff(outcarfile)))
	
	re_energy = re.compile("free  energy")
	re_iteration = re.compile("Iteration")
	re_force = re.compile("TOTAL-FORCE")
	re_mag = re.compile("number of electron")
	re_volume = re.compile("volume of cell")
	
	lastenergy = 0.0
	energy = 0.0
	steps = 1
	iterations = 0
	cputime = 0.0
	totaltime = 0.0
	dE = 0.0
	magmom = 0.0
	spinpolarized = False
	volume = 0.0
	#average = 0.0
	#maxforce = 0.0
	
	i = 0
	for line in outcarlines:
				
		if re_iteration.search(line):
			iterations = iterations + 1
	
		if re_force.search(line):
			# Calculate forces here...
			forces = []
			magnitudes = []
			for j in range(natoms):
				parts = outcarlines[i+j+2].split()
				x = float(parts[3])
				y = float(parts[4])
				z = float(parts[5])
				forces.append([x,y,z])
				magnitudes.append(np.sqrt(x**2 + y**2 + z**2))
	
			average = sum(magnitudes)/natoms
			maxforce = max(magnitudes)
	
		if re_mag.search(line):
			parts = line.split()
			if len(parts) > 5 and parts[0].strip() != "NELECT":
				spinpolarized = True
				magmom = float(parts[5])
	
		if re_volume.search(line):
			parts = line.split()
			if len(parts) > 4:
				volume = float(parts[4])
	
		if re_energy.search(line):
			lastenergy = energy
			energy = float(line.split()[4])
			dE = np.log10(abs(energy-lastenergy+1.0E-12))
	
			# CONSTRUCT OUTPUT STRING
			try:
				stepstr   = str(steps).rjust(4)
				energystr = f'Energy: {energy:3.6f}'.rjust(12)
				logdestr  = f'Log|dE|: {dE:1.3f}'.rjust(6)		
				iterstr   = f'SCF: {iterations:3d}'
				avgfstr   = f'Avg|F|: {average:2.3f}'.rjust(6)
				maxfstr   = f'Max|F|: {maxforce:2.3f}'.rjust(6)
				volstr    = f'Vol.: {volume:3.1f}'.rjust(5)
			except NameError:
				print(f"Cannot understand this OUTCAR file...try to read ahead")
				continue
	
			if iterations == nelmax:
				sys.stdout.write(FAIL)
				#print "         ^--- SCF cycle reached NELMAX. Check convergence!"
	
			if (dE < ediff):
				sys.stdout.write(OKGREEN)
	
			if spinpolarized:
				magstr=f'Mag: {magmom:2.2f}'.rjust(6)
				print (f'{stepstr},{energystr},{logdestr},{iterstr},{avgfstr},{maxfstr},{volstr},{magstr}')
			else:
				print (f'{stepstr},{energystr},{logdestr},{iterstr},{avgfstr},{maxfstr},{volstr}')
	
			sys.stdout.write(ENDC)
	
			steps = steps + 1
			iterations = 0
			totaltime = totaltime + cputime
			cputime = 0.0
	
		i = i + 1
	
