#!/usr/bin/env python3

from lxml import etree
import xml.etree.ElementTree as xml
from math import *
import sys, getopt, os
from sys import stdout
from array import *
import numpy as np

from collections import Counter
from termcolor import colored

#*******************************DOCUMENTATION***********************************
#
#AUTHOR 	 : ASIF IQBAL BHATTI
#CREATED ON	 : 29/04/2020
#USAGE	   	 : To parse an EXCITING xml input file to VASP format.
#Output file has already been defined in the code (POSCAR.vasp).
#CAUTION: Use at your own risk (NOTEVEN IMPLIED GUARANTEED, WHATSOEVER),
#the code has been tested but the user in the end will have to verify the ouput 
#geometry and be vary of the order of atoms in list.
#
#************************END OF DOCUMENTATION***********************************

Bohr2angstrom = 0.529177249
yellow        = '\033[01;33m'
red           = '\033[01;31m'
blue          = '\x1b[01;34m'
cyan          = '\x1b[01;36m'
green         = yellow

#*********************DEFINITION FOR INPUTFILE**********************************

if len(sys.argv) < 2:
	sys.exit(f'Usage: {sys.argv[0]} <Exciting xml inputfile>')

tree = xml.parse(sys.argv[1])
root = tree.getroot()

struct = root.getchildren()[1] # pivot for controling element
lattice = struct.getchildren()


with open('POSCAR.vasp', 'w+') as f:
	print (colored("\nExtracting information from input file... ", 'yellow'))
	print ('-'*80)
	print (colored(" | Title is:  ", 'yellow'), root.findtext('title'))
	f.write(root.findtext('title'))

	#*************************DETECTING # OF SPECIES******************************

	i = 0
	for _ in struct:
		i = i+1
	print (colored(" | Species detected:  ", 'yellow'), i-1	)
	f.write(" |\t Number of species detected:  {:d}\n".format(i-1))

	#***************************SCALING DEFINITION********************************
	#reading from a file
	for sca in struct.findall('crystal'):
		for key in sca.attrib:
			sc = sca.get('scale') if key == 'scale' else 1.0
			print (colored(" | Scaling factor is:  ", 'yellow'), float(sc))
	print (colored(" | Converting bohr to angstrom ... ", 'yellow'))
	print (colored(" | Writing to a File in VASP format ... ", 'yellow'))
	f.write("%f\n" % 1.0)
	print ('-'*80)

	#***********************LATTICE VECTORS SECTION************************
	#reading from a file
	lv = []
	for vect in struct.findall('crystal'):
			#lat = vect.find('basevect')
		lv.extend(k.text.split() for k in vect.findall('basevect'))
	###

	convertion = Bohr2angstrom*float(sc)
	for list in lv:
		print ("\t{:10.8f} {:10.8f} {:10.8f}".format(float(list[0])*convertion, 
		float(list[1])*convertion, float(list[2])*convertion ) )
		f.write ("\t{:10.8f} {:10.8f} {:10.8f}\n".format(float(list[0])*convertion, 
		float(list[1])*convertion, float(list[2])*convertion ) )	

	#******************SPECIES CONTROL SECTION*********************************

	lab = []
	for coord in struct.findall('species'):
		spe = coord.get('speciesfile')
		jj = os.path.splitext(spe)
		lab.append(jj[0])

	sp = []
	for x in range(i-1):
		print("{:s}".format(lab[x]))
		f.write("\t%s" % lab[x])
		sp.extend(lab[x] for _ in struct[x+1].iter('atom'))
	f.write("\n")

	for h in range(i-1):
		print("\t{:d}".format(sp.count(lab[h])) )
		f.write("\t%d" % sp.count(lab[h]))
	f.write("\n")

	#*******************CONTROL SECTION FOR CARTESIAN COORDINATE***************

	lll = []
	ccl = []
	checker = False
	for cart in root.findall('structure'):
		ca = cart.get('cartesian')
		for key in cart.attrib:
			if key == 'cartesian':
				checker = True	
		if ca == "true":
			print ("Cartesian")
			f.write("Cartesian\n")
			for i in range(i-1):
				for s in struct[i+1].iter('atom'):
					atom = s.get('coord')
					ccl.append(atom.split())
					lll.append(lab[i])
#				print (atom, lab[i])

			for list in ccl:
				print ("{:5.8f} {:5.8f} {:5.8f}".format(float(list[0])* Bohr2angstrom, 
				float(list[1])* Bohr2angstrom, float(list[2])* Bohr2angstrom ) )
				f.write ("{:5.8f} {:5.8f} {:5.8f}\n".format(float(list[0])* Bohr2angstrom, 
				float(list[1])* Bohr2angstrom, float(list[2])* Bohr2angstrom ) )
#**************************DIRECT COORDINATE*********************

		else:
			print ("Direct")
			f.write  ("Direct\n")
			for i in range(i-1):
				for s in struct[i+1].iter('atom'):
					atom = s.get('coord')
					ccl.append(atom.split())
					lll.append(lab[i])
#				print (atom, lab[i])
#				print (ccl)

			for list in ccl:
				print ("{:5.8f} {:5.8f} {:5.8f}".format(float(list[0]), 
				float(list[1]), float(list[2]) ) )
				f.write ("{:5.8f} {:5.8f} {:5.8f}\n".format(float(list[0]), 
				float(list[1]), float(list[2]) ) )

	print (colored("File generated ... ", 'green'))
