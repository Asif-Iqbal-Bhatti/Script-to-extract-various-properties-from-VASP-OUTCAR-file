#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   numpy import *
import subprocess
import os.path
import shutil
import numpy, scipy
import math
import sys, os
import matplotlib.pyplot as plt
from colorama import Fore, Back, Style, init
init(autoreset=True)

#_______________________________________________________________________________

def Elastic_strain():
	print("*"*80)
	print ("________| http://exciting-code.org/nitrogen-energy-vs-strain-calculations")
	print ("________| \u03B7  = \u03B5 + 1/2*\u03B5**2")
	print ("________| r' = (I + \u03B5) * r ")
	print ("________| Deformation uses Voigt notation.")
	print("*"*80)
	
	if (str(os.path.exists('CONTCAR'))=='False'): 
		sys.exit("ERROR: Input file CONTCAR not found!\n")
	
	maximum_strain = float(input("\nEnter maximum Lagrangian strain [smax] >>>> "))
	if (maximum_strain == 0): 
		work_directory = 'workdir'
		if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
		os.mkdir(work_directory)
		
		os.system("cp CONTCAR workdir/CONTCAR-01")
		output_info = open('workdir/INFO-elastic-constants',"w")
		output_info.write("\nMaximum Lagrangian strain       = 0"         )
		output_info.write("\nNumber of strain values         = 1",        )
		output_info.write("\nVolume of equilibrium unit cell = 1.0 [A]^3",)
		output_info.write("\nDeformation code                = 000000",   )
		output_info.write("\nDeformation label               = single\n"  )
		output_info.close()
		
		output_eta = open('workdir/strain-01',"w")
		output_info.write("0.00")
		output_eta.close()
		print
		sys.exit("Single unstrained calculation\n")
	
	if (1 < maximum_strain or maximum_strain < 0): 
		sys.exit("ERROR: Maximum Lagrangian strain is out of range [0-1]!\n")
	strain_points = int(input("\nEnter the number of strain values in [-smax,smax] >>>> "))
	strain_points = int(abs(strain_points))
	if (3 > strain_points or strain_points > 99): 
		sys.exit("ERROR: Number of strain values is out of range [3-99]!\n")
	
	print(Style.RESET_ALL)	
	print (Back.GREEN + "------------------------------------------------------------------------" )
	print (Back.GREEN + " List of deformation codes for strains in Voigt notation"                 )
	print (Back.GREEN + "------------------------------------------------------------------------" )
	print (Back.YELLOW + "  0 =>  ( eta,  eta,  eta,    0,    0,    0)  | volume strain "           )
	print (Back.YELLOW + "  1 =>  ( eta,    0,    0,    0,    0,    0)  | linear strain along x "   )
	print (Back.GREEN + "  2 =>  (   0,  eta,    0,    0,    0,    0)  | linear strain along y "   )
	print (Back.GREEN + "  3 =>  (   0,    0,  eta,    0,    0,    0)  | linear strain along z "   )
	print (Back.GREEN + "  4 =>  (   0,    0,    0,  eta,    0,    0)  | yz shear strain"          )
	print (Back.GREEN + "  5 =>  (   0,    0,    0,    0,  eta,    0)  | xz shear strain"          )
	print (Back.GREEN + "  6 =>  (   0,    0,    0,    0,    0,  eta)  | xy shear strain"          )
	print (Back.YELLOW + "  7 =>  (   0,    0,    0,  eta,  eta,  eta)  | shear strain along (111)" )
	print (Back.GREEN + "  8 =>  ( eta,  eta,    0,    0,    0,    0)  | xy in-plane strain "      )
	print (Back.GREEN + "  9 =>  ( eta, -eta,    0,    0,    0,    0)  | xy in-plane shear strain" )
	print (Back.GREEN + " 10 =>  ( eta,  eta,  eta,  eta,  eta,  eta)  | global strain"            )
	print (Back.GREEN + " 11 =>  ( eta,    0,    0,  eta,    0,    0)  | mixed strain"             )
	print (Back.GREEN + " 12 =>  ( eta,    0,    0,    0,  eta,    0)  | mixed strain"             )
	print (Back.GREEN + " 13 =>  ( eta,    0,    0,    0,    0,  eta)  | mixed strain"             )
	print (Back.GREEN + " 14 =>  ( eta,  eta,    0,  eta,    0,    0)  | mixed strain"             )
	print (Back.GREEN + "------------------------------------------------------------------------" )
	print(Style.RESET_ALL)	
		
	deformation_code = int(input("\nEnter deformation code >>>> "))
	if (0 > deformation_code or deformation_code > 14): 
		sys.exit("ERROR: Deformation code is out of range [0-14]!\n")
	
	if (deformation_code == 0 ): dc='EEE000'
	if (deformation_code == 1 ): dc='E00000'
	if (deformation_code == 2 ): dc='0E0000'
	if (deformation_code == 3 ): dc='00E000'
	if (deformation_code == 4 ): dc='000E00'
	if (deformation_code == 5 ): dc='0000E0'
	if (deformation_code == 6 ): dc='00000E'
	if (deformation_code == 7 ): dc='000EEE'
	if (deformation_code == 8 ): dc='EE0000'
	if (deformation_code == 9 ): dc='Ee0000'
	if (deformation_code == 10): dc='EEEEEE'
	if (deformation_code == 11): dc='E00E00'
	if (deformation_code == 12): dc='E000E0'
	if (deformation_code == 13): dc='E0000E'
	if (deformation_code == 14): dc='EE0E00'
	
	#-------------------------------------------------------------------------------
	file1 = open("CONTCAR",'r')
	line1 = file1.readlines()		
	file1.close()
	for i in line1:
		if ("Direct" or "direct" or "d" or "D") in i:
			PP=line1.index(i)
	#-------------------------------------------------------------------------------
	
	input_obj 	= open("CONTCAR","r")
	
	firstline   = input_obj.readline() # IGNORE first line comment
	xml_scale 	= float(input_obj.readline())  # scale
	Latvec1 	= input_obj.readline()
	Latvec2 	= input_obj.readline()            
	Latvec3 	= input_obj.readline()            
	elementtype	= input_obj.readline().split()
	if (str.isdigit(elementtype[0])):
		sys.exit("VASP 4.X POSCAR detected. Please add the atom types")
	atom_number = input_obj.readline()
	Coordtype	= input_obj.readline()
	nat 		= atom_number.split()
	nat 		= [int(i) for i in nat]
	print ("Number of atoms in the cell {} ".format( sum(nat)) )
	
	input_obj.close()
	#-------------------------------------------------------------------------------
	a=[]; b=[]; c=[];
	Latvec1		= Latvec1.split()
	Latvec2		= Latvec2.split()
	Latvec3		= Latvec3.split()	
	for ai in Latvec1: 	a.append(float(ai))
	for bi in Latvec2: 	b.append(float(bi))
	for ci in Latvec3: 	c.append(float(ci))
	
	xml_basevect = numpy.array([a] + [b] + [c])		
	print ("{}".format(xml_basevect),end="\n" )
	axis_matrix = numpy.array(xml_basevect) 
	determinant = numpy.linalg.det(axis_matrix)
	volume = numpy.abs(determinant*xml_scale**3)
	print("Equilibrium volume of cell {} ".format(volume) )
	#-------------------------------------------------------------------------------
	work_directory = 'workdir'
	if (len(sys.argv) > 1): work_directory = sys.argv[1]
	if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
	os.mkdir(work_directory)
	os.chdir(work_directory)
	
	output_info = open('INFO-elastic-constants',"w")
	
	output_info.write("\nMaximum Lagrangian strain       = {}".format( maximum_strain ))
	output_info.write("\nNumber of strain values         = {}".format(strain_points))
	output_info.write("\nVolume of equilibrium unit cell = {} [A]^3".format(volume))
	output_info.write("\nDeformation code                = {}".format(deformation_code))
	output_info.write("\nDeformation label               = {}".format(dc, "\n"))
	
	output_info.close()
	
	#-------------------------------------------------------------------------------
	
	delta=strain_points-1 ;# print (delta)
	convert=1
	if (strain_points <= 1):
		strain_points=1
		convert=-1
		delta=1
	
	eta_step=2*maximum_strain/delta
	#print(eta_step)
	#-------------------------------------------------------------------------------
	
	for i in range(0,strain_points):
	
	#-------------------------------------------------------------------------------
	
		eta=i*eta_step-maximum_strain*convert
		#print (eta)
		if (i+1 < 10): strainfile = 'strain-0'+str(i+1)
		else: strainfile = 'strain-'+str(i+1)
		output_str = open(strainfile,"w")
		output_str.write( "{:11.8f}".format(eta) )
		output_str.close()
	
		if (abs(eta) < 0): eta=0
		ep=eta 
		if (eta < 0.0): em=abs(eta)
		else: em=-eta
		
	#-------------------------------------------------------------------------------
	
		e=[]
		for j in range(6):
			ev=0
			if  (dc[j:j+1] == 'E' ): ev=ep; #print (ev)
			elif(dc[j:j+1] == 'e' ): ev=em
			elif(dc[j:j+1] == '0' ): ev=0
			else: print ("==> "), dc; sys.exit("ERROR: deformation code not allowed!") 
			e.append(ev) 
		e = numpy.array(e)
		#print(numpy.array(e))	
	#-------------------------------------------------------------------------------
	
		eta_matrix=numpy.mat( [
		[ e[0]  , e[5]/2, e[4]/2], 
		[ e[5]/2, e[1]  , e[3]/2], 
		[ e[4]/2, e[3]/2, e[2]  ] ], dtype=float32 )
		if i == 6:	 # for test purposes
			print (eta_matrix)
		one_matrix=numpy.identity(3)
	
	#-------------------------------------------------------------------------------
					
		norma=1
		inorma=0
		eps_matrix=eta_matrix
	
		if (numpy.linalg.norm(eta_matrix) > 0.7):sys.exit("ERROR: too large deformation!") 
	
		while ( norma > 0.0 ):
			x=eta_matrix - (1/2) * numpy.dot(eps_matrix,eps_matrix)
			norma=numpy.linalg.norm(x-eps_matrix)     
			#print(norma)		
			eps_matrix=x
			inorma=inorma+1
	
		def_matrix=one_matrix+eps_matrix
		#print (eps_matrix)
		#if i == 6:
		#     print ( numpy.dot(def_matrix,axis_matrix) )	
		#     print ( "{} {}" .format(def_matrix, axis_matrix )	)
		#     print ( numpy.transpose(numpy.dot(def_matrix,numpy.transpose(axis_matrix)) ) )	
		new_axis_matrix=numpy.transpose(numpy.dot(def_matrix,numpy.transpose(axis_matrix)))
		nam=numpy.mat( new_axis_matrix, dtype=float32  )
	
	#-------------------------------------------------------------------------------
	
		if (i+1 < 10): 
			outputfile = 'POSCAR-0'+str(i+1)
		else: 
			outputfile = 'POSCAR-'+str(i+1)
		output_obj = open(outputfile,"w")
		output_obj.write(firstline)
		output_obj.write("{:10.8f}".format((xml_scale))+"\n")
	
		for j in range(3):
			output_obj.write("{:22.16f} {:22.16f} {:22.16f}\n".format( (nam[j,0]), (nam[j,1]), (nam[j,2]) )  ) 
			#output_obj.write( str(fmt%nam[j,0])+str(fmt%nam[j,1])+str(fmt%nam[j,2])+"\n" )
			
		for j in elementtype:
			output_obj.write("\t" +  j)
		output_obj.write("\n" )
		
		for j in atom_number:
			output_obj.write( "{}".format(j) )	
	
		for i in range(len(line1)-PP):
			output_obj.write(line1[PP+i] )	
		output_obj.close()

#-------------------------------------------------------------------------------
	os.chdir('../')
	print 

#plt.imshow(numpy.log(numpy.abs(numpy.fft.fftn(nam))**2))
#plt.show()
