#!/usr/bin/env python3

'''
# | Elastic Properties
# | Code for calculating mechanical properties 
#   from Energy vs Strain relationship
# 
# Equations can be found at Golesorkhtabar, R., Pavone, P., Spitaler, J., 
# Puschnig, P. & Draxl, C. ElaStic: A tool for calculating second-order elastic 
# constants from first principles. Comput. Phys. Commun. 184, 1861–1873 (2013).
  OR http://wolf.ifj.edu.pl/elastic/index.html 
	https://www.materialsproject.org/wiki/index.php/Elasticity_calculations
'''

import numpy as np
import math, scipy, os, sys
from numpy import linalg as LA
import statistics as st

np.set_printoptions(precision=3)
CRED = '\033[91m';CEND = '\033[0m'
CYEL = '\033[33m'; CEND = '\033[0m'
CPIN = '\033[46m';

def mechanical_properties():
	# Elastic constants have been obtained from SOEC approach mentioned in the above 
	# paper. For three distortions of the crystal three constants are extracted: 
	# C11, C44, C12. For cubic system the matrix is symmetric.

	print (CRED + "Values has to be supplied in the script::" + CEND)

	########################## INPUT PARAMETERS assuming cubic #######################
	
	c11=c22=c33=154 ; 
	c44=c55=c66= 50
	B = 123
	c12=c21=c13=c31=c23=c32=(1/6) * (9 * B - 3 * c11)
	#c12=c21=c13=c31=c23=c32=110
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

def elastic_anisotropy(S):
	#Elastic anisotropy of crystals is important since it correlates with the 
	#possibility to induce micro-cracks in materials 
	v=[]
	f = open('CONTCAR','r')
	lines = f.readlines()
	f.close()
	
	s = float(lines[1])
	for i in range(2,5,1):
		l = s*float(lines[i].split()[0]), s*float(lines[i].split()[1]), s*float(lines[i].split()[2])
		v.append( l )
	v = np.mat(v)
	print("Reading lattice vectors into matrix form")
	print(v[:])

  # l1, l2, l3 are the direction cosines
	n1 = np.linalg.norm(v[0])
	n2 = np.linalg.norm(v[1])
	n3 = np.linalg.norm(v[2])
	l1 = math.degrees(math.acos(np.dot(v[0],np.transpose(v[1]))/(n1*n2)))
	l2 = math.degrees(math.acos(np.dot(v[1],np.transpose(v[2]))/(n2*n3)))
	l3 = math.degrees(math.acos(np.dot(v[2],np.transpose(v[0]))/(n3*n1)))
	print ("l's are direction Cosines:: l1={:6.6f} l2={:6.6f} l3={:6.6f}".format(l1, l2, l3) )	
	B = 1/(3*(S[0,0] + 2*S[0,1] ))
	G = 1/( S[3,3] +4*(S[0,0] - S[0,1] -1/(2*S[3,3]) ) * (l1**2 * l2**2 + l2**2 * l3**2 + l1**2 * l3**2) )
	E = 1/( S[0,0] -2*(S[0,0] - S[0,1] -1/(2*S[3,3]) ) * (l1**2 * l2**2 + l2**2 * l3**2 + l1**2 * l3**2) )
	print ("B={:6.6f} G={:6.6f} E={:6.6f}".format(B, G, E) )
	
############	
def local_lattice_distortion_DEF1():
	#print ("The lattice distortion in paracrystals is measured by the lattice distortion parameter g")
	#print (Back.YELLOW + "Wang, S. Atomic structure modeling of multi-principal-element alloys by the principle")
	#print (Back.YELLOW + "of maximum entropy. Entropy 15, 5536–5548 (2013).")
	print (">"*10,"HUME ROTHERY RULE")
	C_i=C=0.2 ; r_avg = 0.0; del_sum=0.0
	elements = ["Nb", "Hf", "Ta", "Ti", "Zr"]
	eta = {
	"Nb" : 1.98,
	"Hf" : 2.08,
	"Ta" : 2.00,
	"Ti" : 1.76,
	"Zr" : 2.06, }
	
	print ("	{element: atomic radius}")
	print (eta)
	
	for i in elements: 
		r_avg = r_avg + C * eta[i] 
	
	for j in elements:
		del_sum = del_sum + C * ( 1 - float(eta[j]) / r_avg )**2
	del_sum = 100 * np.sqrt(del_sum) 	
	print("HEA_atomic_size_mismatch: \u03B4={}".format(del_sum))

############
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
		xp, yp, zp = lines_contcar[lp+1+i].split()
		x = float(x); y = float(y); z = float(z)
		xp = float(xp); yp = float(yp); zp = float(zp)
		temp = temp + np.sqrt( (x-xp)**2 + (y-yp)**2 + (z-zp)**2 )
	temp = temp/sum_atoms
	print("local lattice distortion: \u0394d={}".format(temp))		

if __name__ == "__main__":	

	S = mechanical_properties()
	print ("_"*30,"Elastic Anisotropy Analysis","_"*30)	
	elastic_anisotropy(S)
	
	print ("")
	print ("_"*30,"Lattice Distortion Analysis","_"*30)
	#local_lattice_distortion_DEF1()
	print ("")
	local_lattice_distortion_DEF2()
	


