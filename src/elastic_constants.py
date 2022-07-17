import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path

ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
ang2bohr   = 6.7483330371   # 1 A^3 = 6.7483330371 [a.u]^3

class Elastic_Matrix:
	def print_Cij_Matrix(): ###EXERCISE
		Bij = []
		C = "C"
		for i in range(6): 
			Bij.append([])
			for j in range(1, 7): 
				Bij[i].append((C + str(i+1) + str(j)))
		return np.matrix(Bij)

	def langragian_strain():
		kk = ('x', 'y', 'z')
		C = "E"
		Eij=[]
		eta = {
		"xx" : 11,
		"yy" : 22,
		"zz" : 33,
		"yz" : 23,
		"xz" : 13,
		"xy" : 12  }
		print ( "The Langragian strain in Voigt notation: ")
		print ( eta.items())

		for n in kk: 
			Eij.extend(C + str(n) + str(m) for m in kk)
		return np.matrix(Eij).reshape(3,3)


def elastic_matrix_VASP_STRESS():
	while True:	
		if not os.path.exists('OUTCAR'):
			print (' ERROR: OUTCAR does not exist here.')
			sys.exit(0)

		s=np.zeros((6,6))
		c=np.zeros((6,6))
		with open("OUTCAR",'r') as file:
			lines = file.readlines()
		for i in lines:
			if "TOTAL ELASTIC MODULI (kBar)" in i:
				ll=lines.index(i)
			if "LATTYP" in i:
				crystaltype=str(i.split()[3])
		print ("DETECTED CRYSTAL FROM OUTCAR:", crystaltype)
		print (" ")
		for i in range(6):
			l=lines[ll+3+i] # indexing a line in huge file
			word = l.split()
			s[i][:] = word[1:7]
			for j in range(6):
				c[i][j] = float(s[i][j])/10.0
##########-------------------stress tensor------------------	

		Cij = np.matrix(c)
		np.set_printoptions(precision=4, suppress=True)
		print (Cij)
		print(Elastic_Matrix.print_Cij_Matrix()  )
		print(Elastic_Matrix.langragian_strain() )
		print ("\nEigen Values of the matrix Cij:")
		evals = np.linalg.eigvals(Cij)
		if evals.all() > 0:
			print(evals)
			print("All Eigen values are positive")
##########-------------------Compliance tensor------------------
##########-------------------  s_{ij} = C_{ij}^{-1}

		Sij = np.linalg.inv(Cij)
#---------------------------- ELASTIC PROPERTIES -----------------------------------

		stability_test(Cij, crystaltype)
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
	print ("{:_^80}".format("GPa"))
	print ("{:25.8s}  {:15.8s} {:15.8s} {:15.8s}".format(" ","Voigt", "Reuss ", "Hill") )
	print ("{:16.20s} {:15.3f} {:15.3f} {:15.3f}".format("Bulk Modulus",Bv, Br, B_H) )
	print ("{:16.20s} {:15.3f} {:15.3f} {:15.3f}".format("Shear Modulus",Gv, Gr, G_H) )
	print ("{:16.20s} {:15.3f} {:15.3f} {:15.3f}".format("Young Modulus",Ev, Er, E_H) )
	print ("{:-^80}".format("-"))
	print ("{:16.20s} {:15.3f} {:15.3f} {:15.3f}".format("Poisson ratio ", NuV, NuR, nu_H) )
	print ("{:16.20s} {:15.3f} {:15.3f} {:15.3f}({:5.3f})".format("B/G ratio ",Bv/Gv,Br/Gr, B_H/G_H, G_H/B_H) )
	print ("{:16.20s} {:15.3s} {:15.3s} {:15.3f}".format("Zener ratio Az",'','', A) )
	print ("{:16.20s} {:15.3s} {:15.3s} {:15.3f}".format("Avr ratio ",'','', (Gv-Gr)/(Gv+Gr)) )
	print ("{:16.20s} {:15.3s} {:15.3s} {:15.3f}".format("AU ",'','', AU) )
	print ("{:16.20s} {:15.3s} {:15.3s} {:15.3f}".format("Cauchy pressure ",'','', (c12-c44)) )
	print ("{:16.20s} {:15.3s} {:15.3s} {:15.3f}".format("C'tetra Shear ",'','',  C) )

	print ("-"*80)
	return Sij
	
		
		
def ductile_test(ratio):
	return "ductile" if (ratio > 1.75) else "brittle"
###
	
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
###

def born_stability_criterion():
	print ("Born stability criteria for the stability of following systems \n")
	print ("Cubic crystal system.... \n")	
	print ("(i) C11 - C12 > 0;    (ii) C11 + 2C12 > 0;   (iii) C44 > 0 \n ")
	print ("Hexagonal crystal system.... \n")	
	print ("(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0 \n ")
	print ("Tetragonal crystal system.... \n")	
	print ("(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0;   (iv) C66 > 0;    (v) 2C16^2 < C66*(C11-C12) \n ")
	print ("rhombohedral crystal system.... \n")	
	print ("(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0; \n ")
	print ("orthorhombic crystal system.... \n")	
	print ("(i) C11 > 0;   (ii) C11*C22 > C12^2;   (iii) C11*C22*C33 + 2C12*C13*C23 - C11*C23^2 - C22*C13^2 - C33*C12^2 > 0;   (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0 \n")
	print ("Monoclinic crystal system.... \n")
	print ("[Ref- Mouhat and Coudert, PRB 90, 224104 (2014), and Wu et al. PRB 76, 054115 (2007)]  \n")
	print ("(i) C11 > 0;  (ii)  C22 > 0; (iii)  C33 > 0; (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0  ")
	print ("(vii) [C11 + C22 + C33 + 2*(C12 + C13 + C23)] > 0;    (viii)  C33*C55 - C35^2 > 0;   (ix)  C44*C66 - C46^2 > 0;   (x) C22 + C33 - 2*C23  > 0 ")
	print ("(xi) C22*(C33*C55 - C35^2) + 2*C23*C25*C35 - (C23^2)*C55 - (C25^2)*C33   > 0  ")
	print ("(xii)  2*[C15*C25*(C33*C12 - C13*C23) + C15*C35*(C22*C13 - C12*C23) + C25*C35*(C11*C23 - C12*C13)] - [C15*C15*(C22*C33 - C23^2) + C25*C25*(C11*C33 - C13^2) + C35*C35*(C11*C22 - C12^2)] + C55*g > 0  ")
	print (" where, g = [C11*C22*C33 - C11*C23*C23 - C22*C13*C13 - C33*C12*C12 + 2*C12*C13*C23 ] "	)	
###
