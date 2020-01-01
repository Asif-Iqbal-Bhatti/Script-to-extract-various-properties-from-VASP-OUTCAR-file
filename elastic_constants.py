import os, sys, spglib
import math, glob
import numpy as np
import subprocess
from os import listdir
from os.path import isfile, join
from pathlib import Path

ang2atomic = 1.889725988579 # 1 A = 1.889725988579 [a.u]
ang2bohr   = 6.7483330371   # 1 A^3 = 6.7483330371 [a.u]^3

def print_Cij_Matrix():
	Bij = []
	C = "C"
	for i in range(0,6):
		Bij.append([])
		for j in range(0,6):
			Bij[i].append((C + str(i) + str(j)))
	l = np.matrix(Bij)
	return l
	
def elastic_matrix():
	while True:
			s=np.zeros((6,6))
			c=np.zeros((6,6))
			
			file = open("OUTCAR",'r')
			lines = file.readlines()		
			file.close()
			
			for i in lines:
				if "TOTAL ELASTIC MODULI (kBar)" in i:
					ll=lines.index(i)
				if "LATTYP" in i:
					crystaltype=str(i.split()[3])
			print ("DETECTED CRYSTAL FROM OUTCAR:", crystaltype)			
			print (" ")			
			for i in range(0,6):
				l=lines[ll+3+i] # indexing a line in huge file
				word = l.split()
				s[i][:] = word[1:7]
				for j in range(0,6):
					c[i][j] = float(s[i][j])/10.0
					
##########--------------------stress tensor------------------	
				
			Cij = np.matrix(c)
			np.set_printoptions(precision=4, suppress=True)
			print (Cij)
			print(print_Cij_Matrix())
			print ("\nEigen Values of the matrix Cij:")
			evals = np.linalg.eigvals(Cij)
			if evals.all() > 0:
				print(evals)
				print("All Eigen values are positive")
			
##########--------------------Compliance tensor------------------
##########--------------------  s_{ij} = C_{ij}^{-1}

			Sij = np.linalg.inv(Cij)
			
#----------------------------- ELASTIC PROPERTIES -----------------------------------

			stability_test(Cij, crystaltype)
			
#-------------------------------- Voigt bulk modulus  K_v  $(GPa)$---------------
#------------------- 9K_v = (C_{11}+C_{22}+C_{33}) + 2(C_{12} + C_{23} + C_{31}) 
			Kv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) + 2 * (Cij[0,1] + Cij[1,2] + Cij[2,0]))/9.0
#-------------------------------- Reuss shear modulus  G_v  $(GPa)$------------------
#------------------- 15/G_R = 4(s_{11}+s_{22}+s_{33}) - 4(s_{12} + s_{23} + s_{31}) + 3(s_{44} + s_{55} + s_{66})$
			Gv = ((Cij[0,0] + Cij[1,1] + Cij[2,2]) - (Cij[0,1] + Cij[1,2] + Cij[2,0]) + 3 * (Cij[3,3] + Cij[4,4] + Cij[5,5]))/15.0
#-------------------------------- Reuss bulk modulus  K_r  $(GPa)$----------------
#-------------------  1/K_R = (s_{11}+s_{22}+s_{33}) + 2(s_{12} + s_{23} + s_{31})$
			Kr = 1/((Sij[0,0] + Sij[1,1] + Sij[2,2]) + 2 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) )
#-------------------------------- Reuss shear modulus  G_r  $(GPa)$------------------
			Gr = 15/(4 * (Sij[0,0] + Sij[1,1] + Sij[2,2]) - 4 * (Sij[0,1] + Sij[1,2] + Sij[2,0]) + 3 * (Sij[3,3] + Sij[4,4] + Sij[5,5]))
		
#------------------------------------------------------------------------------------

			## Young's Modulus "E": Voigt
			Ev = (9*Kv*Gv)/(3*Kv + Gv)				
			## Young's: Reuss
			Er = (9*Kr*Gr)/(3*Kr + Gr)		
			
			## Poisson's ratio: Voigt
			Nu_V = (3*Kv - Ev)/(6*Kv)				
	        ## Poisson's ratio: Reuss
			Nu_R = (3*Kr - Er)/(6*Kr)
			
			## P-wave modulus, M: Voigt
			MV = Kv + (4*Gv/3.0)	        
	        ## P-wave modulus, M: Reuss
			MR = Kr + (4*Gr/3.0)	
				
#------------------------------------------------------------------------------------
#--------------- Voigt-Reuss-Hill Approximation: average of both methods
			Kvrh = (Kv + Kr)/2.0
			Gvrh = (Gv + Gr)/2.0
			Mvrh = (MV + MR)/2.0
			Evrh = (Ev + Er)/2.0
			Nu_vrh = (Nu_V + Nu_R)/2.0
			KG_ratio_V = Kv/Gv
			KG_ratio_R = Kr/Gr
			KG_ratio_vrh = Kvrh/Gvrh			
#--------------- Isotropic Poisson ratio $\mu 
#--------------- $\mu = (3K_{vrh} - 2G_{vrh})/(6K_{vrh} + 2G_{vrh})$
			mu = (3 * Kvrh - 2 * Gvrh) / (6 * Kvrh + 2 * Gvrh )
			
#-----------------------------------------------------------------------
#------------------------------------------------------------------------------------
		
			print ("\n \n                         Voigt     Reuss    Average")
			print ("-------------------------------------------------------")
			print ("Bulk modulus   (GPa)  %9.3f %9.3f %9.3f " % (Kv, Kr, Kvrh))
			print ("Shear modulus  (GPa)  %9.3f %9.3f %9.3f " % (Gv, Gr, Gvrh))
			print ("Young modulus  (GPa)  %9.3f %9.3f %9.3f " % (Ev, Er, Evrh))
			print ("Poisson ratio         %9.3f %9.3f %9.3f " % (Nu_V, Nu_R, Nu_vrh))
			print ("P-wave modulus  (GPa) %9.3f %9.3f %9.3f " % (MV, MR, Mvrh))
			print ("Bulk/Shear ratio      %9.3f %9.3f %9.3f (%s) " %(KG_ratio_V, KG_ratio_R, KG_ratio_vrh,  ductile_test(KG_ratio_vrh) ))
			print ("-------------------------------------------------------")
			print("Isotropic Poisson ratio: ", mu)			
			break
	
			
def ductile_test(ratio):
	if(ratio > 1.75):
		return "ductile"
	else:
		return "brittle"
		
	
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