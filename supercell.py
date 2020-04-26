#!/usr/bin/env python3
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
import numpy as np
import sys, random
from numpy import linalg as LA

#------------------------------INPUT PARAMETRS------------------------------------#
lattice_parameter = 3.405 
supercellx = 5
supercelly = 5
supercellz = 5
#----------------------------------------------------------------------------------#
def HEAs_supercell():
	cartesian_units=[]; reduced_units=[]; randStruct=[]; count=0
	lattice_vector = np.array([[1,0,0],
														[0,1,0],
														[0,0,1]])*lattice_parameter
	
	lattice_bcc = np.array([[0,0,0],
												[0.5,0.5,0.5]])*lattice_parameter
						
	lattice_fcc = np.array([[0,0,0],
												[0.0,0.5,0.5],											
												[0.5,0.0,0.5],											
												[0.5,0.5,0.0]])*lattice_parameter										
	fcc = lattice_fcc
	bcc = lattice_bcc		
		
	for i in range(supercellx):										
		for j in range(supercelly):
			for k in range(supercellz):
				atom_position = np.array([i,j,k])
				cartesian_basis = np.inner(lattice_vector.T, atom_position)
				#print(cartesian_basis)
				if (sys.argv[1] == 'bcc'):
					for atom in lattice_bcc:
						cartesian_units.append(cartesian_basis + atom)
						count+=1
				elif (sys.argv[1] == 'fcc'):		
					for atom in lattice_fcc:
						cartesian_units.append(cartesian_basis + atom)
						count+=1			
	
	with open("POSCAR","w") as POSCAR:
		POSCAR.write('#{}\n'.format(sys.argv[1]))
		POSCAR.write('{:6.6f}\n'.format(1.0))
		POSCAR.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[0][0]*supercellx,lattice_vector[0][1]*supercellx,lattice_vector[0][2]*supercellx ))
		POSCAR.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[1][0]*supercelly,lattice_vector[1][1]*supercelly,lattice_vector[1][2]*supercelly ))
		POSCAR.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[2][0]*supercellz,lattice_vector[2][1]*supercellz,lattice_vector[2][2]*supercellz ))		
		POSCAR.write('Ta\n')
		POSCAR.write('{}\n'.format((count)))
	
		if (sys.argv[2] == "c" or sys.argv[2] == "C"):
			POSCAR.write("Cartesian\n")	
			for i in range(len(cartesian_units) ):
			#print("{:12.9f} {:12.9f} {:12.9f}".format(cartesian_units[i][0], cartesian_units[i][1], cartesian_units[i][2]))
				POSCAR.write("{:12.9f} {:12.9f} {:12.9f}\n".format(cartesian_units[i][0], cartesian_units[i][1], cartesian_units[i][2]))
	
	#------------------------- In fractional/reduced UNITS -------------------------#		
		u = np.cross(lattice_vector[1]*supercelly, lattice_vector[2]*supercellz)
		v = np.cross(lattice_vector[0]*supercellx, lattice_vector[2]*supercellz)
		w = np.cross(lattice_vector[0]*supercellx, lattice_vector[1]*supercelly)
		#V = np.array([ lattice_vector[0]*supercellx,lattice_vector[1]*supercelly,lattice_vector[2]*supercellz ] )
		#print ("Volume of the cell::", LA.det(V) )
		Vx = np.inner(lattice_vector[0]*supercellx,u)
		Vy = np.inner(lattice_vector[1]*supercelly,v)
		Vz = np.inner(lattice_vector[2]*supercellz,w)
		
		if (sys.argv[2] == "d" or sys.argv[2] == "D"):
			POSCAR.write("Direct\n")		
			for i in range(len(cartesian_units) ):
				POSCAR.write("{:12.9f} {:12.9f} {:12.9f}\n".format(np.dot(cartesian_units[i],u)/Vx, np.dot(cartesian_units[i],v)/Vy, np.dot(cartesian_units[i],w)/Vz ) )		
	
	#------------------------- Randomly distribute atoms for HEA -------------------------#
	with open("newPOSCAR","w") as fdata:
		for i in range(int(1e3)):
			randomArrx=random.sample(range(0,count),count)
		if (count%5 == 0): 
			#print(int(count/5) )		
			fdata.write('#{}\n'.format(sys.argv[1]))
			fdata.write('{:6.6f}\n'.format(1.0))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[0][0]*supercellx,lattice_vector[0][1]*supercellx,lattice_vector[0][2]*supercellx ))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[1][0]*supercelly,lattice_vector[1][1]*supercelly,lattice_vector[1][2]*supercelly ))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[2][0]*supercellz,lattice_vector[2][1]*supercellz,lattice_vector[2][2]*supercellz ))		
			fdata.write('Ta Hf Zr Nb Ti\n')
			fdata.write('{0} {0} {0} {0} {0}\n'.format((int(count/5))) )
					
			if (sys.argv[2] == "c" or sys.argv[2] == "C"):
				fdata.write("Cartesian\n")	
				for i in range(len(cartesian_units) ):
					#print("{:12.9f} {:12.9f} {:12.9f}".format(cartesian_units[randomArrx[i]][0],cartesian_units[randomArry[i]][1],cartesian_units[randomArrz[i]][2] ) )
					fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(cartesian_units[randomArrx[i]][0],cartesian_units[randomArrx[i]][1],cartesian_units[randomArrx[i]][2] ) )
	
			if (sys.argv[2] == "d" or sys.argv[2] == "D"):
				fdata.write("Direct\n")		
				for i in range(len(cartesian_units) ):
					fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(np.dot(cartesian_units[randomArrx[i]],u)/Vx, np.dot(cartesian_units[randomArrx[i]],v)/Vy, np.dot(cartesian_units[randomArrx[i]],w)/Vz ) )
				

def help():
	print('A simple script to generate FCC or BCC supercell for HEAs.')
	print('To execute just run python3 <bcc/fcc> <c/d>.')
	print('HEAs consists of five or more elements. The elemnts has already been typed into the')
	print('script just change according to your needs also lattice vectors should be')
	print('equal and atomic composition should corresponds to integer multiple of atoms.')
	
if __name__ == '__main__':
	if len(sys.argv) < 3 or len(sys.argv) > 3:
		help()
	else:
		HEAs_supercell()
			
			
			
			
			
			
			
			
			
