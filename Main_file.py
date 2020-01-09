#!/usr/bin/python3

'''
#####---------------------------------------------------------
#####---------------------------------------------------------
#    Credit	: 	Asif Iqbal BHATTI
#    CODE to: 	convert Cell Matrix to Cell Parameters
#    VERSION: 	This script runs with python3 or later
#    FORMAT	:	POSCAR VASP5 format
#    DATE	: 	28/12/2019
#    USAGE	: 	python3 sys.argv[0]
#####---------------------------------------------------------
#####---------------------------------------------------------
'''

if __name__ == "__main__":

	from termcolor import colored
	import multiprocessing as mp
	from elastic_constants import *
	from energy import *
	from main_contcar_poscar import *
	from poscar import *
	from fitting_energy_vs_volume_curve import *
	
def Introduction():
	global message
	message = "              ____| Python script to process various properties |____"
	print(message)

Introduction()
print("Number of processors Detected: ", mp.cpu_count())
print(Back.MAGENTA + '          NB: POSCAR should be in VASP 5 format', end = '\n', flush=True)
print(Style.RESET_ALL)	
print (colored(' -----------------------------------------------------------','red'), end = '\n', flush=True)
print (colored(' -----------------------------------------------------------','red'), end = '\n', flush=True)
print (colored(' -----------------------------------------------------------','red'), end = '\n', flush=True)
print ('>>> USAGE: execute by typing python3 sys.argv[0]')
print ("**** Following are the options ... ")
print ("(1) To execute only POSCAR file (Convert Lattice Matrix to Lattice parameter)")
print ("(2) To execute POSCAR CELL VOLUME DIFFERENCE with final CONTCAR file")
print ("(3) To extract ENERGY from directories")
print ("(4) To extract ELASTIC CONSTANTS from OUTCAR file (IBRION=6,ISIF=3)")
print ("(5) Fitting energy vs volume curve")

print (colored(' -----------------------------------------------------------','red'), end = '\n', flush=True)
print (colored(' -----------------------------------------------------------','red'), end = '\n', flush=True)
print (colored(' -----------------------------------------------------------','red'), end = '\n', flush=True)
option = input("Enter the option as listed above: ")
option = int(option)
if (option == 1):
	poscar()
elif (option == 2):
	VOL_P = main_poscar()
	VOL_C = main_contcar()
	volume_diff(VOL_P, VOL_C)		
elif (option == 3):
	energy_vs_volume()
elif (option == 4):
	print("Reading OUTCAR. OUTCAR should be in the same directory from which this script is run ")		
	pool = mp.Pool(mp.cpu_count())
	elastic_matrix_VASP_STRESS()
	pool.close()
elif (option == 5):
	fitting_energy_vs_volume_curve()
else:
	print ("INVALID OPTION")
