#!/usr/bin/env python3
#_____________________________________________________________________________
'''
USAGE: 	python3 script to create directories from POSCAR files generated from
		PHONOPY code.
Credit: Asif Iqbal Bhatti
DATED:	19-01-2020

'''
#_____________________________________________________________________________
import 	os, sys, scipy
import  pathlib, subprocess, shutil
from 	os import listdir
from 	os.path import isfile, join
import 	os.path
import 	multiprocessing as mp
import 	fnmatch
import 	time

#_____________________________________________________________________________
os.system('rm -r disp-*')
#os.system('phonopy -d --dim="2 2 2"') ## It will call a phonopy code
#subprocess.call(['phonopy','-v', '-d', '--dim=3 3 3'],shell = False) ## It will call a phonopy code
#_____________________________________________________________________________


def copy_fntn(k):
	for i in range(1, int(k)+1):
		shutil.rmtree(f'disp-{str(i).zfill(3)}', ignore_errors=True)
		os.mkdir(f'disp-{str(i).zfill(3)}')
		#subprocess.check_call(['mkdir', 'disp-'+str(i).zfill(3)])
		#subprocess.check_call('cd disp-'+str(i).zfill(3), shell=True)	
		# print(str(i).zfill(3))
		#
		# The zfill()/rjust() is a function associated with the string object. 
		# 3 is the expected length of the string after padding.
		#
		os.system(f'cp INCAR   disp-{str(i).zfill(3)}')
		os.system(f'cp POTCAR  disp-{str(i).zfill(3)}')
		os.system(f'cp KPOINTS disp-{str(i).zfill(3)}')
		#os.system('cp job.sh  disp-'+str(i).zfill(3) )
		#os.system('cp WAVECAR  disp-'+str(i).zfill(3) )
		subprocess.call(
			['cp', '-r', f'POSCAR-{str(i).zfill(3)}', f'disp-{str(i).zfill(3)}'],
			shell=False,
		)

		os.chdir(f'disp-{str(i).zfill(3)}')
		shutil.copyfile(f"POSCAR-{str(i).zfill(3)}", "POSCAR")
		#os.system('ls')
		os.chdir('../')


def create_phonon_directories():
	print("Number of processors Detected: ", mp.cpu_count())
	pool = mp.Pool(4)
	print("*"*80)
	if not os.path.exists('POSCAR'):
		print (' ERROR: POSCAR does not exist here.')
		sys.exit(0)
	else:
		mypath = os.getcwd()
		counter = sum(
			1
			for entry in os.listdir(mypath)
			if os.path.isfile(os.path.join(mypath, entry))
			and fnmatch.fnmatchcase(entry, 'POSCAR-*')
		)

		print(f"# of POSCAR-* files generated with Phonopy code: ----> {counter}")
		print("*"*80)	

		start_time = time.time()
		#copy_fntn(counter)
		#p = mp.Process(target=copy_fntn, args=('9',))
		#p.start()
		#p.join()
		with pool as p:
			p.map(copy_fntn, [counter])
	end_time = time.time()
	print("total time taken for this loop: ", end_time - start_time)
	
if __name__ == "__main__":	
		create_phonon_directories()
		
		
