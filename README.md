# Python3 Script to convert Cell Matrix to Cell parameter
**_Convert VASP5 POSCAR Lattice Matrix to Lattice parameter_**

If you have performed a number of calculations and you want to obtain the lattice parameters and volume from each directory then run this script in a folder and it will loop over all the directory and find POSCAR files and outut the lattice parameters and volume, moreover it also scans the CONTCAR files and compute the volume and lattice parameters, and on the screen it gives the volume difference upon structure minimization. 
For Elastic constants this gives an estimation how much the cell has deformed and gives us the indication whether we need to increase ENCUT or KPOINTS to minimize the Pulay stress as indicated on VASP website.

 For VASP structural minimization these tag should be used: IBRION = 2; ISIF = 3; EDIFF= 10**-8

```
~/dir/ -->
dir1
dir2
dir3
...
Python3 cellMatrix_cellParameters.py
```

```
Following script calculates Energies, volumetric cell, elastic Matrix from OUTCAR file.
```
```
An excercise to print Cij matrix in a pythonic way

def print_Cij_Matrix():
	Bij = []
	C = "C"
	for i in range(0,6):
		Bij.append([])
		for j in range(0,6):
			Bij[i].append((C + str(i) + str(j)))
	l = np.matrix(Bij)		
	print (l)
 ```




