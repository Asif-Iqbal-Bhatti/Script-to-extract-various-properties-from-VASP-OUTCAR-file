# Code to extract Elastic properties, Energy, lattice parameter from VASP output files 💫
**_Python3 or later: High Entropy Alloys analysis_**

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
![pypi](https://img.shields.io/pypi/v/pybadges.svg)
![versions](https://img.shields.io/pypi/pyversions/pybadges.svg)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)


If you have performed a number of calculations and want to obtain the lattice parameters, energy and volume from each directory (in a given directory) then run this script and it will loop over all the directories and find POSCAR file and output the lattice parameters, energy and volume, moreover it also scans the CONTCAR file and compute the volume and lattice parameters, and print on the screen volume difference upon structure minimization. 
For Elastic constants this gives an indication how much the cell has deformed and gives us the indication whether we need to increase ENCUT or KPOINTS to minimize the Pulay stress as indicated on VASP manual.

```
**USAGE** : Just run the Main_file.py it will call all the modules
The code has been given in one script (**Extract_elastic_energy_VASP.py**) and in modules. 
The module form is much easier to handle than one long complex code.
~/dir/ -->
dir1
dir2
dir3
POSCAR
OUTCAR
Python Main_file.py
```
_______________________
_For VASP structural minimization these tag should be used: IBRION = 2; ISIF = 3; EDIFF= 10**-8. VASP has implemented the stress method (First derivative approach) and with few deformation it gives the Stiffness Constants directly. However as the system size gets larger then it becomes very expensive to compute. The other appraoch is the Energy-strain method where Langragian strain are applied to the cell and the resulting energy is calculated for a number of deformations. After this polynomial fitting is done and second derivative is calculated at equilibrium volume._

```
An excercise to print Cij matrix in a pythonic way: 2D array is created as List of Lists unlike C++ where you define by C[i][j]

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

ref: [exciting](http://exciting-code.org/nitrogen-energy-vs-strain-calculations)
ref: [materialsproject](https://wiki.materialsproject.org/Elasticity_calculations)

[![licensebuttons by-nd](https://licensebuttons.net/l/by-nd/3.0/88x31.png)](https://creativecommons.org/licenses/by-nd/4.0)

