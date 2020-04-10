# Code to extract Elastic properties, Energy, lattice parameter from VASP output files 💫
**_Python3 or later: High Entropy Alloys analysis_**

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
![pypi](https://img.shields.io/pypi/v/pybadges.svg)
![versions](https://img.shields.io/pypi/pyversions/pybadges.svg)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

==> supercell.py script construct bcc/fcc supercell for High Entropy Alloys <==
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

𝗢𝗻𝗲 𝗳𝗶𝗲𝗹𝗱 𝗼𝗳 𝘄𝗼𝗿𝗸 𝗶𝗻 𝘄𝗵𝗶𝗰𝗵 𝘁𝗵𝗲𝗿𝗲 𝗵𝗮𝘀 𝗯𝗲𝗲𝗻 𝘁𝗼𝗼 𝗺𝘂𝗰𝗵 𝘀𝗽𝗲𝗰𝘂𝗹𝗮𝘁𝗶𝗼𝗻 𝗶𝘀 𝗰𝗼𝘀𝗺𝗼𝗹𝗼𝗴𝘆. 𝗧𝗵𝗲𝗿𝗲 𝗮𝗿𝗲 𝘃𝗲𝗿𝘆 𝗳𝗲𝘄 𝗵𝗮𝗿𝗱 𝗳𝗮𝗰𝘁𝘀 𝘁𝗼 𝗴𝗼 𝗼𝗻, 𝗯𝘂𝘁 𝘁𝗵𝗲𝗼𝗿𝗲𝘁𝗶𝗰𝗮𝗹 𝘄𝗼𝗿𝗸𝗲𝗿𝘀 𝗵𝗮𝘃𝗲 𝗯𝗲𝗲𝗻 𝗯𝘂𝘀𝘆 𝗰𝗼𝗻𝘀𝘁𝗿𝘂𝗰𝘁𝗶𝗻𝗴 𝘃𝗮𝗿𝗶𝗼𝘂𝘀 𝗺𝗼𝗱𝗲𝗹𝘀 𝗳𝗼𝗿 𝘁𝗵𝗲 𝗨𝗻𝗶𝘃𝗲𝗿𝘀𝗲, 𝗯𝗮𝘀𝗲𝗱 𝗼𝗻 𝗮𝗻𝘆 𝗮𝘀𝘀𝘂𝗺𝗽𝘁𝗶𝗼𝗻𝘀 𝘁𝗵𝗮𝘁 𝘁𝗵𝗲𝘆 𝗳𝗮𝗻𝗰𝘆. 𝗧𝗵𝗲𝘀𝗲 𝗺𝗼𝗱𝗲𝗹𝘀 𝗮𝗿𝗲 𝗽𝗿𝗼𝗯𝗮𝗯𝗹𝘆 𝗮𝗹𝗹 𝘄𝗿𝗼𝗻𝗴. 𝗜𝘁 𝗶𝘀 𝘂𝘀𝘂𝗮𝗹𝗹𝘆 𝗮𝘀𝘀𝘂𝗺𝗲𝗱 𝘁𝗵𝗮𝘁 𝘁𝗵𝗲 𝗹𝗮𝘄𝘀 𝗼𝗳 𝗻𝗮𝘁𝘂𝗿𝗲 𝗵𝗮𝘃𝗲 𝗮𝗹𝘄𝗮𝘆𝘀 𝗯𝗲𝗲𝗻 𝘁𝗵𝗲 𝘀𝗮𝗺𝗲 𝗮𝘀 𝘁𝗵𝗲𝘆 𝗮𝗿𝗲 𝗻𝗼𝘄. 𝗧𝗵𝗲𝗿𝗲 𝗶𝘀 𝗻𝗼 𝗷𝘂𝘀𝘁𝗶𝗳𝗶𝗰𝗮𝘁𝗶𝗼𝗻 𝗳𝗼𝗿 𝘁𝗵𝗶𝘀. 𝗧𝗵𝗲 𝗹𝗮𝘄𝘀 𝗺𝗮𝘆 𝗯𝗲 𝗰𝗵𝗮𝗻𝗴𝗶𝗻𝗴, 𝗮𝗻𝗱 𝗶𝗻 𝗽𝗮𝗿𝘁𝗶𝗰𝘂𝗹𝗮𝗿, 𝗾𝘂𝗮𝗻𝘁𝗶𝘁𝗶𝗲𝘀 𝘁𝗵𝗮𝘁 𝗮𝗿𝗲 𝗰𝗼𝗻𝘀𝗶𝗱𝗲𝗿𝗲𝗱 𝘁𝗼 𝗯𝗲 𝗰𝗼𝗻𝘀𝘁𝗮𝗻𝘁𝘀 𝗼𝗳 𝗻𝗮𝘁𝘂𝗿𝗲 𝗺𝗮𝘆 𝗯𝗲 𝘃𝗮𝗿𝘆𝗶𝗻𝗴 𝘄𝗶𝘁𝗵 𝗰𝗼𝘀𝗺𝗼𝗹𝗼𝗴𝗶𝗰𝗮𝗹 𝘁𝗶𝗺𝗲. 𝗦𝘂𝗰𝗵 𝘃𝗮𝗿𝗶𝗮𝘁𝗶𝗼𝗻𝘀 𝘄𝗼𝘂𝗹𝗱 𝗰𝗼𝗺𝗽𝗹𝗲𝘁𝗲𝗹𝘆 𝘂𝗽𝘀𝗲𝘁 𝘁𝗵𝗲 𝗺𝗼𝗱𝗲𝗹 𝗺𝗮𝗸𝗲𝗿𝘀." 𝗗𝗶𝗿𝗮𝗰, 𝗣𝗮𝘂𝗹. 𝗢𝗻 𝗺𝗲𝘁𝗵𝗼𝗱𝘀 𝗶𝗻 𝘁𝗵𝗲𝗼𝗿𝗲𝘁𝗶𝗰𝗮𝗹 𝗽𝗵𝘆𝘀𝗶𝗰𝘀. (𝗧𝗿𝗶𝗲𝘀𝘁𝗲. 𝗝𝘂𝗻𝗲 𝟭𝟵𝟲𝟴 .) 

