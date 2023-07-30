# Code to extract Elastic properties, Energy, and lattice parameter from VASP output files

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
![pypi](https://img.shields.io/pypi/v/pybadges.svg)
![versions](https://img.shields.io/pypi/pyversions/pybadges.svg)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

==> supercell.py script constructs a bcc/fcc supercell for high-entropy alloys.

This repository contains a set of Python scripts for analyzing High Entropy Alloys using VASP output files. The scripts can extract lattice parameters, energy, and volume from each directory in a given parent directory. Additionally, it can scan the CONTCAR file, compute the volume and lattice parameters, and print the volume difference upon structure minimization. The extracted data can be used to analyze the Elastic Constants and make decisions on adjusting ENCUT or KPOINTS to minimize Pulay stress, as indicated in the VASP manual.

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
## Structural Minimization in VASP:  
For VASP structural minimization, the following tags should be used: IBRION = 2; ISIF = 3; EDIFF = 10**-8. VASP implements both the stress method (First derivative approach) and the Energy-strain method for Elastic Constants calculation. The Energy-strain method involves applying Lagrangian strains to the cell and calculating the resulting energy for various deformations. Polynomial fitting is then performed, and the second derivative is calculated at equilibrium volume.

```
An exercise to print Cij matrix in a Pythonic way: 2D array is created as a List of Lists unlike C++ where you define by C[i][j]

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

**Parsing Hessian matrix from a file**

Script to parse Hessian Matrix from a file

USAGE: To parse the VASP vasprun.xml input file to extract the Hessian Matrix Output file has already been defined in the code (hessian.dat).

CAUTION: Use at your own risk (NOTEVEN IMPLIED GUARANTEED, WHATSOEVER), the code has been tested but the user in the end will have to verify the ouput.

Xpath is a useful tool for directly accessing the element in a Node

-> CHECK this website for lxml introduccion: https://lxml.de/tutorial.html & -> https://github.com/lxml/lxml

NB: The Matrix is obtained from VASP, vasprun.xml file. It contains the hessian matrix that can be edited with the appropriate script.

The complication that arises by parsing the data from the file is the trailing empty lines and tabs that need to be deleted before it can be read. In the case of a simple file format, there is no need for it. But if the file contains irregular data entry such as empty lines with commas and spaces then it needs to be formatted.

In my case, I have only leading empty lines and spaces and in between columns there are white spaces of different sizes.

**Convert exciting(LMTO) input file to VASP input file (DFT code)**


𝗢𝗻𝗲 𝗳𝗶𝗲𝗹𝗱 𝗼𝗳 𝘄𝗼𝗿𝗸 𝗶𝗻 𝘄𝗵𝗶𝗰𝗵 𝘁𝗵𝗲𝗿𝗲 𝗵𝗮𝘀 𝗯𝗲𝗲𝗻 𝘁𝗼𝗼 𝗺𝘂𝗰𝗵 𝘀𝗽𝗲𝗰𝘂𝗹𝗮𝘁𝗶𝗼𝗻 𝗶𝘀 𝗰𝗼𝘀𝗺𝗼𝗹𝗼𝗴𝘆. 𝗧𝗵𝗲𝗿𝗲 𝗮𝗿𝗲 𝘃𝗲𝗿𝘆 𝗳𝗲𝘄 𝗵𝗮𝗿𝗱 𝗳𝗮𝗰𝘁𝘀 𝘁𝗼 𝗴𝗼 𝗼𝗻, 𝗯𝘂𝘁 𝘁𝗵𝗲𝗼𝗿𝗲𝘁𝗶𝗰𝗮𝗹 𝘄𝗼𝗿𝗸𝗲𝗿𝘀 𝗵𝗮𝘃𝗲 𝗯𝗲𝗲𝗻 𝗯𝘂𝘀𝘆 𝗰𝗼𝗻𝘀𝘁𝗿𝘂𝗰𝘁𝗶𝗻𝗴 𝘃𝗮𝗿𝗶𝗼𝘂𝘀 𝗺𝗼𝗱𝗲𝗹𝘀 𝗳𝗼𝗿 𝘁𝗵𝗲 𝗨𝗻𝗶𝘃𝗲𝗿𝘀𝗲, 𝗯𝗮𝘀𝗲𝗱 𝗼𝗻 𝗮𝗻𝘆 𝗮𝘀𝘀𝘂𝗺𝗽𝘁𝗶𝗼𝗻𝘀 𝘁𝗵𝗮𝘁 𝘁𝗵𝗲𝘆 𝗳𝗮𝗻𝗰𝘆. 𝗧𝗵𝗲𝘀𝗲 𝗺𝗼𝗱𝗲𝗹𝘀 𝗮𝗿𝗲 𝗽𝗿𝗼𝗯𝗮𝗯𝗹𝘆 𝗮𝗹𝗹 𝘄𝗿𝗼𝗻𝗴. 𝗜𝘁 𝗶𝘀 𝘂𝘀𝘂𝗮𝗹𝗹𝘆 𝗮𝘀𝘀𝘂𝗺𝗲𝗱 𝘁𝗵𝗮𝘁 𝘁𝗵𝗲 𝗹𝗮𝘄𝘀 𝗼𝗳 𝗻𝗮𝘁𝘂𝗿𝗲 𝗵𝗮𝘃𝗲 𝗮𝗹𝘄𝗮𝘆𝘀 𝗯𝗲𝗲𝗻 𝘁𝗵𝗲 𝘀𝗮𝗺𝗲 𝗮𝘀 𝘁𝗵𝗲𝘆 𝗮𝗿𝗲 𝗻𝗼𝘄. 𝗧𝗵𝗲𝗿𝗲 𝗶𝘀 𝗻𝗼 𝗷𝘂𝘀𝘁𝗶𝗳𝗶𝗰𝗮𝘁𝗶𝗼𝗻 𝗳𝗼𝗿 𝘁𝗵𝗶𝘀. 𝗧𝗵𝗲 𝗹𝗮𝘄𝘀 𝗺𝗮𝘆 𝗯𝗲 𝗰𝗵𝗮𝗻𝗴𝗶𝗻𝗴, 𝗮𝗻𝗱 𝗶𝗻 𝗽𝗮𝗿𝘁𝗶𝗰𝘂𝗹𝗮𝗿, 𝗾𝘂𝗮𝗻𝘁𝗶𝘁𝗶𝗲𝘀 𝘁𝗵𝗮𝘁 𝗮𝗿𝗲 𝗰𝗼𝗻𝘀𝗶𝗱𝗲𝗿𝗲𝗱 𝘁𝗼 𝗯𝗲 𝗰𝗼𝗻𝘀𝘁𝗮𝗻𝘁𝘀 𝗼𝗳 𝗻𝗮𝘁𝘂𝗿𝗲 𝗺𝗮𝘆 𝗯𝗲 𝘃𝗮𝗿𝘆𝗶𝗻𝗴 𝘄𝗶𝘁𝗵 𝗰𝗼𝘀𝗺𝗼𝗹𝗼𝗴𝗶𝗰𝗮𝗹 𝘁𝗶𝗺𝗲. 𝗦𝘂𝗰𝗵 𝘃𝗮𝗿𝗶𝗮𝘁𝗶𝗼𝗻𝘀 𝘄𝗼𝘂𝗹𝗱 𝗰𝗼𝗺𝗽𝗹𝗲𝘁𝗲𝗹𝘆 𝘂𝗽𝘀𝗲𝘁 𝘁𝗵𝗲 𝗺𝗼𝗱𝗲𝗹 𝗺𝗮𝗸𝗲𝗿𝘀." 𝗗𝗶𝗿𝗮𝗰, 𝗣𝗮𝘂𝗹. 𝗢𝗻 𝗺𝗲𝘁𝗵𝗼𝗱𝘀 𝗶𝗻 𝘁𝗵𝗲𝗼𝗿𝗲𝘁𝗶𝗰𝗮𝗹 𝗽𝗵𝘆𝘀𝗶𝗰𝘀. (𝗧𝗿𝗶𝗲𝘀𝘁𝗲. 𝗝𝘂𝗻𝗲 𝟭𝟵𝟲𝟴 .) 

