# Parsing Hessian matrix from a file
*Script to parse Hessian Matrix from a file*


USAGE	   	 : To parse VASP vasprun.xml input file to extract Hessian Matrix 
Output file has already been defined in the code (hessian.dat).

CAUTION: Use at your own risk (NOTEVEN IMPLIED GUARANTEED, WHATSOEVER),
the code has been tested but the user in the end will have to verify the ouput.

Xpath is useful tool for directly accessing the element in a Node

->   CHECK this website for lxml introduccion: https://lxml.de/tutorial.html &
->   https://github.com/lxml/lxml


NB: The Matrix is obtained from VASP, vasprun.xml file. It contians the hessian matrix that can be edited with appropriate script.

The complication that arises by parsing the data from the file is the trailing empty lines and tabs that need to be deleted before it can be read. In the case of a simple file format, there is no need for it. But if the file contains irregular data entry such as empty lines with commas and spaces then it needs to be formatted. 

In my case, I have only leading empty lines and spaces and in between columns there are white spaces of different sizes.


AUTHOR 	   : ASIF IQBAL BHATTI
CREATED ON	 : 09/02/2020
