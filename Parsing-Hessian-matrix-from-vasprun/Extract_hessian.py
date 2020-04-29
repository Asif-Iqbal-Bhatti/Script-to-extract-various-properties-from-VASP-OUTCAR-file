#!/usr/bin/env python3

import lxml
from yaml import load, dump
import sys, os, yaml
#from   lxml  import etree
import numpy as np
from scipy import linalg as LA

#*******************************DOCUMENTATION***********************************
#
#AUTHOR 	   : ASIF IQBAL BHATTI
#CREATED ON	 : 09/02/2020
#USAGE	   	 : To parse VASP vasprun.xml input file to extract Hessian Matrix
#Output file has already been defined in the code (hessian.dat).
#CAUTION: Use at your own risk (NOTEVEN IMPLIED GUARANTEED, WHATSOEVER),
#the code has been tested but the user in the end will have to verify the ouput. 
#Xpath is useful tool for directly accessing the element in a Node
#
#->   CHECK this website for lxml introduccion: https://lxml.de/tutorial.html &
#->   https://github.com/lxml/lxml
#
#************************END OF DOCUMENTATION***********************************

try:
  from lxml import etree
  print("running with lxml.etree")
except ImportError:
  try:
    # Python 2.5
    import xml.etree.cElementTree as etree
    print("running with cElementTree on Python 2.5+")
  except ImportError:
    try:
      # Python 2.5
      import xml.etree.ElementTree as etree
      print("running with ElementTree on Python 2.5+")
    except ImportError:
      try:
        # normal cElementTree install
        import cElementTree as etree
        print("running with cElementTree")
      except ImportError:
        try:
          # normal ElementTree install
          import elementtree.ElementTree as etree
          print("running with ElementTree")
        except ImportError:
          print("Failed to import ElementTree from any known place")
#-----------------------------------------------------------------------------------
#                             ''' Readng a File '''
#-----------------------------------------------------------------------------------
input_obj = open("vasprun.xml","r")
input_doc = etree.parse(input_obj)
root = input_doc.getroot()

print("*"*80)
print("dynmat")
for type_tag in root.findall('*/dynmat/varray'): # search string in a Node
	value = type_tag.get('name')
	print("   |--> Child detected in a Node", value)
print("*"*80)
#-----------------------------------------------------------------------------------

#*********************************Extracting Hessian********************************
#
#We know the size of Hessian is always 3*N where N is the # of atoms in the unitcell
# H = (3N,3N)
#-----------------------------------------------------------------------------------

hessian = root.xpath('*/dynmat/varray/v/text()')	
hess_v = []

for ind_hess in hessian:
    hess_v.append( ind_hess.split() )
hess_v = np.array(hess_v)

print("*"*80)
#print(hess_v) 								# for debugging ONLY
row = len(hess_v); print("# of rows:", row)
col = len(hess_v[0]); print("# of columns:", col)
print("*"*80)
#-----------------------------------------------------------------------------------

#*********************** writing Hessian to a file ********************************
out = open('hessian.dat','w')

for i in range( int(row/2) ):
	for j in range(col):
		#print(hess_v[i][j], end=' ' )
		out.write( "{:18.11s}".format( hess_v[i][j])  )
	out.write('\n')	
	#print()	
	
out.close()

#******************************End Hessian******************************************


#***************************extracting basevectors**********************************
'''
lst_basevect = root.xpath('structure/crystal/varray/v/text()')	
xml_basevect = []

for ind_basevect in lst_basevect:
    xml_basevect.append( ind_basevect.split() )
print(np.array(xml_basevect))
	
#for i in range( (len(xml_basevect)) ):
#	print(i, xml_basevect[i][0], xml_basevect[i][1], xml_basevect[i][2])

#******************************End basevectors********************************
'''
#---------------------------------'''METHOD 2'''---------------------------------  
	
a = [] ; l=0
f = open('hessian.dat', 'r')

for line in f.readlines():
	a.append([])
	#print (line[0])
	if (line[0] == '#'):
		continue
	for i in line.strip().split():
		a[-1].append(float(i))
	l+=1		
#print (a)

f.close()
#------------------------------------------------------------------------
#a = np.random.random((900, 900))	
s = np.shape(a); print('size of a Matrix', s)
a = np.matrix(a)
#print(a)
print ('--------------------------------EIGVAL AND EIGVECTORS ARE:')

eig, eig_v = LA.eigh(a)

print (eig);
print ( np.transpose(eig_v) )
