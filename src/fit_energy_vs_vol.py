from   sys   import stdin
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk 
import pylab             as pyl
import matplotlib.style

bohr_radius     = 0.529177
bohr32ang3		= 0.14818474347
joule2hartree	= 4.3597482
joule2rydberg   = joule2hartree/2.
unitconv        = joule2hartree/bohr_radius**3*10.**3
	
def energy_vs_volume():
	import fnmatch
	mypath = os.getcwd()
	os.system("rm energy-vs-volume energy-vs-strain")
	eV2Hartree=0.036749309
	Ang32Bohr3=6.74833304162

	E=[]
	dir_list=[]
	count = 0
	dir_E=[];
	vol_cell=[]
	strain_file=[]
	strain_value=[] # strain_value is deformation

	print ("               >>>>> Extracting Energies from directories  <<<<<<")
	for entry in os.listdir(mypath):
		if not os.path.exists('strain-01'):
			print (' ERROR: strain-* files does not exist. create a strain file for each deformation values.')
			sys.exit(0)
		if fnmatch.fnmatchcase(entry,'strain-*'):
			with open(entry,'r') as f:
				lines = f.readline()  #Read  first line only
				strain_value.append( float(lines) )
			if os.path.isfile(os.path.join(mypath, entry)):
				strain_file.append(entry)

		if os.path.isdir(os.path.join(mypath, entry)):
			dir_list.append(entry)

			for file in os.listdir(entry):
				if file == "OUTCAR":
					filepath = os.path.join(entry, file)
					if not os.path.exists(filepath):
						print (' ERROR: OUTCAR does not exist here.')
						sys.exit(0)
					with open(filepath,'r') as f:
						lines = f.readlines()
					for i in lines:
						if "  free  energy   TOTEN  =" in i:
							m=float(i.split()[4])
						if "  volume of cell :"	in i:
							v=float(i.split()[4])
					vol_cell.append(v)
					E.append(m)
					count+=1
	print("# of folders detected: ", count)
	print ("Directory :%10.6s %14s %18s %25.20s " % ("Folder", "Energy(eV)", "Vol_of_cell(A^3)", "strain_deformation" ))

	for i in range(math.floor(count/2)): # 0 to 4
		print ("Folder name: %10.10s %16.8f %16.8f %16.12s %14.4f" %(dir_list[i], E[i], vol_cell[i], strain_file[i], strain_value[i] ))
	if (bool(math.floor(count/2))):
		i = math.floor(count/2)
		print(Back.GREEN + 'Folder name: %10.10s %16.8f %16.8f %16.12s %14.4f <--Ref' % (dir_list[i], E[i], vol_cell[i], strain_file[i], strain_value[i] ))		
		print(Style.RESET_ALL, end="")
	for i in range(math.ceil(count/2), count):
		print ("Folder name: %10.10s %16.8f %16.8f %16.12s %14.4f" %(dir_list[i], E[i], vol_cell[i], strain_file[i], strain_value[i] ))
	#rc = subprocess.Popen(['bash', 'extract_energy.sh'])

	print (colored('ENERGIES & VOLUMES ARE WRITTEN IN ATOMIC UNITS TO A FILE <energy-vs-volume>','yellow'), end = '\n', flush=True)
	print (colored('IT WILL BE READ BY ELASTIC SCRIPTS FOR POSTPROCESSING eV-->Ha; A^3-->Bohr^3','yellow'), end = '\n', flush=True)	

	with open("energy-vs-volume",'w') as file:
		for i in range(count):
			file.write ("%14.6f %14.6f\n" %(vol_cell[i] * Ang32Bohr3, E[i] * eV2Hartree))
	with open("energy-vs-strain",'w') as file:
		for i in range(count):
			file.write ("%12.6f %14.6f\n" %(strain_value[i], E[i] * eV2Hartree))
###

'''
#####---------------------------------------------------------
#            5- FITTING ENERGY vs VOLUME CURVE from VASP OUTCAR file
#####--------------------------------------------------------- 
'''

def fitting_energy_vs_volume_curve():
	from   sys   import stdin
	import matplotlib.pyplot as plt
	import matplotlib.ticker as ptk
	import pylab             as pyl
	import matplotlib.style

	bohr_radius     = 0.529177
	bohr32ang3		= 0.14818474347
	joule2hartree	= 4.3597482
	joule2rydberg   = joule2hartree/2.
	unitconv        = joule2hartree/bohr_radius**3*10.**3

	if (str(os.path.exists('energy-vs-volume'))=='False'): 
		sys.exit("ERROR: file energy-vs-volume not found!\n")
	energy = []
	volume = []
	read_energy = open('energy-vs-volume',"r")

	while True:
		line = read_energy.readline()
		line = line.strip()
		if len(line) == 0: break
		energy.append(float(line.split()[1]))
		volume.append(float(line.split()[0]))
	volume,energy=sortvolume(volume,energy)

	print ("===============================")
	print ("Lattice symmetry codes"         )
	print ("-------------------------------")
	print ("1 --> Simple cubic (sc)"        )
	print ("2 --> Body-centered cubic (bcc)")
	print ("3 --> Face-centered cubic (fcc)")
	print ("-------------------------------")
	print ("0 --> Others")
	print ("===============================\n")
	scheck = input("Enter lattice symmetry code [default 0] >>>> ").replace(" ", "") 

	'''
	These factors are for the conversion from the conventional cell to primitive cell
	BCC: (a^3)/2 primitive cell volume
	FCC: (a^3)/4 primitive cell volume
	
	'''

	isym   = 0
	factor = 1
	if scheck == "1":
		isym = 1
		factor=1
		slabel = "(sc) "
	elif scheck == "2":
		isym = 2
		factor=2
		slabel = "(bcc)"
	elif scheck == "3":
		isym = 3
		factor=4
		slabel = "(fcc)"
	print ("Verification lattice symmetry code      >>>> %d " %(isym) )
#-------------------------------------------------------------------------------
	print(f"{'-' * 105}")
	print ('%20.25s %29.30s %21.30s %11.12s %18.30s' %("Opt_vol Bohr^3 (Ang^3)", "Lattice_const Bohr (A)", "Bulk_modulus [GPa]", "Log(chi)", "Polynomial_order"))
	print(f"{'-' * 105}")
	for order_of_fit in range(2, 11): #order of polynomial fitting
		if order_of_fit % 2 == 0: 
			order_of_fit = int(order_of_fit)
			fitr = np.polyfit(volume,energy,order_of_fit)
			curv = np.poly1d(fitr)
			oned = np.poly1d(np.polyder(fitr,1)) #
			bulk = np.poly1d(np.polyder(fitr,2))
			bpri = np.poly1d(np.polyder(fitr,3)) #
			vmin = np.roots(np.polyder(fitr))
			dmin = [
				vmin[i].real
				for i in range(len(vmin))
				if (abs(vmin[i].imag) < 1.0e-10)
				and volume[0] <= vmin[i] <= volume[-1]
				and (bulk(vmin[i]) > 0)
			]

			xvol = np.linspace(volume[0],volume[-1],100)
			if (len(dmin) > 1): print ("WARNING: Multiple minima are found!\n")
			if not dmin: print ("WARNING: No minimum in the given xrange!\n")

			chi = 0
			for i in range(len(energy)): 
				chi=chi+(energy[i]-curv(volume[i]))**2
			chi=math.sqrt(chi)/len(energy)
#-------------------------------------------------------------------------------			
			for i in range(len(dmin)):
				x0=dmin[len(dmin)-1-i]
				v0=dmin[len(dmin)-1-i]
				a0=(factor*v0)**(0.33333333333)
				b0=bulk(v0)*v0*unitconv
				#if (isym > 0): print ('Lattice constant = %5s %12.6f %12.6f' %(slabel,a0, a0*bohr_radius), '[Bohr, Angstrom]')	
				if (isym > 0): 
					print("%12.6f (%11.6f) %12.6f (%9.6f) %17.6f %13.2f %10d\n" %(v0, v0*bohr32ang3, a0, a0*bohr_radius, b0, math.log10(chi), order_of_fit), end="") 				
				else: 
					print("%12.6f(%12.6f) %12.6f(%12.6f) %17.6f %13.2f %10d\n" %(v0, v0*bohr32ang3, a0, a0*bohr_radius, b0, math.log10(chi), order_of_fit), end="")
	print(f"{'-' * 105}")
#-------------------------------------------------------------------------------

	xlabel = u'Volume [Bohr\u00B3]'
	ylabel = r'Energy [Ha]'

	fontlabel=20
	fonttick=14

	params = {'ytick.minor.size': 6,
			'xtick.major.pad': 8,
			'ytick.major.pad': 4,
			'patch.linewidth': 2.,
			'axes.linewidth': 2.,
			'lines.linewidth': 1.8,
			'lines.markersize': 8.0,
			'axes.formatter.limits': (-4, 6)}

	plt.rcParams.update(params)
	plt.subplots_adjust(left=0.21, right=0.93,
						bottom=0.18, top=0.88,
						wspace=None, hspace=None)

	yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)                           

	figure = plt.figure(1, figsize=(9,9))
	ax     = figure.add_subplot(111)
	ax.text(0.5,-0.18,xlabel,size=fontlabel, transform=ax.transAxes,ha='center',va='center')
	ax.text(-0.25,0.5,ylabel,size=fontlabel, transform=ax.transAxes,ha='center',va='center',rotation=90)
	for line in ax.get_xticklines() + ax.get_yticklines():
		line.set_markersize(6)
		line.set_markeredgewidth(2)
	plt.xticks(size=fonttick)
	plt.yticks(size=fonttick)
	pyl.grid(True)
	plt.plot(xvol, curv(xvol), 'b-', label=f'n={str(order_of_fit)} fit')
	plt.plot(volume,energy,'go',label='calculated')
	plt.plot(dmin,curv(dmin),'ro')
	plt.legend(loc=9,borderaxespad=.8,numpoints=1)

	ymax  = max(max(curv(xvol)),max(energy))
	ymin  = min(min(curv(xvol)),min(energy))
	dxx   = abs(max(xvol)-min(xvol))/18
	dyy   = abs(ymax-ymin)/18
	ax.yaxis.set_major_formatter(yfmt)
	ax.set_xlim(min(xvol)-dxx,max(xvol)+dxx)
	ax.set_ylim(ymin-dyy,ymax+dyy)

	ax.xaxis.set_major_locator(MaxNLocator(7))

	plt.savefig('PLOT.png',orientation='portrait',format='png')
###


def sortvolume(s,e):
	ss=[]
	ee=[]
	ww = [s[i] for i in range(len(s))]
	ww.sort()
	for i in range(len(s)):
	    ss.append(s[s.index(ww[i])])
	    ee.append(e[s.index(ww[i])])
	return ss, ee
###