############################################################
## Make a self consistant model of a dwarf galaxy         ##
## of several populations using DF function and MDF       ##
##                                                        ##  
## By Guillaume F. Thomas                                 ##
## Instituto de Astrofisica de Canarias (IAC) Oct. 2024   ##
##                                                        ##  
## guillaume.thomas.astro at gmail.com                    ##
############################################################

import agama
import numpy as np
import imf
import argparse
import sys
from astropy.table import Table,vstack
import scipy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from scipy.stats import norm


seed=42
np.random.seed(seed)

try:
    from ConfigParser import RawConfigParser  # python 2
except ImportError:
    from configparser import RawConfigParser  # python 3


agama.setUnits(length=1, velocity=1, mass=1) # in Kpc, km/s, Msun

factor_length=agama.getUnits()['length'].value
factor_velocity=agama.getUnits()['velocity'].value
factor_time=agama.getUnits()['time'].value/1000 # then this is in Gyr
factor_mass=agama.getUnits()['mass'].value






debug=False








#########################
## LOAD THE ISOCHRONES ##
#########################
def load_iso(isofileName):
	"""
		Load the isochones used to get the stellar population model
		Inputs:
		-------
		isofileName: Name of the file where is the Padova Isochrone
	"""
	## Load the Isochrone ##
	iso=Table.read(isofileName,format="ascii",names=["Zini","MH","logAge","Mini","int_IMF","Mass","logL","logTe","logg","label","McoreTP","C_O","period0","period1","period2","period3","period4","pmode","Mloss","tau1m","X","Y","Xc","Xn","Xo","Cexcess","Z","mbolmag","Gmag","G_BPmag","G_RPmag"])  	
	
	iso=iso[(iso["label"]<9)] #Remove post AGB-stars

	return iso



#################################################################################
## Extract the isochone of good age and metallicity from the list of isochrone ##
#################################################################################
def extract_iso(iso,MH_sel,logAge_sel):
	"""
		Extract the isochone of good age and metallicity from the list of isochrone
		Inputs:
		-------
		 - iso: [astropy.table] containing the Padova isochrones
		 - MH_sel: selected MH for the interpolation
		 - logAge_sel: selected log(age) for the interpolation
		
		Outputs:
		--------
		 - iso_sel [astropy.table] containing the Padova isochrones of given age and metallicity
	"""
	MH_sel=round(MH_sel,1); 	logAge_sel=round(logAge_sel,5)
	list_metiso=np.array([f"{val:.1f}" for val in np.unique(iso["MH"])]).astype(float)
	list_ageiso=np.array([f"{val:.5f}" for val in np.unique(iso["logAge"])]).astype(float)
	
	# Check that the selected metallicty and age are in the list
	if(MH_sel not in list_metiso):
		sys.exit("ERROR:: MH_sel not in the isochrone file")
	else:# If inside take the real (non approx. value)
		MH_sel=np.unique(iso["MH"])[np.where(list_metiso==MH_sel)[0]]
	if(logAge_sel not in list_ageiso):
		sys.exit("ERROR:: logAge_sel "+str(logAge_sel)+" not in the isochrone file"+ str(list_ageiso))
	else: # If inside take the real (non approx. value)
		logAge_sel=np.unique(iso["logAge"])[np.where(list_ageiso==logAge_sel)[0]]		
		

	iso_sel=iso[((iso["MH"]==MH_sel)&(iso["logAge"]==logAge_sel))]
	
	return iso_sel


#############################
## To get spline isochrone ##
#############################
def get_spline_iso(iso_sel,mag="Gmag"):
	"""
		Interpolate the magnitude as function of the initial mass of a star for a given metallicity and age
		Inputs:
		-------
		iso_sel: [astropy.table] containing the Padova isochrones
		
	"""
	spline=scipy.interpolate.interp1d(iso_sel['Mini'], iso_sel[mag],bounds_error=False,fill_value=-99.9)

	if(debug):
		mplot=np.linspace(np.min(iso_sel['Mini']),np.max(iso_sel['Mini']),num=10000)
		fig=plt.figure(figsize=(10,5))
		ax1=fig.add_subplot(1,1,1)
		ax1.plot(iso_sel['Mini'], iso_sel[mag], "ko",rasterized=True)
		ax1.plot(mplot,spline(mplot), "r-",lw=1,rasterized=True)

		ax1.set_xlabel(r'Mini',fontsize=14)
		ax1.set_ylabel(mag,fontsize=14)

		fig.tight_layout()
		plt.rcParams['pdf.fonttype']=42
		plt.show()
		plt.close() 
		
	return spline
	

#############################################################################
## To know the relative fraction of stars inside a given metallicity range ##
#############################################################################
def get_fmetalbin(mu, sigma, edges_min=-2.15, edges_max=-0.15, interval=0.1):

    mu=float(mu); sigma=float(sigma)
    # Create an array of range edges
    range_edges = np.concatenate([np.array([-np.inf]),np.arange(edges_min, edges_max + interval, interval),np.array([np.inf])])

    # Calculate the fractions and centers for each range
    fractions = []
    centers = []

    for i in range(len(range_edges) - 1):
        lower = range_edges[i]
        upper = range_edges[i + 1]

        fraction = norm.cdf(upper, mu, sigma) - norm.cdf(lower, mu, sigma)
        fractions.append(fraction)
        
        # If it's the first edge, consider -inf
        if i == 0:
        	centers.insert(0, edges_min - 0.5 * interval)  # Center for -inf to first edge
        # If it's the last edge, consider +inf
        elif i == len(range_edges) - 2:
        	centers.append(edges_max - 0.5 * interval)  # Center for last edge to +inf
        else:
            centers.append((lower + upper) / 2)  # Center for middle bins

    # Convert to a NumPy array
    result_array = np.array(list(zip(centers, fractions)))
    
    return result_array




###################
## Make the SSP  ##
###################
def make_SSP(iso,FeH,logAge,Mstar,IMF,list_photoband):
	"""
	Make a SSP populations and get its magnitude in the defined photometric bands
	
	Inputs:
	-------
	 - iso: [astropy.table] containing the Padova isochrones
	 - FeH: selected metallicity
	 - logAge: log10 of the selected Age
	 - Mstar: Mass of the SSP 
	 - IMF: Selected IMF
	 - list_photoband: list of photometric band to include (Should be inside the isochrone provided)
	"""
				
	iso_sel=extract_iso(iso,FeH,logAge) # Get the isochrone corresponding to the metallicity and age
	
	data=Table()		
				
	# Get the mass of each stars assuming an IMF	
	mass = imf.make_cluster(Mstar, massfunc=IMF,mmax=np.max(iso_sel['Mini']),silent=True) 
	min_iso = np.min(iso_sel['Mini'])  # Get the minimum of iso_sel['Mini'] once

	# Update data["Mini"] to retain its values or take the minimum from iso_sel['Mini']
	data["Mini"] = np.maximum(mass, min_iso)
	del(mass)
	for band in list_photoband:
		spline_test=get_spline_iso(iso_sel,mag=band)
		data[band]=spline_test(data["Mini"])
		
	return data
		








#######################################################
## Make the dynamical distribution of the population ##
#######################################################
def make_componant_distribution(ini,popName,Mstellar,data):
	"""
	Inputs:
	-------
	ini: ini file with the different parameters
	popName: Name of the population
	Mstellar: Stellar mass of the pop 
	data: astropy Table contening the photometry and mass of the stars
	"""

	iniPoten = dict(ini.items("Potential"))
	iniDF    = dict(ini.items("DF "+popName))
	iniDF["mass"]=Mstellar # Add the stellar mass
	Nstars=len(data)

	# Define the potential
	pot=agama.Potential(**iniPoten)

	# Define the distribution function
	df = agama.DistributionFunction(potential=pot,**iniDF)

	# Make the model
	xv, mass = agama.GalaxyModel(pot, df).sample(Nstars)

	data['X']=xv[:,0]; data['Y']=xv[:,1];data['Z']=xv[:,2];
	data['Vx']=xv[:,3]; data['Vy']=xv[:,4];data['Vz']=xv[:,5]

	return data		



##########
## Main ##
##########
if __name__=='__main__':

	## Define the input parameters   
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='DESCRIPTION:     Make a realistic initial condition of a stellar componant for a dwarf galaxy with multiple chemokinematical componants',epilog="VERSION:     1.0 \nDATE:        05-Nov-2024 \nAUTHOR:      Guillaume THOMAS")
	parser.add_argument('inifile', help='File containing the parameters of the galaxy to model                     REQUESTED')
	parser.add_argument("--isofileName", default="iso/age12.dat",
                        help="Name of the file where are the isochrones from Parsec   Default: 'iso/age12.dat'")	
	parser.add_argument("--IMF", default="kroupa",
                        help="Initial mass function used   Default: kroupa")
	parser.add_argument("--list_photoband", default=["Gmag","G_BPmag","G_RPmag"],
                        help="List of photometric band to be used. The band should be included in the list of ischrones used Default: =['Gmag','G_BPmag','G_RPmag']")	
	parser.add_argument('--Age', help='log10(Age of the stellar population for the stream [Gyr]                     Default: 12.0 Gyr',type=float,default=12.0)  
	parser.add_argument('--Mmin', help='Minimum mass of a SSP to consider it [Msun]                     Default: 10.0',type=float,default=10.0)  	
	parser.add_argument("--verbose", default=False,type=bool,
                        help="Verbose      Default:False")	
	args = parser.parse_args()
                        
                        
	# Compute the stellar mass of each componant
	MsL_ratio=1.0 # Stellar mass to light ratio
	Mlum=10**6.26 # Luminous mass [Lsun]; From Battaglia and Nipoti 2022
	Mstellar_total=MsL_ratio*Mlum # Total stelalr mass
	
	inifile=args.inifile
	isofileName=args.isofileName
	list_photoband=args.list_photoband
	Age=args.Age
	IMF=args.IMF
	verbose=args.verbose
	Mmin=args.Mmin
	
	# Read the parameters file (.ini)
	ini = RawConfigParser()
	ini.optionxform=str  # do not convert key to lowercase
	ini.read(inifile)
	listPop=[item for item in ini.sections() if 'DF' not in item and 'Potential' not in item] # Get the list of the different population in the .ini file
	
	
	
	# Load the isochrones
	iso=load_iso(isofileName)
	
	output_full=Table()
	# For each population
	for popName in listPop:
		print("Pop.:", popName)
		paramPop    = dict(ini.items(popName))
		Mstellar_pop=float(paramPop["frac"])*Mstellar_total # Compute the stelalr mass in the population
		print("\tMstar= %.2g Msun" %Mstellar_pop)
		
		# Compute the stellar mass in each metallicity range of the population
		fmetal=get_fmetalbin(mu=paramPop["feh_mean"],sigma=paramPop["feh_std"]) 
		fmetal[:,1]*=Mstellar_pop
		
		data=Table()
		
		# For each slice of metallicity of the population range
		for j in range(len(fmetal)):
			Mstar_met=fmetal[j,1]
			FeH_sel=fmetal[j,0]
			if(verbose):
				print("\t\t Age= %.1f Gyr  [Fe/H]= %.1f  Mstar=%.2g Msun" %(logAge,FeH_sel,Mstar_met))
				
			# Make SSP
			if(fmetal[j,1]>Mmin): # Make the SSP  only if its stellar mass is higher than Mmin Msun
				tmp=make_SSP(iso,FeH_sel,np.log10(Age*1e9),Mstar_met,IMF,list_photoband)
				tmp["Age"]=Age
				tmp["FeH"]=FeH_sel

				data=vstack([data,tmp])

		# Make the stellar dynamical distribution 
		data=make_componant_distribution(ini,popName,Mstellar_pop,data)

		if(debug):
			data.write("out/ssp_"+popName+".fits",overwrite=True)
		data["Pop"]=popName
		output_full=vstack([output_full,data])
		del(data)
		
	output_full.write("out/model.fits",overwrite=True)
		
		
