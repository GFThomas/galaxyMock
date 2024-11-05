############################################################
## Project on the sky a model made by makemodel.py and    ##
## add realistic uncertainties to it.                     ##
##                                                        ##  
## By Guillaume F. Thomas                                 ##
## Instituto de Astrofisica de Canarias (IAC) Oct. 2024   ##
##                                                        ##  
## guillaume.thomas.astro at gmail.com                    ##
############################################################

import numpy as np
import argparse
import sys
from astropy.table import Table,vstack
import astropy.units as u
from astropy.table import Table
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import pandas as pd
import scipy.interpolate as interpolate
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty
from pygaia.errors.photometric import magnitude_uncertainty as magunc
from numpy.linalg import inv
from dustmaps.sfd import SFDQuery
import dustmaps.sfd
#dustmaps.sfd.fetch()

seed=42
np.random.seed(seed)

try:
    from ConfigParser import RawConfigParser  # python 2
except ImportError:
    from configparser import RawConfigParser  # python 3







debug=False



def apply_extinction(data,band="Gmag"):
	""" Apply the extinction to the true photometry for a given band (redded the photometry)
		INPUTS:
		-------
		 - data astropy.tabel that contain the photometry
		 - band string. Which band to use
	"""

	if band not in ['Gmag','G_BPmag','G_RPmag']:
		sys.exit("ERROR:: "+str(band)+" not present. You should implement in in the function apply_extinction")

	if band=="Gmag":
		Al=2.664*data["EBV"]
		data[band]+=Al
		
	if band=="G_BPmag":
		Al=3.311*data["EBV"]
		data[band]+=Al
		
		
	if band=="G_RPmag":
		Al=2.021*data["EBV"]
		data[band]+=Al		
		
	return data



def get_photometry(data):
	"""
	Compute the photometry in each band by adding the distance modulus, adding extinction in the direction selected and adding uncertainties
	Inputs:
	-------
	 - Data astropy table
	"""

	# Get the amount of extinction at the position of each stellar particles
	if extinction_map=="SFD": # (Schelgel+1998 Map)
		ext = SFDQuery()
	coords = SkyCoord(data['ra'].data*u.deg, data['dec'].data*u.deg, frame='icrs')
	data['EBV'] = ext(coords)
	
	
	# Add the dm to each photometric band
	for i,band in enumerate(list_photoband):
		# Add the distance modulus
		data[band+"_true"]=data[band].copy() # Copy the true photometry
		data[band]+=float(proj["DM"])
		
		
		# ADD RANDOM FLUX (POISONNIAN NOISE)
		
		# Apply the extinction
		data=apply_extinction(data,band=band)
		
		
		# Compute the uncertainties
		if band=="Gmag":
			data=data[((data[band]<=21.0)&(data[band]>=4.0))]
			tmp=magunc("g", data[band], release=rls)
		elif band=="G_BPmag":
			data=data[((data[band]<=21.0)&(data[band]>=4.0))]
			tmp=magunc("bp", data[band], release=rls)
		elif band=="G_RPmag":
			data=data[((data[band]<=21.0)&(data[band]>=4.0))]
			tmp=magunc("rp", data[band], release=rls)
		data[band+"_err"]=tmp/1000.0 ## mmag to mag

		#Randomize the Observational photometry
		data[band]+=data[band+"_err"]*np.random.normal(0,1,len(data))
		
		del(tmp)
		
		# Remove data brighter than maglim
		data=data[(data[band]<Maglim_photoband[i])]
		
	return data




##################
### Astrometry ###
##################
def get_astrometry(data,rls="dr3"):
	"""
	Compute the astrometry astrometry and get a random value
	
	Inputs:
	-------
	  - data: Astropy.table
	  - rls: Gaia data release ('dr3', 'dr4' or 'dr5')
	"""
	

	## Parallax
	data["parallax_true"]=1.0/data["rh"]
	data["parallax_e"] = parallax_uncertainty(data["Gmag"], release=rls)/1000.0 # mas
	data['parallax']=  data["parallax_true"]+data["parallax_e"]*np.random.normal(0,1,len(data))
	
	# Proper motion
	data["pmra_e"], data["pmdec_e"] = proper_motion_uncertainty(data["Gmag"], release=rls); data["pmra_e"]/=1000.0; data["pmdec_e"]/=1000.0 # mas/yr
	data['pmra_true']=data["pmra"].copy(); data['pmdec_true']=data["pmdec"].copy()
	
	# Add the correlation between PMRA and PMDEC
	cov=np.zeros((len(data),2,2))
	cov[:,0,0]=data["pmra_e"]**2.0
	cov[:,1,1]=data["pmdec_e"]**2.0
	cov[:,0,1]=float(proj["CorrPMRA_PMDEC"])*data["pmra_e"]*data["pmdec_e"]
	cov[:,1,0]=cov[:,0,1]
	L = np.linalg.cholesky(cov)
	Npmra=np.random.normal(0,1,len(data))
	Npmdec=np.random.normal(0,1,len(data))
	data["pmra"]+=Npmra*(L[:,0,0])+Npmdec*(L[:,0,1])
	data["pmdec"]+=Npmra*(L[:,1,0])+Npmdec*(L[:,1,1])
	


	return data
	
	
	
	
	
	
#########################
### Projection on Sky ###
#########################

def do_projection():
	# Compute the ICRS (not the same as Galactocentric! https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html) XYZ Vx Vy Vz of the center of the galaxy
	dist=10**(0.2*(float(proj["DM"])+5)-3) # DM to dist in kpc
	coord_center = SkyCoord(ra=float(proj["ra"])*u.degree, dec=float(proj["dec"])*u.degree,distance=dist*u.kpc,
	pm_ra_cosdec=float(proj["pmra"])*u.mas/u.yr,pm_dec=float(proj["pmdec"])*u.mas/u.yr,radial_velocity=float(proj["Vlos"])*u.km/u.s)



	# Read the data
	data=Table.read(fileName)
	
	

	# Apply the ellipticity
	data['Z']*=(1.0-float(proj["ell"])) # See Arroyo-Polonio et al. 2025. Do it along Z axis for when we will include rotation
	for param in ["X","Y","Z","Vx","Vy","Vz"]:
		data[param+"init"]=data[param].copy()

	
	# Do the rotation to include the Position angle and the inclination. As for Metz, Kroupa, Jerjen et al. 2006 Assume that LOS along X axis. See their appendinx B
	ac=np.deg2rad(float(proj["ra"])); dc=np.deg2rad(float(proj["dec"]));
	PA=np.deg2rad(90.0-float(proj["PA"])); incl=np.deg2rad(90.0-float(proj["incl"]))
	Rrpq=np.array([[np.cos(dc)*np.cos(ac),-np.sin(ac),-np.sin(dc)*np.cos(ac)],[np.cos(dc)*np.sin(ac),np.cos(ac),-np.sin(dc)*np.sin(ac)],[np.sin(dc),0,np.cos(dc)]])
	Rx=np.array([[1.0,0.0,0.0],[0.0,np.cos(PA),np.sin(PA)],[0.0,-np.sin(PA),np.cos(PA)]])
	Ry=np.array([[np.cos(incl),0.0,-np.sin(incl)],[0.0,1.0,0.0],[np.sin(incl),0.0,np.cos(incl)]])
	R=np.matmul(np.matmul(Ry,Rx),inv(Rrpq))


	x_icrs=R[0,0]*data["X"]+R[1,0]*data["Y"]+R[2,0]*data["Z"]
	y_icrs=R[0,1]*data["X"]+R[1,1]*data["Y"]+R[2,1]*data["Z"]
	z_icrs=R[0,2]*data["X"]+R[1,2]*data["Y"]+R[2,2]*data["Z"]
	
	
	vx_icrs=R[0,0]*data["Vx"]+R[1,0]*data["Vy"]+R[2,0]*data["Vz"]
	vy_icrs=R[0,1]*data["Vx"]+R[1,1]*data["Vy"]+R[2,1]*data["Vz"]
	vz_icrs=R[0,2]*data["Vx"]+R[1,2]*data["Vy"]+R[2,2]*data["Vz"]
	
	
	# Move the model to the desired position (Still in ICRS carthesian coordinates)
	x_icrs+=coord_center.cartesian.x.value; y_icrs+=coord_center.cartesian.y.value; z_icrs+=coord_center.cartesian.z.value
	vx_icrs+=coord_center.cartesian.differentials['s'].d_x.value; vy_icrs+=coord_center.cartesian.differentials['s'].d_y.value; vz_icrs+=coord_center.cartesian.differentials['s'].d_z.value
	



	# Project the model at the good position
	coordi=SkyCoord(x=x_icrs.data*u.kpc, y=y_icrs.data*u.kpc,z=z_icrs.data*u.kpc,
	v_x=vx_icrs.data*u.km/u.s, v_y=vy_icrs.data*u.km/u.s, v_z=vz_icrs.data*u.km/u.s,frame="icrs",representation_type="cartesian")
	coordi.representation_type=coord.SphericalRepresentation
	data['ra']=coordi.ra.value
	data['dec']=coordi.dec.value
	data['rh']=coordi.distance.value
	data['pmra']=coordi.pm_ra*np.cos(np.deg2rad(coordi.dec.value)) # Compute PMRA cos dec!
	data['pmdec']=coordi.pm_dec.value
	data['vlos']=coordi.radial_velocity.value
	
	# Get the Galactocentric position
	data["X"]=coordi.galactocentric.x.value
	data["Y"]=coordi.galactocentric.y.value
	data["Z"]=coordi.galactocentric.z.value
	data["Vx"]=coordi.galactocentric.v_x.value
	data["Vy"]=coordi.galactocentric.v_y.value
	data["Vz"]=coordi.galactocentric.v_z.value
	
	
	#######################
	## Compute (xki,eta) ##
	#######################
	def get_xkieta(ra,dec,raZ,decZ):
	  	# ra and dec are coordinates of target in radians
	  	# raZ and decZ are coordinates of tangential point in radians

		ra=ra/180.0*np.pi
		dec=dec/180.0*np.pi
		raZ=raZ/180.0*np.pi
		decZ=decZ/180.0*np.pi
		sdecZ=np.sin(decZ)
		sdec=np.sin(dec)
		cdecZ=np.cos(decZ)
		cdec=np.cos(dec)
		radiff=ra-raZ
		sradiff=np.sin(radiff)
		cradiff=np.cos(radiff)

		denom=sdec*sdecZ+cdec*cdecZ*cradiff

		xki=cdec*sradiff/denom
		eta=(sdec*cdecZ-cdec*sdecZ*cradiff)/denom

		return xki*180.0/np.pi,eta*180.0/np.pi
	data["xki"],data["eta"]=get_xkieta(data["ra"],data["dec"],float(proj["ra"]),float(proj["dec"]))
	
	
	return data



##########
## Main ##
##########
if __name__=='__main__':

	## Define the input parameters   
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='DESCRIPTION:     Make a mock galaxy from an existing model and project it on the sky following the parameters set in the ini file',epilog="VERSION:     1.0 \nDATE:        05-Nov-2024 \nAUTHOR:      Guillaume THOMAS")
	parser.add_argument('fileName', help='File containing the model                     REQUESTED')
	parser.add_argument('inifile', help='File containing the parameters for the projection                     REQUESTED')
	parser.add_argument("--list_photoband", default=["Gmag","G_BPmag","G_RPmag"],
                        help="List of photometric band to be used. The band should be included in the list of ischrones used Default: =['Gmag','G_BPmag','G_RPmag']")	
	parser.add_argument("--Maglim_photoband", default=[21,21,21],
                        help="Magnitude limit of each photometric bands. Default: = [21,21,21]")	
	parser.add_argument("--extinction_map", default="SFD",
                        help="Extinction Map to choose. Default: = SFD")	
	parser.add_argument("--gaia_rls", default="dr3",
                        help="Gaia data release to use for the Astrometry and/or photometry. Default: = 'dr3'")
	parser.add_argument("--with_astro", default=True,type=bool,
                        help="With Astrometric measurement      Default:True")	
	parser.add_argument("-v","--verbose", default=False,type=bool,
                        help="Verbose      Default:False")	



	args = parser.parse_args()
	
		
	
	fileName=args.fileName
	inifile=args.inifile
	list_photoband=args.list_photoband
	Maglim_photoband=args.Maglim_photoband
	extinction_map=args.extinction_map
	rls=args.gaia_rls
	w_astro=args.with_astro
	verbose=args.verbose
	

	# Read the parameters file (.ini)
	ini = RawConfigParser()
	ini.optionxform=str  # do not convert key to lowercase
	ini.read(inifile)
	proj    = dict(ini.items("Projection"))
	


	# Do the projection of the galaxy on the sky
	data=do_projection()
	
	# Include photometry uncertainties
	data=get_photometry(data)
		
	
	# Include astrometric uncertainties
	if w_astro:
		data=get_astrometry(data,rls=rls)
	
		
	data.write("out/mock.fits",overwrite=True)
		
	
	
