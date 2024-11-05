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

seed=42
np.random.seed(seed)

try:
    from ConfigParser import RawConfigParser  # python 2
except ImportError:
    from configparser import RawConfigParser  # python 3







debug=False









data=Table.read("ttmp.fits")
dist=10**(0.2*(19.62+5)-3) # DM to dist in kpc
cart_gal = coord.ICRS(ra=(data["RAdeg"].data)*u.degree, dec=(data["DEdeg"].data)*u.degree,distance=dist*u.kpc).transform_to(coord.Galactocentric)
data["X"]=cart_gal.x.value
data["Y"]=cart_gal.y.value
data["Z"]=cart_gal.z.value



data.write("ftp.fits")
