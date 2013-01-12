""" Some astro related functions """
from __future__ import absolute_import
from .conesearch import sphdist

import math
PI=math.pi
HALFPI = PI/2.0
D2R = PI/180.0
R2D = 1.0/D2R

def hms2deg(self, str, delim=':'):
	""" Convert hex coordinates into degrees """
	if str[0] == '-': 
		neg = -1
		str = str[1:]
	else:
		neg = 1
	str = numpy.str.split(str,delim)
	return neg*((((float(str[-1])/60.+float(str[1]))/60.+float(str[0]))/24.*360.))

def dms2deg(self, str, delim=':'):
	""" Convert hex coordinates into degrees """
	if str[0] == '-': 
		neg = -1
		str = str[1:]
	else:
		neg = 1
	str = numpy.str.split(str,delim)
	return (neg*((float(str[-1])/60.+float(str[1]))/60.+float(str[0])))

def euler(ai_in, bi_in, select, b1950=False, dtype='f8'):
	"""
	Transform between Galactic, celestial, and ecliptic coordinates.

	INPUTS:
	long_in - Input Longitude in DEGREES, scalar or vector.  
	lat_in  - Input Latitude in DEGREES
	select  - Integer (1-6) specifying type of coordinate transformation.  

	select   From          To        |   select      From            To
	1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec    
	2     Galactic       RA-DEC     |     5       Ecliptic      Galactic  
	3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic  

	Celestial coordinates (RA, Dec) should be given in equinox J2000 
	unless the b1950=True keyword is set.

	OUTPUTS:
	long_out - Output Longitude in DEGREES
	lat_out  - Output Latitude in DEGREES

	KEYWORDS:
	b1950 - If this keyword is true then input and output 
	     celestial and ecliptic coordinates should be given in equinox 
	     B1950.

	REVISION HISTORY:
	Written W. Landsman,  February 1987
	Adapted from Fortran by Daryl Yentis NRL
	Converted to IDL V5.0   W. Landsman   September 1997
	Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
	Add option to specify SELECT as a keyword W. Landsman March 2003

	Converted from IDL to numerical Python: Erin Sheldon, NYU, 2008-07-02

	"""

	# Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
	ai = numpy.array(ai_in, ndmin=1, copy=True, dtype=dtype)
	bi = numpy.array(bi_in, ndmin=1, copy=True, dtype=dtype)

	twopi   =   2.0*PI
	fourpi  =   4.0*PI

	#   J2000 coordinate conversions are based on the following constants
	#   (see the Hipparcos explanatory supplement).
	#  eps = 23.4392911111d           Obliquity of the ecliptic
	#  alphaG = 192.85948d            Right Ascension of Galactic North Pole
	#  deltaG = 27.12825d             Declination of Galactic North Pole
	#  lomega = 32.93192d             Galactic longitude of celestial equator  
	#  alphaE = 180.02322d            Ecliptic longitude of Galactic North Pole
	#  deltaE = 29.811438523d         Ecliptic latitude of Galactic North Pole
	#  Eomega  = 6.3839743d           Galactic longitude of ecliptic equator              
	# Parameters for all the different conversions
	if b1950:

		equinox = '(B1950)' 
		psi    = numpy.array([ 0.57595865315, 4.9261918136,  
				      0.00000000000, 0.0000000000,  
				      0.11129056012, 4.7005372834], dtype=dtype)
		stheta = numpy.array([ 0.88781538514,-0.88781538514, 
				      0.39788119938,-0.39788119938, 
				      0.86766174755,-0.86766174755], dtype=dtype)
		ctheta = numpy.array([ 0.46019978478, 0.46019978478, 
				      0.91743694670, 0.91743694670, 
				      0.49715499774, 0.49715499774], dtype=dtype) 
		phi  =   numpy.array([ 4.9261918136,  0.57595865315, 
				      0.0000000000, 0.00000000000, 
				      4.7005372834, 0.11129056012], dtype=dtype)
	else:

		equinox = '(J2000)'

		psi    = numpy.array([ 0.57477043300, 4.9368292465,  
				      0.00000000000, 0.0000000000,    
				      0.11142137093, 4.71279419371], dtype=dtype)
		stheta = numpy.array([ 0.88998808748,-0.88998808748, 
				      0.39777715593,-0.39777715593, 
				      0.86766622025,-0.86766622025], dtype=dtype)
		ctheta = numpy.array([ 0.45598377618, 0.45598377618, 
				      0.91748206207, 0.91748206207, 
				      0.49714719172, 0.49714719172], dtype=dtype)
		phi    = numpy.array([ 4.9368292465,  0.57477043300, 
				      0.0000000000, 0.00000000000, 
				      4.71279419371, 0.11142137093], dtype=dtype)

	# zero offset
	i  = select - 1 
	a  = ai*D2R - phi[i]

	b = bi*D2R
	sb = sin(b)
	cb = cos(b)
	cbsa = cb * sin(a)
	b  = -stheta[i] * cbsa + ctheta[i] * sb
	w, = numpy.where(b > 1.0)
	if w.size > 0:
		b[w] = 1.0
	bo = arcsin(b)*R2D
	a =  arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a) )
	ao = ( (a+psi[i]+fourpi) % twopi) * R2D
	return ao, bo



def conesearch(t, ra_name, dec_name, ra, dec, r, delim=':'):
	""" Perform a cone search on a table
	INPUTS:
		t		table		table containing the data
		ra_name		str		column name to use as RA source
		dec_name 	str		column name to use as DEC source
		ra		float/str	ra to look for (degree of sexa)
		dec		float/str	ra to look for (degree of sexa)
		r		float		distance in degrees
	
	KEYWORDS:
		delim		str		delimiter if position is a string
	"""

	x = t[ra_name]
	y = t[dec_name]
	if x.dtype.kind == 'S':
		_x = 0 
