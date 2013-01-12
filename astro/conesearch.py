"""
Some astronomy related operations
"""

from numpy import deg2rad,rad2deg,sin,cos,sqrt,arcsin

def sphdist (ra1, dec1, ra2, dec2):
	"""measures the spherical distance between 2 points 
	Inputs:
		(ra1, dec1)	in degrees
		(ra2, dec2)	in degrees
	Outputs:
		returns a distance in degrees
	"""
	dec1_r = deg2rad(dec1)
	dec2_r = deg2rad(dec2)
	return 2 * rad2deg ( arcsin ( sqrt ( ( sin((dec1_r - dec2_r) / 2))**2 +  cos(dec1_r) * cos
(dec2_r) * ( sin((deg2rad(ra1 - ra2)) / 2))**2)))	

