	#
	#	Script for calculation of far-field pattern in the liquid case (sphere)
	#

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# User configurable parameters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	# Physical parameters
	rho10 = 5.00 ; # Density ratio
	k_0 = 4 ; # Wave number in media 0
	k_1 = 6 ; # Wave number in media 1	
	a = 1 ; # Sphere radius
	theta_inc = pi/4 ; # Incidence angle in radians
	
	# Software parameters
	grid_size = 200 ; # Angular grid size
	N = 150 ; # Number of terms in the series solution
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculation
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	# Far-Field pattern calculation
	Pattern = Pattern_LiquidSphere( rho10, k_0, k_1, a, grid_size, N, theta_inc ) ;
		
	# Saving to disk
	writedlm( "Out.Pattern.dat", Pattern , '\t' ) ;	