


	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#	Scattering pattern liquid sphere 
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	function Pattern_LiquidSphere( rho_10::Any, k_0::Any, k_1::Any, a::Any, grid_size::Any,
					sum_size::Any, theta_inc::Any )
		#	Calculate the far-field sphere pattern in the liquid case.
		#
		#	input:	rho_10		: density quotient rho_1/rho_0
		#               k_0		: wavenumber in medium 0
		#               K_1		: wavenumber in medium 1
		#               a		: sphere radius
		#               grid_size	: angular grid points 
		#               sum_size	: number of terms in the series solution
		#		theta_inc	: incidence direction 
		#
		#	output:		Matrix of angles (in radians) (column 1) pattern (absolute value) (column 2)
		
		rho_10 = convert( AbstractFloat, rho_10 ) ;
		k_0 = convert( AbstractFloat, k_0 ) ;
		k_1 = convert( AbstractFloat, k_1 ) ;
		a = convert( AbstractFloat, a ) ;
		const Theta = linspace( 0, 2*pi , grid_size ); 
		const Pat = zeros( Complex128, grid_size, 1);	
		const g = 1/rho_10 ; # rho/rho1 ( 'g^(-1)' with the usual 'g' in scattering theory )
		const x = k_0*a ;
		const x1 = k_1*a ; 
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Main loop ( angles )
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		for j = 1 : grid_size 
			acumulator = Complex(0);
			for i = 0 : sum_size
	   			Fn = Coef_LiquidSphere( i, x1 ) ;
				C = -( Sbessj(i,x)*g*Fn - x*Sbessj_p(i,x) )/
						( Sneum(i,x)*g*Fn - x*Sneum_p(i,x) );
				acumulator += (-1)^i*( 2*i + 1 )*( ( C - 1im*C*C )/( 1 + C*C ) )*
						sf_legendre_Pl( i , cos( - Theta[j] + theta_inc + pi ) )  ;
			end
			Pat[ j ] = -1/k_0*acumulator ; # Normalization
		end
		
		return [ Theta  map( Float64, abs( Pat ) ) ] ;
	end 
	
	
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#	Prolate Spheroid coefficients (liquid)		
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	
	function Coeff_Pro( rho_10::Any, k_0::Any, k_1::Any, a::Any, b::Any, M::Any, method::Any, eta_inc::Any )
	
		#	Calculate the (Amn, Bmn) coefficients in the liquid case.
		#
		#	input:	rho_10	: density quotient rho_1/rho_0
		#               k_0	: wavenumber in medium 0
		#               K_1	: wavenumber in medium 1
		#               a	: spheroid major semiaxis
		#               b	: spheroid minor semiaxis
		#               M	: Max size of 'm' index
		#               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
		#               eta_inc	: incidence angle
		# 
		#	output:	Amn	: vector of coefficients
		#               Bmn	: vector of coefficients
		#
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Type conversion
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		
		convert( AbstractFloat, rho_10 ) ;
		convert( AbstractFloat, k_0 ) ;
		convert( AbstractFloat, k_1 ) ;
		convert( AbstractFloat, a ) ;
		convert( AbstractFloat, b ) ;
		convert( Int , M ) ;
		convert( Int , method ) ;
		convert( AbstractFloat, eta_inc ) ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Preliminaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		# Derived parameters
		
		const Size = round( Int64, (M+1)*(M+2)/2 ) ;   
		const d = 2*sqrt( a^2 - b^2 ) ;	
		const xi_0 = 1/sqrt( 1 - (b/a)^2 ) ; # Prolate
		const c_0 = ( d/2 )*k_0 ;
		const c_1 = ( d/2 )*k_1 ;		
		
		# Structure declaration
		
		const Amn = zeros( Complex{BigFloat}, Size, 1) ;
		const Bmn = zeros( Complex{BigFloat}, Size, 1) ;
		const S1 = zeros( BigFloat, Size, 1) ;
		const S1_1 = zeros( BigFloat, Size, 1) ;
		const R1_0 = zeros( BigFloat, Size, 1) ;
		const R1p_0 = zeros( BigFloat, Size, 1) ;
		const R1_1 = zeros( BigFloat, Size, 1) ;
		const R1p_1 = zeros( BigFloat, Size, 1) ;
		const R3_0 = zeros( Complex{BigFloat}, Size, 1) ;	
		const R3p_0 = zeros( Complex{BigFloat}, Size, 1) ;
		const R3_1 = zeros( Complex{BigFloat}, Size, 1) ;	
		const R3p_1 = zeros( Complex{BigFloat}, Size, 1) ;
		const Alpha = zeros( BigFloat, M + 1, M + 1 ) ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Smn, Rmn calculation (Adelman-Gumerov-Duraiswami software)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		for i = 1 : M + 1 # 'm' loop 
			for j = i : M + 1 # 'n' loop
				indice = Index( M, i, j) ;
				pro_lambdamn_approx( c_0, i-1, j-1 ) ;
				run( `./pro_sphwv_S.R.sh $c_0 $(i-1) $(j-1) $eta_inc $xi_0` );
				R_0 = ReadFileToArrayBF( "Out_R.dat", ',', 0 )  ;
				S_0 = ReadFileToArrayBF( "Out_S.dat", ',', 0 )  ;
				pro_lambdamn_approx( c_1, i-1, j-1 ) ;
				run( `./pro_sphwv_S.R.sh $c_1 $(i-1) $(j-1) $eta_inc $xi_0` );
				R_1 = ReadFileToArrayBF( "Out_R.dat", ',', 0 )  ;
				S_1 = ReadFileToArrayBF( "Out_S.dat", ',', 0 )  ;
				S1[indice] = S_0[3] ; 
				S1_1[indice] = S_1[3] ; 
				if (method == 2)# 'Method 2'
					R1_0[indice] = R_0[6] ; 
					R1p_0[indice] = R_0[7] ; 
					R3_0[indice] = complex( R_0[6], R_0[12] ) ; 
					R3p_0[indice] = complex( R_0[7], R_0[13] ) ;   
					R1_1[indice] = R_1[6] ; 
					R1p_1[indice] = R_1[7] ; 
					R3_1[indice] = complex( R_1[6], R_1[12] ) ;  
					R3p_1[indice] = complex( R_1[7], R_1[13] ) ;
				else # 'Method 1'
					R1_0[indice] = R_0[4] ; 
					R1p_0[indice] = R_0[5] ;
					R3_0[indice]= complex( R_0[4], R_0[10] ) ;
					R3p_0[indice]= complex( R_0[5], R_0[11]) ;
					R1_1[indice]= R_1[4] ; 
					R1p_1[indice]= R_1[5] ; 
					R3_1[indice]= complex( R_1[4], R_1[10] ) ;
					R3p_1[indice]= complex( R_1[5], R_1[11] ) ;
				end
			end 
		end
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Main algorithm ( Amn, Bmn calculation)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		
		for m = 0 : M # Start with the M + 1 size matrix
			TAM = M - m + 1; # Structure's size for this 'm'
			# Structure declaration
			Q = zeros( Complex{BigFloat}, TAM, TAM ); 
			F = zeros( Complex{BigFloat}, 1, TAM ); 
			# Alpha matrix calculation
			for i = 1 : TAM # fill in the rows
				sigma = m + i - 1;
				name = @ sprintf("%08d_%03d_%03d", trunc(Int,c_0*1000), m, sigma );
				pro_lambdamn_approx( c_0, m, sigma ) ;
				run(`./pro_sphwv_dr.N.sh $c_0 $m $sigma $name`) ;
				Norm0 = ReadFileToArrayBF( "Out_S_N.dat", '\n', 1) ;
				DR_0 = ReadFileToArrayBF( "Out_dr.dat", '\n', 0) ;
				for j = 1 : TAM # fill in the columns
					n = m + j - 1;
					name = @ sprintf("%08d_%03d_%03d", trunc(Int,c_1*1000), m, n );
					pro_lambdamn_approx( c_1, m, n ) ;
					run(`./pro_sphwv_dr.N.sh $c_1 $m $n $name`) ;
					Norm1 = ReadFileToArrayBF( "Out_S_N.dat", '\n', 1) ;
					DR_1 = ReadFileToArrayBF( "Out_dr.dat", '\n', 0) ;
					size_dr = min( length(DR_1) , length(DR_0) ) ; # Size of dr matrix
					Alpha[i,j] = BigFloat(0) ;	
					for k= 1 : size_dr 
						r = k - 1;
						Alpha[i,j] += DR_0[k]*DR_1[k]*2*factorial( BigInt(2*(m)+r))/
								( ( 2*m + 2*r + 1 )*factorial( BigInt(r) ) );
					end
					Alpha[i,j] = Alpha[i,j]/sqrt(Norm0*Norm1);
				end 
			end	
			# Number of coefficients
			NRO = round(Int64,(M+1)*m -(m-1)*m/2 ) ; 
			# Amn calculation 
			S1_null_flag = 1; # Smn null for default
			for i = 1 : TAM
				if ( S1[NRO+i] != 0 )  
					S1_null_flag = 0 ; # Smn is no null
					break ;	
				end
			end
			if S1_null_flag == 0
				for j = 1 : TAM
					for i = 1 : TAM
						n = m + i - 1 ;
						F[j] += -1im^n*Alpha[i,j]*S1[NRO+i]*( rho_10*R1p_0[NRO+i]*R1_1[NRO+j] -
							R1p_1[NRO+j]*R1_0[NRO+i] );
					end
				end	
				for i = 1 : TAM 
					for j = 1 : TAM 
						n = m + i - 1 ;
						Q[i,j] = 1im^n*Alpha[i,j]*S1[NRO+i]*( rho_10*R3p_0[NRO+i]*R1_1[NRO+j] -
							R1p_1[NRO+j]*R3_0[NRO+i] );
					end 
				end
				if eta_inc == 0 # Q is singular
					for i = 1 : Int( floor(TAM/2) )
						Q[ 2*i, 2*i ] = 1 ;
					end
				end 
 				Amn[ NRO+1 : NRO+TAM ] = \(Q',F') ;
			end
			# Bmn calculation
			for j = 1 : TAM
				if ( S1[ NRO + j ] == 0 ) 
					Bmn[ j + NRO ] = 0.0;
				else
					for i = 1 : TAM
						Bmn[ j + NRO ] += Alpha[i,j]*( 1/(1im) )^(i-j)*( S1[NRO+i]/S1[NRO+j] )*
								((R1_0[NRO+i] + Amn[ i+NRO ]*R3_0[NRO+i]))/R1_1[NRO+j]; 
					end 				
				end		
			end
		end
		
		return Amn, Bmn ;
 	end 		
 	
	
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#	Prolate spheroid far-field angular pattern 
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	function Pattern_LiquidPro( k_0::Any, a::Any , b::Any, M::Any,
					method::Any, eta_inc::Any, Amn::Array, deta::Any )
					
		#	Calculate the far field pattern in the liquid case.
		#
		#	input:	k_0	: wavenumber in medium 0
		#               a	: spheroid major semiaxis
		#               b	: spheroid minor semiaxis
		#               M	: Max size of 'm' index
		#               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
		#		eta_inc	: incidence angle
		#               Amn	: vector of coefficients
		#		deta	: grid size in eta ( 2/deta = # of points in the angular pattern )
		#
		#	output:		the far-field absolute value angular pattern 
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Preliminaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		# Parameter conversion
		convert( AbstractFloat, deta ) ;
		
		# Derived parameters
		const size = map( Int64, (M+1)*(M+2)/2 ); 	
		const xi_0 = 1/sqrt( 1 - (b/a)^2 ) ;		
		const d = 2*sqrt( a^2 - b^2 );
		const c_0 = (d/2)*k_0 ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#	Smn calculations
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		const Smn = cell( 1, size ); # Structure for Smn(eta) indexed by 'idx'
		const Smn_inc = zeros( BigFloat, size, 1 ); # Vector for Smn(eta_inc) 
		
		for i = 1 : M + 1 # 'm' loop
			for j = i : M + 1  # 'n' loop
 				idx = Index( M, i, j) ;
				# Calculate for all the grid and the incidence angle
				pro_lambdamn_approx( c_0, i-1, j-1 ) ;
				run(`./pro_sphwv_Sgrid.Sinc.sh $c_0 $(i-1) $(j-1) $eta_inc $deta`);
				if ( i==1 && j==1 ) # Building the grid (once)
					global ETA = reverse( readdlm( "Out_S.dat",',','\n')[:,2] ) ;
				end
				Smn[idx] = reverse( ReadFileToArrayBF2( "Out_S.dat", ',', :, 5 ) ) ;
				Smn_inc[ idx ] = ReadFileToArrayBF("Out_Sinc.dat", ',', 5 ) ;
			end
		end
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#	Far-field pattern calculation
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		
		# Theta grid
		const ETAsize = length(ETA) ;
		const TETA = map( acos, ETA ); 
		const TETA2 = zeros( Float64, ETAsize - 1, 1 );
		for k = 1 : ETAsize - 1 
			TETA2[k] = pi + TETA[k+1];
		end	
		
		const Pat = zeros( Complex{BigFloat}, 2*ETAsize-1, 1);		
		
		const Phi = [ 0, pi ]; # Only two values for Phi
		const Emn = cell( 1, length(Phi) ); # Emn: neumann_factor*cos( m(pi) )
		
		# Fill in Emn
		for k = 1 : length(Phi)
			temp = zeros( Float64, size, 1 ) ; 
			for i = 1 : (M+1) # 'm' loop
				for j = i : (M+1) # 'n' loop
					idx = Index( M, i, j) ; 
					if ( i==1 ) # m = 0 case
						temp[ idx ]= cos( (i-1)*Phi[k] ) ;
					else
						temp[ idx ]= 2*cos( (i-1)*Phi[k] ) ;
					end 			
				end
			end	
			Emn[k] = d/c_0*temp ; 
		end	
		
		# Fill in the far-field pattern ( 0 - pi )
		for l = 1 : ETAsize 
			for p = 1 : size 
				Pat[ l ] += Emn[1][p]*Amn[p]*Smn[p][l]*Smn_inc[p];
			end
		end
		
		# Fill in the far-field pattern ( pi - 2*pi )
		for l = 1 : ETAsize - 1  
			for p = 1 : size 
				Pat[ ETAsize + l ] += Emn[2][p]*Amn[p]*Smn[p][ETAsize-l]*Smn_inc[p];
			end
		end
		
		# Outputting far-field absolute value
		return [ [TETA ; TETA2] map( Float64, abs( Pat ) ) ] ;
	end  	
	

	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#	Oblate Spheroid coefficients (liquid)		
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function Coeff_Obl( rho_10::Any, k_0::Any, k_1::Any, a::Any, b::Any, M::Any, method::Any, eta_inc::Any )
	
		#	Calculate the (Amn, Bmn) coefficients in the liquid case.
		#
		#	input:	rho_10	: density quotient rho_1/rho_0
		#               k_0	: wavenumber in medium 0
		#               K_1	: wavenumber in medium 1
		#               a	: spheroid major semiaxis
		#               b	: spheroid minor semiaxis
		#               M	: Max size of 'm' index
		#               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
		#               eta_inc	: incidence angle
		# 
		#	output:	Amn	: vector of coefficients
		#               Bmn	: vector of coefficients
		#
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Type conversion
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		
		convert( AbstractFloat, rho_10 ) ;
		convert( AbstractFloat, k_0 ) ;
		convert( AbstractFloat, k_1 ) ;
		convert( AbstractFloat, a ) ;
		convert( AbstractFloat, b ) ;
		convert( Int , M ) ;
		convert( Int , method ) ;
		convert( AbstractFloat, eta_inc ) ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Preliminaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		# Derived parameters
		
		const Size = round( Int64,(M+1)*(M+2)/2 ) ;
		const d = 2*sqrt( a^2 - b^2 ) ;	
		const xi_0 = 1/sqrt( (a/b)^2 - 1 ) ; # Oblate
		const c_0 = (d/2)*k_0 ;
		const c_1 = (d/2)*k_1 ;		
		
		# Structure declaration
		
		const Amn = zeros(Complex{BigFloat}, Size, 1) ;
		const Bmn = zeros(Complex{BigFloat}, Size, 1) ;
		const S1 = zeros(BigFloat, Size, 1) ;
		const S1_1 = zeros(BigFloat, Size, 1) ;
		const R1_0 = zeros(BigFloat, Size, 1) ;
		const R1p_0 = zeros(BigFloat, Size, 1) ;
		const R1_1 = zeros(BigFloat, Size, 1) ;
		const R1p_1 = zeros(BigFloat, Size, 1) ;
		const R3_0 = zeros(Complex{BigFloat}, Size, 1) ;	
		const R3p_0 = zeros(Complex{BigFloat}, Size, 1) ;
		const R3_1 = zeros(Complex{BigFloat}, Size, 1) ;	
		const R3p_1 = zeros(Complex{BigFloat}, Size, 1) ;
		const Alpha = zeros( BigFloat, M + 1, M + 1 ) ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Smn, Rmn calculation (Adelman-Duraiswami-Gumerov software)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		for i = 1 : M + 1 # 'm' loop
			for j = i : M + 1 # 'n' loop
				indice = Index( M, i, j) ;
				obl_lambdamn_approx( c_0, i-1, j-1 ) ;
				run( `./obl_sphwv_S.R.sh $c_0 $(i-1) $(j-1) $eta_inc $xi_0` );
				R_0 = ReadFileToArrayBF( "Out_R.dat", ',', 0)  ;
				S_0 = ReadFileToArrayBF( "Out_S.dat", ',', 0)  ;
				obl_lambdamn_approx( c_1, i-1, j-1 ) ;
				run( `./obl_sphwv_S.R.sh $c_1 $(i-1) $(j-1) $eta_inc $xi_0` );
				R_1 = ReadFileToArrayBF( "Out_R.dat", ',', 0)  ;
				S_1 = ReadFileToArrayBF( "Out_S.dat", ',', 0)  ;
				S1[indice] = S_0[3] ; 
				S1_1[indice] = S_1[3] ; 
				if (method == 2) # 'Method 2'
					R1_0[indice] = R_0[5] ; 
					R1p_0[indice] = R_0[6] ; 
					R3_0[indice] = complex( R_0[5], R_0[11] ) ; 
					R3p_0[indice] = complex( R_0[6], R_0[12] ) ;   
					R1_1[indice] = R_1[5] ; 
					R1p_1[indice] = R_1[6] ; 
					R3_1[indice] = complex( R_1[5], R_1[11] ) ;  
					R3p_1[indice] = complex( R_1[6], R_1[12] ) ;
				else # 'Method 1'
					R1_0[indice] = R_0[3] ; 
					R1p_0[indice] = R_0[4] ;
					R3_0[indice]= complex( R_0[3], R_0[9] ) ;
					R3p_0[indice]= complex( R_0[4], R_0[10]) ;
					R1_1[indice]= R_1[3] ; 
					R1p_1[indice]= R_1[4] ; 
					R3_1[indice]= complex( R_1[3], R_1[9] ) ;
					R3p_1[indice]= complex( R_1[4], R_1[10] ) ;
				end
			end 
		end
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Main algorithm ( Amn, Bmn calculation)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		
		for m = 0 : M # Start with the M + 1 size matrix
			TAM = M - m + 1; # Structure's size for this 'm'
			# Structure declaration
			Q = zeros( Complex{BigFloat}, TAM, TAM ); # Matriz de ceros correspondiente a
			F = zeros( Complex{BigFloat}, 1, TAM ); # Vector correspondiente F	(es fila)	
			# Alpha matrix calculation
			for i = 1 : TAM # fill in the rows
				sigma = m + i - 1;
				name = @ sprintf("%08d_%03d_%03d", trunc(Int,c_0*1000), m, sigma );
				obl_lambdamn_approx( c_0, m, sigma ) ;
				run( `./obl_sphwv_dr.N.sh $c_0 $m $sigma $name` );
				Norm0 = ReadFileToArrayBF( "Out_S_N.dat", '\n', 1) ;
				DR_0 = ReadFileToArrayBF( "Out_dr.dat", '\n', 0) ;
				# Ciclo de 'j' 
				for j = 1 : TAM # fill in the columns
					n = m + j - 1;
					name = @ sprintf("%08d_%03d_%03d", trunc(Int,c_1*1000), m, n );
					obl_lambdamn_approx( c_1, m, n ) ;
					run( `./obl_sphwv_dr.N.sh $c_1 $m $n $name` ) ;
					Norm1 = ReadFileToArrayBF( "Out_S_N.dat", '\n', 1) ;
					DR_1 = ReadFileToArrayBF( "Out_dr.dat", '\n', 0) ;
					size_dr = min( length(DR_1) , length(DR_0) );
					Alpha[i,j] = BigFloat(0) ;	
					for k= 1 : size_dr 
						r = k - 1;
						Alpha[i,j] += DR_0[k]*DR_1[k]*2*factorial( BigInt(2*(m)+r))/
								( ( 2*m + 2*r + 1 )*factorial( BigInt(r) ) );
					end
					Alpha[i,j] = Alpha[i,j]/sqrt(Norm0*Norm1);
				end 
			end	
			# Number of coefficients
			NRO = round(Int64,(M+1)*m -(m-1)*m/2 ) ; 
			# Amn calculation 
			S1_null_flag = 1;  # Smn null for default
			for i = 1 : TAM
				if ( S1[NRO+i] != 0 ) # Smn is not null
					S1_null_flag = 0 ;
					break ;	
				end
			end
			if S1_null_flag == 0
				for j = 1 : TAM
					for i = 1 : TAM
						n = m + i - 1 ;
						F[j] += -1im^n*Alpha[i,j]*S1[NRO+i]*( rho_10*R1p_0[NRO+i]*R1_1[NRO+j] -
							R1p_1[NRO+j]*R1_0[NRO+i] );
					end
				end	
				for i = 1 : TAM 
					for j = 1 : TAM 
						n = m + i - 1 ;
						Q[i,j] = 1im^n*Alpha[i,j]*S1[NRO+i]*( rho_10*R3p_0[NRO+i]*R1_1[NRO+j] -
							R1p_1[NRO+j]*R3_0[NRO+i] );
					end 
				end
				if eta_inc == 0 # Q is singular
					for i = 1 : Int( floor(TAM/2) )
						Q[ 2*i, 2*i ] = 1 ;
					end
				end 
 				Amn[ NRO+1 : NRO+TAM ] = \(Q',F') ;
			end
			# Bmn calculation
			for j = 1 : TAM
				if ( S1[ NRO + j ] == 0 ) 
					Bmn[ j + NRO ] = 0.0;
				else
					for i = 1 : TAM
						Bmn[ j + NRO ] += Alpha[i,j]*( 1/(1im) )^(i-j)*( S1[NRO+i]/S1[NRO+j] )*
								((R1_0[NRO+i] + Amn[ i+NRO ]*R3_0[NRO+i]))/R1_1[NRO+j]; 
					end 				
				end		
			end
		end
		
		return Amn, Bmn ;
 	end 	
 	
 	
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#	Oblate spheroid far-field angular pattern 
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	function Pattern_LiquidObl( k_0::Any, a::Any , b::Any, M::Any,
					method::Any, eta_inc::Any, Amn::Array, deta::Any )
		#	Calculate the far field pattern in the liquid case.
		#
		#	input:	k_0	: wavenumber in medium 0
		#               a	: spheroid major semiaxis
		#               b	: spheroid minor semiaxis
		#               M	: Max size of 'm' index
		#               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
		#		eta_inc	: incidence angle
		#               Amn	: vector of coefficients
		#		deta	: grid size in eta ( 2/deta = # of points in the angular pattern )
		#
		#	output: 	the far-field absolute value angular pattern 	
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Preliminaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		# Parameter conversion
		convert( AbstractFloat, deta ) ;
		
		# Derived parameters
		const size = map( Int64, (M+1)*(M+2)/2 ); 	
		const xi_0 = 1/sqrt( (a/b)^2 - 1 ) ;
		const d = 2*sqrt( a^2 - b^2 );
		const c_0 = (d/2)*k_0 ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#	Smn calculations
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		const Smn = cell( 1, size ); # Structure for Smn(eta) indexed by 'idx'
		const Smn_inc = zeros( BigFloat, size, 1 ); # Vector for Smn(eta_inc) 
		
		for i = 1 : M + 1 # 'm' loop
			for j = i : M + 1  # 'n' loop
 				idx = Index( M, i, j) ;
				# Calculate for all the grid and the incidence angle
				obl_lambdamn_approx( c_0, i-1, j-1 ) ;
				run( `./obl_sphwv_Sgrid.Sinc.sh $c_0 $(i-1) $(j-1) $eta_inc $deta` );
				if ( i==1 && j==1 ) # Building the grid (once)
					global ETA = reverse( readdlm( "Out_S.dat",',','\n')[:,2] ) ; # columna 2
				end
				Smn[idx] = reverse( ReadFileToArrayBF2( "Out_S.dat", ',', :, 5 ) ) ; # Levantamos Smn calculada por el método '2'. Invierto el orden
				Smn_inc[ idx ] = ReadFileToArrayBF("Out_Sinc.dat", ',', 5 ) ;
			end
		end
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#	Far-field pattern calculation
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		
		# Theta grid
		const ETAsize = length(ETA) ;
		const TETA = map( acos, ETA ); 
		const TETA2 = zeros( Float64, ETAsize - 1, 1 );
		for k = 1 : ETAsize - 1 
			TETA2[k] = pi + TETA[k+1];
		end	
		
		const Pat = zeros( Complex{BigFloat}, 2*ETAsize-1, 1);
		
		const Phi = [ 0, pi ];  # Only two values for Phi
		const Emn = cell( 1, length(Phi) ); # Emn: neumann_factor*cos( m(pi) )
		
 		# Fill in Emn
		for k = 1 : length(Phi)
			temp = zeros( Float64, size, 1 ); 
			for i = 1 : (M+1) # 'm' loop
				for j = i : (M+1) # 'n' loop
					idx = Index( M, i, j) ;
					if ( i==1 ) # m = 0 case
						temp[ idx ]= cos( (i-1)*Phi[k] ); 
					else
						temp[ idx ]= 2*cos( (i-1)*Phi[k] );
					end 			
				end
			end	
			Emn[k] = d/c_0*temp; 
		end	
		
		# Fill in the far-field pattern ( 0 - pi )
		for l = 1 : ETAsize 
			for p = 1 : size 
				Pat[ l ] += Emn[1][p]*Amn[p]*Smn[p][l]*Smn_inc[p];
			end
		end
		
		# Fill in the far-field pattern ( pi - 2*pi )
		for l = 1 : ETAsize - 1  
			for p = 1 : size 
				Pat[ ETAsize + l ] += Emn[2][p]*Amn[p]*Smn[p][ETAsize-l]*Smn_inc[p];
			end
		end
		
		# Outputting far-field absolute value
		return [ [TETA ; TETA2] map( Float64, abs( Pat ) ) ] ;
	end  	
	

	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#	Cálculo de coeficientes Prolate Spheroid (Dirichlet-Neumann)		
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 	
 	
	function Coeff_ProDN( k_0::Any, a::Any, b::Any, M::Any, method::Any, eta_inc::Any , flag::AbstractString )
		#	Calculate the Amn coefficients in the Dirichlet-Neumann case.
		#
		#	input:	k_0	: wavenumber in medium 0           
		#               a	: spheroid major semiaxis
		#               b	: spheroid minor semiaxis
		#               M	: Max size of 'm' index
		#               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
		#               eta_inc	: incidence angle
		#		flag	: 'D' or 'N' for Dirichlet o neumann boundary conditions
		# 
		#	output:		vector of coefficients
		#
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Type conversion
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		
		convert( AbstractFloat, k_0 ) ;
		convert( AbstractFloat, a ) ;
		convert( AbstractFloat, b ) ;
		convert( Int , M ) ;
		convert( Int , method ) ;
		convert( AbstractFloat, eta_inc ) ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Preliminaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		# Derived parameters
		const Size = round( Int64, (M+1)*(M+2)/2 );  
		const d = 2*sqrt( a^2 - b^2 ) ;	
		const xi_0 = 1/sqrt( 1 - (b/a)^2 ) ; # Prolate	
		const c_0 = ( d/2 )*k_0 ;
		
		# Structure declaration
		const Amn = zeros( Complex{BigFloat}, Size, 1) ;
		const R1_0 = zeros( BigFloat, Size, 1) ;
		const R1p_0 = zeros( BigFloat, Size, 1) ;
		const R3_0 = zeros( Complex{BigFloat}, Size, 1) ;	
		const R3p_0 = zeros( Complex{BigFloat}, Size, 1) ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Rmn calculation (Adelman-Duraiswami-Gumerov software)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		for i = 1 : M + 1 # 'm' loop
			for j = i : M + 1 # 'n' loop
				indice = Index( M, i, j) ;
				pro_lambdamn_approx( c_0, i-1, j-1 ) ;
				run(`./pro_sphwv_S.R.sh $c_0 $(i-1) $(j-1) $eta_inc $xi_0`);
				R_0 = ReadFileToArrayBF( "Out_R.dat", ',', 0)  ;
				if (method == 2)# 'Method 2'
					R1_0[indice] = R_0[6] ; 
					R1p_0[indice] = R_0[7] ; 
					R3_0[indice] = complex( R_0[6], R_0[12] ) ; 
					R3p_0[indice] = complex( R_0[7], R_0[13] ) ;   
				else # 'Method 1'
					R1_0[indice] = R_0[4] ; 
					R1p_0[indice] = R_0[5] ;
					R3_0[indice]= complex( R_0[4], R_0[10] ) ;
					R3p_0[indice]= complex( R_0[5], R_0[11]) ;
				end
			end 
		end
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Amn calculation
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		if ( flag == "N" ) # Neumann
			for i = 1 : M + 1
				for j = i : M + 1 
					indice = Index( M, i, j )  ;
					if ( R3p_0[indice] != 0 )
						Amn[indice] = - R1p_0[indice]/R3p_0[indice] ;
					end 				
				end		
			end
		elseif ( flag == "D" )	# Dirichlet	
			for i = 1 : M + 1
		 		for j = i : M + 1  
		 			indice = Index( M, i, j ) ;
		 			if ( R3_0[indice] !=0 )
		 				Amn[indice] = - R1_0[indice]/R3_0[indice];
		 			end 				
		 		end		
		 	end		
		else 
			println("ERROR: Boundary conditions not defined.")
		end	
		
		return Amn ;
 	end  	
 	
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#	Oblate Spheroid coefficientes (Dirichlet-Neumann)		
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 	 
 	
	function Coeff_OblDN( k_0::Any, a::Any, b::Any, M::Any, method::Any, eta_inc::Any, flag::AbstractString  )
		#	Calculate the Amn coefficients in the Dirichlet-Neumann case.
		#
		#	input:	k_0	: wavenumber in medium 0
		#               a	: spheroid major semiaxis
		#               b	: spheroid minor semiaxis
		#               M	: Max size of 'm' index
		#               method	: method of calculation. '1' for xi >> 1 and '2' for xi ~ 1
		#               eta_inc	: incidence angle
		#		flag	: 'D' or 'N' for Dirichlet o neumann boundary conditions
		# 
		#	output:	Amn	: Vector of coefficients
		#
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Type conversion
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
		
		convert( AbstractFloat, k_0 ) ;
		convert( AbstractFloat, a ) ;
		convert( AbstractFloat, b ) ;
		convert( Int , M ) ;
		convert( Int , method ) ;
		convert( AbstractFloat, eta_inc ) ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Preliminaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		# Derived parameters
		const Size = round( Int64,(M+1)*(M+2)/2 ); 
		const d = 2*sqrt(a^2 - b^2) ;	# Oblate
		const xi_0 = 1/sqrt((a/b)^2-1) ;
		const c_0 = (d/2)*k_0 ;
		
		# Structure declaration
		const Amn = zeros(Complex{BigFloat}, Size, 1) ;
		const R1_0 = zeros(BigFloat, Size, 1) ;
		const R1p_0 = zeros(BigFloat, Size, 1) ;
		const R3_0 = zeros(Complex{BigFloat}, Size, 1) ;	
		const R3p_0 = zeros(Complex{BigFloat}, Size, 1) ;
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Rmn calculation (Adelman-Duraiswami-Gumerov software)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		
		for i = 1 : M + 1 # 'm' loop
			for j = i : M + 1 # 'n' loop
				indice = Index( M, i, j) ;
				obl_lambdamn_approx( c_0, i-1, j-1 ) ;
				run( `./obl_sphwv_S.R.sh $c_0 $(i-1) $(j-1) $eta_inc $xi_0` );
				R_0 = ReadFileToArrayBF( "Out_R.dat", ',', 0)  ;
				if (method == 2)# 'Method 2'
					R1_0[indice] = R_0[5] ; 
					R1p_0[indice] = R_0[6] ; 
					R3_0[indice] = complex( R_0[5], R_0[11] ) ; 
					R3p_0[indice] = complex( R_0[6], R_0[12] ) ;   
				else # 'Method 1'
					R1_0[indice] = R_0[3] ; 
					R1p_0[indice] = R_0[4] ;
					R3_0[indice]= complex( R_0[3], R_0[9] ) ;
					R3p_0[indice]= complex( R_0[4], R_0[10]) ;
				end
			end 
		end 	
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	Amn Calculation
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
				
		if ( flag == "N" ) # Neumann
			for i = 1 : M + 1
				for j = i : M + 1 
					indice = Index( M, i, j )  ;
					if ( R3p_0[indice] != 0 )
						Amn[indice] = - R1p_0[indice]/R3p_0[indice] ;
					end 				
				end		
			end
		elseif ( flag == "D" ) # Dirichlet		
			for i = 1 : M + 1
		 		for j = i : M + 1  
		 			indice = Index( M, i, j ) ;
		 			if ( R3_0[indice] !=0 )
		 				Amn[indice] = - R1_0[indice]/R3_0[indice];
		 			end 				
		 		end		
		 	end		
		else 
			println("ERROR: Boundary conditions not defined.")
		end 	
		
		return Amn ;
	end  	
 	
	