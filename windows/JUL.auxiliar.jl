

	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	#	Using librarys
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	# GNU Scientific Library
	using GSL ;
	
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	#	Special functions related
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	# Spherical Bessel function first kind j_n(x)
	function Sbessj( i::Int64 , x )
		return sqrt(pi*(1/(2*x)))*besselj(i+1.0/2,x) ;
	end	

	# Spherical Bessel function second kind n_n(x) [Neumann]
	function Sneum( i::Int64 ,x )
		sqrt(pi*(1/(2*x)))*bessely(i+1.0/2,x) ;
	end	
	
	# Spherical Bessel function first kind derivative ( Morse & Feshbach definition )
	function Sbessj_p( i::Int64 , x ) 
 		return (1/(2*i+1))*(i*Sbessj(i-1,x) - (i+1)*Sbessj(i+1,x)) ; 
	end	
	
	# Spherical Bessel second kind ( Neumann ) derivative
	function Sneum_p( i::Int64 , x )
		return Sneum(i-1,x) - ((i+1)/x)*Sneum(i,x) ;
	end	
	
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	#	Liquid Sphere
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	# Liquid sphere coefficient
	function Coef_LiquidSphere( n::Int64, x1::Float64 )
		return x1*Sbessj_p( n, x1 )/Sbessj( n, x1 ) ;  
	end
	
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	#	Prolate/Oblate auxiliar
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	function Index( M::Int64, i::Int64, j::Int64 )
		if j < i # Domain errors
			return	display("ERROR: j < i");
		end
		if i > M + 1 || j > M + 1  
			return	display("ERROR: M < i-1, j-1");
		end
		return map( Int64, j + ( M + 1 )*( i - 1 ) - i*( i - 1 )/2 );
	end

	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	#	File I/O
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	function ReadFileToArrayBF( file::AbstractString, sep::Char, elem::Int )
	#	Read 'file' (as BigFloat) with values 'sep' separated.
	#	A newline in the 'file' means that the retrieved structure will be a matrix
	# 	The 'elem' parameter selects all the entries (elem = 0) or a particular value
	# 	(elem = number)
		if elem == 0 # Return entire vector or matrix
			return map( x -> parse(BigFloat,x), readdlm( file, sep ,ASCIIString) );
		else # Return only the specified element
			return map( x -> parse(BigFloat,x), [readdlm( file, sep ,ASCIIString)[elem]] )[ 1 ] ;
		end	
	end
	
	function ReadFileToArrayBF2( file::AbstractString, sep::Char, row::Any, col::Any )
	#	Read 'file' (as BigFloat) with values 'sep' separated.
	# 	Returned value is a matrix in general.
		return map( x -> parse(BigFloat,x), collect( readdlm( file, sep ,ASCIIString)[row,col] ) ) ;
	end 	
	
	
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	#	Eigenvalues
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	
	function pro_lambdamn_approx( c::Any, m::Int64, n::Int64 )
		# Transcription to Julia of the AGD MATLAB code.
		convert( Float64, c ) ;
		const N = m + n + 200 ; # Size of the matrix
		const A = zeros( Float64, N, N ) ;
		if ( mod( n - m, 2) == 0 )
			r = 0 ; # even
		else
			r = 1 ; # odd
		end
		# Fill in the A matrix
		for i = 1 : N
			if i == 1 
				A[ 1, 1 ] = calculate_betar( c, m, r ) ;
				A[ 1, 2 ] = calculate_alphar( c, m, r ) ;
			elseif  i == N 
				A[ N, N - 1 ] = calculate_gammar( c, m, r ) ;
				A[ N, N ] = calculate_betar( c, m, r ) ;
			else
				A[ i, i - 1 ] = calculate_gammar( c, m, r ) ;
				A[ i, i ] = calculate_betar( c, m, r ) ;
				A[ i, i + 1 ] = calculate_alphar( c, m, r ) ;
			end
			r = r + 2;
		end
		const d = eigvals!(A);
		if ( mod(n - m, 2) == 0 )
			const lambda_approx = real(d[ Int( (n - m + 2) / 2 ) ]) ;
		else
			const lambda_approx = real(d[ Int( (n - m + 1) / 2 ) ]) ;
		end
		filename = @ sprintf("data/pro_%08d_%03d_%03d_lambda_approx.txt", trunc( Int, c*1000 ), m, n ) ;
		writedlm( filename,lambda_approx );
	end
	
	function calculate_alphar( c::Float64, m::Int64, r::Int64 )
		return ((( 2 * m + r + 2 )*( 2 * m + r + 1 ))/
		(( 2 * m + 2 * r + 5 )*( 2 * m + 2 * r + 3 )))*( c*c );
	end
	
	function calculate_betar( c::Float64, m::Int64, r::Int64 )
		return ( m + r )*( m + r + 1 ) + (( 2*( m + r )*( m + r + 1 ) - 2*(m*m) - 1) /
		(( 2*m + 2*r - 1 )*( 2*m + 2 * r + 3 )))*(c*c);
	end
	
	function calculate_gammar( c::Float64, m::Int64, r::Int64 )
		return (( r * (r - 1) )/(( 2*m + 2*r - 3 )*( 2*m + 2*r - 1 )))*(c*c);
	end	
	
	function obl_lambdamn_approx( c::Any, m::Int64, n::Int64 )
		# Transcription to Julia of the AGD MATLAB code.
		convert( Float64, c ) ;
		const N = m + n + 200 ; # Size of the matrix
		const A = zeros( Float64, N, N ) ;
		if ( mod( n - m, 2) == 0 )
			r = 0 ; # even
		else
			r = 1 ; # odd
		end
		# Fill in the A matrix
		for i = 1 : N
			if i == 1 
				A[ 1, 1 ] = calculate_betar_obl( c, m, r ) ;
				A[ 1, 2 ] = calculate_alphar_obl( c, m, r ) ;
			elseif  i == N 
				A[ N, N - 1 ] = calculate_gammar_obl( c, m, r ) ;
				A[ N, N ] = calculate_betar_obl( c, m, r ) ;
			else
				A[ i, i - 1 ] = calculate_gammar_obl( c, m, r ) ;
				A[ i, i ] = calculate_betar_obl( c, m, r ) ;
				A[ i, i + 1 ] = calculate_alphar_obl( c, m, r ) ;
			end
			r = r + 2;
		end
		const d = eigvals!(A);
		if ( mod( n - m, 2 ) == 0 )
			const lambda_approx = real(d[ Int( (n - m + 2) / 2 ) ]) ;
		else
			const lambda_approx = real(d[ Int( (n - m + 1) / 2 ) ]) ;
		end
		filename = @ sprintf("data/obl_%08d_%03d_%03d_lambda_approx.txt", trunc( Int, c*1000 ), m, n ) ;
		writedlm( filename,lambda_approx );
	end	
	
	
	function calculate_alphar_obl( c::Float64, m::Int64, r::Int64 )
		return ((( 2 * m + r + 2 )*( 2 * m + r + 1 ))/
		(( 2 * m + 2 * r + 5 )*( 2 * m + 2 * r + 3 )))*( -c*c );
	end
	
	function calculate_betar_obl( c::Float64, m::Int64, r::Int64 )
		return ( m + r )*( m + r + 1 ) + (( 2*( m + r )*( m + r + 1 ) - 2*(m*m) - 1) /
		(( 2*m + 2*r - 1 )*( 2*m + 2 * r + 3 )))*(-c*c);
	end
	
	function calculate_gammar_obl( c::Float64, m::Int64, r::Int64 )
		return (( r * (r - 1) )/(( 2*m + 2*r - 3 )*( 2*m + 2*r - 1 )))*(-c*c);
	end		
