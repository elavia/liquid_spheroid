	# Liquid Spheroid codes
	
	Julia and C++ code for the acoustic Helmholtz equation applied to liquid spheroids (prolate and oblate).
	This code calculates the scattering far-field acoustic pressure pattern using the software executables 
	for the Spheroidal wave functions developed by Adelman, Gumerov and Duraiswami
	(see http://arxiv.org/abs/1408.0074).
	
	To use our code download to a folder the contents of folder Linux or Windows in accord to your operative
	system. 
	
	In linux start a Julia shell from the command line and load it the file "JUL.main.jl": this loads the core
	routines in "JUL.routines.jl" and auxiliary functions in "JUL.auxiliar.jl". Some example scripts are
	provided for calculation of diverse cases of oblate and prolate far-field pattern in particular cases.
	These examples are based in the material exhibited in a paper submitted to Acta Acustica United with
	Acustica and in the near future probably to arXiv too.
	
	## GSL package on Julia
	
	The Julia codes uses the "GSL" (GNU Scientific Library) package for a marginal routine.
	The basic setup implies:
	
	julia > Pkg.add("GSL")
	julia > Pkg.build("GSL")
	
	but under certain circunstances one command or both fails miserably.
	This is a bug related to a dependency on linux package gsl-devel. A workaround is to install previously
	this package, via
	
	$ sudo dnf install gsl-devel
	
	and then run 
	
	$ sudo ldconfig
	
	to refresh the library database. And of course, update the Julia environment also.
