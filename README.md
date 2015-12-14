# liquid_spheroid
	Julia and C++ code for the acoustic Helmholtz equation applied to liquid spheroids (prolate and oblate).
	This code calculates the scattering far-field acoustic pressure pattern using the software executables 
	for the Spheroidal wave functions developed by Adelman, Gumerov and Duraiswami
	(see http://arxiv.org/abs/1408.0074).
	
	In linux start a Julia shell from the command line and load it the file "JUL.main.jl": this loads the core
	routines in "JUL.routines.jl" and auxiliary functions in "JUL.auxiliar.jl". Some example scripts are
	provided for calculation of diverse cases of oblate and prolate far-field pattern in particular cases.
	These examples are based in the material exhibited in XXXX (journal of computational physics).
