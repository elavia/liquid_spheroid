#Liquid Spheroid codes

	Julia and C++ code for calculations of acoustic scattering by liquid spheroids (prolate and oblate).
	This code calculates the scattering far-field acoustic pressure pattern using the software executables
	(pro_sphwv, obl_sphwv) for the spheroidal wave function calculation developed by Adelman, Gumerov and 
	Duraiswami available at https://github.com/radelman/scattering.git (see http://arxiv.org/abs/1408.0074
	for further details).
	
	To use our code download to a folder the contents of folder Linux or Windows according to your operative
	system. In Linux it can be necessary to give execution rights to the executables and the shell scripts.
	
	In linux start a Julia shell from the command line and load the file "JUL.main.jl": this loads the core
	routines in "JUL.routines.jl" and auxiliary functions in "JUL.auxiliar.jl". Some example scripts are
	provided for calculation of diverse cases of oblate and prolate far-field patterns in particular cases.
	These examples are based on the material exhibited in a paper submitted to Acta Acustica United with
	Acustica and also in the one that will be probably at arXiv in the near future.
	
##License

	The scattering executables (pro_sphwv, obl_sphwv) are Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, 
	and Ramani Duraiswami, and were released under the BSD 2-Clause License 
	(http://opensource.org/licenses/BSD-2-Clause).
	
##GSL package on Julia
	
	The Julia codes uses the "GSL" (GNU Scientific Library) package for a marginal routine.
	The basic setup implies:
```julia
Pkg.add("GSL")
Pkg.build("GSL")
```
	but under certain circumstances one command or both fail miserably.
	This is a bug related to a dependency on linux package gsl-devel. A workaround is to install previously
	this package, via
```console
$ sudo dnf install gsl-devel
```
	and then run 
```console
$ sudo ldconfig
```
	to refresh the library database. And of course, update the Julia environment also.
