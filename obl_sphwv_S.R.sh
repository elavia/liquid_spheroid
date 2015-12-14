	# !/bin/bash

	# Shell script for calculation of S and R in one point. The inputs are 'c', 'm', 'n', 'eta', 'xi' 
	# Sintaxis:
	# 	obl_sphwv_S.R.sh	c	m	n	eta	xi
	#				$1	$2	$3	$4	$5
	
	# Spheroidal oblate function parameters
	source ./obl.parameters
	
	# Preliminaries
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w lambda
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w dr -n_dr $n_dr -dr_min $dr_min
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w dr_neg -n_dr_neg $n_dr_neg -dr_neg_min $dr_neg_min
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w N
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w F
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w k1
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w k2
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w c2k -n_c2k $n_c2k -c2k_min $c2k_min
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w Q
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w B2r -n_B2r $n_B2r -B2r_min $B2r_min
	
	# Spheroidal wave function S
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w S1 -a $4 -b $4 -d 0.01 \
	-arg_type eta -p 34
	cp Salida_S.dat Out_S.dat
	
	# Spheroidal wave function R
	./obl_sphwv -max_memory $max_mem -precision $prec -verbose n -c $1 -m $2 -n $3 -w R -a $5 -b $5 -d 0.1 \
	-arg_type xi -which R1_1,R1_2,R2_1,R2_2 -p 34
	cp Salida_R.dat Out_R.dat