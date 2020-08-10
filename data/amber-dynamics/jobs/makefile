all : dynamics.nc 
	echo DONE

dynamics.out dynamics.info dynamics.en dynamics.rst7 dynamics.nc dynamics.log  : dynamics.in pres.rst7 complex.parm7 
	runcmd -- dynamics.in dynamics.in pres.rst7 complex.parm7  -- dynamics.out dynamics.info dynamics.en dynamics.rst7 dynamics.nc dynamics.log  -- \
	pmemd.cuda -AllowSmallBox -i dynamics.in -c pres.rst7 -ref pres.rst7 -p complex.parm7  \
	  -O -o dynamics.out -inf dynamics.info -e dynamics.en -r dynamics.rst7 -x dynamics.nc -l dynamics.log 

pres.out pres.info pres.en pres.rst7 pres.nc pres.log  : pres.in heat.rst7 complex.parm7 
	runcmd -- pres.in pres.in heat.rst7 complex.parm7  -- pres.out pres.info pres.en pres.rst7 pres.nc pres.log  -- \
	pmemd.cuda -AllowSmallBox -i pres.in -c heat.rst7 -ref heat.rst7 -p complex.parm7  \
	  -O -o pres.out -inf pres.info -e pres.en -r pres.rst7 -x pres.nc -l pres.log 

heat.out heat.info heat.en heat.rst7 heat.nc heat.log  : heat.in min.rst7 complex.parm7 
	runcmd -- heat.in heat.in min.rst7 complex.parm7  -- heat.out heat.info heat.en heat.rst7 heat.nc heat.log  -- \
	pmemd.cuda -AllowSmallBox -i heat.in -c min.rst7 -ref min.rst7 -p complex.parm7  \
	  -O -o heat.out -inf heat.info -e heat.en -r heat.rst7 -x heat.nc -l heat.log 

min.out min.info min.en min.rst7 min.nc min.log  : min.in complex.rst7 complex.parm7 
	runcmd -- min.in min.in complex.rst7 complex.parm7  -- min.out min.info min.en min.rst7 min.nc min.log  -- \
	pmemd -AllowSmallBox -i min.in -c complex.rst7 -ref complex.rst7 -p complex.parm7  \
	  -O -o min.out -inf min.info -e min.en -r min.rst7 -x min.nc -l min.log 


