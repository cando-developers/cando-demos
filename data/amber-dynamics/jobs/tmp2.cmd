cd runjobs || exit
 'pmemd.cuda' '-AllowSmallBox' '-i' 'heat.in' '-c' 'min.rst7' '-ref' 'min.rst7' '-p' 'complex.parm7' '-O' '-o' 'heat.out.tmp' '-inf' 'heat.info.tmp' '-e' 'heat.en.tmp' '-r' 'heat.rst7.tmp' '-x' 'heat.nc.tmp' '-l' 'heat.log'|| exit
