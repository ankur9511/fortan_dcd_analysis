SHELL=/bin/bash

all-stable:
	f2py -c waterdipole.f90 voldensity3.f90 read_header.f90 qfourcalc.f90 pbcwrapxy.f90 pbcwrap.f90 pbcdr.f90 Odensity.f90 findfournn.f90 hbond_shell.f90 densfluc_methsize.f90 angle.f90 t_alldyn_z5.f90 initialize_K.f90 3drdf.f90 2drdf_broadbin_fineres.f90 allframe_alldyn3.f90 allframe_qwwcalc.f90 allframe_rdf.f90  t_allframe_structureandhb.f90 -m fort_analysis2 --opt='-O3 -fopenmp' --fcompiler=intelem --f90flags='-fopenmp -g -traceback -check uninit' --compiler=unix -DF2PY_REPORT_ATEXIT



structure: structure_hbond rdf q4 read wrap hbond_shell

dynamics: twod_msd_rotation_survival 

structure_hbond:
	f2py -c t_allframe_structureandhb.f90 read_header.f90 read.f90 pbcwrapxy.f90 pbcwrap.f90 pbcdr.f90 angle.f90 -m structure_hbond --opt='-O3 -fopenmp' --fcompiler=intelem --f90flags='-fopenmp -g -traceback -check uninit' --compiler=unix -DF2PY_REPORT_ATEXIT

rdf:
	f2py -c read_header.f90 pbcwrapxy.f90 pbcwrap.f90 read.f90 pbcdr.f90 3drdf.f90 2drdf_broadbin_fineres.f90 allframe_rdf.f90 -m fort_analysis2 --opt='-O3 -fopenmp' --fcompiler=intelem --f90flags='-fopenmp -g -traceback -check uninit' --compiler=unix -DF2PY_REPORT_ATEXIT

q4:
	f2py -c  waterdipole.f90 read_header.f90 read.f90 qfourcalc.f90 pbcwrapxy.f90 pbcwrap.f90 pbcdr.f90 findfournn.f90 angle.f90 allframe_qwwcalc.f90 --opt='-O3 -fopenmp' --fcompiler=intelem --f90flags='-fopenmp -g -traceback -check uninit' --compiler=unix -DF2PY_REPORT_ATEXIT



read:
	f2py -c  read_header.f90 read.f90 qfourcalc.f90 pbcwrapxy.f90 pbcwrap.f90 pbcdr.f90 --opt='-O3 -fopenmp' --fcompiler=intelem --f90flags='-fopenmp -g -traceback -check uninit' --compiler=unix -DF2PY_REPORT_ATEXIT



wrap:
	f2py -c  read.f90 read_header.f90 pbcwrapxy.f90 pbcwrap.f90 pbcdr.f90 --opt='-O3 -fopenmp' --fcompiler=intelem --f90flags='-fopenmp -g -traceback -check uninit' --compiler=unix -DF2PY_REPORT_ATEXIT



hbond_shell:
	f2py -c  waterdipole.f90 read.f90 read_header.f90 qfourcalc.f90 pbcwrapxy.f90 pbcwrap.f90 pbcdr.f90 findfournn.f90 hbond_shell.f90 angle.f90 --opt='-O3 -fopenmp' --fcompiler=intelem --f90flags='-fopenmp -g -traceback -check uninit' --compiler=unix -DF2PY_REPORT_ATEXIT



clean:
	rm -f *.o *.m *.so
