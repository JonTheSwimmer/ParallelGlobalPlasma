CC := gfortran -fopenmp -mcmodel=medium

make: _base
	${CC} -o rundata ../rundata.f90 *.o ../dlsode/opk.a

_base:
	${CC} -c -o parameters_base.o ../parameters_base.f90
	${CC} -c -o parameters.o parameters.f90
	${CC} -c -o routines.o ../routines.f90
	${CC} -c -o spectrum.o ../spectrum.f90
	${CC} -c -o grid.o ../grid.f90 
	${CC} -c -o traj.o ../trajectory.f90 
	${CC} -c -o emission.o ../emission.f90
	${CC} -c -o photon.o ../photon.f90
	${CC} -c -o photon_split.o ../photon_split.f90
	${CC} -c -o parallel_prop.o ../parallel_prop.f90

run: _base
	${CC} -o rundata ../rundata.f90 *.o ../dlsode/opk.a
	./rundata 2
	-mv *.dat data

test: _base
	${CC} -o testing testing.f90 *.o ../dlsode/opk.a
	./testing
	-mv *.dat data

process: _base
	${CC} -o process ../process.f90 *.o ../dlsode/opk.a
	./process
	-mv *.dat data

reprocess: _base
	${CC} -o reprocess ../reprocess.f90 *.o ../dlsode/opk.a
	./reprocess
	-mv *.dat data

doubleprocess: _base
	${CC} -o double_process ../double_process.f90 *.o ../dlsode/opk.a
	./double_process
	-mv *.dat data

combine: _base
	${CC} -o combine ../combine.f90 *.o ../dlsode/opk.a
	
#cleanall: 
#	rm -f *.o
#	rm -f *.mod
#	rm -f data/*.dat
#	rm -f testing
#	rm -f rundata
#	rm -f process
#	rm -f double_process
#	rm -f combine
#	rm -f reprocess

clean:
	rm -f *.o
	rm -f *.mod
	rm -f testing
	rm -f rundata
	rm -f process
	rm -f double_process
	rm -f combine
	rm -f reprocess
