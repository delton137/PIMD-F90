geo_tester.x: pflush.f90 fsiesta_pipes.f90 geo_tester.f90 
	gfortran -cpp -O3 pflush.f90 fsiesta_pipes.f90 geo_tester.f90 -o geo_tester.x #--debug --bounds-check --backtrace
clean:
	rm *.o

pristine: 
	rm *.x
