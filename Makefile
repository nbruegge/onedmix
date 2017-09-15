name = 'onedmix'

all: 
	cd cvmix; \
	gfortran -c \
    cvmix_kinds_and_types.f90 \
    cvmix_utils.f90 \
    cvmix_math.f90 \
    cvmix_put_get.f90 \
    cvmix_background.f90 \
    cvmix_convection.f90 \
    cvmix_ddiff.f90 \
    cvmix_idemix.f90 \
    cvmix_kpp.f90 \
    cvmix_shear.f90 \
    cvmix_tidal.f90 \
    cvmix_tke.f90 
	cd src; \
	gfortran -c \
    onedmix_variables.f90 \
    onedmix_eos.f90 \
    onedmix_io.f90 \
    onedmix_vmix_mypp.f90 \
    onedmix_vmix_mytke.f90 \
    onedmix_setup.f90 \
    onedmix_timeloop.f90 \
    onedmix_drive.f90
	gfortran src/*.o cvmix/*.o -o bin/${name}.x
	@echo "All done!"

clean:
	@echo "Cleaning everythin up!"
	cd src; rm *.mod *.o
	cd cvmix; rm *.mod *.o
	cd bin; rm ${name}.x
