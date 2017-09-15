name = 'onedmix'

all: 
	cd src; \
	gfortran \
    onedmix_variables.f90 \
    onedmix_eos.f90 \
    onedmix_io.f90 \
    onedmix_vmix_mypp.f90 \
    onedmix_vmix_mytke.f90 \
    onedmix_setup.f90 \
    onedmix_timeloop.f90 \
    onedmix_drive.f90 \
    -o ../bin/${name}.x
#    #src/cvmix_kinds_and_types.f90 \
#    #src/cvmix_utils.f90 \
#    #src/cvmix_math.f90 \
#    #src/cvmix_put_get.f90 \
#    #src/cvmix_background.f90 \
#    #src/cvmix_convection.f90 \
#    #src/cvmix_ddiff.f90 \
#    #src/cvmix_idemix.f90 \
#    #src/cvmix_kpp.f90 \
#    #src/cvmix_shear.f90 \
#    #src/cvmix_tidal.f90 \
#    #src/cvmix_tke.f90 \
#    -o ${name}.x
	@echo "All done!"

clean:
	@echo "Cleaning everythin up!"
	cd src; rm *.mod #*.o
	cd bin; rm ${name}.x
