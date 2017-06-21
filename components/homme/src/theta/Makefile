FC = mpif90
FFLAGS = -Wall
NP = 4
PLEV = 26

../share/physical_constants.o: ../share/physical_constants.F90
	$(FC) -c $(FFLAGS) ../share/physical_constants.F90 -o ../share/physical_constants.o

../share/kinds.o: ../share/kinds.F90
	$(FC) -c $(FFLAGS) ../share/kinds.F90 -o ../share/kinds.o

../share/dimensions_mod.o: ../share/dimensions_mod.F90
	$(FC) -c $(FFLAGS) -DNP=$(NP) -DPLEV=$(PLEV) ../share/dimensions_mod.F90 -o ../share/dimensions_mod.o

../share/parallel_mod.o: ../share/parallel_mod.F90 ../share/dimensions_mod.o
	$(FC) -c $(FFLAGS) -D_MPI ../share/parallel_mod.F90 -o ../share/parallel_mod.o

../share/coordinate_systems_mod.o: ../share/coordinate_systems_mod.F90 \
	../share/physical_constants.o ../share/parallel_mod.o
	$(FC) -c $(FFLAGS) ../share/coordinate_systems_mod.F90 -o ../share/coordinate_systems_mod.o

../share/edgetype_mod.o: ../share/edgetype_mod.F90
	$(FC) -c $(FFLAGS) ../share/edgetype_mod.F90 -o ../share/edgetype_mod.o

../share/control_mod.o: ../share/control_mod.F90
	$(FC) -c $(FFLAGS) ../share/control_mod.F90 -o ../share/control_mod.o

../share/params_mod.o: ../share/params_mod.F90
	$(FC) -c $(FFLAGS) ../share/params_mod.F90 -o ../share/params_mod.o

../share/gridgraph_mod.o: ../share/gridgraph_mod.F90 ../share/control_mod.o \
	../share/params_mod.o
	$(FC) -c $(FFLAGS) ../share/gridgraph_mod.F90 -o ../share/gridgraph_mod.o

element_state.o: element_state.F90 ../share/kinds.o
	$(FC) -c $(FFLAGS) element_state.F90 -o element_state.o

../share/element_mod.o: ../share/element_mod.F90 ../share/kinds.o \
	../share/coordinate_systems_mod.o ../share/edgetype_mod.o \
	../share/gridgraph_mod.o element_state.o
	$(FC) -c $(FFLAGS) ../share/element_mod.F90 -o ../share/element_mod.o

../share/hybvcoord_mod.o: ../share/hybvcoord_mod.F90
	$(FC) -c $(FFLAGS) ../share/hybvcoord_mod.F90 -o ../share/hybvcoord_mod.o

../share/hybrid_mod.o: ../share/hybrid_mod.F90
	$(FC) -c $(FFLAGS) ../share/hybrid_mod.F90 -o ../share/hybrid_mod.o

../share/quadrature_mod.o: ../share/quadrature_mod.F90
	$(FC) -c $(FFLAGS) ../share/quadrature_mod.F90 -o ../share/quadrature_mod.o

../share/derivative_mod_base.o: ../share/derivative_mod_base.F90 ../share/quadrature_mod.o
	$(FC) -c $(FFLAGS) ../share/derivative_mod_base.F90 -o ../share/derivative_mod_base.o

derivative_mod.o: derivative_mod.F90 ../share/derivative_mod_base.o
	$(FC) -c $(FFLAGS) -ffree-line-length-none derivative_mod.F90 -o derivative_mod.o

homme_nvector.o: homme_nvector.F90 ../share/element_mod.o ../share/hybvcoord_mod.o \
	../share/hybrid_mod.o ../share/parallel_mod.o derivative_mod.o
	$(FC) -c $(FFLAGS) homme_nvector.F90 -o homme_nvector.o

homme_nvector_driver: homme_nvector_driver.F90 homme_nvector.o
	$(FC) -c $(FFLAGS) homme_nvector_driver.F90 -o homme_nvector_driver.o
	$(FC) homme_nvector_driver.o homme_nvector.o ../share/parallel_mod.o \
	  ../share/dimensions_mod.o -o homme_nvector_driver
