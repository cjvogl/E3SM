clean:
	rm -f *.o \
	rm -f ../repro_sum_mod.o \
	rm -f ../share/*.o \
	rm -f $(SHR_PATH)/shr_sys_mod.o $(SHR_PATH)/shr_sys_mod.mod \
	rm -f $(SHR_PATH)/shr_mpi_mod.o $(SHR_PATH)/shr_mpi_mod.mod \
	rm -f homme_nvector_driver

bndry_mod.o: bndry_mod.F90 \
	../share/bndry_mod_base.o
	$(FC) -c $(FFLAGS) -o bndry_mod.o bndry_mod.F90

edge_mod.o: edge_mod.F90 \
	../share/edge_mod_base.o ../share/edgetype_mod.o
	$(FC) -c $(FFLAGS) -o edge_mod.o edge_mod.F90

element_state.o: element_state.F90 \
	../share/kinds.o
	$(FC) -c $(FFLAGS) -o element_state.o element_state.F90

homme_nvector.o: homme_nvector.F90 \
	../share/parallel_mod.o ../share/global_norms_mod.o
	$(FC) -c $(FFLAGS) -o homme_nvector.o homme_nvector.F90

homme_nvector_driver: homme_nvector_driver.F90 homme_nvector.o
	$(FC) -c $(FFLAGS) -o homme_nvector_driver.o homme_nvector_driver.F90
	$(FC) -L$(NETCDF_PATH)/lib -lnetcdff -o homme_nvector_driver \
		../share/parallel_mod.o ../share/thread_mod.o ../share/edgetype_mod.o \
		../share/dimensions_mod.o ../share/global_norms_mod.o ../share/mesh_mod.o \
		../repro_sum_mod.o ../share/reduction_mod.o ../share/physical_constants.o \
		../share/quadrature_mod.o ../share/control_mod.o edge_mod.o \
		../share/edge_mod_base.o bndry_mod.o ../share/bndry_mod_base.o \
		../share/coordinate_systems_mod.o ../share/gridgraph_mod.o \
		../share/schedtype_mod.o ../share/metagraph_mod.o ../share/ll_mod.o \
		../share/element_mod.o ../share/kinds.o element_state.o \
		../share/cube_mod.o ../share/spacecurve_mod.o $(SHR_PATH)/shr_sys_mod.o \
		$(SHR_PATH)/shr_mpi_mod.o \
		homme_nvector.o homme_nvector_driver.o $(GPTL_PATH)/libtiming.a

../repro_sum_mod.o: ../repro_sum_mod.F90 \
	$(SHR_PATH)/shr_sys_mod.o ../share/thread_mod.o
	$(FC) -c $(FFLAGS) ../repro_sum_mod.F90 -o ../repro_sum_mod.o

../share/bndry_mod_base.o: ../share/bndry_mod_base.F90
		$(FC) -c $(FFLAGS) ../share/bndry_mod_base.F90 -o ../share/bndry_mod_base.o

../share/control_mod.o: ../share/control_mod.F90
	$(FC) -c $(FFLAGS) ../share/control_mod.F90 -o ../share/control_mod.o

../share/coordinate_systems_mod.o: ../share/coordinate_systems_mod.F90
	$(FC) -c $(FFLAGS) ../share/coordinate_systems_mod.F90 -o \
	../share/coordinate_systems_mod.o

../share/cube_mod.o: ../share/cube_mod.F90 \
	../share/element_mod.o
	$(FC) -c $(FFLAGS) ../share/cube_mod.F90 -o ../share/cube_mod.o

../share/dimensions_mod.o: ../share/dimensions_mod.F90
	$(FC) -c $(FFLAGS) ../share/dimensions_mod.F90 -o ../share/dimensions_mod.o

../share/edge_mod_base.o: ../share/edge_mod_base.F90 \
	../share/schedtype_mod.o
	$(FC) -c $(FFLAGS) -o ../share/edge_mod_base.o ../share/edge_mod_base.F90

../share/edgetype_mod.o: ../share/edgetype_mod.F90
	$(FC) -c $(FFLAGS) -o ../share/edgetype_mod.o ../share/edgetype_mod.F90

../share/element_mod.o: ../share/element_mod.F90 \
	../share/kinds.o ../share/coordinate_systems_mod.o ../share/edgetype_mod.o \
	../share/gridgraph_mod.o element_state.o
	$(FC) -c $(FFLAGS) -o ../share/element_mod.o ../share/element_mod.F90

../share/global_norms_mod.o: ../share/global_norms_mod.F90 \
	../share/mesh_mod.o ../repro_sum_mod.o ../share/reduction_mod.o \
	../share/physical_constants.o ../share/quadrature_mod.o ../share/control_mod.o \
	edge_mod.o bndry_mod.o
	$(FC) -c $(FFLAGS) -I$(NETCDF_PATH)/include ../share/global_norms_mod.F90 -o \
		../share/global_norms_mod.o

../share/gridgraph_mod.o: ../share/gridgraph_mod.F90
	$(FC) -c $(FFLAGS) ../share/gridgraph_mod.F90 -o ../share/gridgraph_mod.o

../share/kinds.o: ../share/kinds.F90
	$(FC) -c $(FFLAGS) -o ../share/kinds.o ../share/kinds.F90

../share/ll_mod.o: ../share/ll_mod.F90
	$(FC) -c $(FFLAGS) -o ../share/ll_mod.o ../share/ll_mod.F90

../share/metagraph_mod.o: ../share/metagraph_mod.F90 \
	../share/ll_mod.o
	$(FC) -c $(FFLAGS) -o ../share/metagraph_mod.o ../share/metagraph_mod.F90

../share/mesh_mod.o: ../share/mesh_mod.F90 \
	../share/coordinate_systems_mod.o ../share/gridgraph_mod.o \
	../share/cube_mod.o ../share/spacecurve_mod.o
	$(FC) -c $(FFLAGS) -I$(NETCDF_PATH)/include ../share/mesh_mod.F90 -o \
		../share/mesh_mod.o

../share/parallel_mod.o: ../share/parallel_mod.F90 \
		../share/dimensions_mod.o
	$(FC) -c $(FFLAGS) ../share/parallel_mod.F90 -o ../share/parallel_mod.o

../share/physical_constants.o: ../share/physical_constants.F90
	$(FC) -c $(FFLAGS) ../share/physical_constants.F90 -o ../share/physical_constants.o

../share/quadrature_mod.o: ../share/quadrature_mod.F90
	$(FC) -c $(FFLAGS) ../share/quadrature_mod.F90 -o ../share/quadrature_mod.o

../share/reduction_mod.o: ../share/reduction_mod.F90
	$(FC) -c $(FFLAGS) ../share/reduction_mod.F90 -o ../share/reduction_mod.o

../share/schedtype_mod.o: ../share/schedtype_mod.F90 \
	../share/metagraph_mod.o
	$(FC) -c $(FFLAGS) -o ../share/schedtype_mod.o ../share/schedtype_mod.F90

../share/spacecurve_mod.o: ../share/spacecurve_mod.F90
	$(FC) -c $(FFLAGS) ../share/spacecurve_mod.F90 -o ../share/spacecurve_mod.o

../share/thread_mod.o: ../share/thread_mod.F90
	$(FC) -c $(FFLAGS) -o ../share/thread_mod.o ../share/thread_mod.F90

$(SHR_PATH)/shr_mpi_mod.o: $(SHR_PATH)/shr_mpi_mod.F90
	$(FC) -c $(FFLAGS) -o $(SHR_PATH)/shr_mpi_mod.o $(SHR_PATH)/shr_mpi_mod.F90

$(SHR_PATH)/shr_sys_mod.o: $(SHR_PATH)/shr_sys_mod.F90 \
	$(SHR_PATH)/shr_mpi_mod.o
	$(FC) -c $(FFLAGS) -o $(SHR_PATH)/shr_sys_mod.o $(SHR_PATH)/shr_sys_mod.F90
