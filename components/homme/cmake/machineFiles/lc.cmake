# CMake initial cache file for cab
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")

SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")

SET (OPT_FLAGS "-O0" CACHE STRING "")
SET (DEBUG_FLAGS "-g" CACHE STRING "")
#SET (DEBUG_CFLAGS "-g" CACHE STRING "")
#SET (DEBUG_FFLAGS "-g -traceback -ftrapuv -fpe0 -check bounds -check uninit" CACHE STRING "")

SET (NETCDF_DIR $ENV{NETCDF} CACHE FILEPATH "")
SET (WITH_PNETCDF FALSE CACHE BOOL "")
SET (HDF5_DIR $ENV{HDF5} CACHE FILEPATH "")
SET (SUNDIALS_DIR $ENV{SUNDIALS} CACHE FILEPATH "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET (BUILD_HOMME_SWEQX "OFF" CACHE STRING "")
SET (BUILD_HOMME_PREQX "OFF" CACHE STRING "")
SET (BUILD_HOMME_THETA "ON" CACHE STRING "")
SET (BUILD_HOMME_PREQX_ACC "OFF" CACHE STRING "")
SET (BUILD_HOMME_PESE "OFF" CACHE STRING "")
SET (BUILD_HOMME_SWIM "OFF" CACHE STRING "")
SET (BUILD_HOMME_PRIM "OFF" CACHE STRING "")

SET (ENABLE_OPENMP FALSE CACHE BOOL "")

SET (HOMME_PROJID "climalg" CACHE STRING "")
