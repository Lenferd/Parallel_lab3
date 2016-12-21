cmake_minimum_required(VERSION 3.2)
project(Optimization3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(WIN32)
  #This is for finding MS-MPI.
  #set(_MPI_PREFIX_PATH)
  #list(APPEND _MPI_PREFIX_PATH"[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\MPI;Path]/..")
  set(MPIEXEC "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe")
  # For building MPI programs the selected Visual Studio compiler is used,namely cl.exe.
  # So there is no need to set a specific MPI compiler.
  #set(MPI_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
  set(MPI_CXX_INCLUDE_PATH "C:/Program Files (x86)/Microsoft SDKs/MPI/Include")
  # Make sure the correct libraries (64-bit or 32-bit) are selected.
  # Decide between 32-bit and 64-bit libraries for Microsoft's MPI
  if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
    set(MS_MPI_ARCH_DIR x64)
  else()
    set(MS_MPI_ARCH_DIR x86)
  endif()
  set(MPI_CXX_LIBRARIES "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/${MS_MPI_ARCH_DIR}/msmpi.lib")
  set(MPI_C_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}")
  set(MPI_C_LIBRARIES "{${MPI_CXX_LIBRARIES}")
  set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
else()
    find_package(MPI REQUIRED)
endif()

find_package(MPI REQUIRED)
include_directories( ${MPI_INCLUDE_PATH} )

set(SOURCE_FILES src/main.cpp)
set(HEADER mpi.h)
add_executable(Optimization3 ${SOURCE_FILES} ${HEADER})
target_link_libraries( Optimization3 ${MPI_LIBRARIES} )
