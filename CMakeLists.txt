cmake_minimum_required(VERSION 2.6)
project(IdealisedGeometry)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
		"Choose the type of build, options are: None Debug Release."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Location of malloc.h on OS X. 
	include_directories(/usr/include/malloc/)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

enable_language(Fortran)
if(CMAKE_Fortran_COMPILER)
  message(STATUS "Found Fortran compiler: " ${CMAKE_Fortran_COMPILER})
else()
  message(FATAL_ERROR "No Fortran support. Unable to continue")
endif()

add_subdirectory(srcLib/blas)
add_subdirectory(srcLib/linpack)
add_subdirectory(srcLib/dbihar)

add_executable(IdealisedGeometry src/functionality.h src/functions.cpp src/main.cpp)
target_link_libraries(IdealisedGeometry
	dbihar
	linpack
	blas
)

