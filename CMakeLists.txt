PROJECT(klaplace)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

set (CMAKE_CXX_STANDARD 11)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to 'Release' as none was specified.")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
	mark_as_advanced(CMAKE_BUILD_TYPE)
	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

FIND_PACKAGE(VTK REQUIRED)
find_package(OpenMP QUIET)

if (VTK_VERSION VERSION_LESS "8.90.0")
	INCLUDE(${VTK_USE_FILE})
else ()
	if (VTK_VERSION VERSION_LESS "9.0.20210317")
		message (FATAL_ERROR "VTK v9.0.20210317 or higher is required. Your VTK version is ${VTK_VERSION}.")
	endif ()
endif ()
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

ADD_EXECUTABLE(klaplace klaplace.cpp kgeometry.cpp kvolume.cpp vtkio.cpp piOptions.cpp)

TARGET_LINK_LIBRARIES(klaplace ${VTK_LIBRARIES})
