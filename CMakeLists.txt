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

INCLUDE(${VTK_USE_FILE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

ADD_EXECUTABLE(klaplace klaplace.cpp kgeometry.cpp kvolume.cpp vtkio.cpp piOptions.cpp)

TARGET_LINK_LIBRARIES(klaplace ${VTK_LIBRARIES})
