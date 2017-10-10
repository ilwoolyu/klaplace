PROJECT(klaplace)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#cmake_policy(SET CMP0028 OLD)
#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "/Works/Codes/bin")

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(VXL REQUIRED)
INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)

#INCLUDE_DIRECTORIES(/opt/local/include/vxl/core /opt/local/include/vxl/vcl)
#LINK_DIRECTORIES(/opt/local/lib/vxl)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /home-nfs/masi-shared-home/home/local/VANDERBILT/lyui/bin)

ADD_EXECUTABLE(klaplace klaplace.cpp klaplace.h kgeometry.cpp kgeometry.h kstreamtracer.cpp kstreamtracer.h kvolume.cpp kvolume.h kvtkutils.cpp kvtkutils.h vtkio.cpp vtkio.h piOptions.h piOptions.cpp csv_parser.cpp csv_parser.h SimpleOpt.h)

TARGET_LINK_LIBRARIES(klaplace ${VTK_LIBRARIES} vnl)
