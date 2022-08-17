cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(cpp_project LANGUAGES CXX)

set(NBODY_FULL_NAME "Simulation and Visualization for NBody Problem")
set(NBODY_VERSION_MAJOR 1)
set(NBODY_VERSION_MINOR 0)
set(NBODY_VERSION "${NBODY_VERSION_MAJOR}.${NBODY_VERSION_MINOR}")

option(NBODY_ENABLE_DOC "Enable building of the documentation" true)

# C++ 20 is required.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Specify simulation precision, the default precision is Double.
if (SIMULATION_PRECISION)
  string(TOUPPER ${SIMULATION_PRECISION} SIMULATION_PRECISION)
else ()
  set(SIMULATION_PRECISION DOUBLE)
endif ()

# Specify the particular floating-point type for simulator.
if (SIMULATION_PRECISION STREQUAL "SINGLE")
  set(FLOATING_POINT_TYPE "float")
elseif (SIMULATION_PRECISION STREQUAL "DOUBLE")
  set(FLOATING_POINT_TYPE "double")
elseif (SIMULATION_PRECISION STREQUAL "EXTENDED")
  set(FLOATING_POINT_TYPE "long double")
else ()
  message(FATAL_ERROR "Invalid SIMULATION_PRECISION value: ${SIMULATION_PRECISION}")
endif ()

# Set FLOATING_POINT_TYPE marco.
add_compile_definitions(FLOATING_POINT_TYPE=${FLOATING_POINT_TYPE})

# Include the Sanitizers module
include(Sanitizers.cmake)

add_library(NBody INTERFACE)
target_include_directories(NBody INTERFACE include)

find_package(OpenMP REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
# Packages for documentation
find_package(Doxygen)

set(SimLinkLibs NBody OpenMP::OpenMP_CXX)
set(ConfigLinkLibs Boost::program_options)

if (APPLE)
  find_package(OpenGL REQUIRED)
  find_package(GLEW REQUIRED)
  find_package(GLFW3 REQUIRED)
  find_library(CocoaLibrary Cocoa)
  find_library(IOKitLibrary IOKit)
  find_library(CoreVideoLibrary CoreVideo)

  set(AppleLibs ${CocoaLibrary} ${IOKitLibrary} ${CoreVideoLibrary})
  set(VisLinkLibs NBody GLEW::GLEW glfw ${AppleLibs})
else ()
  find_package(OpenGL REQUIRED)
  find_package(PkgConfig REQUIRED)
  pkg_search_module(GLFW REQUIRED glfw3)

  set(VisLinkDirs ${OPENGL_INCLUDE_DIR} ${GLFW_INCLUDE_DIRS})
  set(VisLinkLibs NBody dl ${OPENGL_LIBRARIES} ${GLFW_LIBRARIES})
endif ()

set(VisSrcFiles
  lib/NBody/Visualization/GLLoader.cpp
  lib/NBody/Visualization/Camera.cpp
  lib/NBody/Visualization/Shader.cpp
  lib/NBody/Visualization/VertexArray.cpp
  lib/NBody/Visualization/VertexBuffer.cpp
  lib/NBody/Visualization/VertexBufferLayout.cpp
  lib/NBody/Visualization/IndexBuffer.cpp
  lib/NBody/Visualization/Renderer.cpp)

set(NBodyLinkDirs ${VisLinkDirs})
set(NBodyLinkLibs ${SimLinkLibs} ${VisLinkLibs} ${ConfigLinkLibs})
set(NBodySrcFiles ${VisSrcFiles})

# Executable nbody.
add_executable(nbody nbody/nbody.cpp ${NBodySrcFiles})
target_link_directories(nbody PUBLIC ${NBodyLinkDirs})
target_link_libraries(nbody ${NBodyLinkLibs})

# Install nbody program.
install(TARGETS nbody DESTINATION bin)
# Install demo and the initial data file for demo.
install(PROGRAMS demo DESTINATION bin)
install(FILES input/demo_initial_data.txt DESTINATION bin)

if (NBODY_ENABLE_DOC)
  add_subdirectory(doc)
endif ()