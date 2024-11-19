# CMakeList.txt : CMake project for nummeth2017, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable(nummeth2017 main.cpp
        method.cpp
        mesh.cpp
        MethodGas.cpp
        MethodHeat.cpp)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET nummeth2017 PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
