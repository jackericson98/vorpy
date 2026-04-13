cmake_minimum_required(VERSION 3.18)
project(vorpy_native LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)

# ---- _calc ---------------------------------------------------------
pybind11_add_module(_calc
  vorpy/_native/_calc.cpp
  vorpy/_native/core/calcs.cpp
)

target_include_directories(_calc PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/vorpy/_native
)

target_compile_features(_calc PRIVATE cxx_std_17)

# ---- _vert_cpp -----------------------------------------------------
pybind11_add_module(_vert_cpp
  vorpy/_native/vert_bindings.cpp
  vorpy/_native/core/vert_core.cpp
)

target_include_directories(_vert_cpp PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/vorpy/_native
)

target_compile_features(_vert_cpp PRIVATE cxx_std_17)

# ---- install both into the package --------------------------------
install(TARGETS _calc _vert_cpp
  LIBRARY DESTINATION vorpy/_native
)
