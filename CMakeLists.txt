cmake_minimum_required(VERSION 3.10)
project(TensorCalculus)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Find libraries
find_library(SYMENGINE_LIBRARY symengine REQUIRED)
include_directories(/usr/local/include)

# Add executable
add_executable(tensor_calculus src/main.cpp src/tensor.cpp src/tensor_ops.cpp)

# Link libraries
target_link_libraries(tensor_calculus ${SYMENGINE_LIBRARY} gmp mpfr pthread)
