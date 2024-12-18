cmake_minimum_required(VERSION 3.22)
project(TP_Final)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23) # reducable to 17 or 20 (untested)

# Find CLHEP (required, see ./install-CLHEP.sh)
find_package(CLHEP REQUIRED)

# Add compile options
add_compile_options(-Ofast)
#add_compile_options(-fno-fast-math)  # Comment and Uncomment to test bitwise reproducibility
add_compile_options(-march=native)
add_compile_options(-mtune=native)
add_compile_options(-Wall)
add_compile_options(-Wextra)
add_compile_options(-Wpedantic)

# Add include directories
include_directories(include)

# Add link libraries
link_libraries(CLHEP::Random)

# Add executable target
add_executable(GenStatus src/genmt.cpp)
add_executable(Estimation src/estimation.cpp)
