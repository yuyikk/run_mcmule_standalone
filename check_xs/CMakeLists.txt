# Specify a modern CMake version and project name
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(McMuleUserInterface LANGUAGES CXX)

# Set the C++ standard to C++17 for compatibility
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_language(Fortran)
#----------------------------------------------------------------------------
# Find External Packages
#----------------------------------------------------------------------------

# --- Find the McMule library ---
# It is expected that the MCMULE_PATH environment variable is set
set(MCMULE_LIB_PATH $ENV{MCMULE_PATH})
if(NOT MCMULE_LIB_PATH)
    message(FATAL_ERROR "The MCMULE_PATH environment variable is not set. Please set it to the directory containing libmcmule.so.")
endif()
find_library(MCMULE_LIBRARIES
    NAMES mcmule
    PATHS ${MCMULE_LIB_PATH}
    REQUIRED
)

# --- Find the ROOT library ---
# It is expected that the ROOTSYS environment variable is set
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
# find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Find Source and Header Files
#----------------------------------------------------------------------------

file(GLOB USER_HEADERS "*.h")
set (USER_SOURCE "user.cpp" "io.f95")

#----------------------------------------------------------------------------
# Define Shared Library Targets
#----------------------------------------------------------------------------

add_library(user SHARED
    ${USER_SOURCE}
    ${USER_HEADERS}
)

# --- Link each library against McMule and ROOT ---
target_link_libraries(user PRIVATE
    ${MCMULE_LIBRARIES}
)
