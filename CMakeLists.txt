# Specify a modern CMake version and project name
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(McMuleUserInterface LANGUAGES CXX)

# Set the C++ standard to C++17 for compatibility
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
    NO_DEFAULT_PATH
    REQUIRED
)

# --- Find the ROOT library ---
# It is expected that the ROOTSYS environment variable is set
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Find Source and Header Files
#----------------------------------------------------------------------------

file(GLOB USER_HEADERS "*.h")
set (USER_EE_SOURCE "user_ee.cpp")
set (USER_EP_SOURCE "user_ep.cpp")
set (USER_EP_SOURCE_FROM_FORTRAN "user_convert_from_fortran.cpp")

#----------------------------------------------------------------------------
# Define Shared Library Targets
#----------------------------------------------------------------------------

add_library(user_ee SHARED
    ${USER_EE_SOURCE}
    ${USER_HEADERS}
)

add_library(user_ep SHARED
    ${USER_EP_SOURCE}
    ${USER_HEADERS}
)

add_library(user_convert_from_fortran SHARED
    ${USER_EP_SOURCE_FROM_FORTRAN}
    ${USER_HEADERS}
)
#----------------------------------------------------------------------------
# Set Target Properties (Include Directories and Link Libraries)
#----------------------------------------------------------------------------

# --- Add ROOT's include directory to each target ---
target_include_directories(user_ee PRIVATE ${ROOT_INCLUDE_DIRS})
target_include_directories(user_ep PRIVATE ${ROOT_INCLUDE_DIRS})

# --- Link each library against McMule and ROOT ---
target_link_libraries(user_ee PRIVATE
    ${MCMULE_LIBRARIES}
    ${ROOT_LIBRARIES}
)
target_link_libraries(user_ep PRIVATE
    ${MCMULE_LIBRARIES}
    ${ROOT_LIBRARIES}
)
target_link_libraries(user_convert_from_fortran PRIVATE
    ${MCMULE_LIBRARIES}
)
