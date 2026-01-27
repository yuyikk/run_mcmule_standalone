cmake_minimum_required(VERSION 3.10)
project(user_prad_lib Fortran)

# Set Fortran standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Path to any external shared library you want to link
# Replace with the actual path if needed
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
# Include directory if library has headers (optional)
# include_directories(${MYLIB_DIR}/include)

# Add the shared library target
add_library(user_ep SHARED user_ep.f95)

# Link external shared library if needed
target_link_libraries(user_prad PRIVATE
   ${MCMULE_LIBRARIES}
)