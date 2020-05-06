project(Mex)


# Setup MEX compilation, as in gptoolbox
set( CMAKE_SHARED_LIBRARY_SUFFIX ".${Matlab_MEX_EXTENSION}" )  # set suffix to .mexa64
set( CMAKE_SHARED_LIBRARY_PREFIX )      # remove the "lib" prefix
set( CMAKE_SHARED_LIBRARY_CXX_FLAGS )       # remove the -fPIC option. mex does not accept the "-fPIC" option
set( CMAKE_POSITION_INDEPENDENT_CODE ON)
set( CMAKE_INSTALL_RPATH "\$ORIGIN" )

# Pick a name
get_filename_component(MEX_OUTPUT "curved_hessian" NAME_WE)

# Compile
file(GLOB SRCFILES *.cpp)
list(APPEND SRCFILES ${GLOBAL_CPPS})
add_definitions(-DPROJECT_NAME="${PROJECT_NAME}")
add_library(${MEX_OUTPUT} SHARED ${SRCFILES})
target_link_libraries(${MEX_OUTPUT} ${LINKED_LIBS})


configure_file("curved_hessian.m" "${CMAKE_CURRENT_BINARY_DIR}/curved_hessian.m" COPYONLY)
