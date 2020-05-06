cmake_minimum_required(VERSION 3.7)
project(curvature-hessian)
set(CMAKE_CXX_STANDARD 14)

# Max OS X settings
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13)

# Find libraries
list(APPEND CMAKE_MODULE_PATH
    "${PROJECT_SOURCE_DIR}/libigl/cmake")
set(MATLAB_ADDITIONAL_VERSIONS "R2019b=9.7" "R2018b=9.5" "R2020a=9.8")
find_package(MATLAB)

# libigl
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         "${MATLAB_FOUND}")
#option(LIBIGL_WITH_MOSEK             "Use MOSEK"          "${MOSEK_FOUND}")
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_NANOGUI           "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI"          OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML               "Use XML"            OFF)
include(${PROJECT_SOURCE_DIR}/libigl/cmake/libigl.cmake)
set(LINKED_LIBS igl::core)
if(LIBIGL_WITH_OPENGL_GLFW)
    set(LINKED_LIBS ${LINKED_LIBS} igl::opengl_glfw)
endif()
if(LIBIGL_WITH_MATLAB)
    link_directories("${Matlab_MEX_COMPILER}/../")
    add_definitions(-DLIBIGL_MATLAB_ENABLED)
    set(LINKED_LIBS ${LINKED_LIBS} igl::matlab mx mex mat)
endif()
if(LIBIGL_WITH_MOSEK)
    set(LINKED_LIBS ${LINKED_LIBS} igl::mosek)
endif()
if(LIBIGL_WITH_TRIANGLE)
set(LINKED_LIBS ${LINKED_LIBS} igl::triangle)
endif()

# Static library
if(LIBIGL_USE_STATIC_LIBRARY)
    file(GLOB SOURCES_HESSIAN "${PROJECT_SOURCE_DIR}/include/hessian/*.cpp")
    compile_igl_module("hessian" ${SOURCES_HESSIAN})
    set(LINKED_LIBS ${LINKED_LIBS} igl::hessian)
endif()

# Store location of data directory
set(DATA_PATH "${PROJECT_SOURCE_DIR}/data" CACHE PATH "location of data")
add_definitions(-DDATA_PATH="${DATA_PATH}/")

# include code
include_directories("include" ${LIBIGL_INCLUDE_DIRS})

# Add your projects
add_subdirectory("applications")

