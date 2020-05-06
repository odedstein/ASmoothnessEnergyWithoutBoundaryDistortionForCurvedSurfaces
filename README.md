# A Smoothness Energy without Boundary Distortion for Curved Surfaces

This is the code accompanying the paper ["A Smoothness Energy without Boundary Distortion for Curved Surfaces"](http://www.cs.columbia.edu/cg/curved-hessian/). It is released under the MPLv2 license (included in this repository).
Acknowledgement for the meshes used in the examples can be found in the "Acknowledgements" section of the paper.

The following software is needed to run the code:
- CMake with a C++14 compiler
- [gptoolbox](https://github.com/alecjacobson/gptoolbox) version as of January 10, 2020 (make sure to add it to your MATLAB path)
- MATLAB
- OpenGL

The code is intended to be run on macOS X 10.15.4 with the clang 11.0.3 and MATLAB_R2020a. This code relies heavily on the [libigl](https://github.com/libigl/libigl) library as of January 10, 2020. It might run on other configurations too, but it has not been tested on any of them.


In order to install and compile the code, you have to
- Recursively clone this repository (git clone --recursive).
- Go to into the main directory of the cloned repository.
- cd cpp_interface
- mkdir build; cd build; cmake ..
- Use ccmake to choose between Release and Debug compilation modes.
- make

The result of this compilation is a MEX file which can be used in MATLAB to build the curved Hessian matrix in cpp_interface/build/applications/Mex .
The examples from the paper are separate MATLAB scripts in the directory MATLAB_experiments.


## Known issues and troubleshooting

### CMake is unable to find my MATLAB version
The included CMake script is unable to find newer MATLAB versions. Look at its [documentation](https://github.com/libigl/libigl/blob/master/cmake/FindMATLAB.cmake) to see how to include newer MATLAB versions, specifically MATLAB_ADDITONAL_VERSIONS
```
# .. variable:: MATLAB_ADDITIONAL_VERSIONS
#
#   If set, specifies additional versions of Matlab that may be looked for.
#   The variable should be a list of strings, organised by pairs of release
#   name and versions, such as follows::
#
#     set(MATLAB_ADDITIONAL_VERSIONS
#         "release_name1=corresponding_version1"
#         "release_name2=corresponding_version2"
#         ...
#         )
#
#   Example::
#
#     set(MATLAB_ADDITIONAL_VERSIONS
#         "R2013b=8.2"
#         "R2013a=8.1"
#         "R2012b=8.0")
```
