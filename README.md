# fmm2dbie

## FMM-accelerated boundary integral equation solvers in two dimensions

Currently supports 16th order discretization of smooth curves in the
go2 format

Upcoming support for: 
-  Lower order versions of the above routines
-  Support for other types of input format


This repository has an external dependency - [fmm2d](https://github.com/flatironinstitute/fmm2d)

Make sure you have the shared object for the FMM library installed and
located in an appropriate location (`/usr/local/lib` on MacOSX, and
environment variable of LD_LIBRARY_PATH set to location of libfmm3d.so 
on linux machines)

