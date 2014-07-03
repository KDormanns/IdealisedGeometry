Idealised Geometry
==================

Generation of quad surface meshes of idealised vessel geometry using Bloor and
Wilson's PDE method. The code generates a bifurcation with one parent and two
child branches. All three branches will have the same number of grid points in
both longitudinal and circumfirential directions.

This code has a lot of hard-coded parameters and behaviours in it. In the near
future this code will be superceeded by the code with a more generic, extended
functionality in the SurfaceMeshGenerator repository:

[https://github.com/BlueFern/SurfaceMeshGenerator.git](https://github.com/BlueFern/SurfaceMeshGenerator.git)

How to Compile
==============

At the moment this code compiles only on Foster. To compile run the provided
shell script:

    build.sh
    
which will produce an executable file 'test'.

How to Run
==========

Command Line Arguments
----------------------

TODO: Instead of example values, give the default values, as they are specified in the code.

 * `E` - Number of ECs in axial direction, e.g. `4`.
 * `S` - Number of SMCs in circumferential direction, e.g. `4`.
 * `m` - Number of grid points in axial direction, e.g. `33`. The value must be
 odd and >= 3 to make sure we end up with even number of processors. If we need
 `P` cores in longitudinal direction and `m` is the number of points in axial
 direction, then `m + 2` grid points will be generated resulting in `m + 2 - 1`
 surface mesh quads each representing an MPI process/core. The extra two points
 are the boundary points. 
 * `n` - Number of grid points in circumferential direction for a half-circle,
 e.g. `15`. The value must be odd and >= 3 to make sure we end up with even
 number of processors. Every cyllinder will be generated in two parts, one for
 each side, as if the cyllinder was split along its long axis. If we need `Q`
 cores in circumferential direction and `n` is the number of points in the
 hemicyllinder, then `n + 2` points will be generated for a hemicyllinder
 resulting in `n + 2 - 1` (which is `Q / 2`) surface mesh quad each representing
 an MPI process/core.
 * `D` - Downstream skip, e.g. `0`. This parameter is the outlet skip on the
 'daugter' branches specifying the number of additional points in the mesh which
 will only appear in the STL file (triangulated version of the mesh) to avoid
 boundary condition artifacts in the CFD solution.
 * `U` - Upstream skip, e.g. `0`. This parameter is the inlet skip on the
 'parent' branches specifying the number of points in the mesh which will only
 appear in the STL file (triangulated version of the mesh) to aviod boundary
 condition artifacts in the CFD solution.
 * `l` - Length of a branch in mm, e.g. `25`. Must be langer than radius `r` and
 least 10 times greater than `r` to provide a correct CFD solution.
 * `r` - Raduis r in mm, e.g. 2.5.

For example, to produce a dataset that fits into 4,032 cores and three branches,
we end up with 1,344 cores per branch, which meets the conditions on `n` and
`m`.

Start by evaluating the total number of ECs and SCMs in the longitudinal and
circumferential direction respectively. The number of ECs and SMCs per core are
calculated based on the values of `m` and `n`.

 1. `Total number of SCMs = 2 pi R / length of an SMC`. Hard-coded length of an
 SMC is 50e-6; hard-coded width of an SMC is 5e-6.
 2. `Total number of ECs = l / length of an EC`. Hard-coded length of an EC is
 65e-6; hard-coded width is 10e-6.
 3. Adjust `m` and `n` such that the resulting number of quads is 1,344,
 provided the conditions on `m` and `n` are met.
 4. Numer of ECs per core in axial direction is defined as follows:
 `E = Total number of ECs / (m + 2 - 1)`.
 6. Number of SMCs per core in circumferential direction is devfined as follows:
 `S = Total number of SMCs / ((n + 2 - 1) * 2)`.

In this case the values are as follows:

    ./test -m 13 -n 47 -E 27 -S 3 -r 2.5 -l 25


Output
======

For each branch the output includes VTK files, stereo-lithography STL files, and
TXT files. The VTK files are for imediate visualisation/visual
debugging/validation. The STL files are used for CFD solution generation. The
TXT flies are the intput into the simulation code.

TODO: Add a note abot running the code in a separate directory otherwise the
generate flies mess with git's head and make it harder to make changes.



