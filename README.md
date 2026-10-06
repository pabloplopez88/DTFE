
# The DTFE public software

The DTFE public code is a C++ implementation of the **Delaunay Tessellation Field Interpolation (DTFE)** method. Its purpose is to interpolate quantities stored at the location of an unstructured set of points to a regular grid using the maximum of information contained in the input points set. In particular, the code can calculate the following cosmological quantities:
* the density field - this is calculated directly from the point distribution,
* the velocity field and derivatives (e.g. gradient, divergence, vorticity) - uses the velocity at each particle position, and
* general vector quantities and their derivatives - these quantities must be given as input for each point in the set.

The code was written with the purpose of analysing cosmological simulations and galaxy redshift survey. Even though the code was designed with astrophysics in mind, it can be used for problems in a wide range of fields where one needs to interpolate from a discrete set of points to a grid.

The code was designed using a modular philosophy and with a wide set of features that can easily be selected using the different program options. The DTFE code is also written using OpenMP directives which allow it to run in parallel on shared-memory architectures.

The code comes with a complete [documentation](documentation/DTFE_user_guide.pdf) and with a multitude of examples that detail the program features. A test dataset and analysis of the code output is given in the [demo directory](demo).

The public release of the code is summarised in the arxiv publication [Cautun et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011arXiv1105.0370C/abstract) and it is based on the method paper [Schaap and van de Weygaert (2000)](https://ui.adsabs.harvard.edu/abs/2000A%26A...363L..29S/abstract).


## Compiling (new!)

The easiest way to compile the code on any Linux computer or cluster is with [conda](https://docs.conda.io/en/latest/miniconda.html), which installs the compiler and all the required libraries (GSL, Boost, CGAL, GMP, MPFR and HDF5) in your home directory, without needing administrator rights or cluster modules:

```bash
git clone https://github.com/pabloplopez88/DTFE.git
cd DTFE
module purge                          # only on clusters with environment modules (see notes below)
conda env create -f environment.yml   # only the first time on each computer
conda activate dtfe
make
```

This produces the `DTFE` executable. The library paths are stored inside the executable, so it runs without setting `LD_LIBRARY_PATH` (you do not even need to activate the conda environment to run it, e.g. inside a SLURM job).

Notes:
* Compile-time options (number of dimensions, velocity/scalar fields, default input/output formats, etc.) are set at the top of the `Makefile`.
* Memory: by default each particle uses 44 bytes (velocity + 3 extra components used for the magnetic field). If you do not need them, compile with `make VELOCITY=no SCALARS=no` (20 bytes per particle, e.g. ~310 GB instead of ~690 GB for 2500^3 particles) or `make SCALARS=no` (32 bytes, no magnetic field). The program stops with an error if you ask for a field that was not compiled. Run `make clean` before compiling with different options; to keep several versions, rename the executable after each build (e.g. `make clean; make VELOCITY=no SCALARS=no; mv DTFE DTFE_density`).
* Support for HDF5 snapshots is enabled automatically when the HDF5 C++ library is found. Use `make USE_HDF5=no` to disable it.
* To use libraries installed elsewhere instead of conda (e.g. cluster modules), give their location: `make LIB_PREFIX=/path/to/prefix`, or each one separately with `GSL_PATH`, `BOOST_PATH`, `CGAL_PATH`, `MPRF_PATH` (GMP and MPFR) and `HDF5_PATH`. The compiler can be chosen with `make CXX=g++`.
* The library paths are taken from the conda environment (or from the `make` command line), never from shell variables such as `HDF5_PATH` that some clusters define, and the compiled program uses those libraries even if `LD_LIBRARY_PATH` points to other versions of them.
* On clusters, do not mix cluster modules with the conda environment: run `module purge` before `conda activate dtfe` and `make` (a loaded `gcc` module can make the conda compiler fail with `cannot execute 'cc1plus'`). Also use `module purge` in SLURM job scripts before running `DTFE`.
* Input formats (`--input`): `101` reads Gadget binary snapshots (formats 1 and 2) with either the Gadget-1/2 header or the Gadget-4 header (the default of Gadget-4 when compiled without `GADGET2_HEADER`): the program tries the Gadget-2 header first, then the Gadget-4 one, and stops with an error only if none matches. Positions and velocities can be in single or double precision, and files with the opposite endianness are also supported. `105` reads HDF5 snapshots from Gadget-2/3, Gadget-4 and SWIFT (SWIFT lengths are in Mpc without h, so use `--MpcUnit 1`; particle types 6 and higher are ignored). For snapshots split in several files put `%i` in the file name where the file number goes (e.g. `snapshot_%i.bin` or `snap_0010.%i.hdf5`); a SWIFT "virtual" snapshot file can also be given directly. The second value of `--input` selects the data to read (1 = positions, 2 = masses, 4 = velocities; e.g. `--input 105 7`) and the third one the particle types (1 = type 0, 2 = type 1, 4 = type 2, ...; e.g. `--input 105 7 3` for gas + dark matter).
* Bug fix with respect to the original DTFE code: the interpolation of the velocity (and of the scalar fields) to the grid points used only one of the three terms of the linear interpolation (`=` instead of `+=` in `velocityValue` and `scalarValue`), so the `velocity` / `velocity_a` fields (`.vel`, `.a_vel` files) were wrong and depended on the order of the particles. The velocity gradient and the fields derived from it (divergence, shear, vorticity) and the density were not affected.
* Tested with CGAL 5.6 and 6.2, Boost 1.83, 1.90 and 1.92, HDF5 1.10 and 2.2, and GCC 13 and 15.


## Magnetic fields (SWIFT)

For SWIFT HDF5 snapshots with magnetohydrodynamics, DTFE can interpolate the magnetic field of the gas particles (`/PartType0/MagneticFluxDensities`) in the same run as the density and the velocity fields. Add any of the following to `--field` (or `field = ...` lines in the configuration file):

| Field | Output file | Content |
|---|---|---|
| `magnetic` / `magnetic_a` | `.mag` / `.a_mag` | B (3 components per grid cell) |
| `magneticGradient` / `magneticGradient_a` | `.magGrad` / `.a_magGrad` | dB_c/dx_i, stored as index `c*3+i` (9 components) |
| `magneticDivergence` / `magneticDivergence_a` | `.magDiv` / `.a_magDiv` | div B |
| `magneticCurl` / `magneticCurl_a` | `.magCurl` / `.a_magCurl` | curl B = (dBz/dy-dBy/dz, dBx/dz-dBz/dx, dBy/dx-dBx/dy) |

The magnetic field is read automatically when one of these fields is requested. Since it exists only for the gas, select only the gas particles, e.g.

```bash
./DTFE snap_0010.hdf5 gas --input 105 7 1 --MpcUnit 1 --grid 256 --periodic \
     --field density_a velocity_a divergence_a vorticity_a magnetic_a magneticDivergence_a magneticCurl_a
```

(the dark matter fields need a separate run with `--input 105 7 2`). Internally the magnetic field uses the 3 components of the `scalar` data (`NO_SCALARS=3` in the `Makefile`), so it cannot be combined with the `scalar` fields in the same run.

Note on the averaged fields (`*_a`): with the default averaging method (`--method 1`) the averages are computed by Monte Carlo sampling inside the Delaunay cells and have a sampling noise that decreases with `--samples` (default 100). Method 2 (`--method 2`) samples points inside each grid cell: its noise is much smaller (in a test with linear fields and 100 samples the error of the averaged velocity and magnetic field was ~30 times smaller than with method 1) and the derived fields that are constant inside the Delaunay cells, such as the divergence and the curl, come out exact.


## Running the demo

The [demo](demo) directory contains a small Gadget-4 HDF5 snapshot (64<sup>3</sup> dark matter particles in a periodic box of 50 Mpc, positions in Mpc) and a configuration file with all the options needed to compute its density field on a 256<sup>3</sup> grid. From the main directory of the repository run:

```bash
./DTFE --config demo/config_DTFE.cfg
```

This writes `demo/demo_output.den`: 256<sup>3</sup> single-precision floats (no header), which can be read in Python with

```python
import numpy as np
n = 256
rho = np.fromfile('demo/demo_output.den', dtype=np.float32).reshape(n, n, n)
delta = rho / rho.mean()      # density in units of the mean density
```

To make a quick figure of the result (a slab through the middle of the box and a projection along the full box), run, with the `dtfe` environment active (it includes `numpy` and `matplotlib`):

```bash
python demo/plot_density.py
```

The figure is saved to `demo/demo_density.png` (and also shown on screen if there is a display). In a Jupyter notebook you can use `%run demo/plot_density.py`. Use `python demo/plot_density.py --help` for the options.

The configuration file is equivalent to the command line
`./DTFE demo/gadget4_L50_N64_snap001.hdf5 demo/demo_output --input 105 3 --MpcUnit 1 --grid 256 --field density --periodic`.
Note that `--periodic` is important for simulation boxes: without it, the cells close to the box edges are left empty.


## The DTFE method
The Delaunay Tessellation Field Interpolation (DTFE) method represents the natural way of going from discrete samples/measurements to values on a periodic grid and it is especially suitable for astronomical data due to the following reasons:
* Preserves the multi-scale character of the point distribution. This is the case in numerical simulations of large scale structure where the density varies over more than 6 orders of magnitude.
* Preserves the local geometry of the point distribution. This is important in recovering sharp features like the different components of the cosmic web (i.e. clusters, filaments, walls and voids).
* The method does not depend on user defined parameters or choices.
* The interpolated fields are volume weighted (versus mass weighted quantities in most other interpolation schemes). This can have a significant effect especially when comparing with analytical predictions which are volume weighted.

For detailed information about the DTFE method see [Schaap and van de Weygaert (2000)](https://ui.adsabs.harvard.edu/abs/2000A%26A...363L..29S/abstract), [van de Weygaert and Schaap (2009)](https://ui.adsabs.harvard.edu/abs/2009LNP...665..291V/abstract), and [Cautun et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011arXiv1105.0370C/abstract).

| <img src="figures/DTFE_filament.png" width="600" title="An illustration of the adaptive nature of the DTFE method."> |
|:------:|
| Figure 1: *An illustration of the 2D Delaunay tessellation of a set of particles from a cosmological simulation. Courtesy: Willem Schaap.* |

| <img src="figures/DTFE_paper_density.png" width="800" title="Examples of DTFE density fields."> |
|:------:|
| Figure 2: *An example of the DTFE density field form a cosmological simulation. The right panel shows the same result but now using the smoothed particle hydrodynamics (SPH) method.* |

| <img src="figures/DTFE_paper_velocity.png" width="800" title="Illustration of the DTFE velocity field."> |
|:------:|
| Figure 3: *A map of the DTFE computed velocity flow (left panel) and velocity divergence (right panel) corresponding to the density field shown in Figure 2.* |


## Summary of software features

* Works in both 2 and 3 spatial dimensions.
* Interpolates the fields to three different types of grids:
  + Regular rectangular and cuboid grid - useful for cosmological simulation.
  + Redshift cone (spherical coordinates) grid - useful for galaxy redshift survey or for mock observations.
  + User given sampling points - can describe any complex or non-regular sampling geometry
* Returns both the value at the centre of each cell of the interpolation grid as well as the value averaged over each cell.
* Uses the point distribution to compute the density and interpolates the result to grid.
* Each sample point has a weight associated to it to represent multiple resolution N-body simulations and observational biases for galaxy redshift surveys.
* Interpolates the velocity, velocity gradient, velocity divergence, velocity shear and velocity vorticity.
* Interpolates any additional number of fields and their gradients to grid.
* Periodic boundary conditions.
* Zoom in option for regions of interest.
* Splitting the full data in smaller computational chunks when dealing with limited CPU resources.
* The computation can be distributed in parallel on shared-memory architectures.
* For comparison purposes, the software comes also with three other simpler interpolation techniques: nearest grid point (NGP), triangular shape cloud (TSC; ) and smoothed particle hydrodynamics (SPH; )\citep{1992ARA&A..30..543M}.
* Returns the Delaunay tessellation of the given point set.
* Easy change of input/output data format.
* Easy to use as an external library.
* Extensive documentation of each feature.


## Contributors
* **Marius Cautun (Kapteyn Astronomical Institute, Durham University, Leiden University)** - *code and documentation.*
* **Rien van de Weygaert (Kapteyn Astronomical Institute)** - *various discussions about the method and implementation.*


## License

This project is licensed under GNU GENERAL PUBLIC LICENSE Version 3 - see the [LICENSE.md](LICENSE.md) file for details.
