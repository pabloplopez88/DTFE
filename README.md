
# The DTFE public software

The DTFE public code is a C++ implementation of the **Delaunay Tessellation Field Interpolation (DTFE)** method. Its purpose is to interpolate quantities stored at the location of an unstructured set of points to a regular grid using the maximum of information contained in the input points set. In particular, the code can calculate the following cosmological quantities:
* the density field - this is calculated directly from the point distribution,
* the velocity field and derivatives (e.g. gradient, divergence, vorticity) - uses the velocity at each particle position, and
* general vector quantities and their derivatives - these quantities must be given as input for each point in the set.

The code was written with the purpose of analysing cosmological simulations and galaxy redshift survey. Even though the code was designed with astrophysics in mind, it can be used for problems in a wide range of fields where one needs to interpolate from a discrete set of points to a grid.

The code was designed using a modular philosophy and with a wide set of features that can easily be selected using the different program options. The DTFE code is also written using OpenMP directives which allow it to run in parallel on shared-memory architectures.

The code comes with a complete [documentation](documentation/DTFE_user_guide.pdf) and with a multitude of examples that detail the program features. A test dataset and analysis of the code output is given in the [demo directory](demo).

The public release of the code is summarised in the arxiv publication [Cautun et al. (2011)](https://ui.adsabs.harvard.edu/abs/2011arXiv1105.0370C/abstract) and it is based on the method paper [Schaap and van de Weygaert (2000)](https://ui.adsabs.harvard.edu/abs/2000A%26A...363L..29S/abstract).


## Compiling

With [conda](https://github.com/conda-forge/miniforge) (installs the compiler and all the libraries in your home directory):

```bash
git clone https://github.com/pabloplopez88/DTFE.git
cd DTFE
module purge                          # on clusters: do not mix cluster modules with conda
conda env create -f environment.yml   # only the first time
conda activate dtfe
make
```

`make` builds `DTFE`, which computes only the density. The velocity and magnetic fields need more memory per particle, so they are compiled separately:

| Executable | Compile with | Fields | Bytes/particle |
|---|---|---|---|
| `DTFE` | `make` | density | 20 |
| `DTFE_vel` | `make VELOCITY=yes SCALARS=no` | + velocity, divergence, shear, vorticity | 32 |
| `DTFE_mag` | `make VELOCITY=yes SCALARS=yes` | + magnetic field (SWIFT) | 44 |

To build the three:

```bash
make clean && make VELOCITY=yes SCALARS=yes && mv DTFE DTFE_mag && \
make clean && make VELOCITY=yes SCALARS=no  && mv DTFE DTFE_vel && \
make clean && make
```

The executables do not need the conda environment to run. Before reading the data, DTFE prints an estimate of the memory (RAM) the run needs.


## Input and output

* `--input 101`: Gadget-2/3/4 binary snapshots. `--input 105`: Gadget-2/3/4 and SWIFT HDF5 snapshots. For snapshots in several files use `%i` in the file name (e.g. `snap_0010.%i.hdf5`).
* The second value of `--input` selects the data (1 = positions, 2 = masses, 4 = velocities; add them, e.g. `7`) and the third one the particle types (1 = gas, 2 = dark matter, 3 = both).
* `--MpcUnit`: value of 1 Mpc in the units of the snapshot (1 for SWIFT, 1000 for kpc). Use `--periodic` for simulation boxes.
* The density is written as rho/rho_mean (dimensionless); rho_mean is printed as "Average density in the box".
* Bug fix with respect to the original DTFE: the interpolated velocity (`.vel`, `.a_vel`) was wrong (`=` instead of `+=` in `velocityValue`). The density and the velocity divergence, shear and vorticity were not affected.


## Magnetic fields (SWIFT)

With `DTFE_mag` and only the gas particles, the magnetic field of SWIFT snapshots can be computed together with the other fields:

```bash
./DTFE_mag snap_0010.hdf5 gas --input 105 7 1 --MpcUnit 1 --grid 256 --periodic \
     --field density_a velocity_a divergence_a vorticity_a magnetic_a magneticDivergence_a magneticCurl_a
```

| Field | Output | Content |
|---|---|---|
| `magnetic_a` | `.a_mag` | B |
| `magneticGradient_a` | `.a_magGrad` | dB_c/dx_i (index `c*3+i`) |
| `magneticDivergence_a` | `.a_magDiv` | div B |
| `magneticCurl_a` | `.a_magCurl` | curl B |

(also without `_a`, for the values at the grid points). Averaged fields (`_a`) with the default `--method 1` have Monte Carlo noise; `--method 2` is much more accurate.


## Running the demo

```bash
./DTFE --config demo/config_DTFE.cfg
python demo/plot_density.py
```

computes the density of a small Gadget-4 snapshot (64<sup>3</sup> particles, 50 Mpc box) on a 256<sup>3</sup> grid (`demo/demo_output.den`, raw float32) and saves a figure in `demo/demo_density.png`.


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
