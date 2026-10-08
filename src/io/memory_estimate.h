/*
 *  Copyright (c) 2026       Pablo Lopez
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* Estimate of the peak memory (RAM) needed by a DTFE run.

The estimate is printed as soon as the number of particles is known (before reading the particle data), so the user can check that
the run fits in the available memory. It is an empirical model calibrated against the measured peak memory of DTFE runs with
2.6e5 - 2.1e6 particles, grids of 64^3 - 192^3, 1 - 27 threads and 1 - 512 partitions (29 runs, agreement within -15% / +18%):

    memory = base + max( compute, read )
    read    = N * (bytes read per particle + sizeof(Particle_data))               (the input arrays and the particle vector coexist)
    compute = N * sizeof(Particle_data) * m                                        (m=2 for a periodic box without '--partition')
            + Ntri * (c0 + c1 * sizeof(Particle_data))                             (Delaunay triangulations computed at the same time)
            + kt * Ntri * sizeof(Particle_data) + at * threads   [if threads > 1] (particle copies and memory pools of each thread)
            + grid cells * 4 bytes * (number of field components, including the temporary gradients)

Ntri is the number of particles in the triangulations computed at the same time: the particles of one partition split among the
threads, each thread with its own padding of 'paddingParticles' mean inter-particle distances on each side.
*/

#ifndef MEMORY_ESTIMATE_HEADER
#define MEMORY_ESTIMATE_HEADER

#include <cmath>
#include <cstdio>
#include <vector>
#ifdef OPEN_MP
    #include <omp.h>
#endif


namespace MemoryEstimate
{
    // calibration constants (see above)
    static const double c0   = 589.;          // bytes per triangulation vertex (CGAL vertex and cells)
    static const double c1   = 0.13;          // copies of the particle data per triangulation vertex
    static const double kt   = 1.71;          // copies of the particle data per triangulation vertex when using several threads
    static const double at   = 7.2*1048576.;  // memory per thread (OpenMP thread memory pools)
    static const double base = 26.*1048576.;  // fixed memory (program and libraries)

    // split of the threads along the axes (same algorithm as 'parallelGrid' in 'subpartition.h')
    inline void threadGrid(size_t const noThreads, size_t grid[3])
    {
        size_t n = noThreads;
        int x2 = 0, x3 = 0;
        while ( n/2!=0 ) { n /= 2; ++x2; }
        n = noThreads;
        while ( n/3!=0 ) { n /= 3; ++x3; }
        n = noThreads;
        size_t diff = n; int j2 = 0, j3 = 0; size_t nApprox = 1;
        for (int i2=0; i2<=x2; ++i2)
        {
            size_t temp = nApprox;
            for (int i3=0; i3<=x3; ++i3)
            {
                if ( temp<=n and (n-temp)<diff ) { diff = n-temp; j2 = i2; j3 = i3; }
                temp *= 3;
            }
            nApprox *= 2;
        }
        x2 = j2; x3 = j3;
        for (int i=0; i<3; ++i) grid[i] = 1;
        while ( x3-3>=0 ) { for (int i=0; i<3; ++i) grid[i] *= 3; x3 -= 3; }
        for (int i=2; i>=(3-x3); --i) grid[i] *= 3;
        int const i1 = x2>(3-x3) ? (3-x3) : x2;
        for (int i=0; i<i1; ++i) grid[i%(3-x3)] *= 2;
        for (int i=0; i<x2-i1; ++i) grid[i%3] *= 2;
    }

    // number of float values stored per grid cell for the requested fields (including the temporary gradients used to compute derived fields)
    inline size_t fieldComponents(Field &f)
    {
        size_t n = 0;
        if ( f.density ) n += 1;
        if ( f.velocity ) n += noVelComp;
        if ( f.velocity_gradient or f.selectedVelocityDerivatives() ) n += noGradComp;
        if ( f.velocity_divergence ) n += 1;
        if ( f.velocity_shear ) n += noShearComp;
        if ( f.velocity_vorticity ) n += noVortComp;
        if ( f.velocity_std ) n += 1;
        if ( f.scalar or f.magnetic ) n += noScalarComp;
        if ( f.scalar_gradient or f.selectedMagneticDerivatives() ) n += noScalarGradComp;
        if ( f.magnetic_divergence ) n += 1;
        if ( f.magnetic_curl ) n += 3;
        return n;
    }

    // returns the estimated peak memory in bytes
    inline double estimate(size_t const noParticles, User_options &opt, size_t *noThreadsOut, size_t *noPartitionsOut)
    {
        double const N = double(noParticles);
        double const PD = double( sizeof(Particle_data) );

        // bytes read from the input file per particle
        double readBytes = NO_DIM*sizeof(Real);                                   // positions
        if ( opt.readParticleData.size()>1 and opt.readParticleData[1] ) readBytes += sizeof(Real);          // masses
#ifdef VELOCITY
        if ( opt.readParticleData.size()>2 and opt.readParticleData[2] ) readBytes += noVelComp*sizeof(Real); // velocities
#endif
#ifdef SCALAR
        bool scalarRead = opt.uField.selectedMagnetic() or opt.aField.selectedMagnetic();
        for (size_t i=3; i<opt.readParticleData.size(); ++i) if ( opt.readParticleData[i] ) scalarRead = true;
        if ( scalarRead ) readBytes += noScalarComp*sizeof(Real);
#endif

        // partitions and threads
        size_t part[3] = {1,1,1};
        if ( opt.partitionOn and opt.partition.size()==3 )
            for (int i=0; i<3; ++i) part[i] = opt.partition[i];
        size_t noThreads = 1;
#ifdef OPEN_MP
        noThreads = omp_get_max_threads();
#endif
        size_t thr[3] = {1,1,1};
        if ( noThreads>1 ) threadGrid( noThreads, thr );
        *noThreadsOut = noThreads;
        *noPartitionsOut = part[0]*part[1]*part[2];

        // particles in the triangulations computed at the same time
        double const n1 = std::pow( N, 1./3. );     // mean number of particles along each axis
        double const pad = opt.paddingParticles>0. ? double(opt.paddingParticles) : 0.;
        double Ntri = N / double(part[0]*part[1]*part[2]);
        for (int i=0; i<3; ++i)
        {
            double const side = n1 / double(part[i]*thr[i]);   // size of the region of each thread in mean inter-particle distances
            Ntri *= (side + 2.*pad) / side;
        }

        // grid cells
        double cells = 1.;
        for (size_t i=0; i<opt.gridSize.size(); ++i) cells *= double(opt.gridSize[i]);
        if ( opt.gridSize.empty() ) cells = 0.;
        double const gridBytes = cells * sizeof(Real) * double( fieldComponents(opt.uField) + fieldComponents(opt.aField) );

        double const m = ( opt.periodic and not opt.partitionOn ) ? 2. : 1.;
        double compute = N*PD*m + Ntri*(c0 + c1*PD) + gridBytes;
        if ( noThreads>1 ) compute += kt*Ntri*PD + at*double(noThreads);
        double const read = N*(readBytes + PD);
        return base + (compute>read ? compute : read);
    }

    // prints the estimate (only once per run)
    inline void print(size_t const noParticles, User_options &opt)
    {
        static bool printed = false;
        if ( printed or not opt.DTFE or noParticles==0 ) return;
        printed = true;
        size_t noThreads, noPartitions;
        double const bytes = estimate( noParticles, opt, &noThreads, &noPartitions );
        double const GB = bytes / 1.e9;

        char line[300];
        MESSAGE::Message message( opt.verboseLevel );
        message << "\n================================================================================\n";
        if ( GB>=1. ) std::snprintf( line, sizeof(line), "  ESTIMATED MEMORY (RAM) NEEDED:  ~%.1f GB", GB );
        else std::snprintf( line, sizeof(line), "  ESTIMATED MEMORY (RAM) NEEDED:  ~%.2f GB", GB );
        message << line << "\n";
        std::snprintf( line, sizeof(line), "  (%zu particles, %zu bytes/particle, %zu partition(s), %zu thread(s); typically within 20%%)",
                       noParticles, sizeof(Particle_data), noPartitions, noThreads );
        message << line << "\n";
        message << "  To reduce it: use more partitions ('--partition'), fewer threads (OMP_NUM_THREADS)\n"
                << "  or compile without unused fields ('make SCALARS=no' or 'make VELOCITY=no SCALARS=no').\n";
        message << "================================================================================\n\n" << MESSAGE::Flush;
    }
}

#endif
