/*
 *  Copyright (c) 2011       Marius Cautun
 *
 *                           Kapteyn Astronomical Institute
 *                           University of Groningen, the Netherlands
 *
 *
 *  Modified in 2026 by Pablo Lopez (https://github.com/pabloplopez88/DTFE); see the git history for the changes.
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

#include <vector>
#include <cmath>
#ifdef OPEN_MP
    #include <omp.h>
#endif
#include <boost/math/special_functions/fpclassify.hpp>


#include "define.h"
#include "particle_data.h"
#include "user_options.h"
#include "quantities.h"
#include "subpartition.h"
#include "miscellaneous.h"
#include "message.h"


using namespace std;
#include "user_options.cc"
#include "quantities.cc"

#include "NGP_interpolation.cc"
#include "CIC_interpolation.cc"
#include "TSC_interpolation.cc"
#include "SPH_interpolation.cc"
#include "random.cc"




// splits the total data into different computation regions that are done by different threads (if OPEN_MP is enabled)
void DTFE_parallel(vector<Particle_data> *allParticles,
                   vector<Sample_point> &samples,
                   User_options & userOptions,
                   Quantities *uQuantities,
                   Quantities *aQuantities);

// Computes the velocity divergence, shear and vorticity from the velocity gradient
void computeDivergenceShearVorticity(Field &fields,
                                     int const verboseLevel,
                                     Quantities *quantities);
void computeMagneticDivergenceCurl(Field &fields,
                                   Quantities *q);



/* This function interpolates fileds to a grid using the DTFE method.
It takes the following arguments:
        allParticles - stores the particle positions, velocities, scalars (if any) and weights
        samples - stores the position of the grid points where samples are taken (if this grid points don't have a uniform distribution)
        userOptions - options specified by the user plus additional program wise variables
        uQuantities - structure that will store the output quantities - i.e. the fields at the sampling points (NOTE: "uQuantities" stands for "unaveraged quantities" meaning that this fields were NOT averaged over the sampling cell)
        aQuantities - structure that will store the volume averaged output quantities - i.e. the fields volume averaged over the sampling cells (NOTE: "aQuantities" stands for "averaged quantities")

NOTE: This function clears the vector 'allParticles'. */
void DTFE(vector<Particle_data> *allParticles,
          vector<Sample_point> &samples,
          User_options & userOptions,
          Quantities *uQuantities,
          Quantities *aQuantities)
{
    //! If the user requested for a set of random particles, generate them
    if ( userOptions.poisson!=0 )
        randomParticles( allParticles, &userOptions );
    
    
    //! Select a random subset of the data if the user asked so
    vector<Particle_data> *particlePointer = allParticles; // pointer to the vector containing the particle data
    vector<Particle_data> particlesRandomSubsample;  // vector for storing the random subsample particles
    if (userOptions.randomSample>=Real(0.) )
    {
        randomSample( *particlePointer, &particlesRandomSubsample, userOptions );
        
        // update 'particlePointer' -> should point to 'particlesRandomSubsample'
        particlePointer->clear();   //clear the particle positions - don't need them anymore
        particlePointer = &particlesRandomSubsample;
    }
    
    
    //! Do consistency checks and update some entries in the 'userOptions' structure
    userOptions.updateEntries( particlePointer->size(), not samples.empty() );
    if ( userOptions.averageDensity<0. )
        userOptions.averageDensity = averageDensity( *particlePointer, userOptions );   //computes what is the average density - if not supplied by the user
    
    
    // define some variables
    User_options tempOptions = userOptions;
    // the interpolation functions compute only velocity gradient, hence if the user requested the velocity divergence, shear or vorticity, tell the options to compute only the velocity gradient (the rest of the quantities will be computed from the velocity gradient)
    if ( userOptions.uField.selectedVelocityDerivatives() )
    {
        tempOptions.uField.velocity_gradient = true;
        tempOptions.uField.deselectVelocityDerivatives();
    }
    if ( userOptions.aField.selectedVelocityDerivatives() )
    {
        tempOptions.aField.velocity_gradient = true;
        tempOptions.aField.deselectVelocityDerivatives();
    }
    // the magnetic field is interpolated as a 3-component scalar field (its divergence and curl are computed afterwards from the gradient)
    tempOptions.uField.mapMagneticToScalar();
    tempOptions.aField.mapMagneticToScalar();
    
    
    
    
    //! Select only the particles in the region of interest
    vector<Particle_data> particlesRegion; // vector for storing the particles in '--region'
    if ( userOptions.regionOn )     // Select only the particles in the user defined region
    {
        //check that the box boundaries satisfy some conditions
        tempOptions.region.validSubBox( tempOptions.boxCoordinates, tempOptions.periodic );
        
        // find the box values with padding added and find the particles in that extended region
        tempOptions.paddedBox = tempOptions.region;
        tempOptions.paddedBox.addPadding( tempOptions.paddingLength ); // now stores the padded region of interest
        
        findParticlesInBox( *particlePointer, &particlesRegion, tempOptions );
        tempOptions.periodic = false;   // do not use the periodic translation algorithm after this
        tempOptions.updateFullBox( tempOptions.region );    // now the full box is the region of interest
        
        // update 'particlePointer' -> should point to 'particles_1'
        particlePointer->clear();   //clear the particle positions - don't need them anymore
        particlePointer = &particlesRegion;
    }
    
    vector<Particle_data> particlesPartition; // vector for storing the particles in '--partNo' partition
    if ( userOptions.partitionOn and userOptions.partNo>=0 )   // Select only the particles for a user given '--partNo'
    {
        MESSAGE::Message message( userOptions.verboseLevel );
        message << "The program will interpolate the fields in partition number " << userOptions.partNo << " of partition grid [" << MESSAGE::printElements( userOptions.partition, "," ) << "].\n" << MESSAGE::Flush;
        
        // get the optimal division in a subgrid
        std::vector< std::vector<size_t> > subgridList;
        std::vector< Box > subgridCoords;
        optimalPartitionSplit( *particlePointer, tempOptions, tempOptions.partition, &subgridList, &subgridCoords );
        copySubgridInformation( &tempOptions, subgridList, subgridCoords ); // writes the size of the subgrid to tempOptions.gridSize and the boundaries of the subgrid box to tempOptions.region
        userOptions.region = tempOptions.region;
        
        // find the box values with padding added and find the particles in that extended region
        tempOptions.paddedBox = tempOptions.region;
        tempOptions.paddedBox.addPadding( tempOptions.paddingLength ); // now stores the padded region of interest
        
        findParticlesInBox( *particlePointer, &particlesPartition, tempOptions );
        tempOptions.periodic = false;   // do not use the periodic translation algorithm after this
        
        // update 'particlePointer' -> should point to 'particles_2'
        particlePointer->clear();   //clear the particle positions - don't need them anymore
        particlePointer = &particlesPartition;
        
        subgrid( userOptions, subgridList );
    }
    
    
    //! now compute the DTFE interpolation
    if ( userOptions.partitionOn and userOptions.partNo<0 ) // if option '--partition' was defined (without option '--partNo')
    {
#if NO_DIM==2
        int const totalPartitions = tempOptions.partition[0] * tempOptions.partition[1];
#elif NO_DIM==3
        int const totalPartitions = tempOptions.partition[0] * tempOptions.partition[1] * tempOptions.partition[2];
#endif
        MESSAGE::Message message( userOptions.verboseLevel );
        message << "The program will interpolate the fields in the region of interest using " << totalPartitions << " partitions defined via the grid [" << MESSAGE::printElements( userOptions.partition, "," ) << "].\n" << MESSAGE::Flush;
        
        // reserve memory for the quantities of interest
        uQuantities->reserveMemory( &(tempOptions.gridSize[0]), tempOptions.uField );
        aQuantities->reserveMemory( &(tempOptions.gridSize[0]), tempOptions.aField );
        
        
        // get the optimal division in a subgrid
        std::vector< std::vector<size_t> > subgridList;
        std::vector< Box > subgridCoords;
        optimalPartitionSplit( *particlePointer, tempOptions, tempOptions.partition, &subgridList, &subgridCoords );
        
        
        for (int i=0; i<totalPartitions; ++i)
        {
            message << "\n<<< Interpolating the fields for partition " << i+1 << " ...\n";
            
            // define a temporary 'User_options' object to send values to the DTFE computing function
            User_options tempOpt = tempOptions;
            tempOpt.partNo = i;
            copySubgridInformation( &tempOpt, subgridList, subgridCoords ); // writes the size of the subgrid to tempOptions.gridSize and the boundaries of the subgrid box to tempOptions.region
            
            
            // find the box values with padding added and find the particles in that extended region
            tempOpt.paddedBox = tempOpt.region;
            tempOpt.paddedBox.addPadding( tempOpt.paddingLength ); // now stores the padded region of interest
            
            vector<Particle_data> tempPart;
            findParticlesInBox( *particlePointer, &tempPart, tempOpt );
            tempOpt.periodic = false;   // do not use the periodic translation algorithm in the "triangulation_interpolation.cpp" file
            
            
            //compute the DTFE interpolation
            Quantities temp_uQuantities, temp_aQuantities; // temporary variable to store the grid interpolation for each subgrid
            DTFE_parallel( &tempPart, samples, tempOpt, &temp_uQuantities, &temp_aQuantities );
            
            
            // write the fields from the given partition to the main grid results
            copySubgridResultsToMain( temp_uQuantities, tempOptions.gridSize, tempOpt.uField, tempOpt, subgridList, uQuantities );
            copySubgridResultsToMain( temp_aQuantities, tempOptions.gridSize, tempOpt.aField, tempOpt, subgridList, aQuantities );
        }
    }
    else
    {
        vector<Particle_data> tempPart;
        if ( not userOptions.regionOn and not userOptions.partitionOn )   // if none of the options '--region' or '--partition' were selected, pad the box
        {
            // find the box values with padding added and find the particles in that extended region
            tempOptions.region = tempOptions.boxCoordinates;
            tempOptions.paddedBox = tempOptions.region;
            tempOptions.paddedBox.addPadding( tempOptions.paddingLength ); // now stores the padded region of interest
            
            if ( userOptions.periodic ) // if periodic box, translate the box
            {
                findParticlesInBox( *particlePointer, &tempPart, tempOptions );
                
                // update 'particlePointer' -> should point to 'tempPart'
                particlePointer->clear();   //clear the particle positions - don't need them anymore
                particlePointer = &tempPart;
            }
        }
        
        tempOptions.periodic = false;   // do not use the periodic translation algorithm in the "triangulation_interpolation.cpp" file
        
        
        //compute the DTFE interpolation
        DTFE_parallel( particlePointer, samples, tempOptions, uQuantities, aQuantities );
    }
    
    
    // compute the velocity divergence, shear or vorticity, if any
    computeDivergenceShearVorticity( userOptions.uField, userOptions.verboseLevel, uQuantities );
    computeDivergenceShearVorticity( userOptions.aField, userOptions.verboseLevel, aQuantities );
    
    // compute the magnetic field divergence and curl, if any
    computeMagneticDivergenceCurl( userOptions.uField, uQuantities );
    computeMagneticDivergenceCurl( userOptions.aField, aQuantities );
    
}



// computes the actual Delaunay triangulation and the interpolation to grid
extern void DTFE_interpolation(vector<Particle_data> *p,
                               vector<Sample_point> &samples,
                               User_options &userOptions,
                               Quantities *uQuantities,
                               Quantities *aQuantities);

/* Calls the function corresponding to the chosen interpolation method. */
void interpolate(vector<Particle_data> *allParticles,
                 vector<Sample_point> &samples,
                 User_options & userOptions,
                 Quantities *uQuantities,
                 Quantities *aQuantities)
{
    if ( userOptions.DTFE )
        DTFE_interpolation( allParticles, samples, userOptions, uQuantities, aQuantities );
    else if ( userOptions.NGP )
        NGP_interpolation( allParticles, samples, userOptions, aQuantities );
    else if ( userOptions.CIC )
        CIC_interpolation( allParticles, samples, userOptions, aQuantities );
    else if ( userOptions.TSC )
        TSC_interpolation( allParticles, samples, userOptions, aQuantities );
    else if ( userOptions.SPH )
        SPH_interpolation( allParticles, samples, userOptions, aQuantities );
    else
        throwError( "Unknow interpolation method in function 'interpolate'." );
}




/* This function splits the data in several partitions that can be used by shared memory processors to do the computations in parallel. It uses the OpenMP directives for doing so. */
void DTFE_parallel(vector<Particle_data> *allParticles,
                   vector<Sample_point> &samples,
                   User_options & userOptions,
                   Quantities *uQuantities,
                   Quantities *aQuantities)
{
#ifndef OPEN_MP
    //directly compute the DTFE interpolation
    interpolate( allParticles, samples, userOptions, uQuantities, aQuantities );  // this function deletes the vector 'allParticles'
    return;
#else
    int const noAvailableProcessors = omp_get_max_threads();
    
    // if only 1 processor is available, do nothing
    if ( noAvailableProcessors==1 or not samples.empty() or userOptions.redshiftConeOn )
    {
        interpolate( allParticles, samples, userOptions, uQuantities, aQuantities );
        return;
    }
    
    
    // if more than 1 processor is available, split the data in partitions
    std::vector<size_t> pGrid(NO_DIM,0); // the parallel grid
    parallelGrid( noAvailableProcessors, userOptions, &(pGrid[0]) );
    int const noProcessors = NO_DIM==2 ? pGrid[0]*pGrid[1] : pGrid[0]*pGrid[1]*pGrid[2];    //number of processors actually used (may be different from 'noAvailableProcessors')
    
    // reserve memory for the quantities of interest
    uQuantities->reserveMemory( &(userOptions.gridSize[0]), userOptions.uField );
    aQuantities->reserveMemory( &(userOptions.gridSize[0]), userOptions.aField );
    
    // define some variables to keep track of the time and particle numbers associated to each processor
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "From now on only the master thread will show messages on how the computation is going. Not all threads take the same execution time, so there may be a discrepancy between the messages displayed to the user and the computations across all threads.\n\n" << MESSAGE::Flush;
    size_t processorParticles[noProcessors];    // number of particles associated to each processor
    Real processorTime[noProcessors];           // the actual runtime of each thread
    size_t const noTotalParticles = allParticles->size();
    
    
    // determine the optimal grid to split the computational load
    std::vector< std::vector<size_t> > subgridList;
    std::vector< Box > subgridCoords;
    optimalPartitionSplit( *allParticles, userOptions, pGrid, &subgridList, &subgridCoords );
    
    
    // this contains the parallel section
#pragma omp parallel num_threads( noProcessors )
    {
        int const threadNo = omp_get_thread_num(); // the number of the thread - between 0 to noProcessors
        User_options tempOptions = userOptions;
        tempOptions.noProcessors = noProcessors;    // the number of parallel threads
        tempOptions.threadId = threadNo;            // the id of the thread
        tempOptions.verboseLevel = (userOptions.verboseLevel>0) ? 1 : userOptions.verboseLevel;   // the slave processors will not show runtime messages (with exception of errors and warnings)
        
        if ( threadNo==0 )
            tempOptions.verboseLevel = userOptions.verboseLevel;    // only the master will show runtime messages
        tempOptions.partNo = threadNo;
        tempOptions.partition.clear();
        for (int i=0; i<NO_DIM; ++i)
            tempOptions.partition.push_back( pGrid[i] );
        tempOptions.updateFullBox( userOptions.region );    // this is the new full box
        
        
        // compute the region allocated to each processor
        copySubgridInformation( &tempOptions, subgridList, subgridCoords ); // writes the size of the subgrid to tempOptions.gridSize and the boundaries of the subgrid box to tempOptions.region
        tempOptions.paddedBox = tempOptions.region;
        tempOptions.paddedBox.addPadding( tempOptions.paddingLength );
        
        
        // find the particles in the 'paddedBox'
        vector<Particle_data> particles;
        findParticlesInBox( *allParticles, &particles, tempOptions );
        processorParticles[ threadNo ] = particles.size();
        
        
        //compute the DTFE interpolation
        Quantities temp_uQuantities, temp_aQuantities; // temporary variable to store the grid interpolation for each subgrid
        interpolate( &particles, samples, tempOptions, &temp_uQuantities, &temp_aQuantities ); //this function deletes the vector 'particles'
        processorTime[ threadNo ] = tempOptions.totalTime;  // this is the total CPU time for all threads until this thread ended
        
        
        // write the fields from the given partition to the main grid results
        copySubgridResultsToMain( temp_uQuantities, userOptions.gridSize, tempOptions.uField, tempOptions, subgridList, uQuantities );
        copySubgridResultsToMain( temp_aQuantities, userOptions.gridSize, tempOptions.aField, tempOptions, subgridList, aQuantities );  //copy the results from each processor to a main variable result
        if ( threadNo==0 )
            message << "\nWaiting for all threads to finish the computations ... " << MESSAGE::Flush;
    }
    message << "Done.\n";
    allParticles->clear();
    
    
    // show statistics to the user about the execution in parallel
    userOptions.totalTime += maximum( processorTime, noProcessors );    // update total time by the largest value of CPU time
    approximativeThreadTime( processorTime, noProcessors );   // compute approximative CPU times for each thread
    message << "Statistics of the execution across the " << noProcessors << " threads:\n";
    for (int i=0; i<noProcessors; ++i)
        message << "\t Thread " << i << " had " << processorParticles[i] << " particles (which represent " << setprecision(4) <<  Real(processorParticles[i])/noTotalParticles*100. << "\%) and took " << processorTime[i] << " sec. \n";
    message << MESSAGE::Flush;
    
#endif
}






//! Velocity derivative quantities
/* Computes the velocity divergence from the velocity gradient. */
Real velocityDivergence(Pvector<Real,noGradComp> &velGrad)
{
#if NO_DIM==2
    return velGrad[0] + velGrad[3];
#elif NO_DIM==3
    return velGrad[0] + velGrad[4] + velGrad[8];
#endif
}

/* Computes the velocity shear from the velocity gradient. */
Pvector<Real,noShearComp> velocityShear(Pvector<Real,noGradComp> &velGrad)
{
    Pvector<Real,noShearComp> temp;
    size_t index = 0;
    for (int i=0; i<NO_DIM-1; ++i)
        for (int j=i; j<NO_DIM; ++j)
            temp[index++] = (velGrad[i*NO_DIM+j] + velGrad[j*NO_DIM+i]) / 2.;
    
    // still need to substract the trace since the shear matrix is traceless
#if NO_DIM==2
    Real trace = (temp[0]+velGrad[3]) / 2.;
    temp[0] -= trace;
#elif NO_DIM==3
    Real trace = (temp[0]+temp[3]+velGrad[8]) / 3.;
    temp[0] -= trace;
    temp[3] -= trace;
#endif
    return temp;
}

/* Computes the velocity shear from the velocity gradient. */
Pvector<Real,noVortComp> velocityVorticity(Pvector<Real,noGradComp> &velGrad)
{
    Pvector<Real,noVortComp> temp;
    size_t index = 0;
    for (int i=0; i<NO_DIM; ++i)
        for (int j=i+1; j<NO_DIM; ++j)
            temp[index++] = (velGrad[i*NO_DIM+j] - velGrad[j*NO_DIM+i]) / 2.;
    return temp;
}

/* This function computes the velocity divergence, shear and vorticity. */
void computeDivergenceShearVorticity(Field &fields,
                                     int const verboseLevel,
                                     Quantities *q)
{
    if ( q->velocity_gradient.empty() ) return;
    
    // compute the velocity divergence
    if ( fields.velocity_divergence )
    {
        q->velocity_divergence.reserve( q->velocity_gradient.size() );
        for (std::vector< Pvector<Real,noGradComp> >::iterator it=q->velocity_gradient.begin(); it!=q->velocity_gradient.end(); ++it )
            q->velocity_divergence.push_back( velocityDivergence(*it) );
        
//        for (int i=0; i<q->velocity_divergence.size(); ++i)
//            if ( not boost::math::isfinite(q->velocity_divergence[i]) )
//            {
//                q->velocity_divergence[i] = Real(0.);
//                std::cout << "<<< Found non-numerical value in velocity divergence at array index " << i << ". The value of the velocity divergence at this grid point will be set to 0." << std::flush;
//            }
    }
    
    // compute the velocity shear
    if ( fields.velocity_shear )
    {
        q->velocity_shear.reserve( q->velocity_gradient.size() );
        for (std::vector< Pvector<Real,noGradComp> >::iterator it=q->velocity_gradient.begin(); it!=q->velocity_gradient.end(); ++it )
            q->velocity_shear.push_back( velocityShear(*it) );
    }
    
    // compute the velocity vorticity
    if ( fields.velocity_vorticity )
    {
        q->velocity_vorticity.reserve( q->velocity_gradient.size() );
        for (std::vector< Pvector<Real,noGradComp> >::iterator it=q->velocity_gradient.begin(); it!=q->velocity_gradient.end(); ++it )
            q->velocity_vorticity.push_back( velocityVorticity(*it) );
    }
    
    if ( not fields.velocity_gradient )
        q->velocity_gradient.clear();
}


/* Computes the divergence and the curl of the magnetic field from its gradient. The magnetic field is interpolated as the first 3 components
of the scalar data, so its gradient is stored in 'scalar_gradient' with the layout grad[c*NO_DIM+i] = dB_c/dx_i. */
void computeMagneticDivergenceCurl(Field &fields,
                                   Quantities *q)
{
#if NO_DIM==3
    if ( q->scalar_gradient.empty() or not fields.selectedMagnetic() ) return;
    
    if ( fields.magnetic_divergence )   // div B = dBx/dx + dBy/dy + dBz/dz
    {
        q->magnetic_divergence.reserve( q->scalar_gradient.size() );
        for (size_t n=0; n<q->scalar_gradient.size(); ++n)
        {
            Pvector<Real,noScalarGradComp> &g = q->scalar_gradient[n];
            q->magnetic_divergence.push_back( g[0*3+0] + g[1*3+1] + g[2*3+2] );
        }
    }
    
    if ( fields.magnetic_curl )         // curl B = (dBz/dy - dBy/dz, dBx/dz - dBz/dx, dBy/dx - dBx/dy)
    {
        q->magnetic_curl.reserve( q->scalar_gradient.size() );
        for (size_t n=0; n<q->scalar_gradient.size(); ++n)
        {
            Pvector<Real,noScalarGradComp> &g = q->scalar_gradient[n];
            Pvector<Real,3> curl;
            curl[0] = g[2*3+1] - g[1*3+2];
            curl[1] = g[0*3+2] - g[2*3+0];
            curl[2] = g[1*3+0] - g[0*3+1];
            q->magnetic_curl.push_back( curl );
        }
    }
    
    // keep the gradient only if requested
    if ( not fields.magnetic_gradient and not fields.scalar_gradient )
        q->scalar_gradient.clear();
#endif
}





// intialize the static members
#ifndef VELOCITY
Pvector<Real,noVelComp> Data_structure::_velocity = Pvector<Real,noVelComp>::zero();
#endif
#ifndef SCALAR
Pvector<Real,noScalarComp> Data_structure::_scalar = Pvector<Real,noScalarComp>::zero();
#endif






#ifdef TRIANGULATION
#include "DTFE.h"

// computes the actual Delaunay triangulation and the interpolation to grid
extern void DTFE_interpolation(vector<Particle_data> *p,
                               vector<Sample_point> &samples,
                               User_options &userOptions,
                               Quantities *uQuantities,
                               Quantities *aQuantities,
                               DT &delaunay_triangulation);

/* This function interpolates fileds to a grid using the DTFE method.
It takes the following arguments:
        allParticles - stores the particle positions, velocities, scalars (if any) and weights
        samples - stores the position of the grid points where samples are taken (if this grid points don't have a uniform distribution)
        userOptions - options specified by the user plus additional program wise variables
        quantities - structure that stores vectors for all the allowed quantities that can be computed in the function
        delaunay_triangulation - returns the Delaunay triangulation of the point distribution

NOTE: This function clears the vector 'allParticles'. */
void DTFE(vector<Particle_data> *allParticles,
          vector<Sample_point> &samples,
          User_options & userOptions,,
          Quantities *uQuantities,
          Quantities *aQuantities,
          DT &delaunay_triangulation)
{
    userOptions.field.triangulation = true;
    
    //! If the user requested for a set of random particles, generate them
    if ( userOptions.poisson!=0 )
        randomParticles( allParticles, &userOptions );
    
    
    //! Select a random subset of the data if the user asked so
    vector<Particle_data> *particlePointer = allParticles; // pointer to the vector containing the particle data
    vector<Particle_data> particlesRandomSubsample;  // vector for storing the random subsample particles
    if (userOptions.randomSample>=Real(0.) )
    {
        randomSample( *particlePointer, &particlesRandomSubsample, userOptions );
        
        // update 'particlePointer' -> should point to 'particlesRandomSubsample'
        particlePointer->clear();   //clear the particle positions - don't need them anymore
        particlePointer = &particlesRandomSubsample;
    }
    
    
    //! Do consistency checks and update some entries in the 'userOptions' structure
    userOptions.updateEntries( particlePointer->size(), not samples.empty() );
    if ( userOptions.averageDensity<0. )
        userOptions.averageDensity = averageDensity( *particlePointer, userOptions );   //computes what is the average density - if not supplied by the user
    
    
    // define some variables
    User_options tempOptions = userOptions;
    // the interpolation functions compute only velocity gradient, hence if the user requested the velocity divergence, shear or vorticity, tell the options to compute only the velocity gradient (the rest of the quantities will be computed from the velocity gradient)
    if ( userOptions.uField.selectedVelocityDerivatives() )
    {
        tempOptions.uField.velocity_gradient = true;
        tempOptions.uField.deselectVelocityDerivatives();
    }
    if ( userOptions.aField.selectedVelocityDerivatives() )
    {
        tempOptions.aField.velocity_gradient = true;
        tempOptions.aField.deselectVelocityDerivatives();
    }
    // the magnetic field is interpolated as a 3-component scalar field (its divergence and curl are computed afterwards from the gradient)
    tempOptions.uField.mapMagneticToScalar();
    tempOptions.aField.mapMagneticToScalar();
    
    
    
    
    //! Select only the particles in the region of interest
    vector<Particle_data> particlesRegion; // vector for storing the particles in '--region'
    if ( userOptions.regionOn )     // Select only the particles in the user defined region
    {
        //check that the box boundaries satisfy some conditions
        tempOptions.region.validSubBox( tempOptions.boxCoordinates, tempOptions.periodic );
        
        // find the box values with padding added and find the particles in that extended region
        tempOptions.paddedBox = tempOptions.region;
        tempOptions.paddedBox.addPadding( tempOptions.paddingLength ); // now stores the padded region of interest
        
        findParticlesInBox( *particlePointer, &particlesRegion, tempOptions );
        tempOptions.periodic = false;   // do not use the periodic translation algorithm after this
        tempOptions.updateFullBox( tempOptions.region );    // now the full box is the region of interest
        
        // update 'particlePointer' -> should point to 'particles_1'
        particlePointer->clear();   //clear the particle positions - don't need them anymore
        particlePointer = &particlesRegion;
    }
    
    vector<Particle_data> particlesPartition; // vector for storing the particles in '--partNo' partition
    if ( userOptions.partitionOn and userOptions.partNo>=0 )   // Select only the particles for a user given '--partNo'
    {
        MESSAGE::Message message( userOptions.verboseLevel );
        message << "The program will interpolate the fields in partition number " << userOptions.partNo << " of partition grid [" << MESSAGE::printElements( userOptions.partition, "," ) << "].\n" << MESSAGE::Flush;
        
        subgrid( userOptions.gridSize, &tempOptions );  // computes what is the size of the subgrid (will be written to tempOptions.gridSize) and the boundaries of the box bounding the subgrid (saved in tempOptions.region)
        userOptions.region = tempOptions.region;
        
        // find the box values with padding added and find the particles in that extended region
        tempOptions.paddedBox = tempOptions.region;
        tempOptions.paddedBox.addPadding( tempOptions.paddingLength ); // now stores the padded region of interest
        
        findParticlesInBox( *particlePointer, &particlesPartition, tempOptions );
        tempOptions.periodic = false;   // do not use the periodic translation algorithm after this
        
        // update 'particlePointer' -> should point to 'particles_2'
        particlePointer->clear();   //clear the particle positions - don't need them anymore
        particlePointer = &particlesPartition;
    }
    
    
    //! now compute the DTFE interpolation
    if ( userOptions.partitionOn and userOptions.partNo<0 ) // if option '--partition' was defined (without option '--partNo')
    {
        throwError( "There is no implementation of the option '--partition' in the absence of option '--partNo' when requering that the 'DTFE' function returns the Delaunay triangulation." );
    }
    else
    {
        vector<Particle_data> tempPart;
        if ( not userOptions.regionOn and not userOptions.partitionOn )   // if none of the options '--region' or '--partition' were selected, pad the box
        {
            // find the box values with padding added and find the particles in that extended region
            tempOptions.region = tempOptions.boxCoordinates;
            tempOptions.paddedBox = tempOptions.region;
            tempOptions.paddedBox.addPadding( tempOptions.paddingLength ); // now stores the padded region of interest
            
            if ( userOptions.periodic ) // if periodic box, translate the box
            {
                findParticlesInBox( *particlePointer, &tempPart, tempOptions );
                
                // update 'particlePointer' -> should point to 'tempPart'
                particlePointer->clear();   //clear the particle positions - don't need them anymore
                particlePointer = &tempPart;
            }
        }
        
        tempOptions.periodic = false;   // do not use the periodic translation algorithm in the "triangulation_interpolation.cpp" file
        
        
        //compute the DTFE interpolation
        DTFE_interpolation( particlePointer, samples, tempOptions, uQuantities, aQuantities, delaunay_triangulation );
    }
    
    
    // compute the velocity divergence, shear or vorticity, if any
    computeDivergenceShearVorticity( userOptions.uField, userOptions.verboseLevel, uQuantities );
    computeDivergenceShearVorticity( userOptions.aField, userOptions.verboseLevel, aQuantities );
    
    // compute the magnetic field divergence and curl, if any
    computeMagneticDivergenceCurl( userOptions.uField, uQuantities );
    computeMagneticDivergenceCurl( userOptions.aField, aQuantities );
    
    
    // if the computation was done only for a given partition, output to the user the grid indices used for that
    if ( userOptions.partitionOn and userOptions.partNo>=0 )
    {
        subgrid( userOptions );
    }
}
#endif

