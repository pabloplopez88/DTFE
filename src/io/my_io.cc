#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#ifdef HDF5
	#include <hdf5.h>
	#define N_part_types 7
#endif
#ifdef TYPE_FOUR_GADGET
	#define N_part_types 2
#endif
#ifdef TYPE_TWO_GADGET
	#define N_part_types 6
#endif

#ifdef PERCENT
	#define YELLOW(a)   fprintf(stdout,"\033[01;33m%s\033[22;0m",(a))
#endif

struct cosmoparam
{
  double  omegam		    		;  /* Omega Materia                         */
  double  omegal		    		;  /* Omega Lambda                          */
  double  omegak		    		;  /* Omega Curvatura                       */
  double  hparam            			;  /* Parámetro de Hubble adimensional      */
  double  lbox              			;  /* Lado del box [Kpc / h]                */
  double  Mpart[N_part_types]             	;  /* Masa de la partícula [10^10 Msol / h] */
  unsigned long long  npart[N_part_types] 	;  /* Número de partículas                  */
  int     nfiles 				;  /* Número de files                  	    */
  double  redshift      			;  /* Redshift                              */
  double  aexp           			;  /*                                       */
  double  Hubble_a          			;  /*                                       */
};

#ifdef TYPE_FOUR_GADGET
struct io_header
{
  long long npart[N_part_types];      // Number of particles in the file
  long long npartTotal[N_part_types]; // Number of particles of the simulation
  double    mass[N_part_types];       // Mass particle in 10^10 Msun/h
  double    time;          // Scale Factor a
  double    redshift;      // z
  double    BoxSize;       // Boxsize in Mpc/h
  int       num_files;     // Number of files
  long long Ntrees;
  long long NtreesTot;
};
#else
struct io_header
{
  int npart[N_part_types];             /*!< number of particles of each type in this file */
  double mass[N_part_types];           /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[N_part_types];          /*!< total number of particles of each type in this snapshot. This can be
                                            different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  int num_files;                       /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[N_part_types];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];	                     /*!< fills to 256 Bytes */
};
#endif

static struct cosmoparam cp;

#ifdef HDF5

// -------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------
// ***************************************************  HDF5  **************************************************************
// -------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------

#define h5_open_file(file_id, name)         H5Fopen(file_id, name, H5P_DEFAULT)
#define h5_open_group(file_id, name)        H5Gopen(file_id, name, H5P_DEFAULT)
#define h5_open_dataset(file_id, name)      H5Dopen(file_id, name, H5P_DEFAULT)
#define h5_open_attribute(parent_id, name)  H5Aopen_name(parent_id, name)

static unsigned long long NumPart_ThisFile(const char *filename, const int type_particle)
{
  unsigned long long npart[N_part_types];

  // Read Header
  hid_t hdf5_file = h5_open_file(filename, H5F_ACC_RDONLY);
  hid_t hdf5_headergrp = h5_open_group(hdf5_file, "/Header");
  hid_t hdf5_attribute = h5_open_attribute(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_ULLONG, &npart);
  H5Aclose(hdf5_attribute);
  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  return npart[type_particle];
}

static void lee(const char *filename, const int type_particle, Read_data<Real> *readData, User_options *userOptions)
{
  char key[50];
	long long i, j, pc;
  long long elements;
  herr_t status;
  hid_t hdf5_file;
  hid_t hdf5_typeparticles;
  hid_t hdf5_dataset;
  hid_t hdf5_dataspace;
  hid_t data_type; 
  hsize_t dims[20];
  int rank, type;
  size_t size; 
  void *buff;
  long long numToRead;
  hid_t hdf5_type[6] = {H5T_NATIVE_INT, H5T_NATIVE_LONG, \
                        H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE,\
                        H5T_NATIVE_UINT, H5T_NATIVE_ULONG};
  struct particle_data 
  {
    double      x;
    double      y;
    double      z;
    #ifdef WEIGHT
    double      mp;
    #endif
    #ifdef VELOCITY
    double      vx;
    double      vy;
    double      vz;
    #endif
    #ifdef IDS
    unsigned long long id;
    #endif
  } *particles;

  numToRead = (long long)NumPart_ThisFile(filename, type_particle);
  particles = (struct particle_data *) malloc(numToRead*sizeof(struct particle_data));
	
  sprintf(key,"/PartType%d",type_particle);
  hdf5_file = h5_open_file(filename, H5F_ACC_RDONLY);

  for(j=0;j<N_part_types;j++)
  {
    if(j == type_particle)
    { 
      hdf5_typeparticles = h5_open_group(hdf5_file, key);

      hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "Coordinates");
      hdf5_dataspace     = H5Dget_space(hdf5_dataset);
      data_type          = H5Dget_type(hdf5_dataset);
      size               = H5Tget_size(data_type);
      rank               = H5Sget_simple_extent_ndims(hdf5_dataspace);      
      H5Sget_simple_extent_dims(hdf5_dataspace, dims, NULL); 
      
      elements = 1;
      for(i=0;i<rank;i++)
        elements *= dims[i];

      type = size == 4 ? 2 : 3;
      buff = (void *) malloc(elements*size);
      status = H5Dread(hdf5_dataset, hdf5_type[type], H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
      
      if(status < 0)
      {
        fprintf(stdout,"Error H5Dread Coordinates\n"); fflush(stdout);
        exit(EXIT_FAILURE);
      }

      for(i=0, pc=0;i<numToRead;i++)
      {
	      switch (type) {
  	      case 2:
            particles[pc].x = ((double)(*(float  *)(buff+size*(3*i  ))));
            particles[pc].y = ((double)(*(float  *)(buff+size*(3*i+1))));
            particles[pc].z = ((double)(*(float  *)(buff+size*(3*i+2))));
		      	break;
  	      case 3: 
            particles[pc].x = ((double)(*(double *)(buff+size*(3*i  ))));
            particles[pc].y = ((double)(*(double *)(buff+size*(3*i+1))));
            particles[pc].z = ((double)(*(double *)(buff+size*(3*i+2))));
            break;
		      default:
		      	fprintf(stdout,"Error! Type Coordinates not correct\n"); fflush(stdout);
  	        exit(EXIT_FAILURE);     
	      }
        pc++;
      }
      
      free(buff);
      H5Tclose(data_type);
      H5Dclose(hdf5_dataset);

#ifdef VELOCITY
      if ( userOptions->magnetic_field )
      {
	printf("***** ESTAMOS LEYENDO MAGNETIC FLUX *****\n");
        hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "MagneticFluxDensities");
      }
      else
      {
	printf("***** ESTAMOS LEYENDO VELOCIDADES *****\n");
        hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "Velocities");
      }
      hdf5_dataspace     = H5Dget_space(hdf5_dataset);
      data_type          = H5Dget_type(hdf5_dataset);
      size               = H5Tget_size(data_type);
      rank               = H5Sget_simple_extent_ndims(hdf5_dataspace);      
      H5Sget_simple_extent_dims(hdf5_dataspace, dims, NULL); 

      elements = 1;
      for(i=0;i<rank;i++)
        elements *= dims[i];

      buff = (void *) malloc(elements*size);
      type = size == 4 ? 2 : 3;
      status = H5Dread(hdf5_dataset, hdf5_type[type], H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
      if(status < 0)
      {
        fprintf(stdout,"Error H5Dread Velocities\n"); fflush(stdout);
        exit(EXIT_FAILURE);
      }

      for(i=0, pc=0;i<numToRead;i++)
      {
	      switch (type) {
  	      case 2:
            particles[pc].vx = ((double)(*(float  *)(buff+size*(3*i  ))));
            particles[pc].vy = ((double)(*(float  *)(buff+size*(3*i+1))));
            particles[pc].vz = ((double)(*(float  *)(buff+size*(3*i+2))));
		      	break;
  	      case 3:
            particles[pc].vx = ((double)(*(double *)(buff+size*(3*i  ))));
            particles[pc].vy = ((double)(*(double *)(buff+size*(3*i+1))));
            particles[pc].vz = ((double)(*(double *)(buff+size*(3*i+2))));
            break;
		      default:
		      	fprintf(stdout,"Error! Type Velocities not correct\n"); fflush(stdout);
  	        exit(EXIT_FAILURE);     
          	break;
	      }
        pc++;
      }
    
      free(buff);
      H5Tclose(data_type);
      H5Dclose(hdf5_dataset);
#endif
      
#ifdef WEIGHT
	if(type_particle < 7)
	{
	   hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "Masses");
	   hdf5_dataspace     = H5Dget_space(hdf5_dataset);
	   data_type          = H5Dget_type(hdf5_dataset);
	   size               = H5Tget_size(data_type);
	   rank               = H5Sget_simple_extent_ndims(hdf5_dataspace);      
	   H5Sget_simple_extent_dims(hdf5_dataspace, dims, NULL); 
	
	   elements = 1;
	   for(i=0;i<rank;i++)
	     elements *= dims[i];
	
	   buff = (void *) malloc(elements*size);
	   type = size == 4 ? 2 : 3;
	   status = H5Dread(hdf5_dataset, hdf5_type[type], H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
	   if(status < 0)
	   {
	     fprintf(stdout,"Error H5Dread Masses\n"); fflush(stdout);
	     exit(EXIT_FAILURE);
	   }
	
	   for(i=0, pc=0;i<numToRead;i++)
	   {
	     switch (type) {
	       case 2:
	         particles[pc].mp = ((double)(*(float  *)(buff+size*i)));
	       	break;          
	       case 3:           
	         particles[pc].mp = ((double)(*(double *)(buff+size*i)));
	         break;
	       default:
	       	fprintf(stdout,"Error! Type Masses not correct\n"); fflush(stdout);
	         exit(EXIT_FAILURE);     
	       	break;
	     }
	     pc++;
	   }
	 
	   free(buff);
	   H5Tclose(data_type);
	 	 H5Dclose(hdf5_dataset);
	
	}
        else
        {
	  for(i=0, pc=0;i<numToRead;i++)
          {
    	    particles[pc].mp = cp.Mpart[type_particle];
	    pc++;
          }
	}
#endif


#ifdef IDS   
      hdf5_dataset       = h5_open_dataset(hdf5_typeparticles, "ParticleIDs");
      hdf5_dataspace     = H5Dget_space(hdf5_dataset);
      data_type          = H5Dget_type(hdf5_dataset);
      size               = H5Tget_size(data_type);
      rank               = H5Sget_simple_extent_ndims(hdf5_dataspace);      
      H5Sget_simple_extent_dims(hdf5_dataspace, dims, NULL); 

      elements = 1;
      for(i=0;i<rank;i++)
        elements *= dims[i];

      buff = (void *) malloc(elements*size);
      type = size == 4 ? 0 : 1;
      
      // Special case for unsigned integers - treat as unsigned in memory
      H5T_class_t h5class = H5Tget_class(data_type);
      H5T_sign_t  sign  = H5Tget_sign(data_type);
      if(sign==H5T_SGN_NONE && h5class==H5T_INTEGER)
      {
        if(type == 0) 
          type = 4;
        if(type == 1) 
          type = 5;
      }

      // Read dataset
      status = H5Dread(hdf5_dataset, hdf5_type[type], H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
      if(status < 0)
      {
        fprintf(stdout,"Error H5Dread ParticleIDs\n"); fflush(stdout);
        exit(EXIT_FAILURE);
      }

      for(i=0, pc=0;i<numToRead;i++)
      {
	      switch (type) {
  	      case 0:
            particles[pc].id = ((unsigned long long)(*(int  *)(buff+size*i)));
		      	break;       
  	      case 1:        
            particles[pc].id = ((unsigned long long)(*(long *)(buff+size*i)));
		      	break;       
  	      case 4:        
            particles[pc].id = ((unsigned long long)(*(unsigned int *)(buff+size*i)));
		      	break;       
  	      case 5:        
            particles[pc].id = ((unsigned long long)(*(long unsigned *)(buff+size*i)));
		      	break;
		      default:
		      	fprintf(stdout,"Error! Type ParticleIDs not correct\n"); fflush(stdout);
  	        exit(EXIT_FAILURE);     
		      	break;
	      }
        pc++;
      }
 
      free(buff);
      H5Tclose(data_type);
      H5Dclose(hdf5_dataset);
#endif
      H5Gclose(hdf5_typeparticles);

      if(numToRead != pc)
      {
		     fprintf(stdout,"Error! Bad Read Particles %s - %d - %d\n", filename, pc, numToRead); fflush(stdout);
  	     exit(EXIT_FAILURE);     
      }

    }
  }

  H5Fclose(hdf5_file);

	size_t noParticles = numToRead;
#ifdef SCALAR
  int noScalars = 0;
  for (i=3; i<userOptions->readParticleData.size(); ++i)
      if ( userOptions->readParticleData[i] )
          noScalars += 1;
#endif

  readData->position( noParticles );        // particle positions
#ifdef WEIGHT
  readData->weight( noParticles );        // particle weights (weights = particle masses from the snapshot file)
#endif
#ifdef VELOCITY
  readData->velocity( noParticles );        // particle velocities
#endif
#ifdef SCALAR
  if ( noScalars>0 )
  	readData->scalar( noParticles );        // particle scalar quantity
#endif

  // copy the relevant information to the '*readData' structure
  Real *positions = readData->position();  //particle positions
#ifdef WEIGHT
  Real *weights = readData->weight();    //particle weights (e.g. weights = particle masses)
#endif
#ifdef VELOCITY
  Real *velocities = readData->velocity();  //particle velocity
#endif

#ifdef CUT_REGION
	double *limits;
	double *lado;
	const double delta = userOptions->paddingParticles / (Real)pow( (double)cp.npart[1], 1./NO_DIM) * cp.lbox;

	limits = (double *) malloc(2*NO_DIM*sizeof(double));
	lado   = (double *) malloc(  NO_DIM*sizeof(double));
  		

  //fprintf(stdout,"delta = %f, pad = %d, lbox: %f\n",(float)delta,(int)userOptions->paddingParticles,(float)cp.lbox);

  for(i=0; i<NO_DIM; i++)
	{
    		limits[2*i]   = userOptions->region.coords[2*i]*cp.lbox + pow(-1,2*i+1)*delta;
		limits[2*i+1] = userOptions->region.coords[2*i+1]*cp.lbox + pow(-1,2*i)*delta;
		if(!userOptions->periodic)
		{
			limits[2*i]   = limits[2*i]   < 0       ? 0       : limits[2*i];
			limits[2*i+1] = limits[2*i+1] > cp.lbox ? cp.lbox : limits[2*i+1];
		}
		lado[i] = limits[2*i+1] - limits[2*i];
  		//fprintf(stdout,">> i = %d, region.coords: %f %f\n",i,userOptions->region.coords[2*i],userOptions->region.coords[2*i+1]);
	}
#endif

  for(i=0, pc=0; i<numToRead; i++)
  {

#ifdef CUT_REGION
	double aa = fabs(limits[0] - particles[i].x);
	double bb = fabs(limits[1] - particles[i].x);
	double cc = fabs(limits[2] - particles[i].y);
	double dd = fabs(limits[3] - particles[i].y);
	double ee = fabs(limits[4] - particles[i].z);
	double ff = fabs(limits[5] - particles[i].z);

	if(userOptions->periodic)
	{
		aa = aa > 0.5*cp.lbox ? cp.lbox-aa : aa;
		bb = bb > 0.5*cp.lbox ? cp.lbox-bb : bb;
		cc = cc > 0.5*cp.lbox ? cp.lbox-cc : cc;
		dd = dd > 0.5*cp.lbox ? cp.lbox-dd : dd;
		ee = ee > 0.5*cp.lbox ? cp.lbox-ee : ee;
		ff = ff > 0.5*cp.lbox ? cp.lbox-ff : ff;
	}

	if( (aa <= lado[0] && bb <= lado[0]) &&
	    (cc <= lado[1] && dd <= lado[1]) &&
	    (ee <= lado[2] && ff <= lado[2]) )
#endif
#ifdef PERCENT
	if(drand48()<(double)FRACC)
#endif
  	{
	    positions[((long unsigned)NO_DIM *pc)]      = (Real)particles[i].x;
	    positions[((long unsigned)NO_DIM *pc) + 1]  = (Real)particles[i].y;
	    positions[((long unsigned)NO_DIM *pc) + 2]  = (Real)particles[i].z;
#ifdef VELOCITY
	    velocities[((long unsigned)NO_DIM*pc)]     = (Real)particles[i].vx;
	    velocities[((long unsigned)NO_DIM*pc) + 1] = (Real)particles[i].vy;
	    velocities[((long unsigned)NO_DIM*pc) + 2] = (Real)particles[i].vz;
#endif
#ifdef WEIGHT
            weights[pc] = (Real)particles[i].mp;
#endif
	    pc++;
        }
  }

	//fprintf(stdout,">> Guardo %lu particulas (%.2f %% de las particulas del file)\n",(long unsigned)pc,100.0f*(float)pc/(float)numToRead);
	//fprintf(stdout,">> Particulas que cumplen: %d\n>> Particulas totales: %d", pc,noParticles);

#if defined(CUT_REGION) || defined(PERCENT)
	    readData->position(noParticles, pc);  // particle positions
	#ifdef WEIGHT
	    readData->weight(noParticles,   pc);  // particle weights (weights = particle masses from the snapshot file)
	#endif
	#ifdef VELOCITY
	    readData->velocity(noParticles, pc);  // particle velocities
	#endif
	#ifdef SCALAR
	  if ( noScalars>0 )
	    readData->scalar(noParticles,   pc);  // particle scalar quantity
	#endif
#endif

#ifdef CUT_REGION
	free(limits);
	free(lado);
#endif

#ifdef PERCENT
	sprintf(tmp_fname,"Guardo %lu particulas representa %.2f %% de las particulas del snapshot\n",(long unsigned)pc,100.0f*(float)pc/(float)numToRead);
	YELLOW(tmp_fname);
#endif

  free(particles);

  return;
}

void leeheader(const char fname [], User_options *userOptions)
{
  // HEADER
  struct headerparts
  {
    int CanHaveTypes[7],Flag_Entropy_ICs[7],NumPart_Total[7],NumPart_Total_HighWord[7];
    int Dimension,NumFilesPerSnapshot,NumPartTypes,ThisFile,Virtual;
    long int NumPart_ThisFile[7],Shift[3],TotalNumberOfParticles[7];
    char Code[5],OutputType[10],RunName[15],SelectOutput[7],SnapshotDate[23],System[4],PartTypeNames[128];
    double Redshift,Scale_factor,Time,TimeBase_dloga,TimeBase_dt;
    double BoxSize[3],InitialMassTable[7],MassTable[7];
  } header;

  hid_t file, grupo, att, atttype;
  char labels[N_part_types][15] = {"Gas", "DM", "DMBackground", "Sink", "Stars", "BH", "Neutrino"};

  file  = H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT);
  grupo = H5Gopen(file,"/Header",H5P_DEFAULT);

  att = H5Aopen_name(grupo,"BoxSize");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,&header.BoxSize);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"NumFilesPerSnapshot");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,&header.NumFilesPerSnapshot);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"Redshift");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,&header.Redshift);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"Scale-factor");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,&header.Scale_factor);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"NumPart_ThisFile");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,header.NumPart_ThisFile);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"NumPart_Total");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,header.NumPart_Total);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"NumPart_Total_HighWord");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,header.NumPart_Total_HighWord);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"MassTable");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,header.MassTable);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"InitialMassTable");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,header.InitialMassTable);
  H5Tclose(atttype);
  H5Aclose(att);

  H5Gclose(grupo);

  // COSMOLOGIA
  struct cosmology
  {
    int Cosmological_run;
    double Critical_density,Critical_density_0,H,H0,Hubble_time,lb_time;
    double mbd,mmd,N_eff,N_nu,N_ur,Omega_b,Omega_cdm,Omega_g,Omega_k,Omega_lambda,Omega_m,Omega_nu,Omega_nu_0,Omega_r,Omega_ur;
    double Redshift,Scale_factor,T_CMB_0_K,T_CMB;
  } cosmology;
  
  grupo=H5Gopen(file,"/Cosmology",H5P_DEFAULT);

  att = H5Aopen_name(grupo,"H [internal units]");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,&cosmology.H);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"Omega_m");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,&cosmology.Omega_m);
  H5Tclose(atttype);
  H5Aclose(att);

  att = H5Aopen_name(grupo,"Omega_lambda");
  atttype = H5Aget_type(att);
  H5Aread(att,atttype,&cosmology.Omega_lambda);
  H5Tclose(atttype);
  H5Aclose(att);

  H5Gclose(grupo);
  H5Fclose(file);

  int i;
  for (i=0; i<N_part_types; ++i)
  {
    cp.npart[i] = (unsigned long long)header.NumPart_Total[i] | ((unsigned long long)header.NumPart_Total_HighWord[i] << 32);
    cp.Mpart[i] = header.InitialMassTable[i];
  }

  // Definicion estructura cosmoparam
  cp.nfiles    = header.NumFilesPerSnapshot;
  cp.omegam    = cosmology.Omega_m;
  cp.omegal    = cosmology.Omega_lambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = cosmology.H;
  cp.lbox      = header.BoxSize[0];
  //cp.Mpart    = 3.143E-4*cp.hparam;  /*A-5*/
  //cp.Mpart    = 3.929E-5*cp.hparam;  /*A-4*/
  cp.redshift  = header.Redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;
  fprintf(stdout,"\n--- LEEMOS USANDO HDF5 \n");

  for (int i=0; i<N_part_types; ++i)
  {
    if( userOptions->readParticleSpecies[i] )
    {
      fprintf(stdout,"*********************************** \n");
      fprintf(stdout,"*   Parametros de la simulacion   * \n");
      fprintf(stdout,"*********************************** \n");
      fprintf(stdout,"  Tipo de particula %s\n", labels[i]);
      fprintf(stdout,"  Numero de particulas = %llu \n", cp.npart[i]);
      fprintf(stdout,"  Lado del box = %g \n", cp.lbox);
      fprintf(stdout,"  Redshift = %g \n", cp.redshift);
      fprintf(stdout,"  Omega Materia = %g \n", cp.omegam);
      fprintf(stdout,"  Omega Lambda = %g \n", cp.omegal);
      fprintf(stdout,"  Parametro de Hubble = %g \n",cp.hparam);
      fprintf(stdout,"  Masa inicial por particula = %g \n",cp.Mpart[i]);
      fprintf(stdout,"*********************************** \n");
      fprintf(stdout,"*********************************** \n");
      fflush(stdout);
    }
  }
  
  

  return;
}

#else

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ***************************************************  GADGET  **************************************************************
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------

void leeheader(const char fname [], User_options *userOptions)
{
  struct io_header header;
  FILE *pf;
  int i,d1,d2;
#ifdef TYPE_TWO_GADGET
  int blocksize;
  char label[4];
  char labels[N_part_types][15] = {"Gas", "Dark Matter", "Disk", "Bulge", "Stars", "Bndry"};
#endif
#ifdef TYPE_FOUR_GADGET
  fprintf(stdout,"****** TYPE_FOUR_GADGET (SIMUS POLACO) ******\n");
  char labels[N_part_types][15] = {"Gas", "Dark Matter"};
#endif
  char newfname[1024];
  fprintf(stderr,"open file `%s`\n",fname);
  pf = fopen(fname,"r");
  if (!pf) 
  {
    snprintf(newfname, sizeof(newfname), "%s.0", fname);
    pf = fopen(newfname, "rb");
    if (!pf) 
    {
        fprintf(stderr, "Error: no se pudo abrir el archivo \"%s\" ni \"%s\".\n", fname, newfname);
        exit(EXIT_FAILURE);
    }
  }

#ifdef TYPE_TWO_GADGET
  fread(&d1, sizeof(d1), 1, pf);
  fread(&label,sizeof(char), 4, pf);
  fread(&blocksize, sizeof(int), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
#endif

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
#ifndef TYPE_FOUR_GADGET
  assert(d1==256);

  // Definicion estructura cosmoparam
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize;

  for (i=0; i<N_part_types; ++i)
  {
  	cp.npart[i] = (unsigned long long)header.npartTotal[i] | ((unsigned long long)header.npartTotalHighWord[i] << 32);
  	cp.Mpart[i] = header.mass[i];
  }

  //cp.Mpart    = 3.143E-4*cp.hparam;  //A-5
  //cp.Mpart    = 3.929E-5*cp.hparam;  //A-4
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;
#else
  cp.lbox = header.BoxSize;
  for (i=0; i<N_part_types; ++i)
  {
  	cp.npart[i] = (unsigned long long)header.npartTotal[i];
  	cp.Mpart[i] = header.mass[i];
  }
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
#endif

  fclose(pf);
  //fprintf(stdout," LEEMOS USANDO GADGET \n");

  for (i=0; i<N_part_types; ++i)
  {
    if( userOptions->readParticleSpecies[i] )
    {
  	fprintf(stdout,"*********************************** \n");
  	fprintf(stdout,"*   Parametros de la simulacion   * \n");
  	fprintf(stdout,"*********************************** \n");
  	fprintf(stdout,"  Tipo de particula %s\n", labels[i]);
  	fprintf(stdout,"  Numero de particulas = %llu \n", cp.npart[i]);
  	fprintf(stdout,"  Lado del box = %g \n", cp.lbox);
  	fprintf(stdout,"  Redshift = %g \n", cp.redshift);
#ifdef TYPE_FOUR_GADGET
  	fprintf(stdout,"  Omega Materia = ???? \n");
  	fprintf(stdout,"  Omega Lambda = ???? \n");
  	fprintf(stdout,"  Parametro de Hubble = ???? \n");
#else
  	fprintf(stdout,"  Omega Materia = %g \n", cp.omegam);
  	fprintf(stdout,"  Omega Lambda = %g \n", cp.omegal);
  	fprintf(stdout,"  Parametro de Hubble = %g \n",cp.hparam);
#endif
  	fprintf(stdout,"  Masa por particula = %g \n",cp.Mpart[i]);
  	fprintf(stdout,"*********************************** \n");
  	fprintf(stdout,"*********************************** \n");
  	fflush(stdout);
    }
  }
  return;
}

static void lee(const char *filename, const int type_particle, Read_data<Real> *readData, User_options *userOptions)
{
  char key[50];
  long long i, j, pc;
  long long elements;
  unsigned int d1, d2;
  unsigned int k, n;

  FILE *pf;
  #ifdef TYPE_TWO_GADGET
  unsigned int blocksize;
  char label[4];
  #endif

  struct io_header header;
  struct particle_data 
  {
    double      x;
    double      y;
    double      z;
    #ifdef WEIGHT
    double      mp;
    #endif
    #ifdef VELOCITY
    double      vx;
    double      vy;
    double      vz;
    #endif
    #ifdef IDS
    unsigned long long id;
    #endif
  } *particles;
  
  float r[3];
  #ifdef WEIGHT
  float mm;
  #endif
  #ifdef VELOCITY
  float v[3];
  #endif
  #ifdef IDS
  unsigned int id;
  #endif

  /****** CHEQUEAMOS SI ES UN FILE O VARIOS Y LEEMOS EL HEADER *******/
  char newfname[1024];
  pf = fopen(filename,"rb");
  if (!pf) 
  {
    snprintf(newfname, sizeof(newfname), "%s.0", filename);
    pf = fopen(newfname, "rb");
    if (!pf) 
    {
        fprintf(stderr, "Error: no se pudo abrir el archivo \"%s\" ni \"%s\".\n", filename, newfname);
        exit(EXIT_FAILURE);
    }
  }
  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  fclose(pf);
 
  /****** RESERVAMOS MEMORIA PARA LAS PARTICULAS *******/
  long long NumPart = header.npartTotal[type_particle];
  size_t total_memory = (float)NumPart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  printf("Allocating %.5zu Gb for %u particle properties\n",total_memory,NumPart);
  particles = (struct particle_data *) malloc(NumPart*sizeof(struct particle_data));
  assert(particles != NULL);	// Esta función chequea que la reserva de memoria se haya hecho correctamente

  int ifile, ind; 
  for(ifile = 0, ind = 0; ifile < header.num_files; ifile++)
  {
    if(header.num_files>1)
    {
      sprintf(newfname,"%s.%d",filename,ifile);
      fprintf(stdout,"---- MULTIPLE SNAPSHOT FILES ----\n"); fflush(stdout);
    }
    else
    {
      sprintf(newfname,"%s",filename);
      fprintf(stdout,"---- SINGLE SNAPSHOT FILE ----\n"); fflush(stdout);
    }
    pf = fopen(newfname,"r");
    if(pf == NULL)
    {
      fprintf(stderr,"can't open file `%s`\n",newfname);
      exit(EXIT_FAILURE);
    }

    fprintf(stdout,">> Reading file: %s \n",newfname); fflush(stdout);

    #ifdef TYPE_TWO_GADGET
    fread(&d1, sizeof(d1), 1, pf);
    fread(&label,sizeof(char), 4, pf);
    fread(&blocksize, sizeof(int), 1, pf);
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);
    fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
    #endif
 
    fread(&d1, sizeof(d1), 1, pf);
    fread(&header, sizeof(header), 1, pf);
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);
    #ifndef TYPE_FOUR_GADGET
    assert(d1==256);
    #endif

    #ifdef TYPE_TWO_GADGET
    fread(&d1, sizeof(d1), 1, pf);
    fread(&label,sizeof(char), 4, pf);
    fread(&blocksize, sizeof(int), 1, pf);
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);
    fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
    #endif

    fread(&d1, sizeof(d1), 1, pf);
    for(k = 0, pc = 0; k < N_part_types; k++)
    {
      for(n = 0; n < header.npart[k]; n++){
        fread(&r[0], sizeof(float), 3, pf);
        if(k == type_particle)
        {
          particles[pc+ind].x = r[0];
          particles[pc+ind].y = r[1];
          particles[pc+ind].z = r[2];
          pc++;
        }
      }
    }
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);

    #ifdef TYPE_TWO_GADGET
    fread(&d1, sizeof(d1), 1, pf);
    fread(&label,sizeof(char), 4, pf);
    fread(&blocksize, sizeof(int), 1, pf);
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);
    fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
    #endif

    fread(&d1, sizeof(d1), 1, pf);
    #ifdef VELOCITY
    for(k = 0, pc = 0; k < N_part_types; k++){
      for(n = 0; n < header.npart[k]; n++){
        fread(&v[0], sizeof(float), 3, pf);
        if(k == type_particle){
          particles[pc+ind].vx = v[0];
          particles[pc+ind].vy = v[1];
          particles[pc+ind].vz = v[2];
          pc++;
        }
      }
    }
    #else
    fseek(pf,d1,SEEK_CUR);
    #endif
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);

    #ifdef TYPE_TWO_GADGET
    fread(&d1, sizeof(d1), 1, pf);
    fread(&label,sizeof(char), 4, pf);
    fread(&blocksize, sizeof(int), 1, pf);
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);
    fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
    #endif

    fread(&d1, sizeof(d1), 1, pf);
    #ifdef IDS   
    for(k = 0, pc = 0; k < N_part_types; k++){
      for(n = 0; n < header.npart[k]; n++){
        fread(&id, sizeof(unsigned int), 1, pf);
        if(k == type_particle){
          particles[pc+ind].id = id;
          pc++;
        }
      }
    }
    #else
    fseek(pf,d1,SEEK_CUR);
    #endif
    fread(&d2, sizeof(d2), 1, pf);
    assert(d1==d2);

    #ifdef WEIGHT
    if(type_particle != 1)
    {
      #ifdef TYPE_TWO_GADGET
      fread(&d1, sizeof(d1), 1, pf);
      fread(&label,sizeof(char), 4, pf);
      fread(&blocksize, sizeof(int), 1, pf);
      fread(&d2, sizeof(d2), 1, pf);
      assert(d1==d2);
      fprintf(stdout,"In %s Read %c%c%c%c\n",filename,label[0],label[1],label[2],label[3]);fflush(stdout);
      #endif
      fread(&d1, sizeof(d1), 1, pf);
      for(k = 0, pc = 0; k < N_part_types; k++)
      {
        for(n = 0; n < header.npart[k]; n++)
        {
          fread(&mm, sizeof(float), 1, pf);
          if(k == type_particle)
          {
            particles[pc+ind].mp = mm;
            pc++;
          }
        }
      }
      fread(&d2, sizeof(d2), 1, pf);
      assert(d1==d2);
    }
    else
    {
      for(n=0, pc=0;n<header.npart[type_particle];n++)
      {
        particles[pc+ind].mp = cp.Mpart[type_particle];
        pc++;
      }
    }
    #endif
    fclose(pf);
    ind += header.npart[type_particle];
  }

  size_t noParticles = (size_t)header.npartTotal[type_particle];
  #ifdef SCALAR
  int noScalars = 0;
  for (i=3; i<userOptions->readParticleData.size(); ++i)
      if ( userOptions->readParticleData[i] )
          noScalars += 1;
  #endif
	
  readData->position( noParticles );        // particle positions
  #ifdef WEIGHT
  readData->weight( noParticles );        // particle weights (weights = particle masses from the snapshot file)
  #endif
  #ifdef VELOCITY
  readData->velocity( noParticles );        // particle velocities
  #endif
  #ifdef SCALAR
  if ( noScalars>0 )
  	readData->scalar( noParticles );        // particle scalar quantity
  #endif

  // copy the relevant information to the '*readData' structure
  Real *positions = readData->position();  //particle positions
  #ifdef WEIGHT
  Real *weights = readData->weight();    //particle weights (e.g. weights = particle masses)
  #endif
  #ifdef VELOCITY
  Real *velocities = readData->velocity();  //particle velocity
  #endif

  #ifdef CUT_REGION
  double *limits;
  double *lado;
  const double delta = userOptions->paddingParticles / (Real)pow( (double)cp.npart[1], 1./NO_DIM) * cp.lbox;
  
  limits = (double *) malloc(2*NO_DIM*sizeof(double));
  lado   = (double *) malloc(  NO_DIM*sizeof(double));

  for(i=0; i<NO_DIM; i++)
  {
    limits[2*i]   = userOptions->region.coords[2*i]*cp.lbox + pow(-1,2*i+1)*delta;
    limits[2*i+1] = userOptions->region.coords[2*i+1]*cp.lbox + pow(-1,2*i)*delta;
    if(!userOptions->periodic)
    {
	limits[2*i]   = limits[2*i]   < 0       ? 0       : limits[2*i];
	limits[2*i+1] = limits[2*i+1] > cp.lbox ? cp.lbox : limits[2*i+1];
    }
    lado[i] = limits[2*i+1] - limits[2*i];
  }
  #endif

  for(i=0, pc=0; i<header.npartTotal[type_particle]; i++)
  {
    #ifdef CUT_REGION
    double aa = fabs(limits[0] - particles[i].x);
    double bb = fabs(limits[1] - particles[i].x);
    double cc = fabs(limits[2] - particles[i].y);
    double dd = fabs(limits[3] - particles[i].y);
    double ee = fabs(limits[4] - particles[i].z);
    double ff = fabs(limits[5] - particles[i].z);
    
    if(userOptions->periodic)
    {
    	aa = aa > 0.5*cp.lbox ? cp.lbox-aa : aa;
    	bb = bb > 0.5*cp.lbox ? cp.lbox-bb : bb;
    	cc = cc > 0.5*cp.lbox ? cp.lbox-cc : cc;
    	dd = dd > 0.5*cp.lbox ? cp.lbox-dd : dd;
    	ee = ee > 0.5*cp.lbox ? cp.lbox-ee : ee;
    	ff = ff > 0.5*cp.lbox ? cp.lbox-ff : ff;
    }

    if((aa <= lado[0] && bb <= lado[0]) &&
       (cc <= lado[1] && dd <= lado[1]) &&
       (ee <= lado[2] && ff <= lado[2]))
    #endif

    #ifdef PERCENT
    if(drand48()<(double)FRACC)
    #endif
    {
      positions[((long unsigned)NO_DIM *pc)]      = (Real)particles[i].x;
      positions[((long unsigned)NO_DIM *pc) + 1]  = (Real)particles[i].y;
      positions[((long unsigned)NO_DIM *pc) + 2]  = (Real)particles[i].z;
      #ifdef VELOCITY
      velocities[((long unsigned)NO_DIM*pc)]     = (Real)particles[i].vx;
      velocities[((long unsigned)NO_DIM*pc) + 1] = (Real)particles[i].vy;
      velocities[((long unsigned)NO_DIM*pc) + 2] = (Real)particles[i].vz;
      #endif
      #ifdef WEIGHT
      weights[pc] = (Real)particles[i].mp;
      #endif
      pc++;
    }
  }

  fprintf(stdout,"%d\n", pc);

  #if defined(CUT_REGION) || defined(PERCENT)
    readData->position(noParticles, pc);  	  // particle positions
    #ifdef WEIGHT
    readData->weight(noParticles,   pc);  // particle weights (weights = particle masses from the snapshot file)
    #endif
    #ifdef VELOCITY
    readData->velocity(noParticles, pc);  // particle velocities
    #endif
    #ifdef SCALAR
    if ( noScalars>0 )
      readData->scalar(noParticles,   pc);  // particle scalar quantity
    #endif
  #endif

  #ifdef CUT_REGION
    free(limits);
    free(lado);
  #endif

  #ifdef PERCENT
    sprintf(tmp_fname,"Guardo %lu particulas representa %.2f %% de las particulas del snapshot\n",(long unsigned)pc,100.0f*(float)pc/(float)header.npartTotal[type_particle]);
    YELLOW(tmp_fname);
  #endif

  free(particles);
  return;
}
#endif

void readMyFile(std::string filename, Read_data<Real> *readData, User_options *userOptions)
{
  char fname[200];
  #ifdef PERCENT
  srand48(80);
  #endif
  MESSAGE::Message message( userOptions->verboseLevel );

  strcpy(fname, filename.c_str()); 
 
  leeheader(fname, userOptions);

  if(not userOptions->userGivenBoxCoordinates )
  {
    for (size_t i=0; i<NO_DIM; ++i)
    {
      userOptions->boxCoordinates[2*i]   = 0.;      // this is the left extension of the full box
      userOptions->boxCoordinates[2*i+1] = cp.lbox; // right extension of the full box
    }
  }
  else
    message << "The box coordinates were set by the user using the program options. The program will keep this values and will NOT use the box length information from the Gadget file!" << MESSAGE::Flush;

  #ifdef CUT_REGION		
  assert(userOptions->regionOn);
  #endif

  for (int i=0; i<N_part_types; ++i)
  {
    if( userOptions->readParticleSpecies[i] )
    {
      lee(fname, i, readData, userOptions);
    }
  }
  message << "---- Done ----\n";
  return;
}

/* This function writes the data to a binary file.*/
template <typename T>
void writeMyFile(T const &dataToWrite,
                     std::string filename,
                     std::string variableName,
                     User_options const &userOptions)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    message << "Writing the " << variableName << " to the binary file '" << filename << "' ...  " << MESSAGE::Flush;
    
    // open the file
    std::fstream outputFile;
    openOutputBinaryFile( outputFile, filename );
    
    // write the data to file
    // it is very simple to write data to a binary file, just do: outputFile.write( reinterpret_cast<char const *>(&(dataToWrite[0])), dataSize ). The following implements a recursive method that deals with possible issuess in writing very large data sets.
    size_t dataSize = dataToWrite.size()*sizeof(dataToWrite[0]);    //total number of data bytes to be written
    size_t maxSize = 256*256*256;   // write at most 256^3 elements at once, otherwise the write function might fail
    size_t noRepeats = size_t( dataToWrite.size() / maxSize ), currentPosition = 0;
    size_t tempBuffer = maxSize * sizeof(dataToWrite[0]);
    for (size_t i=0; i<noRepeats; ++i)  //write in blocks of 256^3 elements
    {
        outputFile.write( reinterpret_cast<char const *>(&(dataToWrite[currentPosition])), tempBuffer );
        currentPosition += maxSize;
    }
    tempBuffer = (dataToWrite.size() - currentPosition) * sizeof(dataToWrite[0]);    // write everything else that is left
    outputFile.write( reinterpret_cast<char const *>(&(dataToWrite[currentPosition])), tempBuffer );
    
    checkFileOperations( outputFile, "write to" );   // check that the data writing was succesful
    outputFile.close();
    
    message << "Done.\n";
}
