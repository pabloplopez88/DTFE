# Makefile for compiling the DTFE code on Linux systems
#
# QUICK START (any cluster/computer, see README.md):
#     conda env create -f environment.yml    # only the first time
#     conda activate dtfe
#     make
#
# By default all the libraries (GSL, Boost, CGAL, GMP, MPFR and HDF5) and the C++ compiler are taken from
# the active conda environment, so you do not need to edit any path in this file.


# Directory where the libraries are installed (the program adds the '/lib' and '/include' parts automatically).
# Default: the active conda environment. If you want to use libraries installed somewhere else, set LIB_PREFIX
# (one directory for all the libraries) or the individual paths below, e.g.:  make BOOST_PATH=/opt/boost
# NOTE: these variables are taken only from this file or from the 'make' command line, NOT from environment variables with the
# same name (e.g. many clusters define HDF5_PATH in the shell, which would silently mix the cluster libraries with the conda ones).
LIB_PREFIX = $(CONDA_PREFIX)
GSL_PATH   = $(LIB_PREFIX)
BOOST_PATH = $(LIB_PREFIX)
CGAL_PATH  = $(LIB_PREFIX)
# path to the GMP and MPFR libraries (needed by CGAL)
MPRF_PATH  = $(LIB_PREFIX)
# path to the HDF5 library (it must include the C++ interface, i.e. 'H5Cpp.h' and 'libhdf5_cpp')
HDF5_PATH  = $(LIB_PREFIX)
# support for reading HDF5 gadget files: 'auto' (enabled if the HDF5 C++ library is found), 'yes' or 'no'
USE_HDF5   = auto

# C++ compiler - preferably a version that supports OpenMP. Inside a conda environment with the 'cxx-compiler'
# package, the variable CXX already points to the conda compiler; otherwise 'g++' is used.
CC = $(CXX)


# paths to where to put the object files and the executables files. If you build the DTFE library than you also need to specify the directory where to put the library and the directory where to copy the header files needed by the library (choose an empty directory for the header files).
OBJ_DIR = ./o
BIN_DIR = ./
LIB_DIR = ./
INC_DIR = ./DTFE_include



############################# Choose the compiler directives ##################################

############################# Overall options ##################################
OPTIONS = 
#------------------------ set the number of spatial dimensions (2 or 3 dimensions)
OPTIONS += -DNO_DIM=3 
#------------------------ set type of variables - float (comment the next line) or double (uncomment the next line)
# OPTIONS += -DDOUBLE 

############################# Quantities to be computed ##################################
#------------------------ set which quantities can be computed (can save memory by leaving some out)
# These two options can be switched off to save memory (they can also be given on the command line, e.g. 'make VELOCITY=no SCALARS=no'):
#   VELOCITY = yes : velocity and its derived fields (gradient, divergence, shear, vorticity); 12 bytes per particle
#   SCALARS  = yes : additional fields stored in the 'scalar' variable with 3 components, used for the magnetic field; 12 bytes per particle
# Memory per particle: 44 bytes with both (default), 32 bytes with SCALARS=no, 20 bytes with VELOCITY=no SCALARS=no.
# E.g. for a density map of 2500^3 particles compile with 'make VELOCITY=no SCALARS=no' (~310 GB instead of ~690 GB for the particles).
VELOCITY = yes
SCALARS  = yes
ifeq ($(VELOCITY),yes)
	OPTIONS += -DVELOCITY
endif
ifeq ($(SCALARS),yes)
	OPTIONS += -DSCALAR -DNO_SCALARS=3
endif

############################# Input and output operations default settings ##################################
#------------------------ set which are the default input and output functions for doing data io
# default function to read the input data (101-multiple gadget file, 102-single gadget file, 105-HDF5 gadget file, 111-text file, ... see documentation for more options). The input file type can be set during runtime using the option '--input'. This makefile option only sets a default input file in the case none is given via the program options.
OPTIONS += -DINPUT_FILE_DEFAULT=101 
# default value for the units of the input data (value=what is 1 Mpc in the units of the data - in this example the data is in kpc). You can change this also during runtime using the program option '--MpcUnit'.
OPTIONS += -DMPC_UNIT=1000. 
# default function to write the output data (101-binary file, 111-text file, ... see documentation for more options). The output file type can be set during runtime using the option '--output'. This makefile option only sets a default output file in the case none is given via the program options.
OPTIONS += -DOUTPUT_FILE_DEFAULT=101  
#101 for binary file, 100 my density file

############################# additional compiler options ##################################
# enable this option if to use OpenMP (share the workload between CPU cores sharing the same RAM)
OPTIONS += -DOPEN_MP 
# enable to check if the padding gives a complete Delaunay Tesselation of the region of interest
# OPTIONS += -DTEST_PADDING 
# enable this option to shift from position space to redshift space; You also need to activate this option during run-time using '--redshift-space arguments'
OPTIONS += -DREDSHIFT_SPACE

#------------------------ options usefull when using DTFE as a library
# uncomment the line to get access to a function that returns the Delaunay triangulation of the point set
# OPTIONS += -DTRIANGULATION 


############################# Help menu messages options ##################################
#------------------------ compiler directive that affect only the help messages when using the '-h / --help' option (it does not affect the program in any other way)- if the option is uncommented, than it will show that set of options in the help menu
OPTIONS += -DFIELD_OPTIONS 
OPTIONS += -DREGION_OPTIONS 
# OPTIONS += -DPARTITION_OPTIONS 
# OPTIONS += -DPADDING_OPTIONS 
OPTIONS += -DAVERAGING_OPTIONS 
# OPTIONS += -DREDSHIFT_CONE_OPTIONS 
OPTIONS += -DADDITIONAL_OPTIONS 








###############  DO NOT MODIFY BELOW THIS LINE  ###########################
# do not modify below this line
SRC = ./src

# list of the library directories given by the user (duplicates removed)
LIB_DIRS := $(sort $(strip $(GSL_PATH) $(BOOST_PATH) $(CGAL_PATH) $(MPRF_PATH)))
ifeq ($(USE_HDF5),auto)
	USE_HDF5 := $(if $(wildcard $(strip $(HDF5_PATH))/include/H5Cpp.h),yes,no)
endif
ifeq ($(USE_HDF5),yes)
	LIB_DIRS := $(sort $(LIB_DIRS) $(strip $(HDF5_PATH)))
	OPTIONS += -DHDF5
	HDF5_LIB = -lhdf5_cpp -lhdf5
endif

INCLUDES  = $(foreach dir,$(LIB_DIRS),-I$(dir)/include)
# '-rpath' stores the library paths in the executable, so it runs without having to set LD_LIBRARY_PATH
# '--disable-new-dtags' makes the stored paths take precedence over LD_LIBRARY_PATH, so the program always uses the libraries it
# was compiled with, even if LD_LIBRARY_PATH points to other versions of them (e.g. cluster modules or libraries in the home directory)
LIBRARIES = $(foreach dir,$(LIB_DIRS),-L$(dir)/lib -Wl,-rpath,$(dir)/lib) -Wl,--disable-new-dtags

COMPILE_FLAGS = -frounding-math -O3 -fopenmp -DNDEBUG $(OPTIONS)
DTFE_INC = $(INCLUDES)
# the following libraries should work in most cases (CGAL >= 5 is header-only, so there is no '-lCGAL')
DTFE_LIB = $(LIBRARIES) $(HDF5_LIB) -lboost_filesystem -lboost_program_options -lgsl -lgslcblas -lmpfr -lgmp -lm



IO_SOURCES = $(addprefix io/, input_output.h memory_estimate.h gadget_reader_header.cc gadget_reader_binary.cc gadget_reader_HDF5.cc gadget_reader_HDF5_Cristian.cc gadget_reader_MOG.cc hdf5_input_my_DESI.cc text_io.cc binary_io.cc my_io.cc)
MAIN_SOURCES = main.cpp DTFE.h message.h user_options.h input_output.cc $(IO_SOURCES)
DTFE_SOURCES = DTFE.cpp define.h particle_data.h user_options.h box.h quantities.h user_options.cc quantities.cc subpartition.h random.cc CIC_interpolation.cc TSC_interpolation.cc SPH_interpolation.cc kdtree/kdtree2.hpp Pvector.h message.h miscellaneous.h
TRIANG_SOURCES = $(addprefix CGAL_triangulation/, triangulation.cpp triangulation_miscellaneous.cc unaveraged_interpolation.cc averaged_interpolation_1.cc averaged_interpolation_2.cc padding_test.cc CGAL_include_2D.h CGAL_include_3D.h vertexData.h particle_data_traits.h) define.h particle_data.h user_options.h box.h quantities.h Pvector.h message.h math_functions.h

ALL_FILES = $(DTFE_SOURCES) $(TRIANG_SOURCES) $(MAIN_SOURCES) kdtree/kdtree2.hpp kdtree/kdtree2.cpp
LIB_FILES = $(DTFE_SOURCES) $(TRIANG_SOURCES)

HEADERS_1 = DTFE.h define.h user_options.h particle_data.h quantities.h Pvector.h math_functions.h  message.h box.h miscellaneous.h
HEADERS_2 = $(addprefix CGAL_triangulation/, CGAL_include_2D.h CGAL_include_3D.h vertexData.h particle_data_traits.h)



DTFE: set_directories $(OBJ_DIR)/DTFE.o $(OBJ_DIR)/triangulation.o $(OBJ_DIR)/main.o $(OBJ_DIR)/kdtree2.o Makefile
	$(CC) $(COMPILE_FLAGS) $(OBJ_DIR)/DTFE.o $(OBJ_DIR)/triangulation.o $(OBJ_DIR)/main.o $(OBJ_DIR)/kdtree2.o $(DTFE_LIB) -o $(BIN_DIR)/DTFE


$(OBJ_DIR)/main.o: $(addprefix $(SRC)/, $(MAIN_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/main.o -c $(SRC)/main.cpp

$(OBJ_DIR)/DTFE.o: $(addprefix $(SRC)/, $(DTFE_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/DTFE.o -c $(SRC)/DTFE.cpp

$(OBJ_DIR)/kdtree2.o: $(SRC)/kdtree/kdtree2.hpp $(SRC)/kdtree/kdtree2.cpp Makefile
	$(CC) -O3 -ffast-math -fomit-frame-pointer $(DTFE_INC) -o $(OBJ_DIR)/kdtree2.o -c $(SRC)/kdtree/kdtree2.cpp

$(OBJ_DIR)/triangulation.o: $(addprefix $(SRC)/, $(TRIANG_SOURCES)) Makefile
	$(CC) $(COMPILE_FLAGS) $(DTFE_INC) -o $(OBJ_DIR)/triangulation.o -c $(SRC)/CGAL_triangulation/triangulation.cpp


library: set_directories set_directories_2 $(addprefix $(SRC)/, $(LIB_FILES) ) copy_headers Makefile
	$(CC) $(COMPILE_FLAGS) -fPIC $(DTFE_INC) -o $(OBJ_DIR)/DTFE_l.o -c $(SRC)/DTFE.cpp
	$(CC) -O3 -ffast-math -fomit-frame-pointer -fPIC $(DTFE_INC) -o $(OBJ_DIR)/kdtree2_l.o -c $(SRC)/kdtree/kdtree2.cpp
	$(CC) $(COMPILE_FLAGS) -fPIC $(DTFE_INC) -o $(OBJ_DIR)/triangulation_l.o -c $(SRC)/CGAL_triangulation/triangulation.cpp
	$(CC) $(COMPILE_FLAGS) -shared $(OBJ_DIR)/DTFE_l.o $(OBJ_DIR)/triangulation_l.o $(OBJ_DIR)/kdtree2_l.o $(DTFE_LIB) -o $(LIB_DIR)/libDTFE.so


clean:
	rm -f $(BIN_DIR)/DTFE $(OBJ_DIR)/*.o

copy_headers:
	cp $(addprefix $(SRC)/, $(HEADERS_1)) $(INC_DIR)
	cp $(addprefix $(SRC)/, $(HEADERS_2)) $(INC_DIR)/CGAL_triangulation

set_directories:
	@ if !( test -d $(OBJ_DIR) ); \
	then mkdir $(OBJ_DIR); \
	fi
	@ if !( test -d $(BIN_DIR) ); \
	then mkdir $(BIN_DIR); \
	fi

set_directories_2:
	@ if !( test -d $(LIB_DIR) ); \
	then mkdir $(LIB_DIR); \
	fi
	@ if !( test -d $(INC_DIR) ); \
	then mkdir $(INC_DIR); \
	fi
	@ if !( test -d $(INC_DIR)/CGAL_triangulation ); \
	then mkdir $(INC_DIR)/CGAL_triangulation; \
	fi
