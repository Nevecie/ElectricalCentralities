#COMPILER
FC = gfortran #mpif90 #ifort #
FC_MPI = mpif90 #ifort #
#FLAGS
CHECKFLAG = #-fbacktrace -fbounds-check  #-fcheck=all   #-CB 
PARALLELFLAG = -L/usr/local/lib/# -mkl #-fopenmp
LIBFLAGS =  -lpmi
LONGLINE = -ffree-line-length-none
COMPFLAGS = $(LONGLINE) $(CHECKFLAG) $(PARALLELFLAG) #-heap-arrays
ALLFLAGS = $(COMPFLAGS) #-heap-arrays
#OTHER DETAILS
NAME = bulat.gmres
NAME_MPI = bulat.gmres.mpi
NCPUs = 2
CPUs = -n $(NCPUs)
NODEPAT = node-
NODELIST = --nodelist=$(NODEPAT)
NODEIDS = 0[3-9]
NODES = #$(NODELIST)$(NODEIDS)
SCRIPT = bulat.run
LAUNCHER = sbatch
MKLLIB = #/share/home/golovin/intel_11.3/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64_lin
MKLINCLUDE = #/share/home/golovin/intel_11.3/compilers_and_libraries_2016.3.210/linux/mkl/include
COMPILERLIB = #/share/home/golovin/intel_11.3/compilers_and_libraries_2016.3.210/linux/compiler/lib/intel64_lin

#SOURCE = bulat.hi.f90
SOURCE = bulat.gmres.f90
SOURCE_MPI = bulat.gmres.mpi.f90
#SOURCE = bulat.f90

all: $(NAME)

$(NAME): blas.f90 ilut.f90 $(SOURCE)
#	export LD_LIBRARY_PATH="/share/home/golovin/intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64":/share/mpi/mvapich-1.1.0-3143-intel/lib/shared:/opt/intel/mkl/10.0.011/lib/em64t:/opt/intel/itac/7.1/itac/slib_impi3:/opt/intel//fce/10.1.008/lib:/opt/intel//cce/10.1.008/lib:/share/home/golovin/intel/composer_xe_2011_sp1.7.256/mpirt/lib/intel64:/share/home/golovin/intel/composer_xe_2011_sp1.7.256/debugger/lib/intel64:/usr/lib
#	$(FC)  -I $(MKLINCLUDE) -w $^ -Wl,--start-group $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_scalapack_lp64.a $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a $(COMPILERLIB)/libiomp5.a -Wl,--end-group -lpthread -L $(ALLFLAGS) -g -o $(NAME) $(LIBFLAGS)
	$(FC) $(LONGLINE) -o $(NAME) $^

$(NAME_MPI): blas.f90 ilut.f90 $(SOURCE_MPI)
#	export LD_LIBRARY_PATH="/share/home/golovin/intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64":/share/mpi/mvapich-1.1.0-3143-intel/lib/shared:/opt/intel/mkl/10.0.011/lib/em64t:/opt/intel/itac/7.1/itac/slib_impi3:/opt/intel//fce/10.1.008/lib:/opt/intel//cce/10.1.008/lib:/share/home/golovin/intel/composer_xe_2011_sp1.7.256/mpirt/lib/intel64:/share/home/golovin/intel/composer_xe_2011_sp1.7.256/debugger/lib/intel64:/usr/lib
#	$(FC)  -I $(MKLINCLUDE) -w $^ -Wl,--start-group $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_scalapack_lp64.a $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a $(COMPILERLIB)/libiomp5.a -Wl,--end-group -lpthread -L $(ALLFLAGS) -g -o $(NAME) $(LIBFLAGS)
	$(FC_MPI) $(ALLFLAGS) -g -o $(NAME_MPI) $^ $(LIBFLAGS)

go: $(NAME_MPI)
	$(LAUNCHER) $(CPUs) $(NODES) $(SCRIPT) 

help:
	#echo "make builds bulat.exe; bulat.run starts it with srun. Parameters are in params.nml; data are in files, their names are in params.nml"
	echo "make builds the executable bulat.gmres. It uses gmres.nml list of parameters. It can also accept an integer command-line argument for the position of 1 in the right-hand side vector. Result is written to the results.dat file (or other, see param.nml). You can compose a script to solve problems one by one."

