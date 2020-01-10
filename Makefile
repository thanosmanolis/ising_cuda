#############################################################
#    C/C++ Makefile											#
#    Author: Thanos Manolis <thanos.m14@gmail.com			#
#############################################################
#															#
#   'make'  		  build all executable files			#
#   'make exec_name'  build executable file 'test_*'		#
#   'make clean'  	  removes .o .a and executable files    #
#															#
#############################################################

# define the C/C++ compiler to use, default here is gcc-7
CC = gcc-7

# define the CUDA compiler to use
NVCC = nvcc

# all the executables
EXECS = sequential v1 v2

# define flags
CFLAGS =
CUDAFLAGS =

# define command to remove files
RM = rm -rf

# always build those, even if "up-to-date"
.PHONY: $(EXECS)

all: $(EXECS)

sequential:
	cd ising; make lib; cd ..
	cd ising; cp lib/*.a inc/*.h ../; cd ..
	$(CC) main.c ising_sequential.a -o $@ $(CFLAGS)

v1:
	cd ising; make lib; cd ..
	cd ising; cp lib/*.a inc/*.h ../; cd ..
	$(NVCC) main_cuda.cu ising_v1.a -o $@ $(CUDAFLAGS)

v2:
	cd ising; make lib; cd ..
	cd ising; cp lib/*.a inc/*.h ../; cd ..
	$(NVCC) main_cuda.cu ising_v2.a -o $@ $(CUDAFLAGS)

clean:
	$(RM) *.h *.a ising/src/*.o ising/lib/*.a $(EXECS)
