#############################################################
#    C/C++ Makefile											#
#    Author: Thanos Manolis <thanos.m14@gmail.com			#
#############################################################
#															#
#   'make lib'	  	  build the libraries .a				#
#   'make clean'  	  removes .o .a files				    #
#															#
#############################################################

# define the C/C++ compiler to use, default here is gcc-7
CC = gcc-7

# define the CUDA compiler to use
NVCC = nvcc

# define flags
CFLAGS =
CUDAFLAGS = -arch=sm_35

# define command to remove files
RM = rm -rf

# always build those, even if "up-to-date"
.PHONY: lib

lib:
	cd src; $(CC) -c ising_sequential.c $(CFLAGS); cd ..
	cd src; $(NVCC) -c ising_v1.cu $(CUDAFLAGS); cd ..
	cd src; $(NVCC) -c ising_v2.cu $(CUDAFLAGS); cd ..
	cd src; ar rcs ../lib/ising_sequential.a ising_sequential.o; cd ..
	cd src; ar rcs ../lib/ising_v1.a ising_v1.o; cd ..
	cd src; ar rcs ../lib/ising_v2.a ising_v2.o; cd ..

clean:
	$(RM) src/*.o lib/*.a
