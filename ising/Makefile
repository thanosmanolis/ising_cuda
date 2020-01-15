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

# all the libraries
LIBS = lib_seq lib_v1 lib_v2 lib_v3

# always build those, even if "up-to-date"
.PHONY: $(LIBS)

lib: $(LIBS)

lib_seq:
	cd src; $(CC) -c ising_sequential.c $(CFLAGS); cd ..
	cd src; ar rcs ../lib/lib_seq.a ising_sequential.o; cd ..

lib_v1:
	cd src; $(NVCC) -c ising_v1.cu $(CUDAFLAGS); cd ..
	cd src; ar rcs ../lib/lib_v1.a ising_v1.o; cd ..

lib_v2:
	cd src; $(NVCC) -c ising_v2.cu $(CUDAFLAGS); cd ..
	cd src; ar rcs ../lib/lib_v2.a ising_v2.o; cd ..

lib_v3:
	cd src; $(NVCC) -c ising_v3.cu $(CUDAFLAGS); cd ..
	cd src; ar rcs ../lib/lib_v3.a ising_v3.o; cd ..

clean:
	$(RM) src/*.o lib/*.a
