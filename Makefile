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
EXECS = test_sequential test_v1

# define flags
CFLAGS =
LDFLAGS =

# define command to remove files
RM = rm -rf

# always build those, even if "up-to-date"
.PHONY: $(EXECS)

all: $(EXECS)

test_sequential:
	$(CC) src/ising_sequential.c -o $@ $(CFLAGS) $(LDFLAGS)

test_v1:
	$(NVCC) src/ising_v1.cu -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	$(RM) *.h *.a *.o $(EXECS)
