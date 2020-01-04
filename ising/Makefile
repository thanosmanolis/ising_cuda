#############################################################
#    C/C++ Makefile											#
#    Author: Thanos Manolis <thanos.m14@gmail.com			#
#############################################################
#															#
#   'make lib'	  	  build the libraries .a				#
#   'make clean'  	  removes .o .a and executable files    #
#															#
#############################################################

# define the C/C++ compiler to use,default here is gcc-7
CC = gcc-7

# define flags
CFLAGS = -O3
LDFLAGS =

# define command to remove files
RM = rm -rf

# always build those, even if "up-to-date"
.PHONY: lib

lib:
	cd src; $(CC) -c ising_sequential.c $(CFLAGS); cd ..
	cd src; ar rcs ../lib/ising_sequential.a ising_sequential.o; cd ..

clean:
	$(RM) src/*.o lib/*.a $(EXECS)
