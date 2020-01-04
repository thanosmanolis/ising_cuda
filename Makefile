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

# define the C/C++ compiler to use,default here is gcc-7
CC = gcc-7

# all the executables
EXECS = test_sequential

# define flags
CFLAGS =
LDFLAGS =

# define command to remove files
RM = rm -rf

# always build those, even if "up-to-date"
.PHONY: $(EXECS)

all: $(EXECS)

test_sequential:
	cd ising; make lib; cd ..
	cd ising; cp lib/*.a inc/*.h ../; cd ..
	$(CC) tester.c ising_sequential.a -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	$(RM) *.h *.a $(EXECS)
