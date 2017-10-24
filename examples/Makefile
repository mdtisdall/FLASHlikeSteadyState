
MKLROOT = /opt/intel/mkl


CPPFLAGS += -DMKL_ILP64 -m64 
CPPFLAGS += -I${MKLROOT}/include 

CXXFLAGS += -g
#CXXFLAGS += -O3

# Flags passed to the C++ compiler.
CXXFLAGS += -Wall -Wextra -pthread -std=c++14
CXXFLAGS += -DMKL_DIRECT_CALL_SEQ
CXXFLAGS += -I../

#LDFLAGS += -L ~/lib
LDLIBS += ${MKLROOT}/lib/libmkl_intel_ilp64.a
LDLIBS += ${MKLROOT}/lib/libmkl_sequential.a
LDLIBS += ${MKLROOT}/lib/libmkl_core.a
LDFLAGS += -lpthread -lm -ldl

all : test 

test : ../blas_local_MKL.o ../mkl_ops_mkl.o

clean:
	rm -rf *.o
	rm -rf test
