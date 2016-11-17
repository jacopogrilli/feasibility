COMP = gcc
LIB = -L/opt/local/lib/ -lgsl -lgslcblas -lm
OPT = -Wall -O3 -DHAVE_INLINE -I/opt/local/include/
#OPT = -Wall -O0 -g -DHAVE_INLINE -I/opt/local/include/ # needed for profiling using valgrind

OBJECTS = StructStability.o SpectralStats.o main.o

StructuralStab: ${OBJECTS}
	${COMP} ${OBJECTS} -o StructuralStab ${OPT} ${LIB}

SpectralStats.o: SpectralStats.c
	${COMP} SpectralStats.c -c ${OPT}

StructStability.o: StructStability.c SpectralStats.h
	${COMP} StructStability.c -c ${OPT}

main.o: main.c SpectralStats.h StructStability.h
	${COMP} main.c -c ${OPT}

clean: 
	rm -rf *.o *.*~ *.log StructuralStab
