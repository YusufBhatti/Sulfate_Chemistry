F90=$(fortran)
#OPT = -O

OBJS1 = main.lltoeq.o lltoeq.o
OBJS2 = main.eqtoll.o eqtoll.o

lltoeq: ${OBJS1}
	${F90} -o lltoeq.exec ${OBJS1}; \
	mv lltoeq.exec ../bin/.
main.lltoeq.o: main.lltoeq.F90
	${F90} ${OPT} -c main.lltoeq.F90
lltoeq.o: lltoeq.F90
	${F90} ${OPT} -c lltoeq.F90

eqtoll: ${OBJS2}
	${F90} -o eqtoll.exec ${OBJS2}; \
	mv eqtoll.exec ../bin/.
main.eqtoll.o: main.eqtoll.F90
	${F90} -c main.eqtoll.F90
eqtoll.o: eqtoll.F90
	${F90} -c eqtoll.F90

coasts: coasts.o
	${F90} -o coasts.exec coasts.o; \
	mv coasts.exec ../bin/.; \
	cd ../bin; \
	coasts.exec; \
        cd ../source
coasts.o: coasts.F90
	${F90} -c coasts.F90
