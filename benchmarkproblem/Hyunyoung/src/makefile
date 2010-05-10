CC = g++
DEBUG = -g
INCLUDEDIR = ../include/
CFLAGS = -Wall -c $(DEBUG) -I$(INCLUDEDIR)
LFLAGS = -Wall $(DEBUG)
OFLAGS = -lm -pthread

all: sim_run pool 

sim_run: sim_run.o
	$(CC) $(LFLAGS)$  sim_run.o -o sim_run

pool: pool.o random.o
	${CC} $(LFLAGS) $(OFLAGS) pool.o random.o -o pool

$(INCLUDEDIR)pool.h: $(INCLUDEDIR)defines.h

sim_run.o: sim_run.cpp $(INCLUDEDIR)pool.h
	$(CC) $(CFLAGS)$  sim_run.cpp -o sim_run.o
	
pool.o: pool.cpp $(INCLUDEDIR)pool.h
	$(CC) $(CFLAGS) pool.cpp -o pool.o

random.o: random.cpp
	$(CC) ${CFLAGS} random.cpp -o random.o

clean:
	\rm *.o pool sim_run

# DO NOT DELETE

pool.o: ../include/defines.h 
random.o: ../include/defines.h
sim_run.o: ../include/defines.h
