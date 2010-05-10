all: sim_run pool 
sim_run: sim_run.c
	gcc sim_run.c -o sim_run
pool: pool.o random.o objective.o
	g++ -lm -pthread pool.o random.o objective.o -o pool
pool.o: pool.cpp pool.h
	g++ -c pool.cpp -o pool.o
random.o: random.cpp
	g++ -c random.cpp -o random.o
objective.o:objective.cpp
	g++ -c objective.cpp -o objective.o
clean:
	\rm *.o pool sim_run
