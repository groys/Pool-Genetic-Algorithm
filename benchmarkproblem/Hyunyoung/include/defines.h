#ifndef __DEFINES__
#define __DEFINES__

/***** General ****/
#define BITS_PER_BYTE 8
#define UINTSIZE (BITS_PER_BYTE*sizeof(unsigned))
#define INFINITY1 1e7
#define EPSILON  1e-6
#define PI 3.1415927
#define MAXVECSIZE 100	/* The maximum number of real variables in the objective function */
#define MAXPOPSIZE 1000
#define TRUE 1
#define FALSE 0
#define BLX 0
#define SBX 1
#define ONESITE 1
#define UNIF 2
#define ONLINE 3
#define square(x)  ((x)*(x))

/***** Objective Function ******/
//#define sync 	//(sin x)/ x
//#define bin1
#define f4
#define MINMAX 1 //1 means minimization, -1 means maximization, needed by sim_run

/*Related to Chromosome Length*/
//#define NUM_LIFECYCLE 3
#define BITS_PER_LIFECYCLE 24 // 24 bits + 1 floating point number
#define MAX_LIFECYCLE 8
#define LC_ENCODE_LENGTH 4


/* Regarding Failures */
//#define CRASHFAILURE
//#define BYZANTINE_RANDOM	//Inserts random individuals into the pool
#define BYZANTINE_ANTIELITISM	//Writes individual back if it is worse than previous
#define BYZANTINE_FRACTION (0.40)//This always needs to be defined, sim_run needs it
/*HLEE*/
#define MAX_GEN  1000

/*For getting snapshot of the pool*/
//#define POOL_SNAPSHOT

/* Synchronous Operation */
//#define BARRIER_SYNC



/********* For Sim Run ****************/

#define VARY_NUM_THREADS  6 /* 1, 2, 4, 8, 16, 32 */
#define MAX_GEN  1000
#define RUN 10
#define FILENAME_LEN  80
#define LINE_LEN  80
#define KEY_STR  "BestEverFitness="
#define KEY_LEN   17 //21
#define KEY_STR1 "Best in each gen = "
#define KEY_LEN1  20
#define VAL_LEN  10
#define GEN_LEN  3
#define EOS  '\0'

#define INTERVALS 10
#define INPUTDIR "../inputs/"
#define OUTPUTDIR "../outputs/"


#endif
