/*
 * For the Benchmark Problems
 * f1 minimization, minima = 0 at all xi = 0;
 * f2 minimization, minima = 0 at x1 = x2 = 1
 * f3 minimization, minima = 0 at all xi = 0;
 * f4 minimization, minima = -4189.83 at all xi = 420.9687
 * f4 has second minimum at xj = -302.5253, and all other xi = 420.9657
 */


#ifndef __POOL_
#define __POOL_

#include<ctime> //HLEE
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include<limits>
#include<ctime>
using namespace std;

#define BITS_PER_BYTE 8
#define UINTSIZE (BITS_PER_BYTE*sizeof(unsigned))
#define INFINITY1 1e7
#define EPSILON  1e-6
#define PI 3.1415927
#define MAXVECSIZE 100
#define MAXPOPSIZE 1000
#define TRUE 1
#define FALSE 0
#define BLX 0
#define SBX 1
#define ONESITE 1
#define UNIF 2
#define ONLINE 3
#define square(x)  ((x)*(x))

/***** Current Objective Function ******/
//#define prob1  /* define your problem at the end in objfunc() */

//#define sync
//#define bin1
//#define BARRIER_SYNC
#define f2
#define VECTOR_SIZE     150
#define MINMAX 1 //1 means minimization, -1 means maximization, needed by sim_run

//Related to Chromosome Length
//#define NUM_LIFECYCLE 3
#define BITS_PER_LIFECYCLE 24 // 24 bits + 1 floating point number
#define MAX_LIFECYCLE 8
#define LC_ENCODE_LENGTH 4


// Regarding Failures
//#define CRASHFAILURE
//#define BYZANTINE_RANDOM	//Inserts random individuals into the pool
//#define BYZANTINE_ANTIELITISM	//Writes individual back if it is worse than previous
#define BYZANTINE_FRACTION (0.0)
/*HLEE*/
#define MAX_GEN  1000

//For getting snapshot of the pool
//#define POOL_SNAPSHOT

/*=================
TYPE DEFINTIONS :
=================*/
struct indiv
{  double x[MAXVECSIZE];     /* variables            */
	double obj;               /* objective fn. value  */
	double mod_obj;           /* modified objective   */
	unsigned *chrom;         /* chrosome string 24 bits   */
	double* lastbit;           /* cont var between 0 1 */
	int parent1;
	int parent2;             /* s.no. of parents     */
	double cross_var;         /* cross over variable  */
};
typedef struct indiv INDIVIDUAL ;
typedef INDIVIDUAL *POPULATION ;        /* array of individuals */

/*====================
  FUNCTION PROTOTYPES :
  ====================*/
//double objectivenew(unsigned* px2, double* x3);
double   degraded_value();
double   randomperc1();
double   get_beta();
double   get_delta();
double  noise();
double   rndreal();
void    input_parameters();
void    ignore_comment();
void    initreport(FILE*);
void    nomemory(const char*);
void    error_ptr_null(const char*);
void    initpool();
void copy_individual(INDIVIDUAL* source , INDIVIDUAL* target, int);
int flip1(double);

void writechrom(unsigned *chrom);
void writechrom1(unsigned *chrom,double* lb);

double objectivecopy(double *x);

/*
 * As the name suggests represents one thread of the pool GA
 */
struct GAThread
{
	GAThread();
	void copy_global();
	static void* ThreadFunction(void*);

	//GA Algorithm
	void allocnewold();
	void statistics(int gen);
	void decode_string(INDIVIDUAL*);
	void decodevalue(unsigned *chrom,double value[]);
	double objective(double *x);
	double objectivenew(unsigned* px2, double* x3);
	void generate_new_pop();
	void preselect_tour();
	void reset1();
	int tour_select();
	void cross_over_unif(int first,int second,int childno1,int childno2);
	void binary_xover (unsigned *parent1,double* plb1 ,unsigned *parent2,double* plb2, unsigned *child1,double* clb1, unsigned *child2, double* clb2);
	void create_children(double p1,double p2,double *c1,double *c2,double low,double high,double *rand_var);
	double get_beta(double u);
	void mutation(INDIVIDUAL  *indiv);
	double get_delta(double u);
	//void binmutation(unsigned *child);
	void binmutation(INDIVIDUAL *indiv);
	void free_all();


	//Random Number Generation
	void   initrandomnormaldeviate();
	double noise(double,double);
	double randomnormaldeviate();
	void   advance_random();  
	int    flip(double);
	void   randomize();
	double  randomperc();
	int    rnd(int,int);
	double  rndreal(double, double);
	void   warmup_random(double);



	private:
	int     pop_size,               /* Population Size                      */
		gen_no,                 /* Current generation number            */
		max_gen,                /* Maximum no. of generations           */
		no_xover,               /* No. of cross overs done              */
		no_mutation,            /* No. of mutations done                */
		best_ever_gen,          /* Generation no. of best ever indiv.   */
		num_var,                /* Number of total design variables     */
		num_discr_var,          /* Number of discrete variables         */
		lchrom,                 /* Length of chromosome                 */
		chromsize,              /* Number of bytes needed to store
					   lchrom strings          */
		cross_type,             /* Cross over type ( SBX / BLX )        */
		x_strategy,s_strategy,  /* Cross-over strategy UNIF,ONLINE etc. */
		maxrun,                 /* Maxm no. of GA runs for each set of
					   parameter values          */
		run,                    /* Actual run no.                       */
		SHARING,                /* Flag for Sharing ( True / False)     */
		REPORT,                 /* Flag for Full reports (True/False)   */
		RIGID,                  /* Flag for rigid boundaries (T/F)      */
		BINGA,                  /* Flag for binary GA (T/F)             */
		REALGA,                 /* Flag for real-GA (T/F)               */
		READFILE,               /* Flag for reading input from file     */
		tourneylist[MAXPOPSIZE],/* List of indices of individuals for
					   tournament selection routine    */
		tourneypos,             /* Current position of tournament       */
		tourneysize,            /* Tournament size ( = 2 for binary )   */
		MINM;
	public:
	double   seed,                  /* Random seed number                   */
		basic_seed;             /* Basic seed number                    */
	private:
	double  	n_distribution_c, n_distribution_m,
		p_xover,                /* Cross over probability               */
		p_mutation,             /* Mutation probability                 */
		sum_obj,                /* Sum of objective fn. values          */
		avg_obj,                /* Average of objective fn. values      */
		max_obj,                /* Maximum objective fn. value          */
		min_obj,                /* Minimum objective fn. value          */
		minx[MAXVECSIZE],       /* Minimum and maximum values of design */
		maxx[MAXVECSIZE],       /*        variables in a population     */
		x_lower[MAXVECSIZE],    /* Lower and Upper bounds on each       */
		x_upper[MAXVECSIZE],    /*        design variable               */
		sigma_share;            /* Sharing distance                     */

	POPULATION oldpop, newpop;      /* Old and New populations      */
	INDIVIDUAL best_ever;   /* Best fit individual till current gen.*/
	/*HLEE: to record the best fit individual in each generation */
	double   best_in_each_gen[MAX_GEN];

	double 	failure_prob;

	public:
	pthread_t thread;
	int id;
	int byzthreshold;
	private: 
	//For random number generators
	double oldrand[55];   /* Array of 55 random numbers */
	int jrand;                 /* current random number */
	double rndx1,rndx2;    /* used with random normal deviate */
	int rndcalcflag; /* used with random normal deviate */


};

#endif
