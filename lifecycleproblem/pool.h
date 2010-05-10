#ifndef __POOL_
#define __POOL_

#include<ctime> //HLEE
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<vector>
#include<string>

using namespace std;

#define BITS_PER_BYTE 8
#define UINTSIZE (BITS_PER_BYTE*sizeof(unsigned))
#define INFINITY1 1e7
#define EPSILON  1e-6
#define PI 3.1415927
#define MAXVECSIZE 30
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
#define bin1
#define BARRIER_SYNC

#define VECTOR_SIZE     150

//Related to Chromosome Length
//#define NUM_LIFECYCLE 3
#define BITS_PER_LIFECYCLE 24 // 24 bits + 1 floating point number
#define MAX_LIFECYCLE 8
#define LC_ENCODE_LENGTH 4
#define NUM_CUSTOMER_GROUP 1
#define NUM_VAR_PER_LIFECYCLE 12

#define CONSTANT_POOL
#define POOL_SIZE 640
/*HLEE*/
#define MAX_GEN  1000


/*=================
TYPE DEFINTIONS :
=================*/
struct indiv
            {  float x[MAXVECSIZE];     /* variables            */
               float obj;               /* objective fn. value  */
               float mod_obj;           /* modified objective   */
               unsigned *chrom;         /* chromosome		  */
               float* lastbit;           /* cont var between 0 1 */
	       int parent1;
               int parent2;             /* s.no. of parents     */
               float cross_var;         /* cross over variable  */
            };
typedef struct indiv INDIVIDUAL ;
typedef INDIVIDUAL *POPULATION ;        /* array of individuals */

/*====================
FUNCTION PROTOTYPES :
====================*/
//float objectivenew(unsigned* px2, float* x3);
float   degraded_value();
float   randomperc1();
float   get_beta();
float   get_delta();
double  noise();
float   rndreal();
void    input_parameters();
void    ignore_comment();
void    initreport(FILE*);
void    nomemory(const char*);
void    error_ptr_null(const char*);
void    initpool();
void copy_individual(INDIVIDUAL* source , INDIVIDUAL* target, int);
int flip1(float);

void writechrom(unsigned *chrom);
void writechrom1(unsigned *chrom,float* lb);



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
	float objective(float *x);
	float objectivenew(unsigned* px2, float* x3);
	void generate_new_pop();
	void preselect_tour();
	void reset1();
	int tour_select();
	void cross_over_unif(int first,int second,int childno1,int childno2);
	void binary_xover (unsigned *parent1,float* plb1 ,unsigned *parent2,float* plb2, unsigned *child1,float* clb1, unsigned *child2, float* clb2);
void create_children(float p1,float p2,float *c1,float *c2,float low,float high,float *rand_var);
	float get_beta(float u);
	void mutation(INDIVIDUAL  *indiv);
	float get_delta(float u);
//void binmutation(unsigned *child);
	void binmutation(INDIVIDUAL *indiv);
	void free_all();


//Random Number Generation
	void   initrandomnormaldeviate();
	double noise(double,double);
	double randomnormaldeviate();
	void   advance_random();  
	int    flip(float);
	void   randomize();
	float  randomperc();
	int    rnd(int,int);
	float  rndreal(float, float);
	void   warmup_random(float);



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
float   seed,                  /* Random seed number                   */
        basic_seed;             /* Basic seed number                    */
	private:
float  	n_distribution_c, n_distribution_m,
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
float   best_in_each_gen[MAX_GEN];


	public:
	pthread_t thread;
	int id;
	private: 
//For random number generators
	double oldrand[55];   /* Array of 55 random numbers */
	int jrand;                 /* current random number */
	double rndx1,rndx2;    /* used with random normal deviate */
	int rndcalcflag; /* used with random normal deviate */


};

#endif
