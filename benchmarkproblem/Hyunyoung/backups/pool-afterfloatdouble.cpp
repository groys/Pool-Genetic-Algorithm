/* Instructions for Running 
 * 1. Always set gMINM in the main function = 1 for minimization and = -1 for maximization
 */


#include "pool.h"
#include <pthread.h>
#include <vector>
using namespace std;

/*==================
GLOBAL VARIABLES  :
==================*/

int     gpop_size,               /* Population Size                      */
        gmax_gen,                /* Maximum no. of generations           */
        gbest_ever_gen,          /* Generation no. of best ever indiv.   */
        gnum_var,                /* Number of total design variables     */
        gnum_discr_var,          /* Number of discrete variables     
    				    Means how many variables are embedded
			     	    in a binary GA String                 */
        glchrom,                 /* Length of chromosome in bits          */
        gchromsize,              /* Number of unsigned ints needed to store
                                    strings of length glchrom            */
        gcross_type,             /* Cross over type ( SBX / BLX )        */
        gx_strategy,gs_strategy,  /* Cross-over strategy UNIF,ONLINE etc. */
        gmaxrun,                 /* Maxm no. of GA runs for each set of
                                              parameter values          */
        gSHARING,                /* Flag for Sharing ( True / False)     */
        gREPORT,                 /* Flag for Full reports (True/False)   */
        gRIGID,                  /* Flag for rigid boundaries (T/F)      */
        gBINGA,                  /* Flag for binary GA (T/F)             */
        gREALGA,                 /* Flag for real-GA (T/F)               */
        gREADFILE,               /* Flag for reading input from file     */
        gtourneysize,            /* Tournament size ( = 2 for binary )   */
        gMINM,
	gnum_threads,		 /* The number of threads to spawn	 */
	gnum_failed;		 /* The number of threads that failed	 */

double   gseed,                   /* Random seed number                   */
        gbasic_seed,             /* Basic seed number                    */
  	gn_distribution_c, gn_distribution_m,
        gp_xover,                /* Cross over probability               */
        gp_mutation,             /* Mutation probability                 */
        gsum_obj,                /* Sum of objective fn. values          */
        gavg_obj,                /* Average of objective fn. values      */
        gmax_obj,                /* Maximum objective fn. value          */
        gmin_obj,                /* Minimum objective fn. value          */
        gminx[MAXVECSIZE],       /* Minimum and maximum values of design */
        gmaxx[MAXVECSIZE],       /*        variables in a population     */
        gx_lower[MAXVECSIZE],    /* Lower and Upper bounds on each       */
        gx_upper[MAXVECSIZE],    /*        design variable               */
        gsigma_share,            /* Sharing distance                     */
	gfailure_prob;		 /* probability that a thread fails	 */

POPULATION goldpop, gnewpop;     /* Old and New populations              */
POPULATION pool;		/* The pool of individuals		 */
INDIVIDUAL gbest_ever;           /* Best fit individual till current gen.*/

int pool_size;

/* Array of Mutexes */
pthread_mutex_t *vmutex;
pthread_mutex_t malloc_mutex;
pthread_mutex_t	failure_mutex;

/*====================================================================
SUBROUTINE FOR INPUTTING GLOBAL PARAMETERS :
====================================================================*/
void input_parameters()
{
	int k;
	char ans;

	/*printf("       ");
	  puts("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::");
	  printf("       ");
	  puts("::::::::       REAL-CODED GENETIC ALGORITHM        :::::::");
	  printf("       ");
	  puts("::::::::       ============================        :::::::");
	  printf("       ");
	  puts("::::::::     (c) R.B.Agrawal and K.Deb, 1995;      :::::::");
	  printf("       ");
	  puts("::::::::            All rights reserved.           :::::::");
	  printf("       ");
	  puts("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::");

	  printf("\n ARE YOU READING IT THROUGH A COMMENTED FILE (y/n) ?");
	 */   
	do { ans = getchar(); } while (ans!= 'y' && ans !='n');
	if (ans == 'y')      gREADFILE = TRUE;
	else                 gREADFILE = FALSE;

	if (gREADFILE) printf("\n Reading data from file ............");
	if (!gREADFILE) printf("\nHow many generations ? ------------- : ");
	ignore_comment();
	scanf("%d",&gmax_gen);
	if (!gREADFILE)  printf("\nPopulation Size ? ------------------ : ");
	ignore_comment();
	scanf("%d", &gpop_size );
	if (gpop_size > MAXPOPSIZE)
	{
		printf("\n Increase the value of MAXPOPSIZE in program");
		printf("  and re-run the program");
		exit(-1);
	}
	if (!gREADFILE) printf("\nCross Over Probability ? ( 0 to 1 )  : ");
	ignore_comment();
	scanf("%lf",&gp_xover);
	if (!gREADFILE) printf("\nMutation Probability ? ( 0 to 1 ) -- : ");
	ignore_comment();
	scanf("%lf",&gp_mutation);
	if (!gREADFILE)
		printf("\nNumber of variables (Maximum %d) ---- : ",MAXVECSIZE);
	ignore_comment();
	scanf("%d",&gnum_var);
	gBINGA = gREALGA = FALSE;
	if (!gREADFILE) printf("\n Binary or Real-coded parameters? (b for binary, r for real-coded) ");
	ignore_comment();
	do { ans = getchar(); } while (ans!= 'b' && ans !='r');
	if (ans == 'b') gBINGA = TRUE;
	else            gREALGA = TRUE;
	if (gBINGA) gnum_discr_var = gnum_var;
	if(gREALGA) gnum_discr_var = 0;
	ignore_comment(); 
	for (k=0; k<= gnum_var-1; k++)
	{
		if (!gREADFILE) printf("\nLower and Upper bounds of x[%d] ----- : ",k+1);
		scanf("%lf %lf",&gx_lower[k],&gx_upper[k]);
	}
	if (gREALGA) {
		if (!gREADFILE) printf("\n Are these bounds rigid ? (y/n)");
		ignore_comment();
		do { ans = getchar(); } while (ans!= 'y' && ans !='n');
		if (ans == 'y')      gRIGID = TRUE;
		else  gRIGID = FALSE;
	} else {
		if (gREADFILE) {
			ignore_comment();
			ans = getchar();
			gRIGID = TRUE;
		}
	}
	if (gBINGA)
	{
		if (!gREADFILE) printf("\n Total string length (each variable has equal string length)?");
		ignore_comment();
		scanf("%d",&glchrom);
	}
	else 
		if (gREADFILE) {
			ignore_comment();
			scanf("%d",&glchrom);
			glchrom=0;
		}
	if (!gREADFILE) printf("\nSharing to be done ? (y/n) --------- :");
	ignore_comment();
	do { ans = getchar(); } while (ans!= 'y' && ans !='n');
	if (ans == 'y')
	{ gSHARING = TRUE;
		if (!gREADFILE) printf("\nSigma share value  ? --------------- :");
		scanf("%lf",&gsigma_share);
	}
	else  gSHARING = FALSE;
	if (!gREADFILE) printf ("\n Reports to be printed ? (y/n) ");
	ignore_comment();
	do { ans = getchar(); } while (ans!= 'y' && ans !='n');
	if (ans == 'y') gREPORT = TRUE;
	else            gREPORT = FALSE;
	if (!gREADFILE) printf("\n How many runs ?");
	ignore_comment();
	scanf("%d",&gmaxrun);
	if (!gREADFILE) 
	{
		printf("\n Enter selection operator --> ");
		printf("\n  1 : Tournament selection (min or max set by MINM in the code)");
		printf("\n  2 : Roulette wheel selection (always max)");
		printf("\n  3 : Stochastic remainder roulette wheel selection (always max)");
		printf("\n  Give your choice :");
	}
	ignore_comment();
	scanf("%d",&gs_strategy);
	if (gs_strategy == 1) {
		if (!gREADFILE) printf("\n Enter tournament size ");
		ignore_comment();
		scanf("%d", &gtourneysize);
	} else 
		if (gREADFILE) {
			ignore_comment();
			scanf("%d",&gtourneysize);
			gtourneysize=0;
		}
	if (gSHARING) gs_strategy = 3; /* Stoch. Rem. RW is default */
	printf("TTTTT  %d \n",gtourneysize);
	if (gREALGA)
	{
		if (!gREADFILE)
		{
			printf("\n  Give the strategy for X-over");
			printf("\n  1 : Polynomial distribution in one variable");
			printf("\n  2 : Polynomial distribution in all variables");
			printf("\n  3 : Polynomial distribution on a straight line");
			printf("\n  Give your choice :");
		}
		ignore_comment();
		scanf("%d",&gx_strategy);

		if (!gREADFILE) printf("\n Type of cross over ? ( s for SBX, b for BLX) ");
		ignore_comment();
		do { ans = getchar(); } while (ans!= 's' && ans !='b');
		if (ans == 's') gcross_type = SBX;
		else            gcross_type = BLX;
		if (gcross_type == SBX)
		{
			if (!gREADFILE) printf("\n Give eta for SBX and mutation?");
			ignore_comment();
			scanf("%lf %lf",&gn_distribution_c,&gn_distribution_m);
		}
	}
	else {
		ignore_comment();
		scanf("%d",&gx_strategy);
		gx_strategy = 1; /* binary crossover */ //@gautam whats this hack!
		ignore_comment();
		ans = getchar();
		ignore_comment();
		scanf("%lf %lf",&gn_distribution_c,&gn_distribution_m);
		gn_distribution_c = 0; gn_distribution_m = 0;
	}
	if (!gREADFILE) printf("\n Give random seed (0 to 1.0)");
	ignore_comment();
	scanf("%lf",&gbasic_seed);
	//HLEE
	gbasic_seed = (double) (time(0) % 10)/ (double) (time(0) % 10 + 1);
	printf("\ngbasic_seed =  %lf\n", gbasic_seed);


	//Calculate the length of the chromosome to be used
	//the + 1 is now added in the LC_ENCODE_LENGTH
	glchrom = LC_ENCODE_LENGTH + BITS_PER_LIFECYCLE*MAX_LIFECYCLE; //+ 1; 
	//   input_app_parameters();
}
/*====================================================================
  Ignores the comment from input ( ended by a ':')
  ====================================================================*/
void ignore_comment()
{
	if (gREADFILE == FALSE) return;

	do
	{}
	while (getchar() != ':');
}


GAThread::GAThread()
{
}


void GAThread::copy_global()
{
        pop_size = gpop_size;
	gen_no = 0;              		/* Current generation number            */
        max_gen = gmax_gen;             	/* Maximum no. of generations           */
        no_xover = 0;     			/* No. of cross overs done              */
        no_mutation = 0;           		/* No. of mutations done                */
        best_ever_gen = 0;         		/* Generation no. of best ever indiv.   */
        num_var = gnum_var;                	/* Number of total design variables     */
        num_discr_var = gnum_discr_var;         /* Number of discrete variables         */
        lchrom = glchrom;                 	/* Length of chromosome                 */
        chromsize = gchromsize;              	/* Number of bytes needed to store
                                                		lchrom strings          */
        cross_type = gcross_type;             	/* Cross over type ( SBX / BLX )        */
        x_strategy = gx_strategy;
	s_strategy = gs_strategy; 		/* Cross-over strategy UNIF,ONLINE etc. */
        maxrun = gmaxrun;                 	/* Maxm no. of GA runs for each set of
                                              			parameter values        */
        run = 0;       				/* Actual run no.                       */
        SHARING = gSHARING;                	/* Flag for Sharing ( True / False)     */
        REPORT = gREPORT;     			/* Flag for Full reports (True/False)   */
        RIGID = gRIGID;                  	/* Flag for rigid boundaries (T/F)      */
        BINGA = gBINGA;                  	/* Flag for binary GA (T/F)             */
        REALGA = gREALGA;                 	/* Flag for real-GA (T/F)               */
        READFILE = gREADFILE;               	/* Flag for reading input from file     */
        tourneysize = gtourneysize;            	/* Tournament size ( = 2 for binary )   */
        MINM = gMINM;
	seed = gseed;                   	/* Random seed number                   */
        basic_seed = gbasic_seed;           	/* Basic seed number                    */
  	n_distribution_c = gn_distribution_c;
	n_distribution_m = gn_distribution_m;
        p_xover = gp_xover;                	/* Cross over probability               */
        p_mutation = gp_mutation;             	/* Mutation probability                 */
        sum_obj = 0;                		/* Sum of objective fn. values          */
        avg_obj = 0;            		/* Average of objective fn. values      */
        max_obj = 0;                		/* Maximum objective fn. value          */
        min_obj = 0;            		/* Minimum objective fn. value          */
        sigma_share = gsigma_share;            	/* Sharing distance                     */
	for (int k=0; k<= num_var-1; k++)
     	{
		x_upper[k] = gx_upper[k];
		x_lower[k] = gx_lower[k];
//cout << x_lower[k] << " " << x_upper[k] << endl;
	}
	failure_prob = gfailure_prob;

	byzthreshold = gnum_threads*BYZANTINE_FRACTION;
//	cout << "Byzantine Threshold = " << byzthreshold << '\n';
}



/*====================================================================
Allocates memory for oldpop and newpop
====================================================================*/
void GAThread::allocnewold()
{
   double u;
   int k,k1,i,j,j1,stop;
   double temp[MAXVECSIZE],coef;	//stack
   unsigned mask=1,nbytes;

   randomize();					//in random.cpp note that seed must be set

   //lock as malloc is not thread safe
//   pthread_mutex_lock(&malloc_mutex);
   oldpop = (INDIVIDUAL *)malloc(pop_size*sizeof(INDIVIDUAL));
   if (oldpop == NULL) nomemory("oldpop in allocnewold()");
   newpop = (INDIVIDUAL *)malloc(pop_size*sizeof(INDIVIDUAL));
   if (newpop == NULL) nomemory("newpop in allocnewold()");

   nbytes = chromsize*sizeof(unsigned);
   for(j = 0; j < pop_size; j++)
   {
      if((oldpop[j].chrom = (unsigned *) malloc(nbytes)) == NULL)
        nomemory("oldpop chromosomes");

      if((newpop[j].chrom = (unsigned *) malloc(nbytes)) == NULL)
        nomemory("newpop chromosomes");
      
      if((oldpop[j].lastbit = (double *) malloc(MAX_LIFECYCLE*sizeof(double))) == NULL)
        nomemory("oldpop lastbit");

      if((newpop[j].lastbit = (double *) malloc(MAX_LIFECYCLE*sizeof(double))) == NULL)
        nomemory("newpop lastbit");
   }
   if((best_ever.chrom = (unsigned *) malloc(nbytes)) == NULL)
      nomemory("best_ever chromosomes");
//   pthread_mutex_unlock(&malloc_mutex);
   if((best_ever.lastbit = (double *) malloc(MAX_LIFECYCLE*sizeof(double))) == NULL)
      nomemory("best_ever lastbit");
   
   no_xover = no_mutation = 0;
}

/*====================================================================
Releases the memory for all mallocs
====================================================================*/
void GAThread::free_all()
{
  int i;

  for(i = 0; i < pop_size; i++)
  {
      free(oldpop[i].chrom);
      free(newpop[i].chrom);
      free(oldpop[i].lastbit);
      free(newpop[i].lastbit);
  }
  free(oldpop);
  free(newpop);
  free(best_ever.chrom);
  free(best_ever.lastbit);
  //app_free();
}




/*====================================================================
Calculates statistics of current generation :
====================================================================*/
void GAThread::statistics(int gen)
{
   INDIVIDUAL current_best;
   int k,j;
   double f;
   double bitpow,coef;

   //Calculate all the objectives
   for (k=0; k<=pop_size-1; k++)
   {
     newpop[k].obj = objective(newpop[k].x);
     newpop[k].mod_obj = newpop[k].obj;
   }
   if ( gen == 0)
   {
	     best_ever.mod_obj = best_ever.obj;
       /*HLEE*/
	     //@gautam
	     //commented out to eliminate negative values that were output in some generations. Now minimum value output is 0 based on initialization of best_in_each_gen
	     //
	//		 if (best_ever.obj > best_in_each_gen[gen])
	//		     best_in_each_gen[gen] = best_ever.obj;
   }

   current_best = newpop[0];
   sum_obj = avg_obj = newpop[0].obj;
   max_obj = min_obj = newpop[0].obj;
   for (k=0;  k<= num_var-1 ; k++) 
	   maxx[k] = minx[k] = newpop[0].x[k];

   //Find the current best in this generation
   for(k=1; k<= pop_size-1; k++)
   {
     if(MINM * current_best.obj  >  MINM * newpop[k].obj)
                   current_best = newpop[k];		//note this is bitwise copy not a copy constructor
     if(MINM * max_obj < MINM * newpop[k].obj)
                   max_obj = newpop[k].obj;
     if(MINM * min_obj > MINM * newpop[k].obj)
                   min_obj = newpop[k].obj;
     sum_obj += newpop[k].obj;
     for (j=0; j<= num_var-1; j++)
     {
        if (MINM * maxx[j] < MINM * newpop[k].x[j]) 
		maxx[j] = newpop[k].x[j];
        if (MINM * minx[j] > MINM * newpop[k].x[j]) 
		minx[j] = newpop[k].x[j];
     }
   }
   avg_obj = sum_obj/pop_size;
   if (MINM * current_best.obj < MINM * best_ever.obj)
   { 
	   copy_individual(&current_best, &best_ever,chromsize);
     	   best_ever_gen = gen;
   }
   //cout << gen << "statistics function done " << '\n';
  /*HLEE :keep the current best in each gen. Other statistics done in sim_run*/
 //@gautam best_in_each_gen[gen] = better(current_best.obj, best_in_each_gen[gen-1] 
 //this gives it the monotonically better property
 //18March 2009 this seems same as best_in_each_gen[gen] = best_ever.obj
   if (current_best.obj < MINM*best_in_each_gen[gen-1])
       best_in_each_gen[gen] = current_best.obj;
   else /* keep the previous best */
       best_in_each_gen[gen] = best_in_each_gen[gen-1];
}
/*====================================================================
GENERATION OF NEW POPULATION through SELECTION, XOVER & MUTATION :
====================================================================*/
void GAThread::generate_new_pop()
{
   int k,mate1,mate2;

   preselect_tour();
   for (k=0; k<= pop_size-1; k +=2)
   {
	mate1 = tour_select();
	mate2 = tour_select();
	
     cross_over_unif(mate1,mate2,k,k+1);
     mutation(&newpop[k]);
     mutation(&newpop[k+1]);
     newpop[k].parent1 = newpop[k+1].parent1 = mate1+1;
     newpop[k].parent2 = newpop[k+1].parent2 = mate2+1;
   }
}


void GAThread::preselect_tour()
{
    reset1();
    tourneypos = 0;
}

void GAThread::reset1()
/* Name changed from reset because of clash with lib. function - RBA */
/* Shuffles the tourneylist at random */
{
    int i, rand1, rand2, temp_site;

    for(i=0; i<pop_size; i++) tourneylist[i] = i;

    for(i=0; i < pop_size; i++)
    {
        rand1= rnd(0,pop_size-1);
        rand2=  rnd(0,pop_size-1);
        temp_site = tourneylist[rand1];
        tourneylist[rand1]=tourneylist[rand2];
        tourneylist[rand2]=temp_site;
    }
}

int GAThread::tour_select()
{
    int pick, winner, i;

    /* If remaining members not enough for a tournament, then reset list */
start_select :
    if((pop_size - tourneypos) < tourneysize)
    {
        reset1();
        tourneypos = 0;
    }

    /* Select tourneysize structures at random and conduct a tournament */
    winner=tourneylist[tourneypos];
/* Added by RBA */
    if( winner < 0 || winner > pop_size-1) 
    {
                                             printf("\n Warning !! ERROR1");
                                             printf(" tourpos = %d",tourneypos);
                                             printf(" winner = %d",winner);
                                             preselect_tour();
                                             goto start_select;
    }
    for(i=1; i<tourneysize; i++)
    {
        pick=tourneylist[i+tourneypos];
/* Added by RBA */
        if (pick < 0 || pick > pop_size-1) 
	{ 
					     preselect_tour();
                                             printf("\n Warning !! ERROR2");
                                             goto start_select;
       	}
        if(MINM * oldpop[pick].mod_obj < MINM * oldpop[winner].mod_obj) 
		winner=pick;
    }

    /* Update tourneypos */
    tourneypos += tourneysize;
    return(winner);
}


/*====================================================================
CROSS - OVER  USING strategy of uniform 50% variables
  For one variable problem, it is crossed over as usual.
  For multivariables, each variable is crossed over with a probability
  of 50 % , each time generating a new random beta.
====================================================================*/
void GAThread::cross_over_unif(int first1,int second1,int childno1,int childno2)
{
    double difference,x_mean,beta;
    double u = 0.0;
    int site,k;

    if (flip(p_xover))   /* Cross over has to be done */
    {
     no_xover++;
     if (BINGA)
     binary_xover(oldpop[first1].chrom,oldpop[first1].lastbit, oldpop[second1].chrom,oldpop[first1].lastbit, newpop[childno1].chrom,newpop[childno1].lastbit, newpop[childno2].chrom, newpop[childno2].lastbit);
     if (REALGA)
     {
      for (site = num_discr_var; site<=num_var-1; site++)
      {
       if(flip(0.5) || (num_var==1))
       {
         create_children(oldpop[first1].x[site],oldpop[second1].x[site],
                      &(newpop[childno1].x[site]),&(newpop[childno2].x[site]),
                      x_lower[site],x_upper[site],&u);
       }
       else
       {
        newpop[childno1].x[site] = oldpop[first1].x[site];
        newpop[childno2].x[site] = oldpop[second1].x[site];
       }
      }               /* for loop */
     newpop[childno1].cross_var = newpop[childno2].cross_var = u;
     }                /* if REALGA      */
    }                 /* Cross over done */

    else              /* Passing x-values straight */
    {
      if (BINGA)
      {
	      for (k=0; k< chromsize; k++)
	      {
	        newpop[childno1].chrom[k] = oldpop[first1].chrom[k];
	        newpop[childno2].chrom[k] = oldpop[second1].chrom[k];
	      }
	      for(k = 0; k < MAX_LIFECYCLE ; k++)
	      {
		newpop[childno1].lastbit[k]  = oldpop[first1].lastbit[k];
		newpop[childno2].lastbit[k]  = oldpop[second1].lastbit[k];
	      }
      }
      for (site=0; site<=num_var-1; site++)
      {
        newpop[childno1].x[site] = oldpop[first1].x[site];
        newpop[childno2].x[site] = oldpop[second1].x[site];
      }
    }
}

/*====================================================================
Binary cross over routine.
====================================================================*/
/* edited aug 1 
 * adding support for floating point number as last bit
 * note still considering on 24 bit positions must add for the last bit
 */
void GAThread::binary_xover (unsigned *parent1,double* plastbit1, unsigned *parent2, double* plastbit2, unsigned *child1,double* clastbit1, unsigned *child2,double* clastbit2)
/* Cross 2 parent strings, place in 2 child strings */
{
    int j, jcross, k,nfloat;
    unsigned int mask, temp;
    unsigned int nlc1, nlc2,nlc;
    
    mask = 1;
    
    if (BINGA == FALSE) return;
    if (parent1 == NULL) error_ptr_null("parent1 in binary_xover");
    if (parent2 == NULL) error_ptr_null("parent2 in binary_xover");
    if (child1== NULL) error_ptr_null("child1 in binary_xover");
    if (child2== NULL) error_ptr_null("child2 in binary_xover");

    nlc1 = (parent1[0] & mask) + (parent1[0] & (mask << 1)) + (parent1[0] & (mask << 2)) + 1;
    nlc2 = (parent2[0] & mask) + (parent2[0] & (mask << 1)) + (parent2[0] & (mask << 2)) + 1;
    nlc = nlc1 < nlc2 ? nlc1 : nlc2;

//Note there are 24 binary bits per lifecycle and the 25th bit is a floating point number
//jcross = rnd(1, (lchrom + MAX_LIFECYCLE));
 
    jcross = rnd(1, ((BITS_PER_LIFECYCLE + 1)*nlc + LC_ENCODE_LENGTH ));
    //find how many floats will be crossed over = MAX_LIFECYCLE - nfloat
    nfloat = (jcross-LC_ENCODE_LENGTH)/((int) (BITS_PER_LIFECYCLE + 1));
    
    //The parts of the chromosome that get exchanged
    for(k = nfloat; k < MAX_LIFECYCLE; k++)
    {
	clastbit2[k] = plastbit1[k];
	clastbit1[k] = plastbit2[k];
    }
    //parts of chromosome that go through unchanged
    for(k = 0; k < nfloat; k++)
    {
	clastbit2[k] = plastbit2[k];
	clastbit1[k] = plastbit1[k];
    }
    //Now do crossover on the bits
    jcross -= nfloat;
    for(k = 1; k <= chromsize; k++)
    {
            if(jcross >= (k*UINTSIZE))
            {
                child1[k-1] = parent1[k-1];
                child2[k-1] = parent2[k-1];
            }
            else if((jcross < (k*UINTSIZE)) && (jcross > ((k-1)*UINTSIZE)))
            {
                mask = 1;
                for(j = 1; j <= (jcross-1-((k-1)*UINTSIZE)); j++)
                {
                    temp = 1;
                    mask = mask<<1;
                    mask = mask|temp;
                }
                child1[k-1] = (parent1[k-1]&mask)|(parent2[k-1]&(~mask));
                child2[k-1] = (parent1[k-1]&(~mask))|(parent2[k-1]&mask);
            }
            else
            {
                child1[k-1] = parent2[k-1];
                child2[k-1] = parent1[k-1];
            }
    }
    
}

/*====================================================================
Creates two children from parents p1 and p2, stores them in addresses
pointed by c1 and c2.  low and high are the limits for x values and
rand_var is the random variable used to create children points.
====================================================================*/
void GAThread::create_children(double p1,double p2,double *c1,double *c2,double low,double high,double *rand_var)
{
    double difference,x_mean,beta;
    double u,distance,umax,temp,alpha;
    int flag;

    if (c1 == NULL) error_ptr_null("c1 in create_children");
    if (c2 == NULL) error_ptr_null("c2 in create_children");
    if (rand_var == NULL) error_ptr_null("rand_var in create_children");
    flag = 0;
    if ( p1 > p2) { temp = p1; p1 = p2; p2 = temp; flag = 1; }
    x_mean = ( p1 + p2) * 0.5;
    difference = p2 - p1;
    if ( (p1-low) < (high-p2) ) distance = p1-low;
    else                        distance = high-p2;
    if (distance < 0.0) distance = 0.0;
    if (RIGID && (difference > EPSILON))
    {
      alpha = 1.0 + (2.0*distance/difference);
      umax = 1.0 - (0.5 / pow((double)alpha,(double)(n_distribution_c+1.0)));
      *rand_var = umax * randomperc();
    }
    else *rand_var = randomperc();
    beta = get_beta(*rand_var);
    if (fabs(difference*beta) > INFINITY1) beta = INFINITY1/difference;
    *c2 = x_mean + beta * 0.5 * difference;
    *c1 = x_mean - beta * 0.5 * difference;
    if (flag == 1) { temp = *c1; *c1 = *c2; *c2 = temp; }
}

/*===================================================================
Calculates beta value for given random number u (from 0 to 1)
If input random numbers (u) are uniformly distributed for a set of
inputs, this results in uniform distribution of beta values in case
of BLX , and Binary Probability distribution simulation in case of
SBX.
====================================================================*/
double GAThread::get_beta(double u)
{
   double beta;

   if (cross_type == BLX) return(2.0*u);
   if (1.0-u < EPSILON ) u = 1.0 - EPSILON;
   if ( u < 0.0) u = 0.0;
   if (u < 0.5) beta = pow(2.0*u,(1.0/(n_distribution_c+1.0)));
   else beta = pow( (0.5/(1.0-u)),(1.0/(n_distribution_c+1.0)));
   return beta;
}

/*===================================================================
Mutation Using polynomial probability distribution. Picks up a random
site and generates a random number u between -1 to 1, ( or between
minu to maxu in case of rigid boudaries) and calls the routine
get_delta() to calculate the actual shift of the value.
====================================================================*/
void GAThread::mutation(INDIVIDUAL  *indiv)
{
   double distance1,distance2,x,delta,minu,maxu,u;
   int k, site;

   if (indiv == NULL) error_ptr_null("indiv in mutation");
   if(flip (p_mutation) && REALGA)
   {
        site = rnd(0,num_var - 1);
        no_mutation++;
        if(fabs(x_upper[site] -x_lower[site]) < EPSILON) return;

        /* calculation of bounds on delta */
        if(RIGID)
        { x = indiv->x[site];
          distance1 = x - x_lower[site];
          distance2 = x_upper[site] - x;

          delta = 2.0 * distance1 / (x_upper[site] - x_lower[site]);
          if (delta > 1.0)   delta = 1.0;
          minu = -1.0 + pow((1.0 - delta),(n_distribution_m + 1.0));

          delta = 2.0 * distance2 / (x_upper[site] - x_lower[site]);
          if (delta > 1.0)   delta = 1.0;
          maxu = 1.0 - pow((1.0 - delta),(n_distribution_m + 1.0));
          u = rndreal(minu,maxu);
        }
        else u = rndreal(-1.0,1.0);

        /* calculation of actual delta value */
        delta = get_delta(u) *  0.5 * (x_upper[site] - x_lower[site]);
        indiv->x[site] += delta;
  }    /* if flip() */
  if (BINGA) 
	  binmutation(indiv);
}

/*==================================================================
For given u value such that   -1 <= u <= 1, this routine returns a
value of delta from -1 to 1. Exact value of delta depends on specified
n_distribution. This is called by mutation().
====================================================================*/
double GAThread::get_delta(double u)
{
   double delta;
   int negative = FALSE;   /* Flag for negativeness of delta */

   if (cross_type == BLX) return(u);
   if(u <= -1.0) u = -1.0;
   if(u >1.0)  u =  1.0;
   if(u < 0.0)  {  u = -u;
                   negative = TRUE;
                }
   delta = 1.0 - pow((1.0 - u),(1.0 / (n_distribution_m + 1.0)));
   if(negative)  return (-delta);
   else          return delta;

}
/*==================================================================
Binary mutation routine ( borrowed from sga.c )
====================================================================*/
/* edited aug 1 
 * adding support for floating point number as last bit*/
void GAThread::binmutation(INDIVIDUAL *ind)
/* Mutate an allele w/ pmutation, count # of mutations */
{
    int j, k, stop;
    unsigned mask, temp = 1,*child;
    child = ind->chrom;

    if (BINGA == FALSE) return;
    if (child== NULL) error_ptr_null(" child in binmutation");
    for(k = 0; k < chromsize; k++)
    {
        mask = 0;
        if(k == (chromsize-1))
            stop = lchrom - ((k-1)*UINTSIZE);
        else
            stop = UINTSIZE;
        for(j = 0; j < stop; j++)
        {
            if(flip(p_mutation))
            {
		mask = mask|(temp<<j);
        	no_mutation++;
            }
        }
        child[k] = child[k]^mask;
    }
    for(k = 0; k < MAX_LIFECYCLE ; k++)
    {
    	if(flip(p_mutation))
    	{
	    ind->lastbit[k] = rand()/(double) RAND_MAX;
	    no_mutation++;
    	}
    }
}



/*
 * Remember When copying in from pool copy to oldpop
 * When copying back to pool, copy from newpop
 */

void* GAThread::ThreadFunction(void* p)
{
	/*HLEE*/
	int g; /* for generation number when output */

	int i,j,ind;
	double u;
	int mask = 1;
	GAThread* pg = (GAThread*) p;
	//	cout << "In thread " << pg->id << endl;


	pg->gen_no = 0;
	pg->allocnewold();

	//Initialize the best_in_each_gen array to zeros
	for(i = 0; i < pg->max_gen; i++)
		pg->best_in_each_gen[i] = pg->MINM*1e10;	
	//Initialize the Best ever individual
	pthread_mutex_lock(&malloc_mutex);	
	copy_individual(&pool[0],&pg->best_ever,pg->chromsize);
	pg->best_ever.obj = pg->objective(pg->best_ever.x);
	pthread_mutex_unlock(&malloc_mutex);	
	pg->best_ever.mod_obj = pg->best_ever.obj;


	//Read into oldpop
	//Copy out from newpop
	for(pg->gen_no = 1; pg->gen_no<= pg->max_gen; pg->gen_no++)
	{
		//Stop Thread if it has to crash
#ifdef CRASHFAILURE
		//Should the thread fail?
		if(pg->flip(pg->failure_prob))
		{
			pthread_mutex_lock(&failure_mutex);
			gnum_failed++;
			cout << "Thread " << pg->id << " failed in generation " << pg->gen_no << '\n';
			pthread_mutex_unlock(&failure_mutex);
			return NULL;
		}
#endif

		//Generate the new set of individuals
		//
		for(i = 0; i < pg->pop_size; i++)
		{
			//Generate a random index and read from that index in the pool
			ind = pg->rnd(0,pool_size-1);
			pthread_mutex_lock(&vmutex[ind]);
			copy_individual(&pool[ind],&pg->oldpop[i],pg->chromsize);
			pthread_mutex_unlock(&vmutex[ind]);
		}
		pg->generate_new_pop();
		pg->statistics(pg->gen_no);

		//Writing Back to the Pool	
		for(i = 0; i < pg->pop_size; i++)
		{
			ind = pg->id*pg->pop_size + i;
			pthread_mutex_lock(&vmutex[ind]);
			//If no Byzantine Failute writeback with Elitism :)
			//The first byzthreshold proccessors are byzantine
			//for BAE writeback to pool only if worse than existing indiv
			//for BR writeback the randomly generated individual to pool

#ifdef BYZANTINE_ANTIELITISM

			//HLEE-Apr-18-2009: changed <= to < since this is more intuitive:
			//   consider, for example, 2 threads with 0.5 byz fraction, then one 
			//   would expect to have 1 thread faulty and the other thread non-faulty.
			if(pg->id < pg->byzthreshold )
			{	
				if(pg->MINM*pool[ind].obj <  pg->MINM*pg->newpop[i].obj)
					copy_individual(&pg->newpop[i],&pool[ind],pg->chromsize);
			}
			else
			{
				if(pg->MINM*pool[ind].obj >  pg->MINM*pg->newpop[i].obj)
					copy_individual(&pg->newpop[i],&pool[ind],pg->chromsize);
			}
#else
			if(pg->MINM*pool[ind].obj >  pg->MINM*pg->newpop[i].obj)
				copy_individual(&pg->newpop[i],&pool[ind],pg->chromsize);
#endif

			pthread_mutex_unlock(&vmutex[ind]);
		}
		//force context switch
		struct timespec interval, remainder; 
		interval.tv_sec = 0; interval.tv_nsec = 50000;
		nanosleep(&interval,&remainder);

	}       /* One GA run is over  */


	pthread_mutex_lock(&malloc_mutex);	

	cout << "In thread " << pg->id << '\n';
	printf("\n===================================================");
	//      printf("\nMax = %8.5f  Min = %8.5f   Avg = %8.5f",
	//	      pg->max_obj,pg->min_obj,pg->avg_obj);
	printf("\nNo. of mutations = %d ;  No. of x-overs = %d",
			pg->no_mutation,pg->no_xover);

	printf("\nFinestIndividual ");
	for (int j=0; j< pg->num_var; j++)
		printf(" %lf",pg->best_ever.x[j]);
	printf("\n");
	//HLEE-4-18-2009: added pg->id printed for Byz-faulty simulation
	printf("\n BestEverFitness= %d %lf from generation = %d",
			pg->id, pg->best_ever.obj, pg->best_ever_gen);

	//If this is a binary GA i.e. UIUC problem 
	if (pg->num_discr_var) 
	{
		unsigned u = pg->best_ever.chrom[0];
		int nlc;  
		nlc = (u & mask) + (u & (mask << 1)) + (u & (mask << 2)) + 1;

		//cout << '\n' << (u & mask) << " " << (u & (mask << 1)) << " " << (u & (mask << 2)) << '\n';
		printf("\nNumber of Lifecyles %d \n", nlc); 
		printf("\nString =\n");
		writechrom1(pg->best_ever.chrom,pg->best_ever.lastbit);
		//HLEE-4-18-2009: added pg->id printed for Byz-faulty simulation
		printf("\n BestEverFitness= %d %lf from generation = %d",
				pg->id, pg->best_ever.obj,pg->best_ever_gen);
	}
	/*HLEE*/
	printf("\n Best in each gen = %d  %d ", pg->max_gen, pg->id);
	for (int g = 0; g <= pg->max_gen; ++g)
		printf(" %lf ", pg->best_in_each_gen[g]);

	printf("\n===================================================");
	printf("\n\n");
	pg->free_all();

	pthread_mutex_unlock(&malloc_mutex);


}

//Initialize the pool with random variables
void initpool()
{
	double u;
	int k,k1,i,j,j1,stop;
	unsigned mask=1,nbytes;

	//This is used for binary GA
	gchromsize = (glchrom/UINTSIZE);
	if(glchrom%UINTSIZE) gchromsize++;
	nbytes = gchromsize*sizeof(unsigned);

	for(j = 0; j < pool_size; j++)
	{
		if((pool[j].chrom = (unsigned *) malloc(nbytes)) == NULL)
			nomemory("pool chromosomes");
	}
	if((gbest_ever.chrom = (unsigned *) malloc(nbytes)) == NULL)
		nomemory("best_ever chromosomes");


	for(j = 0; j < pool_size; j++)
	{
		if((pool[j].lastbit = (double *) malloc(MAX_LIFECYCLE*sizeof(double))) == NULL)
			nomemory("pool lastbit");
	}
	if((gbest_ever.lastbit = (double *) malloc(MAX_LIFECYCLE*sizeof(double))) == NULL)
		nomemory("best_ever chromosomes");


	for (k=0; k< pool_size; k++)
	{
		pool[k].obj = pool[k].mod_obj = 0;//This is why the snapshot showed 0
		pool[k].parent1 = pool[k].parent2 = 0;
		//implicitly doing this only in case of real GA
		for (j=gnum_discr_var; j<=gnum_var-1; j++)
		{
			u = rand()/(double) RAND_MAX;
			pool[k].x[j] = gx_lower[j] * (1-u) + gx_upper[j] * u;
		}
		for(k1 = 0; k1 < MAX_LIFECYCLE;k1++)
			pool[k].lastbit[k1] = rand()/(double) RAND_MAX;

		//This is executed because an inconsequential value of 72 is 
		//provided in the input file	
		for(k1 = 0; k1 < gchromsize; k1++)
		{
			pool[k].chrom[k1] = 0;
			if(k1 == (gchromsize-1))
				stop = glchrom - (k1*UINTSIZE);
			else
				stop = UINTSIZE;
			// A fair coin toss /
			for(j1 = 1; j1 <= stop; j1++)
			{
				if(flip1(0.5))
					pool[k].chrom[k1] = pool[k].chrom[k1]|mask;
				if (j1 != stop) pool[k].chrom[k1] = pool[k].chrom[k1]<<1;
			}

		}
		pool[k].obj = objectivecopy(pool[k].x);
		pool[k].mod_obj = pool[k].obj;

		//writechrom(pool[k].chrom);
	}
	copy_individual(&pool[0],&gbest_ever,gchromsize);
#if 0 
	//keep array of 10 and classify stuff
	//open a file
	//write this data to file
	//close file
	//carry on
	double distribution_array[10];
	double minrange = std::numeric_limits<double>::max(),maxrange = std::numeric_limits<double>::min(),intlength;

	//find minima maxima and poollength
	for(i = 0; i < pool_size ;i++)
	{
		pool[i].obj = objectivecopy(pool[i].x);
		if(pool[i].obj < minrange)
			minrange = pool[i].obj;
		if(pool[i].obj > maxrange)
			maxrange = pool[i].obj;
	}
	intlength = (maxrange - minrange)/10;
	cout << "minrange = " << minrange << " maxrange = " << maxrange <<  " intlength = " << intlength << '\n';

#endif
}

void free_pool()
{
	int j;
	for(j = 0; j < pool_size; j++)
	{
		free(pool[j].chrom);
		free(pool[j].lastbit);
	}
	free(gbest_ever.chrom);
	free(gbest_ever.lastbit);
}
/*====================================================================
  Writes a given string of 0's and 1's
  puts a `-` between each substring (one substring for one variable)
  Leftmost bit is most significant bit.
  ====================================================================*/
void writechrom(unsigned *chrom)
{
	int j, k, stop,bits_per_var,count=0;
	unsigned mask = 3, tmp;

	if (gBINGA == FALSE) return;
	if (chrom == NULL) error_ptr_null("chrom in writechrom");
	//bits_per_var = glchrom/gnum_discr_var;
	for(k = 0; k < gchromsize; k++)
	{
		tmp = chrom[k];
		if(k == (gchromsize-1))
			stop = glchrom - (k*UINTSIZE);
		else
			stop = UINTSIZE;

		for(j = 0; j < stop; j+=2)
		{
			int a = tmp & mask;
			switch(a)
			{
				case 3:
					cout << "11";
					break;
				case 2:
					cout << "10";
					break;
				case 1: 
					cout << "01";
					break;
				case 0: 
					cout << "00";
					break;
			}
			tmp = tmp >> 2;
		}
	}
	printf("\n");
}

#if 1 
/* Write Chrom Formatted to send to UIUC*/
void writechrom1(unsigned *chrom,double* lb)
{
	int i,i1,i2,ind,stage,op;
	int j, k, stop,bits_per_var,count=0;
	unsigned mask = 1, tmp;

	if (gBINGA == FALSE) return;
	if (chrom == NULL) error_ptr_null("chrom in writechrom");

	int nlc = (chrom[0] & 1) + (chrom[0] & 2) + (chrom[0] & 4) + 1;

	for(stage = 0; stage < nlc ; stage++)
	{
		cout << "Lifecycle " << stage + 1 << '\n';
		for(i = 0; i < 12; i++)
		{
			i1  = (stage)*BITS_PER_LIFECYCLE + 2*i + LC_ENCODE_LENGTH; 
			ind = i1 / (8*sizeof(unsigned));  //The index to use for the unsigned array
			i2  = i1 - ind*8*sizeof(unsigned);//The offset to use for the element of the unsigned array
			op = (((chrom[ind])>>(i2+1)) & mask )*2 + (((chrom[ind]) >> (i2)) & mask);
			cout << op << " "; 
		}
		cout << lb[stage];	
		cout << "\n\n";
	}

}
#endif

/*
 * construct array of individuals
 *1. construct all the ga classes
 *2. start running the threads
 *
 *
 *
 *
 *
 */
int main(int argc, char* argv[])
{
	int i,j,t,rc;
	void* status;
	GAThread *ga;	//should these be created on the stack/heap/readwrite mem?
	if(argc != 2)
	{
		cout << "Usage: ./pool NUM_THREADS < inputfile\n";
		return -1;
	}
	else
	{
		gnum_threads = atoi(argv[1]);
		cout << "Number of Threads = " << gnum_threads << '\n';
		cout << "ByzPercentage=" << BYZANTINE_FRACTION << '\n'; 
	}
	gnum_failed = 0;
	pthread_attr_t attr;
	ga = new GAThread[gnum_threads];
	vmutex = new pthread_mutex_t[gnum_threads*MAXPOPSIZE];
	input_parameters();		//get input parameters
	gMINM = MINMAX;			//Indicates that we must minimize objective function 

	//initreport(fp_out);
	
	//Initialise the mutexes
	for(t = 0; t < gnum_threads*gpop_size; t++)
	{
 	       pthread_mutex_init(&vmutex[t], NULL);
	}
	pthread_mutex_init(&malloc_mutex,NULL);
#ifdef CRASHFAILURE
	//gfailure_prob = FAILURE_PROBABILITY;
	gfailure_prob = 1.0/(2.0*gmax_gen); 	//aim for num_threads/2 to fail in num_thread*gmax_gen tries
	pthread_mutex_init(&failure_mutex,NULL);
#endif
	//Create the pool of Individuals
	pool_size = gnum_threads*gpop_size;
	pool = new INDIVIDUAL[pool_size];	//this size is necessary for the strategy we choose to write back individuals to the pool

	//Initialize the pool of Individuals
	initpool();
//return -100;
	//Set thread attribute to joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//Spawn the threads
	for(i = 0; i < gnum_threads; i++)
	{
		ga[i].copy_global();
		ga[i].seed = gbasic_seed + (1.0-gbasic_seed)*(double)(i)/(double)gnum_threads;	//setting the seeds	
		if (ga[i].seed > 1.0) printf("\n Warning !!! seed number exceeds 1.0");

		ga[i].id = i;	
		rc = pthread_create(&(ga[i].thread), &attr,GAThread::ThreadFunction , (void *)&ga[i]);
		if (rc)
		{
	        	printf("ERROR; return code from pthread_create() is %d\n", rc);
		}
	}
	//attribute
	pthread_attr_destroy(&attr);
	//wait for the threads to end
	for(i=0; i<gnum_threads; i++)
        {
		rc = pthread_join( ga[i].thread, &status);
		if (rc)
		{
        	printf("ERROR; return code from pthread_join() is %d\n", rc);
        	exit(-1);
      		}
		
	}
#ifdef POOL_SNAPSHOT
	//keep array of 10 and classify stuff
	//open a file
	//write this data to file
	//close file
	//carry on
	double distribution_array[10];
	double minrange = std::numeric_limits<double>::max(),maxrange = std::numeric_limits<double>::min(),intlength;
	
	//find minima maxima and poollength
	for(i = 0; i < pool_size ;i++)
	{
		pool[i].obj = objectivecopy(pool[i].x);
		if(pool[i].obj < minrange)
			minrange = pool[i].obj;
		if(pool[i].obj > maxrange)
			maxrange = pool[i].obj;
	}
	intlength = (maxrange - minrange)/10;
cout << "minrange = " << minrange << " maxrange = " << maxrange <<  " intlength = " << intlength << '\n';
		

#endif
	for(t = 0; t < gnum_threads*gpop_size; t++)
	{
 	       pthread_mutex_destroy(&vmutex[t]);
	}
	pthread_mutex_destroy(&malloc_mutex);
#ifdef CRASHFAILURE
	pthread_mutex_destroy(&failure_mutex);
	cout << "Number of threads that failed = " << gnum_failed << '\n';
#endif
	//fclose(fp_out);
	free_pool();
	delete [] pool;
	delete [] ga;
	delete [] vmutex;
	cout << '\n' << " O.K Good bye !!!" << '\n';

	pthread_exit(NULL);

	

}


/*====================================================================
Reporting the user-specified parameters :
fp is the file pointer to output file.
====================================================================*/
void initreport(FILE* fp)
{
   int k;

   if (fp == NULL) error_ptr_null(" File fp in initreport");
   fprintf(fp,"\n\n=============================================");
   fprintf(fp,"\n             INITIAL REPORT                  ");
   fprintf(fp,"\n=============================================");
   if (gBINGA) fprintf(fp,"\n BINARY-CODED GA ");
   if (gREALGA)
   {
      if (gcross_type ==SBX) fprintf(fp,"\n REAL-CODED GA (SBX)");
      else                  fprintf(fp,"\n REAL-CODED GA (BLX)");
   }
   switch (gs_strategy) 
     {
     case 1 : fprintf(fp,"\n Tournament Selection Used (Size = %d)",gtourneysize); break;
     case 2 : fprintf(fp,"\n Roulette Wheel Selection Used"); break;
     case 3 : fprintf(fp,"\n Stochastic Remainder RW Selection Used"); break;
     }
   switch (gx_strategy)
   {

    case ONESITE : fprintf(fp,"\n Crossover Strategy : 1 xsite with swapping"); break;
    case UNIF : fprintf(fp,"\n Crossover Strategy : Uniformly all variables 50 %%");
                 break;
    case ONLINE : fprintf(fp,"\n Crossover Strategy : On a straight line"); break;
    default  : fprintf(fp,"\n CROSSOVER NOT SET CORRECTLY "); break;
   }
   if (gBINGA) 
     fprintf(fp,"\n Mutation Strategy: Bit-wise Mutation");
   else
     fprintf(fp,"\n Mutation Strategy: Polynomial Mutation");
   fprintf(fp,"\n Variable Boundaries : ");
   if (gRIGID) fprintf(fp," Rigid");
   else       fprintf(fp," Flexible");
   fprintf(fp,"\n Population size            : %d",gpop_size);
   fprintf(fp,"\n Total no. of generations   : %d",gmax_gen);
   fprintf(fp,"\n Cross over probability     : %6.4f",gp_xover);
   fprintf(fp,"\n Mutation probability       : %6.4f",gp_mutation);
   if (gSHARING)
   {
        fprintf(fp,"\n Sharing to be done :");
        fprintf(fp,"\n Sigma-share value          : %6.4f",gsigma_share);
   }
   if (gBINGA)
      fprintf(fp,"\n String length              : %d",glchrom);
   fprintf(fp,"\n Number of variables        : %d",gnum_var);
   fprintf(fp,"\n Total Runs to be performed : %d",gmaxrun);
   if ((gREALGA) && (gcross_type == SBX)) {
     fprintf(fp,"\n Exponent (n for SBX)       : %7.2f",gn_distribution_c);
     fprintf(fp,"\n Exponent (n for Mutation)  : %7.2f",gn_distribution_m);
   }
   if (gs_strategy == 1)
   fprintf(fp,"\n Lower and Upper bounds     :");
   for (k=0; k<=gnum_var-1; k++)
     fprintf(fp,"\n   %8.4f   <=   x%d   <= %8.4f",gx_lower[k],k+1,gx_upper[k]);
   fprintf(fp,"\n=================================================\n");
//   app_initreport();
}


/*====================================================================
Prints an error message and terminates the program
====================================================================*/
void nomemory(const char* string)
{
   printf("\nmalloc: out of memory making %s!!\n",string);
   printf("\n Program is halting .....");
   exit(-1);
}
/*==============================================================
Gives error message of null pointer  and terminates the program.
==============================================================*/
void error_ptr_null(const char *string)
{
   printf("\n Error !! Pointer %s found Null !",string);
   printf("\n Program is halting .....");
   exit(-1);
}

/*====================================================================
Copys contents of one individual into another.
====================================================================*/
void copy_individual(INDIVIDUAL *source,INDIVIDUAL *target,int chromsz)
{
   int k;

   if (source==NULL) error_ptr_null("source in copy_individual");
   if (target==NULL) error_ptr_null("target in copy_individual");
   for (k=0; k< MAXVECSIZE; k++)
      target->x[k] = source->x[k];
   target->mod_obj = source->mod_obj;
   target->obj = source->obj;
   target->parent1 = source->parent1;
   target->parent2 = source->parent2;
   target->cross_var = source->cross_var;
   for (k=0; k< chromsz; k++)
        target->chrom[k] = source->chrom[k];
  
   for (k=0; k< MAX_LIFECYCLE; k++)
   	target->lastbit[k] = source->lastbit[k];
}


/*====================================================================
Decodes the string of the individual (if any) and puts the values in
the array of floats.
====================================================================*/
//@gautam allocated temp on stack
void GAThread::decode_string(INDIVIDUAL* ptr_indiv)
{
   double coef,*temp;
   
   int j;

   if (ptr_indiv == NULL) error_ptr_null("ptr_indiv in decode_string");
   if (BINGA)
   {
	   //note that only one double per discrete variable
	   //entire string gets decoded to one double
     temp = (double *) malloc(num_discr_var * sizeof(double));
     for(j=0; j<=num_discr_var - 1; j++) temp[j] = 0.0;
     decodevalue(ptr_indiv->chrom,temp);
     coef = pow(2.0,(double)(lchrom/num_discr_var)) - 1.0;
     for(j=0; j<=num_discr_var - 1; j++)
     {
        temp[j] = temp[j]/coef;
        ptr_indiv->x[j] = temp[j]*x_upper[j] + (1.0 - temp[j])*x_lower[j];
     }
     free(temp);
   }
}



/*====================================================================
Decodes the value of a group of binary strings and puts the decoded
values into an array 'value'.
====================================================================*/
void GAThread::decodevalue(unsigned *chrom,double value[])
{
    int k,j,stop,tp,bitpos,mask=1,position,bits_per_var;
    double bitpow;

    if (BINGA == FALSE) return;
    if (chrom == NULL) error_ptr_null("chrom in decodevalue");


    bits_per_var = lchrom/num_discr_var; //good insight...seems the discrete variables are concatenated
    for(k = 0; k < chromsize; k++)
    {
        if(k == (chromsize-1))
            stop = lchrom-(k*UINTSIZE);
        else
            stop = UINTSIZE;
        /* loop thru bits in current byte */
        tp = chrom[k];
        for(j = 0; j < stop; j++) {
            bitpos = j + UINTSIZE*k;
            /* test for current bit 0 or 1 */
            if((tp&mask) == 1) {
                position = bitpos / bits_per_var;
                bitpos -= position * bits_per_var;
                bitpow = pow(2.0,(double)(bits_per_var- bitpos-1));
                value[position] += bitpow;
            }
            tp = tp>>1;
        }
    }
}


/*
 * For the Benchmark Problems
 * f1 minimization, minima = 0 at all xi = 0;
 * f2 minimization, minima = 0 at x1 = x2 = 1
 * f3 minimization, minima = 0 at all xi = 0;
 * f4 minimization, minima = -4189.83 at all xi = 420.9687
 */
/*====================================================================
OBJECTIVE FUNCTION  ( Supposed to be minimized) :
Change it for different applications
====================================================================*/
double objectivecopy(double *x)
{

  double term1,term2, term3, pi;
  double g, penalty_coef;
  double fit;
   if (x == NULL) error_ptr_null("x in objective()");

//-10 < xi < 10
//i going from 0 to 4
#ifdef f1
//	gMINM = 1;
	fit = 0;
	fit = x[0]*x[0] + 10*x[1]*x[1] + 100*x[2]*x[2] + 1000*x[3]*x[3] + 10000*x[4]*x[4];
	return fit;
#endif

//only two variables
//1.5 <= xi <= 1.5 
#ifdef f2
//	gMINM = 1;
	fit = 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1-x[0])*(1-x[0]);
	return fit;
#endif

//currently considering 20 variables as given in the paper
//-5.12 <= xi <= 5.12
#ifdef f3
//	gMINM = 1;
	fit = 20;
	for(int i = 0; i < 20; i++)
	{
		fit += (x[i]*x[i] - cos(2*pi*x[i]));
	}
	return fit;
#endif

//using 10 variables as in the paper
//-500 <= xi <= 500
#ifdef f4
//	gMINM = 1;
	fit = 0;
	for(int i = 0; i < 10 ; i++)
	{
//		cout << x[i] << " " ;
		fit += (-x[i]*sin(sqrt(abs(x[i]))));
	}
//	cout << endl;
	return fit;	
#endif

#ifdef sync
   gMINM = 1;
   if( *x == 0 ) return .5;
   else return sin(*x)/(2*(*x));
#endif

#ifdef prob1
   gMINM  = -1;
   term1 = (x[0]*x[0]+x[1]-11.0)*(x[0]*x[0]+x[1]-11.0);
   term2 = (x[0]+x[1]*x[1]- 7.0)*(x[0]+x[1]*x[1]- 7.0);
   term3 = term1+term2;
   
   penalty_coef = 0.0;
   g = (square(x[0]-5.0) + square(x[1]))/26.0 - 1.0;
   if (g < 0.0) term3 = term3 + penalty_coef * g * g;
   return(1.0/(1.+term3));
#endif

#ifdef can
   gMINM = 1;
   pi = 4.0 * atan(1.0);
   term3 = pi * x[0] * x[0]/2.0 + pi * x[0] * x[1];

   penalty_coef = 1.0e4;
   g = (pi * x[0] * x[0] * x[1]/4.0 - 400.0)/400.0;
   if (g < 0.0) term3 = term3 + penalty_coef * g * g;
   return(term3);
#endif

#ifdef bin1
   gMINM = -1;
   return (*x) < 0 ? -(*x): (*x);
#endif
}

/*
 * For the Benchmark Problems
 * f1 minimization, minima = 0 at all xi = 0;
 * f2 minimization, minima = 0 at x1 = x2 = 1
 * f3 minimization, minima = 0 at all xi = 0;
 * f4 minimization, minima = -4189.83 at all xi = 420.9687
 */
/*====================================================================
OBJECTIVE FUNCTION  ( Supposed to be minimized) :
Change it for different applications
====================================================================*/
double GAThread::objective(double *x)
{

  double term1,term2, term3, pi;
  double g, penalty_coef;
  double fit;
   if (x == NULL) error_ptr_null("x in objective()");

//-10 < xi < 10
//i going from 0 to 4
#ifdef f1
	MINM = 1;
	fit = 0;
	fit = x[0]*x[0] + 10*x[1]*x[1] + 100*x[2]*x[2] + 1000*x[3]*x[3] + 10000*x[4]*x[4];
	return fit;
#endif

//only two variables
//1.5 <= xi <= 1.5 
#ifdef f2
	MINM = 1;
	fit = 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1-x[0])*(1-x[0]);
	return fit;
#endif

//currently considering 20 variables as given in the paper
//-5.12 <= xi <= 5.12
#ifdef f3
	MINM = 1;
	fit = 20;
	for(int i = 0; i < 20; i++)
	{
		fit += (x[i]*x[i] - cos(2*pi*x[i]));
	}
	return fit;
#endif

//using 10 variables as in the paper
//-500 <= xi <= 500
#ifdef f4
	fit = 0;
	for(int i = 0; i < 10 ; i++)
	{
//		cout << x[i] << " " ;
		fit += (-x[i]*sin(sqrt(abs(x[i]))));
	}
//	cout << endl;
	return fit;	
#endif

#ifdef sync
   MINM = 1;
   if( *x == 0 ) return .5;
   else return sin(*x)/(2*(*x));
#endif

#ifdef prob1
   MINM  = -1;
   term1 = (x[0]*x[0]+x[1]-11.0)*(x[0]*x[0]+x[1]-11.0);
   term2 = (x[0]+x[1]*x[1]- 7.0)*(x[0]+x[1]*x[1]- 7.0);
   term3 = term1+term2;
   
   penalty_coef = 0.0;
   g = (square(x[0]-5.0) + square(x[1]))/26.0 - 1.0;
   if (g < 0.0) term3 = term3 + penalty_coef * g * g;
   return(1.0/(1.+term3));
#endif

#ifdef can
   MINM = 1;
   pi = 4.0 * atan(1.0);
   term3 = pi * x[0] * x[0]/2.0 + pi * x[0] * x[1];

   penalty_coef = 1.0e4;
   g = (pi * x[0] * x[0] * x[1]/4.0 - 400.0)/400.0;
   if (g < 0.0) term3 = term3 + penalty_coef * g * g;
   return(term3);
#endif

#ifdef bin1
   MINM = -1;
   return (*x) < 0 ? -(*x): (*x);
#endif
}

