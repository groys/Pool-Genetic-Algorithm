//in mono_over_gens for each different number of threads the value corresponding to generation g is the best amonst all threads at generation g

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#define VARY_NUM_THREADS  6 /* 1, 2, 4, 8, 16, 32 */
#define MAX_GEN  1000 
#define RUN  10
#define FILENAME_LEN  25
#define NUM_FILES_TO_GO_THRU  60 /* 60 = 6 x 10 */
#define LINE_LEN  80
#define KEY_STR  "Best ever fitness = "
#define KEY_LEN   21 
#define KEY_STR1 "Best in each gen = "
#define KEY_LEN1  20 
#define VAL_LEN  10
#define GEN_LEN  3 
#define EOS  '\0'

#define MINM -1	//1 for minimization and -1 for maximization

/* function prototypes */
int execute(void);
int gather_data(void);


int main(void)
{

/* This program is to perform the two steps automatically in the GA simulation.

   1. Execute five times each with varying numbers of threads:
      1, 2, 4, 8, 16, 32, and 64. 
      The results are written in new files named as t_n1_n2, 
      where n1 is the number of threads and n2 is the number of run 
      (varying from 1 to 5).

   2. On the data produced by the step 1, analyze the followings:

    2a. Collect the data "Best ever" and compute an average over all 
        the threads in one run and all the five runs.

    2b. Collect the data "from generation" and compute an average over all 
        the threads in one run and all the five runs.

        The results of steps 2a and 2b will be written in a file called 
        "best_val_gen".
        The file has three columns: the first column is the number of threads
        and will be plotted x-axis. 
        The second colum is the average best value for the corresponding 
        number of threads and will be plotted y-axis in the first graph 
        (use gnuplot).
        The third colum is the average number of generation in which best 
        value was seen for the corresponding number of threads and will be 
        plotted y-axis in the second graph.
 
    2c. Collect the data "Best in each gen" in two two-dimensional float array 
		    "mono_best_vals" of VARY_NUM_THREADS x MAX_GEN, and
		    "each_gen_best_avg" of VARY_NUM_THREADS x MAX_GEN.
	In mono_best_vals, keep the best among all the threads.
	This will be written in a file called "mono_over_gens":
		col1: gen#, col2: row1 of "mono_best_vals", 
		col3: row2 of "mono_best_vals", and so on. 
	In each_gen_best_avg, keep the average of the best values seen by 
	all threads in each gen. This will be written in "avg_over_gens":
		col1: gen#, col2: row1 of "each_gen_best_avg", 
		col3: row2 of "each_gen_best_avg", and so on. 

    2d. (Will be used for the other graphs with the synchronous version...)

    Questions?  E-mail to hlee@cs.tamu.edu
*/


    int exec_code = execute();
    if (exec_code != 0) {
       printf("Error in execute() ....\n");
       return 0;
    }


    int gather_code = gather_data();

/*    if (gather_code != NUM_FILES_TO_GO_THRU) {
       printf("Error in gather_data() ....\n");
       return 0;
    } 
*/

    printf("Successful! The output file name is \"best_val_gen\".\n");
    return 0;
}


int execute(void)
{
    int i, t, run_num;
    char cmd_str[LINE_LEN], int_to_str[4];
    char ofstring[FILENAME_LEN], ofname[FILENAME_LEN]; 

    system("rm ./rslt/*.rslt");
    system("rmdir ./rslt");
    system("mkdir ./rslt");

    system("rm best_val_gen");
    system("rm mono_over_gens");
    system("rm avg_over_gens"); 

    i = 1; /* i is used as the number of threads */
    for (t = 1; t <= VARY_NUM_THREADS; ++t) { 
     /* construct the file name for output */
        ofstring[0] = EOS;
        int_to_str[0] = EOS;
        sprintf(ofstring, "./rslt/t_"); 
        sprintf(int_to_str, "%d", i);
        strcat(ofstring, int_to_str);

        for (run_num = 1; run_num <= RUN; run_num++) {
            int_to_str[0] = EOS;
            strcpy(ofname, ofstring);
            sprintf(int_to_str, "_%d", run_num);
            strcat(ofname, int_to_str);
            printf("%s\n", ofname);
            sprintf(cmd_str, "./pool %d < input1 > %s.rslt", i, ofname);
            printf("%s\n", cmd_str);
            system(cmd_str);
        } /* end inner for */
        i *= 2;
    }
            
    return 0;
}

int gather_data(void)
{
    FILE *ifp, *ofp, *ofp1, *ofp2;
    char ifstring[FILENAME_LEN], ifname[FILENAME_LEN]; 
    char c, int_str[4];
    char *a_line;
    char key[KEY_LEN], key1[KEY_LEN1], val_str[VAL_LEN], gen_str[GEN_LEN];
    int  z, g, i, j, k, t, val_num, gen_num;
    int  gen, max_gen;
    int  num_threads, run_num, num_files = 0;
    double val, val_sum, val_avg, gen_sum, gen_avg;
    double each_gen_val_sum[MAX_GEN] = {0};
    double mono_best_val_sum[MAX_GEN] = {0};
    double each_gen_best_avg[VARY_NUM_THREADS][MAX_GEN] = {0}; 
    float  bval;  //, best_vals[VARY_NUM_THREADS][MAX_GEN] = {0};
    float  mono_best_vals[VARY_NUM_THREADS][MAX_GEN];
	       	
		for(i = 0; i < VARY_NUM_THREADS; i++)
			for(j = 0; j < MAX_GEN ; j++)
				mono_best_vals[i][j] = MINM*1e10;

    ofp = fopen("best_val_gen", "w");  /* file open for writing */
    ofp1 = fopen("mono_over_gens", "w"); /*output file for monotone incr */
    ofp2 = fopen("avg_over_gens", "w");  /* output for best in each gen */

    num_threads = 1;
    for (t = 0; t < VARY_NUM_THREADS; ++t) {  // outer for varying threads
	/* initialize variables for each thread number */
        val_sum = 0.0;
        val_num = 0;
        gen_sum = 0.0; 
        gen_num = 0;
        /* construct the file name for input */
				ifstring[0] = EOS;
				int_str[0] = EOS;
        sprintf(ifstring, "./rslt/t_");
        sprintf(int_str, "%d", num_threads);
        strcat(ifstring, int_str);

        for (run_num = 1; run_num <= RUN; run_num++) 
	{  //inner for multiple runs
            int_str[0] = EOS;
            memset(ifname, EOS, FILENAME_LEN);
            strcpy(ifname, ifstring);
            printf("ifname %s ifstring %s \n",ifname,ifstring);
            sprintf(int_str, "_%d", run_num);
            strcat(ifname, int_str);
            strcat(ifname, ".rslt");
            printf("%s\n", ifname);
            num_files++;  /* increment the number of files that are handled */

            ifp = fopen(ifname, "r");  /* open for reading */
            if(ifp == NULL) {
               printf("Error opening file %s \n",ifname);
               return -2;
            } 
            while ((c = getc(ifp)) != EOF) {
               if (c == 'B') { //Get to the point of Best ever fitness
                   key[0] = c;
                   for (i = 1; i < KEY_LEN - 1; i++) 
                       key[i] = getc(ifp);
                   key[i] = EOS;

                   printf("key = %s\n", key);
                   if (!(strncmp(key, KEY_STR, KEY_LEN-1))) {/* key matched */ 
                      j = 0; 
                      while ((!isspace(c = getc(ifp))) && (j < VAL_LEN)) {
                          val_str[j] = c;
                          j++;
                      }
                      val_str[j] = EOS;
                      printf("val_str = %s\n", val_str);
                      val = atof(val_str); 
                      val_num++;
                      val_sum += val;
                      printf("val_sum = %f  val_num = %d\n", val_sum, val_num);

                      /* move the ifp to the gen number */
                     //   while (getc(ifp) != '(');
                      while (getc(ifp) != '=');
                      getc(ifp); /* to ignore the space before the gen no. */
                      k = 0; 
                      while (!isspace((c = getc(ifp))) && (k < GEN_LEN)) {
                           gen_str[k] = c;  
                           k++;
                      }
                      gen_str[k] = EOS;
                      printf("gen_str = %s\n", gen_str);
                      gen = atoi(gen_str); 
                      gen_num++;
                      gen_sum += gen;
                      printf("gen_sum = %f  gen_num = %d\n", gen_sum, gen_num);
                   } /* endif */
                   while ((c = getc(ifp)) != 'B'); 
                   //Got to the point of "Best in each gen" 
                   key1[0] = c;
                   for (i = 1; i < KEY_LEN1 - 1; i++) 
                       key1[i] = getc(ifp);
                   key1[i] = EOS;

                   printf("key1 = %s\n", key1);
                   if (!(strncmp(key1, KEY_STR1, KEY_LEN1-1))) {// key matched 
                       fscanf(ifp, "%d", &max_gen);
		      printf("max_gen = %d\n", max_gen); 
                       fscanf(ifp, "%f", &bval); // to skip gen 0's value
                       for (g = 1; g <= max_gen; ++g) 
		       {  // for each gen
                          fscanf(ifp, "%f", &bval);
                          each_gen_val_sum[g] += bval;
			if (MINM*mono_best_vals[t][g] > MINM*bval)
                             mono_best_vals[t][g] = bval;
                        if (MINM*mono_best_vals[t][g-1] < MINM*mono_best_vals[t][g])
                             mono_best_vals[t][g] = mono_best_vals[t][g-1];
                       } /* endfor each gen */
                   } /* endif */
               } /* endif */
            } /* endwhile not EOF */
            printf("End of file found \n" );
            fclose(ifp);
            for (g = 1; g <= max_gen; ++g)
               mono_best_val_sum[g] += mono_best_vals[t][g]; 
        } /* end inner for multiple run*/
        val_avg = val_sum / val_num;
        gen_avg = gen_sum / gen_num;
        fprintf(ofp, "%d  %f  %f\n", num_threads, val_avg, gen_avg); 

/* HLEE-Nov-9: 
   Each_gen_val_sum and mono_best_val_sum have the values added over 
   all threads over all runs per gen, thus divide it by (num_threads * RUN) */

        z = num_threads * RUN;
	for (g = 1; g <= max_gen; ++g) { 
            each_gen_best_avg[t][g] = each_gen_val_sum[g] / z;
            mono_best_vals[t][g] = mono_best_val_sum[g] / RUN;
        }

        for (g = 1; g <= max_gen; ++g)
            each_gen_val_sum[g] = 0.0;

        for (g = 1; g <= max_gen; ++g)
            mono_best_val_sum[g] = 0.0;

	num_threads *= 2;
    }/* end outer for */ 

    for (g = 1; g <= max_gen; ++g) {
        fprintf(ofp1, "%d", g);
        for (t = 0; t < VARY_NUM_THREADS; ++t)
            fprintf(ofp1, " %f", mono_best_vals[t][g]);
        fprintf(ofp1, "\n");
    }

    for (g = 1; g <= max_gen; ++g) {
        fprintf(ofp2, "%d", g);
        for (t = 0; t < VARY_NUM_THREADS; ++t)
            fprintf(ofp2, " %f", each_gen_best_avg[t][g]);
        fprintf(ofp2, "\n");
    }

    fclose(ofp);
    fclose(ofp1);
    fclose(ofp2);
    return num_files;
}
