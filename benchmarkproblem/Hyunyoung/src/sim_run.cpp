
#include "pool.h"



/* function prototypes */
int execute(char*, char*);
int gather_data(void);

double comp;
bool pred(double d){
	return d <= comp;
}

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
int main(int argc, char* argv[])
{

	if(argc != 3)
	{
		printf("Usage: ./sim_run benchmarkfunction byzantinefraction\n");
		printf("specify a dummy byzantine fraction if BYZANTINE_ANTIELITISM is not defined\n");
		return -1;
	}
	int exec_code = execute(argv[1], argv[2]);
	if (exec_code != 0) {
		printf("Error in execute() ....\n");
		return 0;
	}

	int gather_code = gather_data();
	if (gather_code != VARY_NUM_THREADS * RUN) {
		printf("Error in gather_data() ....\n");
		return 0;
	}

	printf("Successful! The output file name is \"best_val_gen\".\n");
	return 0;
}

int execute(char* function, char* byzfraction) {
	int i, t, run_num;
	char cmd_str[LINE_LEN], temp_str[LINE_LEN], int_to_str[4];
	char ofstring[FILENAME_LEN], ofname[FILENAME_LEN];


	sprintf(cmd_str,"rm %srslt/*.rslt",OUTPUTDIR);
	system(cmd_str);
	sprintf(cmd_str,"rmdir %srslt",OUTPUTDIR);
	system(cmd_str);
	sprintf(cmd_str,"mkdir %srslt",OUTPUTDIR);
	system(cmd_str);

	sprintf(cmd_str,"rm %sbest_val_gen",OUTPUTDIR);
	system(cmd_str);
	sprintf(cmd_str,"rm %smono_over_gens",OUTPUTDIR);
	system(cmd_str);
	sprintf(cmd_str,"rm %savg_over_gens",OUTPUTDIR);
	system(cmd_str);

	i = 1; /* i is used as the number of threads */
	for (t=1; t <= VARY_NUM_THREADS; ++t) { /* construct file name for output */
		ofstring[0] = EOS;
		int_to_str[0] = EOS;
		sprintf(ofstring, OUTPUTDIR);
		strcat(ofstring, "rslt/t_");
		sprintf(int_to_str, "%d", i);
		strcat(ofstring, int_to_str);

		for (run_num=1; run_num <= RUN; run_num++) { /* inner for each run */

			strcpy(ofname, ofstring);
			int_to_str[0] = EOS;
			sprintf(int_to_str, "_%d", run_num);
			strcat(ofname, int_to_str);
			strcat(ofname, ".rslt");
			printf("%s\n", ofname);

			sprintf(cmd_str, "./pool %s %d %s < ",function, i, byzfraction);
			printf("%s\n", cmd_str);
			strcat(cmd_str,INPUTDIR);
			sprintf(temp_str,"input%s > %s", function, ofname);
			strcat(cmd_str,temp_str);
			printf("%s\n", cmd_str);
			system(cmd_str);
		} /* end inner for each run */
		i *= 2;
	}

	return 0;
}


/*
Update for pool snapshot
Aug 1: Gautam
	0. There will be two sets of plots one for the initial pool and one for the final pool
	1. data for all runs of a particular number of threads, will be put in a multiset.
	2. Using some stl algorithms we can determine max,min, int-length = max-min/10 and number of points in each interval
	3. there will be one output file for each number of threads(have to choose how to represent the data for initial and final pool
*/


int gather_data(void) {
	FILE *ifp, *ofp, *ofp1, *ofp2;
	char ifstring[FILENAME_LEN], ifname[FILENAME_LEN],tempfilename[FILENAME_LEN];
	char c, int_str[4];
	char key[KEY_LEN], key1[KEY_LEN1];//, val_str[VAL_LEN], gen_str[GEN_LEN];
	int  g, i, j, t, val_num, gen_num, byz_threshold;
	int  gen, max_gen, th_id;
	int  num_threads, run_num, num_files = 0;
	double byz_frac; // to read in ByzPercentage
	double val, val_sum, val_avg, gen_sum, gen_avg;
	double each_gen_val_sum[MAX_GEN] = {0};
	double mono_best_val_sum[MAX_GEN] = {0};
	double each_gen_best_avg[VARY_NUM_THREADS][MAX_GEN] = {{0}};
	double  bval;  //, best_vals[VARY_NUM_THREADS][MAX_GEN] = {0};
	double  mono_best_vals[VARY_NUM_THREADS][MAX_GEN];
	int threadcount;	//counts only the threads that executed to termination, crashed/byzantine thread


	for (i = 0; i < VARY_NUM_THREADS; i++)
		for (j = 0; j < MAX_GEN ; j++)
			mono_best_vals[i][j] = MINMAX*1e10;

	sprintf(tempfilename,OUTPUTDIR);
	strcat(tempfilename,"best_val_gen");
	ofp = fopen(tempfilename, "w");    /* file open for writing */

	sprintf(tempfilename,OUTPUTDIR);
	strcat(tempfilename,"mono_over_gens");
	ofp1 = fopen(tempfilename, "w");    /* file open for writing */

	sprintf(tempfilename,OUTPUTDIR);
	strcat(tempfilename,"avg_over_gens");
	ofp2 = fopen(tempfilename, "w");    /* file open for writing */

	//Create file name
	//snapshot_init_thread
	//snapshot_final_thread
#ifdef POOL_SNAPSHOT
	multiset <double> snapshot_init, snapshot_final;
	char ssfilei[FILENAME_LEN], ssfilef[FILENAME_LEN];
	char int_str1[4];
	double d;
	int minl,maxl,it;
	int poolsize;
	FILE *sifp, *sffp;
	double min_range, max_range, int_length;
	unsigned int cnt,prev;
        char cmd_str[FILENAME_LEN];

	sprintf(cmd_str,"rm %ssnapshot/*",OUTPUTDIR);
	system(cmd_str);
	sprintf(cmd_str,"rmdir %ssnapshot",OUTPUTDIR);
	system(cmd_str);
	sprintf(cmd_str,"mkdir %ssnapshot",OUTPUTDIR);
	system(cmd_str);
#endif

	num_threads = 1;
	for (t = 0; t < VARY_NUM_THREADS; ++t) {  // for1 varying threads
		/* initialize variables for each thread number */
		val_sum = 0.0;
		val_num = 0;
		gen_sum = 0.0;
		gen_num = 0;
		/* construct the file name for input */
		ifstring[0] = EOS;
		int_str[0] = EOS;
		sprintf(ifstring,OUTPUTDIR);
		strcat(ifstring, "rslt/t_");
		sprintf(int_str, "%d", num_threads);
		strcat(ifstring, int_str);
		threadcount = 0;

#ifdef POOL_SNAPSHOT
		snapshot_init.clear();
		snapshot_final.clear();
		ssfilei[0] = ssfilef[0] = EOS;
		sprintf(ssfilei,OUTPUTDIR);
		strcat(ssfilei,"snapshot/snapshot_init");
		sprintf(ssfilef,OUTPUTDIR);
		strcat(ssfilef,"snapshot/snapshot_final");
		int_str1[0] = EOS;
		sprintf(int_str1,"_%d",1<<t);
		strcat(ssfilei,int_str1);
		strcat(ssfilef,int_str1);
		cout <<"ssfilei " << ssfilei << '\n';
		cout <<"ssfilef " << ssfilef << '\n';
		sifp = fopen(ssfilei,"w");
		sffp = fopen(ssfilef,"w");
#endif
		for (run_num = 1; run_num <= RUN; run_num++) {// for2 multiple runs

			int_str[0] = EOS;
			memset(ifname, EOS, FILENAME_LEN);
			strcpy(ifname, ifstring);
			//printf("ifname %s ifstring %s \n", ifname, ifstring);
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
			else
				cout << "file opened" << '\n';

			// proceed until the first 'B' that is for ByzPercentage
			while (((c = getc(ifp)) != EOF) && c != 'B') {}
			while ((c = getc(ifp)) != '=') {} // skip till '=' sign
			// read in ByzPercentage
			fscanf(ifp, "%lf", &byz_frac);
			printf("byz_frac=%lf\n", byz_frac);
			byz_threshold = byz_frac * num_threads;
			printf("byz_threshold=%d\n", byz_threshold);

#ifdef POOL_SNAPSHOT
			//Get to the point Initial Pool
			while((c = getc(ifp)) != 'I');
			fscanf(ifp,"nitial Pool %u",&poolsize);
			printf("Initial Pool Size is %u \n", poolsize);
			for(it = 0; it < poolsize; it++)
			{
				fscanf(ifp,"%lf",&d);
				snapshot_init.insert(d);
			}
			printf("Initial Multiset size is %u\n", snapshot_init.size());
#endif

			// now read the data corresponding to all threads in this file
			while ((c = getc(ifp)) != EOF) {
				while (c != 'B' && c != 'P' && c != EOF) {
					c = getc(ifp);
				} // get to the point of BestEverFitness or EOF
				if (c == EOF) break;  // to exit the while loop
#ifdef POOL_SNAPSHOT
				if(c == 'P')
				{
					fscanf(ifp,"ool %u",&poolsize);
					printf("Final Pool Size is %u \n", poolsize);
					for(it = 0; it < poolsize; it++)
					{
						fscanf(ifp,"%lf",&d);
						snapshot_final.insert(d);
					}
					printf("Final Multiset size is %u\n", snapshot_final.size());
				}
				else
				{
#endif
				key[0] = c;
				for (i = 1; i < KEY_LEN - 1; i++)
					key[i] = getc(ifp);
				key[i] = EOS;

				printf("key = %s\n", key);
				if (!(strncmp(key, KEY_STR, KEY_LEN-1))) {/* key matched */
					// read in the value of BestEverFitness
					fscanf(ifp, "%d", &th_id);

					if (th_id >= byz_threshold) { // gather data only from  so for threshold 0 all threads are getting read
						// non-Byz-faulty processor
						threadcount++;

						fscanf(ifp, "%lf", &val);
						val_num++;
						val_sum += val;
						printf("val_sum = %lf  val_num = %d\n",val_sum,val_num);

						/* move the ifp to the gen number */
						while (getc(ifp) != '=') {}
						getc(ifp); /* to ignore the space before the gen no. */
						//read in the generation in which the best fitness was seen
						fscanf(ifp, "%d", &gen);
						gen_num++;
						gen_sum += gen;
						printf("gen_sum = %lf  gen_num = %d\n",gen_sum,gen_num);
					} /* endif non-Byz-faulty */
				} /* endif key matched */
				while ((c = getc(ifp)) != 'B') {}
				// got to the point of "Best in each gen"
				key1[0] = c;
				for (i = 1; i < KEY_LEN1 - 1; i++) key1[i] = getc(ifp);
				key1[i] = EOS;

				printf("key1 = %s\n", key1);
				if (!(strncmp(key1, KEY_STR1, KEY_LEN1-1))) { // key matched
					fscanf(ifp, "%d", &max_gen); // read the # of generations
					fscanf(ifp, "%d", &th_id);   // read the processor id
					printf("max_gen = %d  th_id = %d\n", max_gen, th_id);

					if (th_id >= byz_threshold) { // gather data only from
						// non-Byz-faulty processor
						printf("proc %d: non-faulty...\n", th_id);
						fscanf(ifp, "%lf", &bval); // to skip gen 0's value
						for (g = 1; g <= max_gen; ++g) { // for3 each gen
							//read the value for this gen for this thread
							fscanf(ifp, "%lf", &bval);
							each_gen_val_sum[g] += bval;

							//note what this step means...
							//collecting the best value for this run at gen g
							//check for(g = 1; g <=max_gen ... below which averages this
							if (MINMAX*mono_best_vals[t][g] > MINMAX*bval)
								mono_best_vals[t][g] = bval;
							//@g this step seems unnecessary because everything is
							//monotonic anyway
							if (MINMAX*mono_best_vals[t][g-1]<MINMAX*mono_best_vals[t][g])
								mono_best_vals[t][g] = mono_best_vals[t][g-1];
						} /* end for3 each gen */
					} /* endif non-Byz-faulty */
					printf("continue with next thread data ...\n");
				} /* endif key matched */
#ifdef POOL_SNAPSHOT
				}
#endif
			} /* end while not EOF */
			printf("End of file found \n" );
			fclose(ifp);
			printf("File closed\n");
			for (g = 1; g <= max_gen; ++g)
				mono_best_val_sum[g] += mono_best_vals[t][g];
			printf("mono best done\n");
		} /* end for2 multiple run*/

#ifdef POOL_SNAPSHOT
//do some analysis and write stuff out to files
//need two modes based on the separation between minimum and maximum
//one normal mode interval are just (max-min)/numint
//one log mode where intervals are logarithmic
		min_range = *(snapshot_init.begin());
		max_range = *(--snapshot_init.end());
cout << "Init min_range " << min_range << endl;
cout << "Init max_range " << max_range << endl;
//fprintf(sifp,"Init min_range %lf\n",min_range);
//fprintf(sifp,"Init max_range %lf\n",max_range);
		minl = log10(min_range);
		maxl = log10(max_range);
//fprintf(sifp,"minl %d maxl %d\n",minl,maxl);
		//log mode
		if(!(min_range == 0))// && (maxl - minl > 3))
		{
			prev =0;
			for(it = minl ; it <= maxl +1; it++)
			{
				comp = pow(10.0,(double)it);
				cnt =  count_if(snapshot_init.begin(),snapshot_init.end(),pred);
				fprintf(sifp,"%1.0e %d\n",comp,cnt-prev);
				cout << "Init Interval " <<it << " " << cnt - prev<< endl;
				prev = cnt;
			}
		}
		else
		{
			int_length = (max_range - min_range)/INTERVALS;
			prev = 0;
			for(it = 0; it <= INTERVALS ; it++)
			{
				comp = min_range + it*int_length;
				cnt =  count_if(snapshot_init.begin(),snapshot_init.end(),pred);
				fprintf(sifp,"%1.3e %d\n",comp,cnt-prev);
				cout << "Init Interval " <<it << " " << cnt - prev<< endl;
				prev = cnt;
			}
		}
		min_range = *(snapshot_final.begin());
		max_range = *(--snapshot_final.end());

cout << "Final min_range " << min_range << endl;
cout << "Final max_range " << max_range << endl;
//fprintf(sffp,"final min_range %lf\n",min_range);
//fprintf(sffp,"final max_range %lf\n",max_range);
		minl = log10(min_range);
		maxl = log10(max_range);
//fprintf(sffp,"minl %d maxl %d\n",minl,maxl);
		//log mode

		if(!(min_range == 0))// && (maxl - minl > 3))
		{
		printf("(maxl - minl > 3)");
			prev =0;
			for(it = minl ; it <= maxl +1; it++)
			{
//printf("it = %d\n",it);
				comp = pow(10.0,(double)it);
				cnt =  count_if(snapshot_final.begin(),snapshot_final.end(),pred);
				fprintf(sffp,"%1.0e %d\n",comp,cnt-prev);
				cout << "Final Interval " <<it << " " << cnt - prev<< endl;
				prev = cnt;
			}
		}
		else
		{
			prev = 0;
			int_length = (max_range - min_range)/INTERVALS;
			for(it = 0; it <= INTERVALS ; it++)
			{
				comp = min_range + it*int_length;
				cnt =  count_if(snapshot_final.begin(),snapshot_final.end(),pred);
				fprintf(sffp,"%1.3e %d\n",comp,cnt-prev);
				cout << "final Interval " <<it << " " << cnt - prev<< endl;
				prev = cnt;
			}
		}
		fclose(sifp);
		fclose(sffp);
#endif
		if (val_num > 0)
			val_avg = val_sum / val_num;
		else val_avg = -1;
		if (gen_num > 0)
			gen_avg = gen_sum / gen_num;
		else gen_avg = -1;
		fprintf(ofp, "%d  %lf  %lf\n", num_threads, val_avg, gen_avg);

		/* HLEE-Nov-9:
		   each_gen_val_sum and mono_best_val_sum have the values added over
		   all threads over all runs per gen, thus
		   divide it by (num_threads * RUN)

		   HLEE-Apr-18-2009:
		   for Byzantine failure simulation, the values of non-faulty proc.s
		   are added, thus divide it by ((num_threads - byz_threshold) * R
		 */

//		z = (num_threads - byz_threshold) * RUN;
		for (g = 1; g <= max_gen; ++g) {
			each_gen_best_avg[t][g] = each_gen_val_sum[g] / threadcount;

			mono_best_vals[t][g] = mono_best_val_sum[g] / RUN;
		}
cout << "ThreadCount = " << threadcount << endl;

		for (g = 1; g <= max_gen; ++g) each_gen_val_sum[g] = 0.0;

		for (g = 1; g <= max_gen; ++g) mono_best_val_sum[g] = 0.0;

		num_threads *= 2;
	}/* end for1 varying threads */
cout << "out of loops "  << endl;
	for (g = 1; g <= max_gen; ++g) {
		fprintf(ofp1, "%d", g);
		for (t = 0; t < VARY_NUM_THREADS; ++t)
			fprintf(ofp1, " %lf", mono_best_vals[t][g]);
		fprintf(ofp1, "\n");
	}
cout << "finished writing mono over " << endl;
	for (g = 1; g <= max_gen; ++g) {
		fprintf(ofp2, "%d", g);
		for (t = 0; t < VARY_NUM_THREADS; ++t)
			fprintf(ofp2, " %lf", each_gen_best_avg[t][g]);
		fprintf(ofp2, "\n");
	}

cout << "finished writing avg over " << endl;
	fclose(ofp);
	fclose(ofp1);
	fclose(ofp2);
	return num_files;
}
