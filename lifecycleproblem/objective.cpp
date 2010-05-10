/*customer group and optimum values seen
 * Neutral  	nlc = 4 fitness = 0.83
 * Utilitarian	nlc = 3 fitness = 0.82
 * Technophile 	nlc = 8 fitness = 0.65
 * Greens	nlc = 2 fitness = 0.87

*/
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "pool.h"
//#include <vector>
using namespace std;

unsigned int planhorizon = 10;
//the operation design impact on the age of each component
float agecoefficient[] = {0,1,0.5,0.1};

//get customer attribute interval and scaling constant

//Neutral
//float customer_info[] = {130, 10680, 1050 ,15790, 0.64, 0.95, 0.5, 0.5, 0.5};

//Utilitarian
//float customer_info[] = {  500,  6000, 2100,  9000, 0.65, 0.85, 0.75, 0.35, 0.45 };

//Technophile
float customer_info[] = { 6000, 10000, 4200, 12800, 0.84 ,   1 , 0.3,  0.1,  0.8 };

//Greens
//float customer_info[] = {  500 ,10000, 1050,  7000,  0.5,  0.8, 0.15 ,0.85, 0.15 };


//get cost and environment data
float costdatamatrix[][4] = { {156.2539,  26.8318,     93.00390367, 104.851643},
                  	 {37.1714,   5.144327856, 28.0465226,  29.84743834},
                  {19.979,    2.371710153, 13.6944418,  14.90640359},
                  {88.2269,   8.270501073, 58.68280845, 64.30979485},
                  {27.0815,   1.787185936, 16.79436684, 18.72745049},
                  {260.9658,  4.497161839, 139.4638802, 161.8898537},
                  {45.0689,   5.111606414, 30.18686851, 33.04458381},
                  {35.519,    0.556967249, 20.50344498, 23.27418079},
                  {113.519,   1.780071655, 65.52916947, 74.38446264},
                  {42.519,    0.666733029, 24.54421513, 27.86100095},
                  {15.2294,   2.092057037, 10.83673327, 11.69432708},
                  {23.4708,   3.7094,      17.1605,     18.4097} };

float envdatamatrix[][4] = {286.945,    76.891,  82.718,  161.297,
                 81.545,     52.800,  53.812,  61.107,
                 70.058,     45.786,  46.148,  55.068,
                 85.165,     49.308,  51.169,  61.391,
                 74.496,     49.308,  49.678,  56.563,
                 211.625,    28.176,  29.107,  54.214,
                 231.034,    27.138,  45.509,  109.013,
                 39.144,     28.176,  28.183,  31.559,
                 39.144,     28.176,  28.183,  31.559,
                 39.144,     28.176,  28.183,  31.559,
                 83.477,     2.422,   3.942 ,  11.079,
                 36.383,     4.523,   8.992 ,  20.917};

  

float  GAThread::objectivenew(unsigned* px2, float* x3)
{
	unsigned int i,ind, i1,i2;
	int num_lifecycle;
	int operationdecision;
	float costvector[NUM_VAR_PER_LIFECYCLE];
	float envvector[NUM_VAR_PER_LIFECYCLE];
	unsigned int ageinput[NUM_VAR_PER_LIFECYCLE];
	float agestart[NUM_VAR_PER_LIFECYCLE];
	float usagetime[MAX_LIFECYCLE];
	float sumusagetime;
	float ageoutput[NUM_VAR_PER_LIFECYCLE];
	float product_cost[MAX_LIFECYCLE];
	float product_env[MAX_LIFECYCLE];
	float product_rel[MAX_LIFECYCLE];
	unsigned int mask = 1;
	float SumCost, SumEnv, MinRel;
	double tb[NUM_VAR_PER_LIFECYCLE];
	float v[3];
	float total_violation;
    	float RC = 0, P0 = 0, P1 = 0, P2 = 0;
	double theta1 = 66.68;
	double theta2 = 25.95;
	double theta3 = 130.01;
	double theta4 = 35.54;
	double theta5 = 56.64;
	double theta6 = 53.58;
	double theta7 = 52.27;
	double theta8 = 50.49;
	double theta9 = 50.49;
	double theta10 = 50.49;
	double theta11 = 100;
	double theta12 = 128.49;

	double b1 = 1;
	double b2 = 1;
	double b3 = 1;
	double b4 = 1;
	double b5 = 1;
	double b6 = 1;
	double b7 = 1;
	double b8 = 1;
	double b9 = 1;
	double b10 = 1;
	double b11 = 1;
	double b12 = 1;
	

	SumCost = 0; 
	SumEnv = 0;
	MinRel = 0;//1e100;

	num_lifecycle = (px2[0] & mask) + (px2[0] & (mask << 1)) + (px2[0] & (mask << 2)) + 1;
//cout << "num_lifecycle = " << num_lifecycle << '\n';	

	
	sumusagetime = 0;
	for(i = 0; i < num_lifecycle ; i++)
	{
		sumusagetime += x3[i];
		product_rel[i] = 1;
	}
	for(i = 0; i < num_lifecycle ; i++)
	{
		usagetime[i] = x3[i]*planhorizon/sumusagetime;
//cout << "usagetime "<< usagetime[i] << '\n';
	}
	for(i = 0; i < NUM_VAR_PER_LIFECYCLE ; i++)
	{
		ageoutput[i] = 0;
	}
for(int stage = 0; stage < num_lifecycle ; stage++)
{
	product_cost[stage] = 100;
	product_env[stage] = 0;
	for(i = 0; i < NUM_VAR_PER_LIFECYCLE; i++)
	{
		//    % get the operation decision for each component from chromosome x2
		i1  = (stage)*BITS_PER_LIFECYCLE + 2*i + LC_ENCODE_LENGTH ;		//The overall index into the chromosome 
		//note that the + 1 at the end is make the 2bit operation decisions end at word boundary
		//essentially we are letting the first 4 bits represent the lifecycle length but ignore the 4th bit

		ind = i1 / (8*sizeof(unsigned));	//The index to use for the unsigned array
		i2  = i1 - ind*8*sizeof(unsigned);//The offset to use for the element of the unsigned array
		operationdecision = (((px2[ind])>>(i2+1)) & mask )*2 + (((px2[ind]) >> (i2)) & mask);
		if(stage == 0)
			operationdecision = 0;
    		//% Based on each component operation decison, we can get the cost and
    		//environment impact
//		cout << ind << " " << i1 << " "  << i2 << " " << operationdecision << '\n';
    		costvector[i] = costdatamatrix[i][operationdecision];
    		envvector[i]  = envdatamatrix[i][operationdecision];
    
    		//Get the component age information based on the operation decision and age coefficient 
    		agestart[i]  = agecoefficient[operationdecision] * ageoutput[i];
		ageoutput[i] = agestart[i] + usagetime[stage];
		product_cost[stage] += costvector[i];
    		product_env[stage]  += envvector[i];
//cout << costvector[i] << '\n';
	}
//cout << stage << " pc = " << product_cost[stage] << " pe = " << product_env[stage] << '\n';
	//% thetaX and bX are the parameters to calculate the component reliability
	tb[0] = pow(M_E,-(pow((ageoutput[0]/theta1),b1)));
      	tb[1] = pow(M_E,-(pow((ageoutput[1]/theta2),b2)));
      	tb[2] = pow(M_E,-(pow((ageoutput[2]/theta3),b3)));
        tb[3] = pow(M_E,-(pow((ageoutput[3]/theta4),b4)));
      	tb[4] = pow(M_E,-(pow((ageoutput[4]/theta5),b5)));
      	tb[5] = pow(M_E,-(pow((ageoutput[5]/theta6),b6)));
      	tb[6] = pow(M_E,-(pow((ageoutput[6]/theta7),b7)));
      	tb[7] = pow(M_E,-(pow((ageoutput[7]/theta8),b8)));
      	tb[8] = pow(M_E,-(pow((ageoutput[8]/theta9),b9)));
      	tb[9] = pow(M_E,-(pow((ageoutput[9]/theta10),b10)));
      	tb[10] =pow(M_E,-(pow((ageoutput[10]/theta11),b11)));
      	tb[11] =pow(M_E,-(pow((ageoutput[11]/theta12),b12)));

	float R10 = tb[3];       //%Hard Drive
	float R11 = tb[5];       //%Motherboard
	float R12 = tb[8];       //%Video Card

	float R1 = tb[0];
	float R2 = tb[1];
	float R3 = tb[2];
	float R4 = tb[4];
	float R5 = tb[6];
	float R6 = tb[7];
	float R7 = tb[9];
	float R8 = tb[10];
	float R9 = tb[11];

	P0 = R1*R2*R3*R4*R5*R6*R7*R8*R9;

	P1 = (1-R1)*(R2*R3*R4*R5*R6*R7*R8*R9) + (1-R2)*(R1*R3*R4*R5*R6*R7*R8*R9)+ \
     		(1-R3)*(R1*R2*R4*R5*R6*R7*R8*R9) + (1-R4)*(R1*R2*R3*R5*R6*R7*R8*R9)+ \
     		(1-R5)*(R1*R2*R3*R4*R6*R7*R8*R9) + (1-R6)*(R1*R2*R3*R4*R5*R7*R8*R9)+ \
     		(1-R7)*(R1*R2*R3*R4*R5*R6*R8*R9) + (1-R8)*(R1*R2*R3*R4*R5*R6*R7*R9)+ \
     		(1-R9)*(R1*R2*R3*R4*R5*R6*R7*R8);

	P2 = (1-R1)*(1-R2)*(R3*R4*R5*R6*R7*R8*R9) + (1-R1)*(1-R3)*(R2*R4*R5*R6*R7*R8*R9)+ \
	     (1-R1)*(1-R4)*(R2*R3*R5*R6*R7*R8*R9) + (1-R1)*(1-R5)*(R2*R3*R4*R6*R7*R8*R9)+ \
	     (1-R1)*(1-R6)*(R2*R3*R4*R5*R7*R8*R9) + (1-R1)*(1-R7)*(R2*R3*R4*R5*R6*R8*R9)+ \
	     (1-R1)*(1-R8)*(R2*R3*R4*R5*R6*R7*R9) + (1-R1)*(1-R9)*(R2*R3*R4*R5*R6*R7*R8)+ \
	     (1-R2)*(1-R3)*(R1*R4*R5*R6*R7*R8*R9) + (1-R2)*(1-R4)*(R1*R3*R5*R6*R7*R8*R9)+ \
	     (1-R2)*(1-R5)*(R1*R3*R4*R6*R7*R8*R9) + (1-R2)*(1-R6)*(R1*R3*R4*R5*R7*R8*R9)+ \
	     (1-R2)*(1-R7)*(R1*R3*R4*R5*R6*R8*R9) + (1-R2)*(1-R8)*(R1*R3*R4*R5*R6*R7*R9)+ \
	     (1-R2)*(1-R9)*(R1*R3*R4*R5*R6*R7*R8) + (1-R3)*(1-R4)*(R1*R2*R5*R6*R7*R8*R9)+ \
	     (1-R3)*(1-R5)*(R1*R2*R4*R6*R7*R8*R9) + (1-R3)*(1-R6)*(R1*R2*R4*R5*R7*R8*R9)+ \
	     (1-R3)*(1-R7)*(R1*R2*R4*R5*R6*R8*R9) + (1-R3)*(1-R8)*(R1*R2*R4*R5*R6*R7*R9)+ \
	     (1-R3)*(1-R9)*(R1*R2*R4*R5*R6*R7*R8) + (1-R4)*(1-R5)*(R1*R2*R3*R6*R7*R8*R9)+ \
	     (1-R4)*(1-R6)*(R1*R2*R3*R5*R7*R8*R9) + (1-R4)*(1-R7)*(R1*R2*R3*R5*R6*R8*R9)+ \
	     (1-R4)*(1-R8)*(R1*R2*R3*R5*R6*R7*R9) + (1-R4)*(1-R9)*(R1*R2*R3*R5*R6*R7*R8)+ \
	     (1-R5)*(1-R6)*(R1*R2*R3*R4*R7*R8*R9) + (1-R5)*(1-R7)*(R1*R2*R3*R4*R6*R8*R9)+ \
	     (1-R5)*(1-R8)*(R1*R2*R3*R4*R6*R7*R9) + (1-R5)*(1-R9)*(R1*R2*R3*R4*R6*R7*R8)+ \
	     (1-R6)*(1-R7)*(R1*R2*R3*R4*R5*R8*R9) + (1-R6)*(1-R8)*(R1*R2*R3*R4*R5*R7*R9)+ \
	     (1-R6)*(1-R9)*(R1*R2*R3*R4*R5*R7*R8) + (1-R7)*(1-R8)*(R1*R2*R3*R4*R5*R6*R9)+ \
	     (1-R7)*(1-R9)*(R1*R2*R3*R4*R5*R6*R8) + (1-R8)*(1-R9)*(R1*R2*R3*R4*R5*R6*R7);
	RC = R10*R11*R12;
    
	
//% The above is to calculate the product reliability
    	product_rel[stage]  = RC*(P0 + P1 + P2);


}

MinRel = 1e100;
for(i =0; i < num_lifecycle ; i++)
{
	SumCost += product_cost[i];
	SumEnv  += product_env[i];
	if(MinRel > product_rel[i]  )
		MinRel = product_rel[i];
}
//cout << "SC =" << SumCost << " SE = " << SumEnv << " M= "<< MinRel << '\n';
//% finding constraint violation that makes the solution very bad if any of
//% the constraints is violated.
v[0] = -(SumCost - customer_info[1])/(customer_info[1] - customer_info[0]);
v[1] = -(SumEnv - customer_info[3])/(customer_info[3] - customer_info[2]);
v[2] = -(customer_info[4]- MinRel)/(customer_info[5] - customer_info[4]);


//Total violation added to the fitness to punish the bad solutions and reduce their fitness
total_violation = 0;
for(i = 0; i < 3;i++)
{
//cout << "v " << i << " " << v[i] << '\n';
	if(v[i] < 0) 
		total_violation += v[i];
}
total_violation *= 10000;//suggested by yuan
//cout <<"total voilation " << total_violation << '\n';
//% get customer attribute scaling constant
float k1 = customer_info[6];
float k2 = customer_info[7];
float k3 = customer_info[8];
float a = k1*k2*k3;
float b = k1*k2 + k2*k3 + k3*k1;
float c = k1 + k2 + k3 -1;
float K = (-b+sqrt(b*b - 4*a*c))/(2*a);
float u1 = (customer_info[1] - SumCost)/(customer_info[1] - customer_info[0]);
float u2 = (customer_info[3] - SumEnv)/(customer_info[3] - customer_info[2]);
float u3 = (MinRel - customer_info[4])/(customer_info[5] - customer_info[4]);


if (u1 < 0 || u2 < 0 || u3 < 0 ) 
	return 0;
if(u1 > 1) u1 = 1;
if(u2 > 1) u2 = 1;
if(u3 > 1) u3 = 1;

// fitness - total_violation  = utility.
//
float fitness = ((k1*K*u1 + 1)*(k2*K*u2 + 1)*(k3*K*u3 + 1)-1)/K ;//+ total_violation;  	
//cout << "rem fitness = " << fitness - total_violation << '\n';
//cout <<"fitness = " << fitness << '\n';
return fitness;
}

#if 0 
int main()
{
//	float x3[MAX_LIFECYCLE] = {1.39,1.30,1.28,1.23,1.16,1.22,1.15,1.28};
	float x3[MAX_LIFECYCLE] = {2.28,1.92,1.87,1.90,2.03,0,0,0};
	float fit;
	int i = 0;
	
//	for(int i = 0; i < MAX_LIFECYCLE ; i++)	
//		x3[i] = 0.5;//rand()/RAND_MAX;
	
	unsigned int x22[7];
	x22[0] = 0x50000000;//0xD0000007;
	x22[1] = 0xefdd49dd;//0xEF068620;
	x22[2] = 0xb3cd57f9;//0x01E27744;
        x22[3] = 0x03b9D2e5;//0xBB2DC1DB;
	x22[4] = 0;//0xB33DF4EE; 
	x22[5] = 0;//0xB0315679; 
	x22[6] = 0;//0x0000000A;
/*	
	for(i = 0; i < LENGTH*NUM_LIFECYCLE; i++)
	{
		if(!(i%2))
			x22 |= (1 << i);
	}
*/	 	
	fit= objectivenew(x22,x3);	
	cout << "fitness " << fit << '\n';	
	return 0;
}
#endif
