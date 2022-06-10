/* This code was written by Loïc Marrec */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Init(int *X1, int X10, double *t, double *tlist, double deltat, double tmax, int *ind)
{
	int i;
	double q = 0;
	
	*X1 = X10;
	*t = 0;
	*ind = 0;
	
	for (i = 0; i <= tmax/deltat; i++)
	{
		tlist[i] = q;
		q = q+deltat;
	}
}

double DivisionRate(double F, int X1, int K)
{
	return (double)F*log((double)K/X1)*X1;
}

double DeathRate(double G, int X)
{
	return G*X;
}

void TotalTransitionRate(double FA, double GA, int XA, int K, double *array)
{
	array[0] = DivisionRate(FA, XA, K);
	array[1] = DeathRate(GA, XA);
}

void CumSum(double *array, int index) 
{
    if(index <= 0) return;
    CumSum(array, index-1);
    array[index] += array[index-1];
}

double rnd()
{
    return (double)rand()/(double)RAND_MAX;
}

int SamplingTowerLinear(double *array, int index) 
{
	int i = 0;
	double R1, R2;
	
	R1 = rnd();
	R2 = R1*array[index-1];

	while (array[i] < R2)
	{
		i++;
	}

	return i;
}

double UpdateTime(double *t, double T)
{
	double r, x;
	r = rnd();
	x = 1/r;
	*t += 1/T*log(x);
}

void Reaction(int *X1, int ir1)
{
	if (ir1 == 0)
	{
		*X1 += 1;
	}
	else
	{
		*X1 -= 1;
	}
}

int main()
{
	/* DEFINE YOUR PARAMETERS BELOW */
	const double FA = 1.0, GA = 0; /* FA: division rate; GA: death rate */
	const int K = 100000, XA0 = 100, Nit = 10000; /* K: carrying capacity; XA0: initial population size; Nit: number of stochastic realizations */
	double t = 0, Tvect[2], tlist[301], deltat = .1, tmax = 30; /* deltat: time step; tmax: maximum time */
	int XA, k, i, j, ir1, ind;
	auto XAlist = new int[301][10000];

	for (k = 0; k < Nit; k++)
	{
		Init(&XA, XA0, &t, tlist, deltat, tmax, &ind);
		XAlist[ind][k] = XA;
		ind = ind+1;

		while (XA != 0 && t < tlist[300])
		{
			TotalTransitionRate(FA, GA, XA, K, Tvect);
			CumSum(Tvect, 1);
			UpdateTime(&t, Tvect[1]);
			while (t > tlist[ind] && ind < tmax/deltat)
			{
		            XAlist[ind][k] = XA;
		            ind = ind+1;
			}
			ir1 = SamplingTowerLinear(Tvect, 2); 
			Reaction(&XA, ir1);
		}
		
		XAlist[ind][k] = XA;
	}

	FILE *fp1;
	fp1 = fopen("Gom.txt", "w");
	for (i = 0; i < tmax/deltat+1; i++)
	{
		for (j = 0; j < Nit; j++)
		{
			fprintf(fp1, "%d", XAlist[i][j]);
			fprintf(fp1, " ");
		}
		fprintf(fp1, "\n");
	}
	fclose(fp1); 

	return 0;
} 
