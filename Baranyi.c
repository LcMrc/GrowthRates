/* This code was written by Loïc Marrec */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Init(double *X1, double X10, double *t, double *tlist, double deltat, double tmax, int *ind)
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
	return F*(1-(double)X1/K)*X1;
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
	const double FA = 1.0, GA = 0, h0 = 2, K = 100000, XA0 = 100; /* FA: division rate; GA: death rate; h0: lag phase; K: carrying capacity; XA0: initial population size */
	const int Nit = 10000; /* Number of stochastic realizations */
	double t = 0, S[2], T[2], Tvect[2], Tvectbis[2], tlist[301], deltat = .1, tmax = 30, Dt[2], XA; /* deltat: time step; tmax: maximum time; tlist[tmax/deltat+1] */
	int k, i, j, ir1, ind;
	auto XAlist = new double[301][10000]; /* double[tmax/deltat+1][Nit] */

	for (k = 0; k < Nit; k++)
	{
		Init(&XA, XA0, &t, tlist, deltat, tmax, &ind);
		XAlist[ind][k] = XA;
		ind = ind+1;
		S[0] = log(1/rnd());
		S[1] = log(1/rnd());
		T[0] = 0;
		T[1] = 0;

		while (XA > .5 && t < tlist[300]) /* tlist[tmax/deltat] */
		{
			Dt[0] = -t+log(1-exp(h0)+exp((-S[0]+T[0])*K/XA/(XA-K))*(-1+exp(h0)+exp(FA*t)))/FA;
			Dt[1] = (S[1]-T[1])/GA/XA;

			if (Dt[0] < Dt[1])
			{
				T[0] += XA*(XA-K)*log((1-exp(h0)-exp(FA*t))/(1-exp(h0)-exp(FA*(t+Dt[0]))))/K;
				T[1] += GA*XA*Dt[0];
				t += Dt[0];
				while (t > tlist[ind] && ind < tmax/deltat)
				{
				    XAlist[ind][k] = XA;
				    ind = ind+1;
				}
				XA += 1;
				S[0] += log(1/rnd());
			}
			else
			{
				T[0] += XA*(XA-K)*log((1-exp(h0)-exp(FA*t))/(1-exp(h0)-exp(FA*(t+Dt[1]))))/K;
				T[1] += GA*XA*Dt[1];
				t += Dt[1];
				while (t > tlist[ind] && ind < tmax/deltat)
				{
				    XAlist[ind][k] = XA;
				    ind = ind+1;
				}
				XA -= 1;
				S[1] += log(1/rnd());
			}
			TotalTransitionRate(FA, GA, XA, K, Tvect);
		}
		
		XAlist[ind][k] = XA;
	}

	FILE *fp1;
	fp1 = fopen("Bar.txt", "w");
	for (i = 0; i < tmax/deltat+1; i++)
	{
		for (j = 0; j < Nit; j++)
		{
			fprintf(fp1, "%f", XAlist[i][j]);
			fprintf(fp1, " ");
		}
		fprintf(fp1, "\n");
	}
	fclose(fp1);

	return 0;
} 
