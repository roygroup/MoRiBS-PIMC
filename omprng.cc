#include <omp.h>
#include "rngstream.h"
#include <iostream>
#include "sys/time.h"
#include <cmath>
#include <stdio.h>

using namespace std;


const double pi = std::acos(-1.0);


void fixedseed() {
	//you can change seed changing one or more of the integers below
	unsigned long mySeed[6] = {12345, 12345, 12345, 12345, 12345, 12345};
	RngStream::SetPackageSeed(mySeed);
}


void randomseed() {
	timeval tim;
	gettimeofday(&tim, NULL);
	unsigned long seed[6] = {tim.tv_sec, tim.tv_usec, 2, 3, 4, 5};
	RngStream::SetPackageSeed(seed);
}


void noparallel() {
	omp_set_num_threads(1);
}


double runif(RngStream* myRng) 
{
	//generate random number between 0 and 1, i.e. from Unif(0,1)

	return(myRng[omp_get_thread_num()].RandU01());
}


double runifab(RngStream* myRng, double a, double b) 
{
	// generate random number from X~Unif(a,b)
	// f(x) = 1/(b-a)
	// a < x < b, -infty < a < b < infty
	// E(X) = (a+b)/2, Var(X) = (b-a)^2/12

	return(a + (b-a)*runif(myRng));
}


double rnorm(RngStream* myRng, double mu, double sigma) 
{
	// generate random number from X~N(mu,sigma)
	// f(x) = 1/sqrt(2*pi*si^2) * exp(-0.5*((x-mu)/si)^2)
	// -infty < x < infty, -infty < mu < infty, sigma > 0
	// E(X) = mu, Var(X) = sigma^2
	
	double u, v, z, x;

	u = runif(myRng);
	v = runif(myRng);

	z = sqrt(-2*log(u)) * cos(2*pi*v);
	x = mu + z*sigma;

	return(x);
}


double rexp(RngStream* myRng, double theta) 
{
	// generate random number from X~exp(theta)
	// f(x) = 1 /theta * exp(-x/theta)
	// x > 0, theta > 0
	// E(X) = theta, Var(X) = theta^2

	double u, x;

	u = runif(myRng);
	x = -log(u)*theta;

	return(x);
}


double rgamma(RngStream* myRng, double alpha, double beta) 
{
	// generate random number from X~gamma(alpha,beta)
	// f(x) = 1/(gamma(alpha)*beta^alpha) * x^(alpha-1) * exp(-x/beta)
	// x > 0, alpha > 0, beta > 0
	// E(X) = alpha*beta, Var(X) = alpha*beta^2

	double x=0.0;
	double delta, v0, v1, v2, v3, psi=1.0, nu=1.0e10;

	for(int i=0; i<floor(alpha); i++)
		x = x + rexp(myRng, 1.0);

	if(alpha > floor(alpha)) { // alpha has fractional part
		delta = alpha - floor(alpha);
		while(nu > pow(psi,delta-1.0)*exp(-psi)) {
			v1 = runif(myRng); v2 = runif(myRng); v3 = runif(myRng);
			v0 = exp(1.0) / (exp(1.0) + delta);
			if(v1 <= v0) {
				psi = pow(v2, 1.0/delta);
				nu = v3 * pow(psi,delta-1.0);
			}
			else {
				psi = 1.0 - log(v2);
				nu = v3 * exp(-psi);		
			}
		}
		x = x + psi;
	}

	return(beta*x);
}


double rchisq(RngStream* myRng, int df) 
{
	//generate random number from X~chisq(df)
	// f(x) = 1/(gamma(df/2)*2^(df/2)) * x^(df/2-1) * exp(-x/2)
	// x > 0, df = 1,2,3,...
	// E(X) = df, Var(X) = 2*df

	return(rgamma(myRng, double(df)/2.0, 2.0));
}


double rbeta(RngStream* myRng, double alpha, double beta)
{
	//generate random number from X~beta(alpha,beta)
	// f(x) = gamma(alpha+beta)/(gamma(alpha)*gamma(beta)) * x^(alpha-1) * (1-x)^(beta-1)
	// 0 < x < 1, alpha > 0, beta > 0
	// E(X) = alpha/(alpha+beta), Var(X) = alpha*beta / ((alpha+beta+1)*(alpha+beta)^2)

	double x, y;
	
	x = rgamma(myRng, alpha, 1.0);	
	y = rgamma(myRng, beta, 1.0);

	return(x/(x+y));
}

