/************************************************************************************/
/* Copyright (C) 2006  Sebastian E. Ramos-Onsins                                    */
/*                                                                                  */ 
/* This program is free software; you can redistribute it and/or                    */
/* modify it under the terms of the GNU General Public License                      */
/* as published by the Free Software Foundation; either version 2                   */
/* of the License, or (at your option) any later version.                           */
/*                                                                                  */
/* This program is distributed in the hope that it will be useful,                  */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of                   */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    */
/* GNU General Public License for more details.                                     */
/*                                                                                  */
/* You should have received a copy of the GNU General Public License                */
/* along with this program; if not, write to the Free Software                      */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.  */
/************************************************************************************/

#include <math.h>
#define PI 3.14159265359

#include "mlsp_sm.h"
#include <stdio.h>
#include <stdlib.h>

/*Routines based on Numerical Recipes in C*/

double gammalogn(double zz)
{
	/*Based on Numerical Recipes in C, Press et al. 1992. p. 213. and on 
	Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*gamma distribution for a z integer*/
	double loggammaz;
	double z,logg,h,sumc;
	static double gamma = 5.0;
	static double c0 =  1.000000000178;
	static double c1 = 76.180091729406;
	static double c2 = 86.505320327112;
	static double c3 = 24.014098222230;
	static double c4 =  1.231739516140;
	static double c5 =  0.001208580030;
	static double c6 =  0.000005363820;
	
	if(zz <= 0.) {
		puts("Error gamma");
		return (double)-10000.;
	}
	
	z = (double)zz;
	h = sqrt(2. * PI);
	sumc = c0 + c1/(z+1.) - c2/(z+2.) + c3/(z+3.)  - c4/(z+4.) + c5/(z+5.) - c6/(z+6.);
	logg = (z + 0.5)*log((double)(z + gamma + 0.5)) - (z + gamma + 0.5);
	loggammaz = log((double)h);
	loggammaz += logg + log((double)sumc);
	loggammaz -= log((double)z);
	
	return loggammaz;
}

double factln(long int x)
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*G(n+1) = n!*/
	/*do the log(n!)*/
	double gammalogn(double);
	static double factlog[120];

	if(x == 0) return 0.;
	
	if(x < 120) {
		if(factlog[x] == 0.0) {
			factlog[x] = gammalogn((double)x+1.0);
			return factlog[x];
		}
		else return factlog[x];
	}
	return (gammalogn((double)x+1.0));
}

double gammadist(double alfa) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992, on
	Cheng anf Feast 1979, Appl. Statist. 28, No. 3, pp. 290-295, and on
	PSeq-Gen v1.1. Nicholas C. Grassly, Jun Adachi and Andrew Rambaut*/
	
	double ran1(void);
	double rndgamma1 (double);
	double rand0,rand1,rand2;
	double a,b,c,d,f,W;
	
	if(alfa <= 0.0) {
		return (double)-10000.0;
	}
	if(alfa < 1.0) {
		a = rndgamma1(alfa);
		return a;
	}
	if(alfa == 1.0) return (-log(ran1()));
	
	a = alfa - 1.0;
	b = (alfa - 1.0/(6.0*alfa))/a;
	c = 2.0/a;
	d = c + 2.0;
	f = sqrt((double)alfa);
		
	do {
		if(alfa < 3.0) {
			rand1 = ran1();
			rand2 = ran1();
		}
		else {
			do {
				rand0 = ran1();
				rand1 = ran1();
				rand2 = rand1 + 1.0/f * (1.0-1.86*rand0);
			}while(rand2 < 0.0 || rand2 > 1.0);
		}
		W = b * rand1/rand2;
		if(c*rand2-d+W+1.0/W <= 0.0) break;
	}while(c*log(rand2)-log(W)+W-1.0 >= 0.0);
	
	return a*W;
}

/*From PSeq-Gen v1.1. Nicholas C. Grassly, Jun Adachi and Andrew Rambaut */
double rndgamma1 (double s)
{

	double ran1(void);
	double			r, x=0.0, small=1e-37, w;
	static double	a, p, uf, ss=10.0, d;
	
	if (s!=ss) {
		a  = 1.0-s;
		p  = a/(a+s*exp(-a));
		uf = p*pow(small/a,s);
		d  = a*log(a);
		ss = s;
	}
	for (;;) {
		r = (double)ran1();
		if (r > p) {
			x = a-log((1.0-r)/(1.0-p));
			w=a*log(x)-d;
		}
		else {
			if (r>uf) {
				x = a*pow(r/p,1/s);
				w=x;
			}
			else return (0.0);
		}
		r = (double)ran1();
		if (1.0-r <= w && r > 0.0)
			if (r*(w+1.0) >= 1.0 || -log(r) <= w)  
				continue;
		break;
	}
	return ((double)x);
}

double poissondist(double lambda) 
{
	/*Based on Atkinson 1979 Appl. Statist. 28: No. 1, pp, 29-35.*/
	double ran1(void);	
	double r,s;
	int N;
	
	double factln(long int);
	double alfa,beta,k,X;
	double rand1,rand2;
	static double c = 0.6;
	
	if(lambda < 0.0) {
		puts("Error poissondist");
		return (double)-10000.;
	}
	
	if(lambda == 0.0) return 0.0;
	if(lambda <= 20.0) {
		/*included for having not biased small values with mhits=0...*/
		r = exp(-lambda);
        N = 0;
        s = 1.0;
        do {
            s *= ran1();
            if(s >= r) N += 1;
            else break;
        }while(1);
	}
	else {
		beta = PI * 1.0/sqrt(3.0*lambda);
		alfa = beta * lambda;
		k = log(c) - lambda - log(beta);
		do{
			rand1 = ran1();
			X = (alfa-log((1.0-rand1)/rand1))/beta;
			if(X >= -0.5) {
				N = (int)(X + 0.5);
				rand2 = ran1();
				if(alfa - beta*X +
				  log((double)(rand2/((1.0+exp((alfa-beta*X)))*(1.0+exp((alfa-beta*X)))))) 
				  <= k + (double)N*log(lambda) - factln((long int)N)) 
					break;
			}
		}while(1);
	}
	return (double)N;
}

double binomialdist(double pp, int n) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	Fishman 1979 J. American Statistical Association Vol. 74, No. 366, pp 418-423*/
	
	double ran1(void);	
	double p,np;
	int N;
	int nn;
	double r;
	double A,B,C,D,V,s;
	double m,mu;	
	static double *f=0;
	double poissondist(double);
	static int max = 200;
	
	if(f == 0) {
		if((f=(double *)calloc(max+1,sizeof(double))) == 0) {
			printf("Error allocation memory");
			return -10000;
		}
		f[1] = 0.0;
		for(N=1;N<max;N++)
			f[N+1] = f[N] + log((double)N);
	}
	if(n > max) {
		if((f=realloc(f,n*sizeof(double))) == 0) {
			printf("Error allocation memory");
			return -10000;
		}
		for(N=max;N<n;N++)
			f[N+1] = f[N] + log((double)N);
		max = n;
	}
	
	if(pp > 0.5) p = 1.0-pp;
	else p = pp;
	
	np = n * p;
	
	if(n==0) {
		puts("Error bindist");
		return (double)-10000.;
	}
	if(p==0.0) {
		if(pp > 0.5) return (double)n;
		return 0.0;
	}
	
	if(n < 20) {
		/*Bernouilli Method*/
		nn = n;
		N=0;
		while(nn--) 
			if(ran1()<p) N++;		
	}
	else {
		if(np < 10.0) {
			/*Rejection Method: BI Algorithm*/
			s = 1.0- p;
			A = 1.0;
			B = p/s;
			C = ((double)n+1.0)*B;
			D = A;
			N = 0;
			V = ran1()/pow(s,(double)n);
			while(V > A) {
				N++;
				D *= (C/(double)N - B);
				A += D;
				if(N > n) break;
			}
		}
		else {
			/*Poisson method: BP Algorithm*/
			mu = n - (double)floor((double)(n*(1.0 - p)));
			if(n*(1.0-p) - (double)floor((double)(n*(1.0-p))) > p)
				mu = p*((double)floor((double)(n*(1.0-p))) + 1.0) / (1.0-p);
			r = (1.0/p - 1.0) * mu;
			s = log((double)r);
			m = (double)floor((double)(r));
			do {
				do {
					N = (int)poissondist(mu);
				}while((int)N > n);
				V = -log((double)ran1());
			}while(V < (m-(double)(n - N))*s - f[(int)m+1] + f[(int)(n-N)+1]);
		}
	}
	if(pp > 0.5) N = n - N;
	return (double)N;
}

double largebinomialdist(double pp, double n) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	 Fishman 1979 J. American Statistical Association Vol. 74, No. 366, pp 418-423*/
	
	double ran1(void);	
	double p,np;
	int N;
	double A,B,C,D,V,s;
	double g,plog,pclog,sq,angle,y,em,tt;
	double gammln(double);	
	
	if(pp > 0.5) p = 1.0-pp;
	else p = pp;
	
	np = n * p;
	
	if(n==0) {
		puts("Error bindist");
		return -10000.0;
	}
	if(p==0.0) {
		if(pp > 0.5) return n;
		return 0.0;
	}
	
	if(np < 10.0) {
		/*Rejection Method: BI Algorithm*/
		s = 1.0- p;
		A = 1.0;
		B = p/s;
		C = (n+1.0)*B;
		D = A;
		N = 0;
		V = ran1()/pow(s,n);
		while(V > A) {
			N++;
			D *= (C/(double)N - B);
			A += D;
			if(N > n) break;
		}
	}
	else { /*Rejection method with a Lorentzian comparison distribution*/
		g = gammln(n+1.);
		plog  = log(p);
		pclog = log((1.0 - p));
		sq = sqrt(2.0*np*(1.0 - p));
		do {
			do {
				angle = PI*ran1();
				y = tan(angle);
				em=sq*y+np;
			} while(em < 0.0 || em >= (n + 1.0));
			em = floor(em);
			tt = 1.2*sq*(1.0+y*y)*exp(g-gammln(em+1.0) - gammln(n-em+1.0)+em*plog+(n-em)*pclog);
		} while(ran1() > tt);
		N = (int)em;
	}

	if(pp > 0.5) N = (int)n - N;
	return (double)N;
}

int zbracn(double (*func)(double,double,double,double,double,double,double,double,double,double),double *x1,double *x2,double a0,double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8)
{
	/* Based on Numerical Recipes in C. Press et al. 1992. 
    We need *x1 and *x2 be the range where a root is within them. We expand geometrically the range until finding
	(one a positive and one a negative value). If not, return 0.
	*/

    double f1,f2;
    int k=60;
    
    if(*x1 == *x2) return 0;

    f1 = (*func)(*x1,a0,a1,a2,a3,a4,a5,a6,a7,a8);
    f2 = (*func)(*x2,a0,a1,a2,a3,a4,a5,a6,a7,a8);
	
	if(f1*f2 < 0.0) return 1;

    while(k--) {
        if(fabs(f1) < fabs(f2)) {
            *x1 += 1.5 * (*x1 - *x2);
			if(*x1 < 0.) 
				*x1 = 0.;
            f1 = (*func)(*x1,a0,a1,a2,a3,a4,a5,a6,a7,a8);
        }
        else {
            *x2 += 1.5 * (*x2 - *x1);
            f2 = (*func)(*x2,a0,a1,a2,a3,a4,a5,a6,a7,a8);
        }
        if(f1*f2 < 0.0) return 1;
    }
    
	return 0;
}

double zriddrn(double (*func)(double,double,double,double,double,double,double,double,double,double),double xlow,double xhigh,double xacc,double a0,double a1,double a2,double a3,double a4,double a5,double a6,double a7, double a8)
{
	/* Based on Numerical Recipes in C. Press et al. 1992., p. 358 an on
	Ridders, 1979, IEEE Transactions on Circuits and systems, Vol. Cas-26, No. 11, pp. 979-980.
	*/
    int k=60;
	double flow,fhigh;
	double f1,f2,f3,f4;
	double x1,x2,x3,x4;
	double den,num,nsign;
	

    flow  = (*func)(xlow ,a0,a1,a2,a3,a4,a5,a6,a7,a8);
    fhigh = (*func)(xhigh,a0,a1,a2,a3,a4,a5,a6,a7,a8);

	if(flow  == 0.0) return xlow;
	if(fhigh == 0.0) return xhigh;
	if(flow*fhigh > 0.0) 
		return (double)-1e32;
	
	x1 = xlow;
	x2 = xhigh;
	f1 = flow;
	f2 = fhigh;
		
	while(k--) {
		x3 = (x1+x2)/2.0;
		f3 = (*func)(x3,a0,a1,a2,a3,a4,a5,a6,a7,a8);
		if(f1 - f2 < 0.0) nsign = -1.0;
		else nsign = 1.0;
		num = (x3-x1) * f3 * nsign;
		den = sqrt(f3*f3 - f1*f2);
		if(den <= xacc && -den <= xacc) return x3;
		x4 = x3 + num/den;
		f4 = (*func)(x4,a0,a1,a2,a3,a4,a5,a6,a7,a8);
		if(f4 <= xacc && -f4 <= xacc) return x4;
		if(f3*f4<0.0) {
			x1 = x3;
			f1 = f3;
			x2 = x4;
			f2 = f4;
		}
		else {
			if(f1*f4<0.0) {
				x2 = x4;
				f2 = f4;
			}
			else {
				if(f2*f4<0.0) {
					x1 = x4;
					f1 = f4;
				}
			}
		}
		if(fabs(x1-x2) <= xacc) return x1;
	}	
	return (double)-1e32;
}

/*CALCULATE PRIORS AND INCLUDE ON THE STRUCT var_priors */
int make_priord(struct var_priors *priorx,long int niter) 
{
	long int i,j;
	double value;
	double arg[6];
	double gammadist(double);
	double betai(double,double,double);
	double ran1(void);
	void init_seed1(long int);
	
	double val1;
	FILE *read_prior;
	int x,c;
	char number[100];
	
	if((*priorx).prior_file[0] != '\0') {
		if (!(read_prior = fopen ((*priorx).prior_file,"r"))) {
			printf("Error reading the prior file.\n");
			exit(1);
		}
		if(!((*priorx).priordist = (double *)calloc((long unsigned)niter,sizeof(double)))) return 0;
		/*read a header*/
		while((c = getc(read_prior))!= 10 && c!= 13 && c!= 0);
		if(c==0) perror("Error reading the prior file.\n");
		while((c = getc(read_prior))== 10 || c== 13);
		if(c==0) perror("Error reading the prior file.\n");
		/*read numbers*/
		for(i=0;i<niter;i++) {
			x = 0;
			while(c!=0 && c!=10 && c!=13) {
				number[x] = c;
				x++;
				c = getc(read_prior);
				if(x==99) break;
			}
			number[x] = '\0';
			(*priorx).priordist[i] = (double)atof(number);
			if(c==0 && i<niter-1) perror("Error reading the prior file: not enough values\n");
			c = getc(read_prior);
		}
		fclose(read_prior);
	}
	else {
		init_seed1((*priorx).seed_prior);
		for(i=0;i<6;i++) arg[i] = (*priorx).dist_par[i];
		if(!((*priorx).priordist = (double *)calloc((long unsigned)niter,sizeof(double)))) return 0;

		switch((*priorx).kind_dist) {
			case 0:
				/*uniform*/
				if(arg[1] == arg[2]) {
					printf("\nWarning: min and max values are equal when calculating prior distribution.\n");
				}
				break;
			case 1:
				/*log10-uniform*/ /*min, max*/
				if(arg[1] == arg[2]) {
					printf("\nWarning: min and max values are equal when calculating prior distribution.\n");
				}
				if(arg[1] == 0.) arg[1] = 1E-100;
				if(arg[2] == 0.) arg[2] = 1E-100;
				break;
			case 2:
				/*gamma*/
				if(arg[1] == 0.) {
					printf("\nError: alpha value must be > 0");	
					return 0;
				}
				if(arg[2] == 0.) {
					printf("\nError: p value must be > 0");	
					return 0;
				}
				if(arg[3] == 0.) arg[3] = 1.;
				break;
			case 3:
				/*beta*/
				if(arg[1] == 0.) {
					printf("\nError: alpha value must be > 0");	
					return 0;
				}
				if(arg[2] == 0.) {
					printf("\nError: beta value must be > 0");	
					return 0;
				}
				if(arg[3] == 0. && arg[4] == 0.) {
					arg[3] = 0.;
					arg[4] = 1.;
				}
				break;
			default:
				break;
		}
		
		for(i=0;i<niter;i++) {
			switch((*priorx).kind_dist) {
				case 0:
					/*uniform*/ /*min, max*/
					value = (arg[2]-arg[1])*(double)ran1() + arg[1]; /*max-min*ran + min*/
					break;
				case 1:
					/*log10-uniform*/ /*min, max*/
					if(arg[1] < 0 && arg[2] < 0) {
						val1 = -1.; 
						value = pow(10.0,log10((double)arg[2]*val1) + (log10((double)arg[1]*val1) - log10((double)arg[2]*val1)) * ran1());
						value *= val1;
					}
					else {
						val1 = +1.;
						value = pow(10.0,log10((double)arg[1]*val1) + (log10((double)arg[2]*val1) - log10((double)arg[1]*val1)) * ran1());
					}
					break;
				case 2:
					/*gamma*/
					value = gammadist(arg[1])/arg[2] * arg[3]; /*gammadist(alpha)/p_gamma * correct*/
					break;
				case 3:
					/*beta*/ /*alpha, beta, min, max*/
					value = betai(arg[1],arg[2],(double)ran1()); /*distribution function*/
					value = arg[3] + (arg[4] - arg[3]) * value;
					break;
				default:
					return 0;
					break;
			}
			if((*priorx).kind_var == 0) {
				j = (long int)value;
				(*priorx).priordist[i] = (double)j;
			}
			if((*priorx).kind_var == 1) (*priorx).priordist[i] = value;
		}	
	}
	return 1;
}

#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define MAXIT 99

double gammln(double zz)
{
	/*Based on Numerical Recipes in C, Press et al. 1992. p. 213. and on 
	 Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*gamma distribution for a z integer*/
	double loggammaz;
	double z,logg,h,sumc;
	static double gamma = 5.0;
	static double c0 =  1.000000000178;
	static double c1 = 76.180091729406;
	static double c2 = 86.505320327112;
	static double c3 = 24.014098222230;
	static double c4 =  1.231739516140;
	static double c5 =  0.001208580030;
	static double c6 =  0.000005363820;
	
	if(zz <= 0.) {
		puts("Error gamma");
		return -10000.0;
	}
	
	z = zz;
	h = sqrt(2. * PI);
	sumc = c0 + c1/(z+1.) - c2/(z+2.) + c3/(z+3.)  - c4/(z+4.) + c5/(z+5.) - c6/(z+6.);
	logg = (z + 0.5)*log((z + gamma + 0.5)) - (z + gamma + 0.5);
	loggammaz = log(h);
	loggammaz += logg + log(sumc);
	loggammaz -= log(z);
	
	return loggammaz;
}

double betacf(double a, double b, double x)
{
    /*
	 used by betai: Evaluates continued fraction for incomplete beta function by modified
	 Lentz's method. numerical recipes in C. 2nd ed. p. 227.
	 */
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;
    d = 1.0 - qab*x/qap;
    if(fabs(d) < FPMIN) d = (double)FPMIN;
    d = 1.0/d;
    h = d;
    for(m = 1;m <= MAXIT; m++) {
        m2 = 2*m;
        aa = m*(b-m)*x/((qam+m2)*(a+m2));
        d = 1.0 + aa*d;
        if(fabs(d) < FPMIN) d = (double)FPMIN;
        c = 1.0 + aa/c;
        if(fabs(c) < FPMIN) c = (double)FPMIN;
        d = 1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d = 1.0 + aa*d;
        if(fabs(d) < FPMIN) d = (double)FPMIN;
        c = 1.0 + aa/c;
        if(fabs(c) < FPMIN) c = (double)FPMIN;
        d = 1.0/d;
        del = d*c;
        h *= del;
        if((double)fabs(del-1.0) < EPS) break;
    }
    if(m > MAXIT)  {
		puts("Error in betacf. MAXIT is too small.");
		return -10000.0;
	}
    return h;
}

double betai(double a, double b,double x)
{
    /*
	 Returns the incomplete beta function Ix(a,b). numerical recipes in C. 2nd ed. p. 227.
	 */
    double betacf(double a,double b,double x);
    double gammln(double xx);
    double bt;
    
    if(x < 0.0 || x > 1.0) {
        puts("Error in betai.");
        return -1.0;
    }
    if(x == 0.0 || x == 1.0) bt = 0.0;
    else bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
    
    if(x < (a + 1.0)/(a + b + 2.0)) return bt*betacf(a,b,x)/a;
    else return 1.0 - bt*betacf(b,a,1.0-x)/b;
}

/*INTEGRATION BIONOMIAL: SEARCH FOR p FOR A GIVEN CUMMULATIVE VALUE*/

#define FUNC(n,k,x) ((*distp_binom)(n,k,x))
#define EPS2 1.0e-6
#define JMAX 20

double distp_binom(int n, int k, double p)
/*p-value for a given n niter,k cases and p probability using a binomial*/
{
	double s1=0.0;
	int n1,nmax;
	
	if(p==0 || p==1) return 0;
	nmax =(k>n-k? k:n-k);
	for(n1=n;n1>nmax;n1--) s1+=log(n1);
	for(n1=1;n1<=n-nmax;n1++) s1-=log(n1);
	return(exp(s1+(double)k*log(p)+(double)(n-k)*log(1.0-p)));
}

double trapzd_inv(double a,double b, int n, int nsam, int k, double *xm, double *psum)
/*nth stage of trapezoidal rule. accuracy adding 2^n-2 interior points*/
{
	double x,tnm,sum,del;
	double s=0.0;
	int it,j;
	
	if(n==1) {
		s=0.5*(b-a)*(FUNC(nsam,k,a)+FUNC(nsam,k,b));
	} else {
		for(it=1,j=1;j<n-1;j++) it <<= 1;/*!!(2^n-2)*/
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for(sum=0.0,j=1;j<=it;j++,x+=del) {
			sum += FUNC(nsam,k,x);
			xm[j-1]=x;
			psum[j-1]=sum;
		}	
		s=0.5*(s+(b-a)*sum/tnm);
		for(j=1;j<=it;j++) {
            psum[j - 1] = psum[j - 1] / sum;
        }
	}
    return s;
}

