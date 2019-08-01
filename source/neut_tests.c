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
#include <stdlib.h>
#include <stdio.h>
#include "zf_log.h"

void init_coef(double *p,int sample_size)
{
/*
 Tajima and Fu coefficients... and Achaz coefficients. Total = 18
*/
	int x,n;
	double an,bn,cn;
		
	an = bn = 0.0;
	n = sample_size;

	for(x=1;x<n;x++) {
            an += 1.0/(double)x;
            bn += 1.0/((double)x*(double)x);
	}
        	
	p[0] = an;
	p[1] = bn;

	/* vt */
	p[2] = (2.0*((double)n*(double)n + (double)n + 3.0)/(9.0*(double)n*((double)n-1.0))
	        -((double)n+2.0)/(an * (double)n) + bn/(an * an)) / (an*an + bn);
	        
	/* ut */
	p[3] = (( ((double)n+1.0)/(3.0*((double)n-1.0)) - 1.0/an)/an) - p[2];
	
	/* vd* */
	p[4] = (bn/(an*an) - 2.0/(double)n *(1.0 + 1.0/an - an + an/(double)n) - 1.0/((double)n*(double)n))
	       / (an*an + bn);
	
	/* ud* */
	p[5] = ((((double)n-1.0)/(double)n - 1.0/an) / an) - p[4];
	
	/* vf* */
	p[6] = (((2.0*(double)n*(double)n*(double)n + 110.0 * (double)n*(double)n - 255.0 * (double)n + 153.0)
	        / (9.0 * (double)n * (double)n * ((double)n-1.0)) + (2.0*((double)n-1.0) *an)/ ((double)n*(double)n) 
	        - (8.0 * bn)/(double)n)) / (an*an + bn);
	
	/* uf* */
	p[7] = (((4.0*(double)n*(double)n + 19.0*(double)n + 3.0 - 12.0 * ((double)n +1.0) * (an + 1.0/(double)n))
	      / (3.0*(double)n *((double)n-1.0))) / an) - p[6];

	/* cn */
	cn = 2.0 * ((double)n*an-2.0*((double)n-1.0)) / (((double)n-1.0)*((double)n-2.0));
	
	/* vd */
	p[8] = 1.0 + (an*an/(bn+an*an)) * (cn - (((double)n+1.0)/((double)n-1.0)));
	
	/* ud* */
	p[9] = an -1.0 - p[8];
	
	/* vf */
	p[10] = (cn + (2.0*((double)n*(double)n+(double)n+3.0))/(9.0*(double)n*((double)n-1.0)) - 2.0
	       /((double)n-1.0)) / (an*an + bn);
	
	/* uf */
	p[11] = (1.0 + ((double)n+1.0)/(3.0*((double)n-1.0)) - 4*((double)n+1.0)
	       /(((double)n-1.0)*((double)n-1.0)) * (an+(1.0/(double)n) - 2.0*(double)n/((double)n+1.0)))
	       / an  -  p[10];
	
	/* Achaz f */
	p[12] = ((double)n - 2.0)/((double)n * (an - 1.0));
	
	/* Achaz f* */
	p[13] = ((double)n - 3.0)/(an * ((double)n - 1.0) - (double)n);
	
	/* alfa_n */
	p[14] = p[12] * p[12] * (an - 1.0) + 
			p[12] * (an * (4.0*((double)n+1.0))/(((double)n-1.0)*((double)n-1.0)) - (2.0 *((double)n + 1.0)*((double)n + 2.0))/((double)n * ((double)n-1.0))) -
			an * (8.0 * ((double)n+1.0)/((double)n*((double)n-1.0)*((double)n-1.0))) +
			((double)n*(double)n*(double)n + (double)n*(double)n + 60.0*(double)n + 12.0)/(3.0*(double)n*(double)n*((double)n-1.0));
	
	/* beta_n */
	p[15] = p[12] * p[12] * (bn + an * (4.0/(((double)n - 1.0) * ((double)n - 2.0))) - 4.0/((double)n-2.0)) +
			p[12] * (-an * (4.0*((double)n+2.0))/((double)n*((double)n-1.0)*((double)n-2.0)) - ((double)n*(double)n*(double)n - 3.0*(double)n*(double)n - 16.0*(double)n + 20.)/((double)n*((double)n-1.0)*((double)n-2.0))) +
			an * (8.0/((double)n*((double)n-1.0)*((double)n-2.0))) +
			2.0 * ((double)n*(double)n*(double)n*(double)n - (double)n*(double)n*(double)n - 17.0*(double)n*(double)n - 42.0*(double)n + 72.0)/(9.0*(double)n*(double)n*((double)n-1.0)*((double)n-2.0));
	
	/* alfa*_n */
	p[16] = p[13] * p[13] * (an - (double)n/((double)n-1.0)) + 
			p[13] * (an * 4.0 * ((double)n+1.0)/(((double)n-1.0)*((double)n-1.0)) - 2.0 * ((double)n+3.0)/((double)n-1.0)) -
			an * (8.0*((double)n+1.0))/((double)n*((double)n-1.0)*((double)n-1.0)) +
			((double)n*(double)n + (double)n + 60.0)/(3.0*(double)n*((double)n-1.0));
	
	/* beta*_n */
	p[17] = p[13] * p[13] * (bn - (2.0*(double)n-1.0)/(((double)n-1.0)*((double)n-1.0))) +
			p[13] * (bn * 8.0/((double)n-1.0) - an * 4.0/((double)n*((double)n-1.0)) - ((double)n*(double)n*(double)n + 12.0*(double)n*(double)n - 35.0*(double)n +18.0)/((double)n*((double)n-1.0)*((double)n-1.0))) -
			bn * 16.0/((double)n*((double)n-1.0)) +
			an * 8.0/((double)n*(double)n*((double)n-1.0)) +
			2.0 * ((double)n*(double)n*(double)n*(double)n + 110.0*(double)n*(double)n - 255.0*(double)n + 126.0)/(9.0*(double)n*(double)n*((double)n-1.0)*((double)n-1.0));
}
/*Tajima's D*/
double tajima_d(double k_, int S_, double *coef_taj)
{
    ZF_LOGD("START");
	double an,ut,vt;
	double S_D = -10000;
        
        if(S_ == 0 || *(coef_taj+0) < 1.51) {
            ZF_LOGD("END");
            return(-10000);
        }
        
	an = *(coef_taj+0);
	ut = *(coef_taj+3);
	vt = *(coef_taj+2);

	S_D = (k_ - ((double)S_/an)) / (sqrt((ut*(double)S_) + (vt*(double)S_*((double)S_))));
	
	if (fabs(S_D) < 1.0E-15)
		S_D = 0.0;

    ZF_LOGD("END");
	return S_D;
}
/*Tajima's D/Dmin*/
double tajima_dvsdmin(double k_, int S_, double *coef_taj,int sample_size)
{
    ZF_LOGD("START");
	double an,kmin;
	double D_Dmin = -10000;
        
        if(S_ == 0 || *(coef_taj+0) < 1.51) {
            ZF_LOGD("END");
            return(-10000);
        }
        
	an = *(coef_taj+0);
        kmin = (double)S_ * ((double)2/(double)sample_size);

	D_Dmin = (k_ - ((double)S_/an)) / (double)fabs(kmin - ((double)S_/an));

    ZF_LOGD("END");
	return D_Dmin;
}
/*Fu and Li's D*/
double fl_d(int sample_size,int fr1,int S, double *coef) /* amb outgroup */
{
    ZF_LOGD("START");
	double an;
	double ud,vd;
	int re,n;
	double D = -10000;
	n = sample_size;
	
	if(S == 0 || *(coef+0) < 1.5) {
        ZF_LOGD("END");
	    return(-10000);
	}
                
	re = fr1;	
	an = *(coef+0);
	
	vd = *(coef+8);
	ud = *(coef+9);
	D  = ((double)S - an*(double)re) / 
	     sqrt(ud*(double)S + vd*(double)S*(double)S);

    ZF_LOGD("END");
	return D;
}
/* Fu and Li's D* */
double fl_d2(int sample_size,int fr1w,int S, double *coef) /* NO outgroup */
{
    ZF_LOGD("START");
	double an;
	int n;
	double ud2,vd2;
	int rs;
	double D2 = -10000;
        
        if(S == 0 || *(coef+0) < 1.51) {
            ZF_LOGD("END");
            return(-10000);
        }
	
	rs = fr1w;

	n = sample_size;
	an = *(coef+0);
	
	vd2 = *(coef+4);
	ud2 = *(coef+5);
	D2  = ((double)S/an - (double)rs*(((double)n-(double)1)/(double)n)) /
	      sqrt(ud2*(double)S + vd2*(double)S*(double)S);

    ZF_LOGD("END");
	return D2;
}
/*Fu and Li's F*/
double fl_f(int sample_size,int fr1, int S, double pi, double *coef) /* amb outgroup */
{
    ZF_LOGD("START");
	double uf,vf;
	int re,n;
	double F;
	n = sample_size;
	
	if(S == 0 || *(coef+0) < 1.5) {
        ZF_LOGD("END");
	    return(-10000);
	}

	re = fr1;		
	vf = *(coef+10);
	uf = *(coef+11);

	F  = (pi - (double)re) / sqrt(uf*(double)S + vf*(double)S*(double)S);

    ZF_LOGD("END");
	return F;
}
/* Fu and Li's F* */
double fl_f2(int sample_size,int fr1w, int S, double pi, double *coef) /* NO outgroup */
{
    ZF_LOGD("START");
	int n;
	double uf2,vf2;
	int rs;
	double F2;
        
        if(S == 0 || *(coef+0) < 1.51) {
            ZF_LOGD("END");
            return(-10000);
        }
	
	rs = fr1w;
	
	n   = sample_size;
	vf2 = *(coef+6);
	uf2 = *(coef+7);
	
	F2  = (pi - ((((double)n-1.0)/(double)n)*(double)rs)) / 
	        sqrt(uf2*(double)S + vf2*(double)S*(double)S);

    ZF_LOGD("END");
	return F2;
}

/*Fay and Wu H*/
double fay_wu(int sample_size,int *fr,double pi) /* nomes outgroup */
{
    ZF_LOGD("START");
    int i;
    double Th,H;
    
    if(sample_size < 2) {
        ZF_LOGD("END");
        return(-10000);
    }
    
    Th = 0.;
    for(i=1;i<sample_size;i++) Th += ((double)*(fr+i))*((double)i*(double)i);
    Th *= (double)2/((double)sample_size*(double)(sample_size-1));
    
    H = pi - Th;

    ZF_LOGD("END");
    return H;
}

/*Fay and Wu H divided by their minimum given S*/
double fay_wuvsminH(int sample_size,int *fr,double pi,int S) /* nomes outgroup */
{
    ZF_LOGD("START");
    int i;
    double Th,Hmin,Thmax,pimin;
    
    if(sample_size < 2) {
        ZF_LOGD("END");
        return(-10000);
    }
    
    Th = 0.;
    for(i=1;i<sample_size;i++) Th += ((double)*(fr+i))*((double)i*(double)i);
    Th *= 2.0/((double)sample_size*(sample_size-1));
    
    Thmax  = (double)S * ((double)(sample_size-1)*(double)(sample_size-1)) * ((double)2/((double)sample_size*(double)(sample_size-1)));
    /*pimin  = (double)S * (2./((double)sample_size*(double)sample_size) * ((double)sample_size-1.);*/
    pimin  = (double)S * ((double)2/((double)sample_size));
    
    if(pimin == Thmax) {
        ZF_LOGD("END");
        return(-10000);
    }
    else Hmin = (pi - Th)/fabs((pimin - Thmax));

    ZF_LOGD("END");
    return Hmin;
}

double fay_wu_normalized(int n,int *fr,double pi) /* Fay and Wu H nomes outgroup NORMALIZED (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    ZF_LOGD("START");
    int i;
    double TL,H,varpiTL,thetaw,an,bn,S;

    if(pi == (double)0.0 || n < 2) {
        ZF_LOGD("END");
        return(-10000);
    }
    
    TL = thetaw = an = bn = (double)0;
	for(i=1;i<n;i++) {
		TL += ((double)*(fr+i))*((double)i);
		thetaw += (double)*(fr+i); 
		an += (double)1/(double)i;
		bn += (double)1.0/((double)i*(double)i);
	}
    TL *= (double)1.0/((double)(n-(double)1));
    S = thetaw;
	thetaw = thetaw/an;
	varpiTL = thetaw * ((double)(n-(double)2))/((double)6*((double)(n-(double)1))) + 
			S*(S-(double)1)/(an*an+bn) / ((double)9*((double)n*(n-(double)1)*(n-(double)1))) * 
			  ((double)18*(double)n*(double)n*((double)3*(double)n+(double)2)*(bn+(double)1.0/((double)n*(double)n)) - 
			  ((double)88*(double)n*(double)n*(double)n + (double)9*(double)n*(double)n - 13*(double)n + (double)6)) ;
	
	H = (pi - TL)/(double)sqrt(varpiTL);

    ZF_LOGD("END");
    return H;
}

double fay_wu_normalized2(int n,double thetaL,double thetaw,double S,double *coef,double pi) /* eq 11 and 12 from Fay and Wu H NORMALIZED (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    ZF_LOGD("START");
    double H,varpiTL,an,bn;
    
    if(pi == (double)0.0 || n < 4) {
        ZF_LOGD("END");
        return(-10000);
    }
    an = coef[0];
	bn = coef[1];

	varpiTL = thetaw * ((double)(n-(double)2))/((double)6*((double)(n-(double)1))) + 
	          S*(S-(double)1)/(an*an+bn) / ((double)9*((double)n*(n-(double)1)*(n-(double)1))) *
				((double)18*(double)n*(double)n*((double)3*(double)n+(double)2)*(bn+(double)1.0/((double)n*(double)n)) - 
				((double)88*(double)n*(double)n*(double)n + (double)9*(double)n*(double)n - 13*(double)n + (double)6));
	
	H = (pi - thetaL)/(double)sqrt(varpiTL);

    ZF_LOGD("END");
    return H;
}


double E_zeng(int n,double thetaL,double thetaw,double S,double *coef) /* (eq 13 and 14 from Zeng et al. Genetics 2006 174: 1431-9)*/
{
    ZF_LOGD("START");
    double E,varLW,an,bn;
    
    if(thetaw == (double)0.0 || n < 4) {
        ZF_LOGD("END");
        return(-10000);
    }
    an = coef[0];
	bn = coef[1];
	varLW = thetaw * ((double)n/((double)2*(double)(n-1)) - (double)1/an) +
			S*(S-(double)1)/(an*an+bn) * 
			(bn/(an*an) + (double)2*bn*((double)n/(double)(n-1))*((double)n/(double)(n-1)) - 
			 (double)2*((double)n*bn-(double)n+(double)1)/((double)(n-1)*an) - 
			 ((double)3*(double)n+(double)1)/((double)(n-1)));
	
	E = (thetaL - thetaw)/(double)sqrt(varLW);

    ZF_LOGD("END");
    return E;
}

double Fst(double piwithin, double pibetween,int ntotpop)
{
    ZF_LOGD("START");
    double fst; /*equation 3 or 6 in Hudson, Slatkin and Maddison 1992*/

    if((pibetween == 0.0 && piwithin == 0.0) || ntotpop < 2){
        ZF_LOGD("END");
        return(-10000.);
    }
	if(pibetween == 0. || pibetween == -10000. || piwithin == -10000.){
        ZF_LOGD("END");
	    return(-10000.);
	}
    fst = (double)1 - (piwithin/pibetween); /* eq. 3 */

    ZF_LOGD("END");
    return(fst);
}

/*Ewens-Watterson test*/
double EWtest(int Nsample, int *Freqhap)
{
    ZF_LOGD("START");
    int i;
    double H;
    
 	if(Nsample == 0) {
        ZF_LOGD("END");
 	    return -10000;
 	}
	if(Nsample < 2) {
        ZF_LOGD("END");
	    return 0.;
	}
	
	H = 0;
    for(i=0;i<Nsample;i++)
        H += Freqhap[i] * Freqhap[i];
    H = (double)H/((double)Nsample*(double)Nsample);

    ZF_LOGD("END");
    return H;
}

double pwh(int n,double *pid)
{/*try to find two  divergent groups in the mismatch distribution with low variation within in at least one*/
    ZF_LOGD("START");
	int ncomp;
	int i;
	double within,between;
	double diff=(double)0;
	int compare_(const void *,const void *);
	
	ncomp = n * (n-1) / 2;
	/*the pairwise vector has to be sorted*/
	qsort(pid,(int)ncomp,sizeof(double),compare_);/*all values are sorted*/
	/*divide sorted vector in two. find the highest difference*/
	for(i=0;i<ncomp-1;i++) {
		between= pid[i+1]-pid[i];
		within = fabs((pid[ncomp-1]-pid[i+1]) - (pid[i]-pid[0]));
		if((between > within) && (between + within > diff)) 
			diff = between + within;
	}

    ZF_LOGD("END");
	return diff;
}

/*Haplotype tests from Depaulis et al.*/
double testHap(int Nsample, int *Freqhap)
{/*Depaulis statistics*/
    ZF_LOGD("START");
    int i;
    double H;
    
	if(Nsample == 0) {
        ZF_LOGD("END");
	    return -10000;
	}
	if(Nsample < 2) {
        ZF_LOGD("END");
	    return 0.;
	}
	
    H = 0;
    for(i=0;i<Nsample;i++)
        H += Freqhap[i] * Freqhap[i];
    H = (double)1 - H/((double)Nsample*(double)Nsample);
    /*and weighted: haplotype diversity*/
    H = H*(double)Nsample/(double)(Nsample-1);

    ZF_LOGD("END");
    return H;
}

/*Fs from Fu */
double Fs(int Nsample, double pi, int NumAlelos)
{
    ZF_LOGD("START");
    /* Rozas program */

    double SumaP;
    double RestaP;
    int AleloI;
    long int i;
    double ValorFs;
    double *qew;
    double est_var;
    double FunEq23Ewens(int, int, double, double *);

    if(pi == 0.0 || Nsample < 2) {
        ZF_LOGD("END");
        return(-10000);
    }
    est_var = pi;
    qew  = (double *)malloc((long int)Nsample*(long int)Nsample*sizeof(double));

    for(i=0;i<(long int)Nsample*(long int)Nsample;i++)
    	qew[i] = -1.0;

    SumaP=RestaP=0.0;
    for (AleloI=1;AleloI<NumAlelos;AleloI++) {
        /* calculo q(n,aleloI)   ecuacion 21 (recurrente con eq. 19 y 20) */
        SumaP += FunEq23Ewens(Nsample, AleloI, est_var, qew);
    }

    if(SumaP > 1.-1E-37) {
    	for (AleloI = NumAlelos;AleloI <= Nsample; AleloI++)
            RestaP += FunEq23Ewens(Nsample, AleloI, est_var, qew);
        if(RestaP < 1E-37) {
            free(qew);
            ZF_LOGD("END");
            return -10000;
        }
        ValorFs = log((double)RestaP) - log((double)1.0-RestaP);
    }
    else {
        if(SumaP < 1E-37) {
            free(qew);
            ZF_LOGD("END");
            return -10000;
        } else
            ValorFs = log((double)1.0-(double)SumaP) - log((double)SumaP);
    }
    if (fabs(ValorFs) < 1.0E-37)
        ValorFs = 0.0;
    if (fabs(ValorFs) < -1.0E37) {
        free(qew);
        ZF_LOGD("END");
		return -10000;
    }
	free(qew);

    ZF_LOGD("END");
    return ValorFs;
}

double FunEq23Ewens(int N,int i,double theta, double *qew_)
{
    ZF_LOGD("START");
    /* Rozas program */

	long int acceso;
    double ValorN;  /* log del numerador */
    double ValorD;  /* log del denominador */

    acceso= (long int)((N-1) * N + i - 1);
    ValorN=0.0;
    ValorD=0.0;        
    if (qew_[acceso] < 0.0) {
        int jj;
        if (i==1) {
            /* calculo de qj,1   (i = 1)   Antigua equacion 19  */
            if(N > 2) {             
                for (jj=2;jj<N;jj++)
                    ValorN = ValorN + log((double)jj);  
            }
            ValorN = ValorN + log(theta);
            for(jj=0;jj<N;jj++)
                ValorD  = ValorD + log(theta + (double)jj);
            qew_[acceso] = exp(ValorN - ValorD);
        }    
        if(i==N) {          
            /* calculo de qj,j   (n = i)   antigua equacion 20 */
            ValorN = log(theta) * (double)N;
            for(jj=0;jj<N;jj++)     
                ValorD  = ValorD + log((double)theta + (double)jj);
            qew_[acceso] = exp((double)(ValorN - ValorD));
	}
	    if(i>1 && i<N) {
	        // TODO: try to explore tail recursion (https://en.m.wikipedia.org/wiki/Tail_call)
            /*  recursividad  */
            qew_[acceso] = FunEq23Ewens(N-1,i,  theta,qew_) * ((N-1)/(theta + (double)N-1.0))
                         + FunEq23Ewens(N-1,i-1,theta,qew_) *         (theta/(theta + (double)N-1.0));
        }    
    }

    ZF_LOGD("END");
    return(qew_[acceso]);
}

/* Programa de Ying */
double estnm(int npop,int nsam,int *config,long int segsites,char **list)
{
    ZF_LOGD("START");
    int i;
    double gst,within,between;
    void diff_within_between(int,long int,int,int *,char **,double *,double *);

    nsam = 0 ;
    for(i=0;i<npop; i++) nsam += config[i] ;

    diff_within_between(nsam,segsites,npop,config,list,&within,&between);
    if(within == 0) gst = 0.;
    else gst = 1. - within/(within/npop + (1.-1./npop)*between);

    ZF_LOGD("END");
    return (gst);
}

void diff_within_between(int nsam,long int ns,int npop,int *config,char **list,double *pwithin,double *pbetween)
{
    ZF_LOGD("START");
    int pop,ind,ind2,indpop,endpop,count ;
    double diff(long int,char *,char *);
	
    *pwithin = *pbetween = 0. ;
    for(pop=indpop=ind=endpop=count=0; pop<npop; pop++,ind++)
        for(endpop += config[pop];ind < endpop-1; ind++) 
            for(indpop = 1; (indpop+ind)<endpop; indpop++) {
                *pwithin += diff(ns,list[ind],list[ind+indpop]);
                count++;
            }
    *pwithin /= count ; /*Es divideix pel conjunt de totes les combinacions. no es separa el valor de cada subpoblacio*/
    for(pop=ind=count=endpop=0;pop<npop-1;pop++)
        for(endpop += config[pop]; ind<endpop; ind++)
            for(ind2 = endpop ; ind2<nsam ; ind2++) {
                *pbetween += diff(ns,list[ind],list[ind2]);
                count++;
            }
    *pbetween /= count;

    ZF_LOGD("END");
}

double diff(long int ns,char *gam1,char *gam2)
{
    ZF_LOGD("START");
    long int  i;
    double count;
	
    count = 0.;
    for(i=0;i<ns;i++) if(gam1[i] != gam2[i]) count += 1.;

    ZF_LOGD("END");
    return(count);
}

/*R2 Ramos & Rozas: "*unic" is the number of singletons in each sequence (in comparison to the sample studied)*/	
double R2(long int *unic,double pi,int sample_size,long int S)
{
    ZF_LOGD("START");
    double sm2 = 0.0;
    int i;
    
    if(S == 0 || sample_size == 0) return(-10000);
    for (i=0;i<sample_size;i++)
            sm2 += ((double)unic[i] - pi/2.0)*((double)unic[i] - pi/2.0);
    
    sm2 = sqrt(sm2/((double)sample_size))/(double)S;
            
    if (fabs(sm2) < 1.0E-15)
            sm2 = 0.0;


    ZF_LOGD("END");
    return (double)sm2;
}

/*Gxi(outgroup) test from Fu(1995) based on frequency of segregating sites, but using theta from watterson*/
double Gxi(int sample_size,int *fr, double thetaw) /* with outgroup */
{
    ZF_LOGD("START");
    int i;
    double G;
    double varxi(int,int,double);
    
    if(thetaw == 0.0 || sample_size < 2) return(-10000);

    G = 0.;
    for(i=1;i<sample_size;i++) 
        G += ((fr[i] - thetaw/(double)i) * (fr[i] - thetaw/(double)i)) / varxi(i,sample_size,thetaw);
    G /= ((double)sample_size - (double)1);

    ZF_LOGD("END");
    return G;
}

double varxi(int i,int sample_size,double theta)
{
    ZF_LOGD("START");
    double sigma(int,int);
    double result = (double)1/((double)i)*theta + sigma(i,sample_size)*theta*theta;

    ZF_LOGD("END");
    return result;
}

double sigma(int i,int sample_size)
{
    ZF_LOGD("START");
    double ai(int);
    double bn(int,int);
    double result = -10000;
    
    if(i <  (double)sample_size/(double)2) {
        result = bn(i+1,sample_size);
    }
    if(i == (double)sample_size/(double)2) {
        result = (double)2*(ai(sample_size) - ai(i))/((double)(sample_size - i)) - (double)1/((double)(i*i));
    }
    if(i >  (double)sample_size/(double)2) {
        result = bn(i,sample_size) - (double)1/((double)(i*i));
    }

    ZF_LOGD("END");
    return result;
}

double ai(int i)
{
    ZF_LOGD("START");
    int j;
    double a = 0;
    for(j=1;j<i;j++) a += (double)1/(double)j;

    ZF_LOGD("END");
    return a;
}

/*Gximod (outgroup) test from Fu(1995) based on frequency of segregating sites, but using theta from watterson,
 and not variance...*/
double Gximod(int sample_size,int *fr, double thetaw) /* with outgroup */
{
    ZF_LOGD("START");
    int i;
    double G;
    
    if(thetaw == 0.0 || sample_size < 2) {
        ZF_LOGD("END");
        return(-10000);
    }

    G = 0.;
    for(i=1;i<sample_size;i++) 
        G += ((fr[i] - thetaw/(double)i) * (fr[i] - thetaw/(double)i))/((double)ceil(thetaw/(double)i));

    ZF_LOGD("END");
    return G;
}

double frabs(int sample_size,int *fr, double thetaw) /* with outgroup */
{
    ZF_LOGD("START");
    int i;
    double G;
    
    if(thetaw == 0.0 || sample_size < 2) {
        ZF_LOGD("END");
        return(-10000);
    }

    G = 0.;
    for(i=1;i<sample_size;i++) 
        G += (double)fabs((double)fr[i] - thetaw/(double)i);

    ZF_LOGD("END");
    return G/(double)sample_size;
}

/*Frequency tests from Achaz Genetics 2009: outgroup*/
double freqtesto_achaz(int sample_size,int *fr,int singleton,double *w1,double *w2) /* nomes outgroup */
{
    ZF_LOGD("START");
    int i,j,ss;
    double Th1,Th2,Test,Thw,Thw2,*ww;
	double sumw1,sumw2,sumww,omi;
	double alfan,betan,alfat,betat;
	double omegai(int,int,double *,double *);
	double sigmaii(int,int),sigmaij(int,int,int);
	double an(int),a2n(int);
    
    if(sample_size < 4) {
        ZF_LOGD("END");
        return(-10000);
    }
	
	ww = (double *)calloc(sample_size,sizeof(double));
    
    Th1 = 0.;
	sumw1 = 0.;
    for(i=1;i<sample_size;i++) {
		Th1 += ((double)*(fr+i))*((double)i*(double)w1[i]);
		sumw1 += w1[i];
	}
    Th1 /= sumw1;
	
	Th2 = 0.;
	sumw2 = 0.;
    for(i=1;i<sample_size;i++) {
		Th2 += ((double)*(fr+i))*((double)i*(double)w2[i]);
		sumw2 += w2[i];
	}
    Th2 /= sumw2;
	
    if(Th1 == 0. && Th2 == 0.) {
        ZF_LOGD("END");
        return(-10000);
    }
	
	Thw = 0.;
	sumww = 0.;
    for(i=1;i<sample_size;i++) {
		if(i==1) ss = singleton;
		else ss = 1;
		ww[i] = (double)1/(double)i * (double)ss;
		Thw += (double)*(fr+i)*(double)i*ww[i];
		sumww += ww[i];
	}
    Thw /= sumww;
	
	alfan = 0.;
	for(i=1;i<sample_size;i++) {
		omi = omegai(sample_size,i,w1,w2);
		alfan += i*(omi*omi);
	}
	
 	betan = 0.;
	for(i=1;i<sample_size;i++) {
		omi = omegai(sample_size,i,w1,w2);
		betan += i*i * (omi*omi) * sigmaii(sample_size,i);
		/*printf("\nbetan=%f\tsigmaii[n=%d,i=%d]=%f",betan,sample_size,i,sigmaii(sample_size,i));*/
		for(j=i+1;j<sample_size;j++) {
			betan += 2.0 * i*j * omegai(sample_size,i,w1,w2) * omegai(sample_size,j,w1,w2) * sigmaij(sample_size,j,i);
			/*printf("\nbetan=%f\tsigmaij[n=%d,i=%d,j=%d]=%f",betan,sample_size,i,j,sigmaij(sample_size,j,i));*/
		}
	}
	
	/*Theta2*/
	alfat = 0.;
	for(i=1;i<sample_size;i++) {
		alfat += (ww[i]/sumww * ww[i]/sumww)*i;
	}	
 	betat = 0.;
	for(i=1;i<sample_size;i++) {
		betat += i*i * (ww[i]/sumww * ww[i]/sumww) * sigmaii(sample_size,i);
		for(j=i+1;j<sample_size;j++) {
			betat += 2.0 * i*j * ww[i]/sumww * ww[j]/sumww * sigmaij(sample_size,j,i);
		}
	}
	Thw2 = (Thw*Thw - alfat*Thw)/(1.0 + betat);	
	
	/*Test*/
	Test = (Th1 - Th2)/(sqrt(alfan*Thw + betan*Thw2));
	
	free(ww);
	if (fabs(Test) < 1.0E-15)
		Test = 0.0;

    ZF_LOGD("END");
    return Test;
}
/*Frequency tests from Achaz Genetics 2009: NO outgroup*/
double freqtestn_achaz(int sample_size,int *fr,int singleton,double *w1,double *w2) /* NO outgroup */
{
    ZF_LOGD("START");
    int i,j,ss;
    double Th1,Th2,Test,Thw,Thw2,*ww;
	double sumw1,sumw2,sumww,omi,omj,psi,psj;
	double psii(int,int),rhoii(int,int),rhoij(int,int,int);
	double alfan,betan,alfat,betat;
	double omegain(int,int,double *,double *);
	double sigmaii(int,int),sigmaij(int,int,int);
	double an(int),a2n(int);
    
    if(sample_size < 4 || (sample_size < 6 && singleton==0)) {
        ZF_LOGD("END");
        return(-10000);
    }
	
	ww = (double *)calloc(sample_size,sizeof(double));
    
    Th1 = 0.;
	sumw1 = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
		Th1 += ((double)*(fr+i))*(double)w1[i]/((double)psii(sample_size,i));
		sumw1 += w1[i];
	}
    Th1 /= sumw1;
	
	Th2 = 0.;
	sumw2 = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
		Th2 += ((double)*(fr+i))*(double)w2[i]/((double)psii(sample_size,i));
		sumw2 += w2[i];
	}
    Th2 /= sumw2;
	
    if(Th1 == 0. && Th2 == 0.) {
        ZF_LOGD("END");
        return(-10000);
    }
	
	Thw = 0.;
	sumww = 0.;
	for(i=1;i<=floor(sample_size/2);i++) {
		if(i==1) ss = singleton;
		else ss = 1;
		if(i == sample_size-i) ww[i] = (double)sample_size/((double)i*(double)(sample_size - i)*2.0)*(double)ss;
		else ww[i] = (double)sample_size/((double)i*(double)(sample_size - i)*1.0)*(double)ss;
		Thw += ((double)*(fr+i))*ww[i]/((double)psii(sample_size,i));
		sumww += ww[i];
	}
    Thw /= sumww;
	
	alfan = 0.;
	for(i=1;i<=floor(sample_size/2);i++) {
		omi = omegain(sample_size,i,w1,w2);
		psi = psii(sample_size,i);
		alfan += (omi*omi)/psi;
	}
	
 	betan = 0.;
	for(i=1;i<=floor(sample_size/2);i++) {
		omi = omegain(sample_size,i,w1,w2);
		psi = psii(sample_size,i);
		betan += omi/psi * omi/psi * rhoii(sample_size,i);
		for(j=i+1;j<=floor(sample_size/2);j++) {
			omj = omegain(sample_size,j,w1,w2);
			psj = psii(sample_size,j);
			betan += 2.0 * omi/psi * omj/psj * rhoij(sample_size,j,i);
		}
	}
	
	/*Theta2*/
	alfat = 0.;
	for(i=1;i<=floor(sample_size/2);i++) {
		psi = psii(sample_size,i);
		alfat += (ww[i]/sumww * ww[i]/sumww)/psi;
	}	
 	betat = 0.;
	for(i=1;i<=floor(sample_size/2);i++) {
		psi = psii(sample_size,i);
		betat += (ww[i]/sumww)/psi * (ww[i]/sumww)/psi * rhoii(sample_size,i);
		for(j=i+1;j<=floor(sample_size/2);j++) {
			psj = psii(sample_size,j);
			betat += 2.0 * (ww[i]/sumww)/psi * (ww[j]/sumww)/psj * rhoij(sample_size,j,i);
		}
	}
	Thw2 = (Thw*Thw - alfat*Thw)/(1.0 + betat);	
	
	/*Test*/
	Test = (Th1 - Th2)/(sqrt(alfan*Thw + betan*Thw2));
	
	free(ww);
	if (fabs(Test) < 1.0E-15)
		Test = 0.0;

    ZF_LOGD("END");
    return Test;
}

double an(int n)
{
    ZF_LOGD("START");
    static double *an_m=0;
    static int nsam=0;
    long int i;
    
    if(nsam == 0) {
        an_m = (double *)calloc(n+2, sizeof(double));
        for(i=1;i<n+1;i++)
            an_m[i+1] = an_m[i] + 1.0/i;
        nsam = n;
    } else {
        if(nsam < n) {
            an_m = (double *)realloc(an_m,(n+2)*sizeof(double));
            for(i=nsam+1;i<n+1;i++)
                an_m[i+1] = an_m[i] + 1.0/i;
            nsam = n;
        }
    }

    ZF_LOGD("END");
    return an_m[n];
}
double a2n(int n)
{
    ZF_LOGD("START");
    double a2ni = 0.;
    int i;
    for(i=1;i<n;i++)
        a2ni += (1.0/((double)i*(double)i));

    ZF_LOGD("END");
    return a2ni;
}
double bn(int n,int i)
{
    ZF_LOGD("START");
    double an(int);
    double bni;
    bni = (2.0*(double)n)/((double)(n-i+1)*(double)(n-i)) * (an(n+1) - an(i)) - 2.0/((double)(n-i));

    ZF_LOGD("END");
    return bni;
}

double sigmaii(int n,int i)
{
    ZF_LOGD("START");
	double an(int),bn(int,int);
	double sii=0;
	if(2*i<n)
		sii = bn(n,i+1);
	if(2*i==n)
		sii = 2.0 * (an(n)-an(i))/((double)(n-i)) - 1.0/((double)i*(double)i);
	if(2*i>n)
		sii = bn(n,i) - 1.0/((double)i*(double)i);

    ZF_LOGD("END");
	return sii;
}
double sigmaij(int n,int i,int j)
{
    ZF_LOGD("START");
	double an(int),bn(int,int);
	double sigmaii(int,int);
	double sij=0;

	if (i==j) {
	    sij = sigmaii(n,i);
	} else {
        int ii,jj;

        if(i < j) {
            ii = j;
            jj = i;
        } else {
            ii = i;
            jj = j;
        }

        if((ii+jj)<n)
            sij = (bn(n,ii+1) - bn(n,ii))/2.0;
        if((ii+jj)==n)
            sij = (an(n)-an(ii))/((double)(n-ii)) + (an(n)-an(jj))/((double)(n-jj)) - (bn(n,ii) + bn(n,jj+1))/2.0 - 1.0/((double)ii*(double)jj);
        if((ii+jj)>n)
            sij = (bn(n,jj) - bn(n,jj+1))/2.0 - 1.0/((double)ii*(double)jj);
	}

    ZF_LOGD("END");
	return sij;
}
double omegai(int n,int i,double *w1,double *w2)
{
    ZF_LOGD("START");
	double omi;
	int x;
	double sumw1,sumw2;
	
	sumw1=0.0;
	for(x=1;x<n;x++) sumw1 += w1[x];
	sumw2=0.0;
	for(x=1;x<n;x++) sumw2 += w2[x];
	omi = w1[i]/sumw1 - w2[i]/sumw2;

    ZF_LOGD("END");
	return omi;
}
double psii(int n,int i)
{
    ZF_LOGD("START");
	double psi;
	int krond;
	
	if(i==n-i) krond = 1;
	else krond = 0;
	psi= (double)n/((double)(1.+krond)*i*(n-i));

    ZF_LOGD("END");
	return psi;
}
double rhoii(int n,int i)
{
    ZF_LOGD("START");
	double sigmaii(int,int);
	double sigmaij(int,int,int);
	double rhoi;
	int krond;
	
	if((int)i==(int)(n-i)) 
		krond = 1;
	else 
		krond = 0;
	
	rhoi  = (sigmaii(n,i)+sigmaii(n,n-i)+2.0*sigmaij(n,i,n-i));
	rhoi /= ((1.0 + krond) * (1.0 + krond));

    ZF_LOGD("END");
	return rhoi;
}
double rhoij(int n,int i,int j)
{
    ZF_LOGD("START");
	double sigmaii(int,int);
	double sigmaij(int,int,int);
	double rhoj;
	int krondi,krondj;
	
	if(i==(n-i)) krondi = 1;
	else krondi = 0;
	if(j==(n-j)) krondj = 1;
	else krondj = 0;
	
	rhoj  = (sigmaij(n,i,j)+sigmaij(n,i,n-j)+sigmaij(n,n-i,j)+sigmaij(n,n-i,n-j));
	rhoj /= ((1.0 + krondi) * (1.0 + krondj));

    ZF_LOGD("END");
	return rhoj;
}
double omegain(int n,int i,double *w1,double *w2)
{
    ZF_LOGD("START");
	double omi;
	int x;
	double sumw1,sumw2;
	
	sumw1=0.0;
	for(x=1;x<=floor(n/2);x++) sumw1 += w1[x];
	sumw2=0.0;
	for(x=1;x<=floor(n/2);x++) sumw2 += w2[x];
	omi = w1[i]/sumw1 - w2[i]/sumw2;

    ZF_LOGD("END");
	return omi;
}


double fl_d_achaz(int nsam,int *freq,long int S)
{
    ZF_LOGD("START");
	double Test;
	double *w1,*w2;
	int x;
	double freqtesto_achaz(int,int *,int,double *,double *);
	
	w1 = (double *)calloc(nsam,sizeof(double));
	w2 = (double *)calloc(nsam,sizeof(double));
	
	for(x=1;x<nsam;x++) {
		w1[x] = 1.0/(double)x;
		if(x==1) w2[x] = 1.0;
		else w2[x] = 0.0;
	}
	Test = freqtesto_achaz(nsam,freq,1,w1,w2);
	
	free(w1);
	free(w2);

    ZF_LOGD("END");
	return Test;
}
double fl_d2_achaz(int nsam,int *freq,long int S)
{
    ZF_LOGD("START");
	double Test;
	double *w1,*w2;
	int x,*freqn;
	double freqtestn_achaz(int,int *,int,double *,double *);
	
	w1 = (double *)calloc(nsam,sizeof(double));
	w2 = (double *)calloc(nsam,sizeof(double));
	freqn = (int *)calloc(nsam,sizeof(double));
	
	for(x=1;x<=floor(nsam/2);x++) {
		if(x == nsam-x) w1[x] = (double)nsam/((double)x*(double)(nsam - x)*2.0);
		else w1[x] = (double)nsam/((double)x*(double)(nsam - x)*1.0);
		if(x==1) w2[x] = nsam;
		else w2[x] = 0.0;
		if(x == nsam-x) freqn[x] = (int)((freq[x] + freq[nsam-x])/2.0);
		else freqn[x] = freq[x] + freq[nsam-x];
	}
	Test = freqtestn_achaz(nsam,freqn,1,w1,w2);
	
	free(w1);
	free(w2);
	free(freqn);

    ZF_LOGD("END");
	return Test;
}
double fl_f_achaz(int nsam,int *freq,long int S)
{
    ZF_LOGD("START");
	double Test;
	double *w1,*w2;
	int x;
	double freqtesto_achaz(int,int *,int,double *,double *);
	
	w1 = (double *)calloc(nsam,sizeof(double));
	w2 = (double *)calloc(nsam,sizeof(double));
	
	for(x=1;x<nsam;x++) {
		w1[x] = (double)(nsam - x);
		if(x==1) w2[x] = 1.0;
		else w2[x] = 0.0;
	}
	Test = freqtesto_achaz(nsam,freq,1,w1,w2);
	
	free(w1);
	free(w2);

    ZF_LOGD("END");
	return Test;
}
double fl_f2_achaz(int nsam,int *freq,long int S)
{
    ZF_LOGD("START");
	double Test;
	double *w1,*w2;
	int x,*freqn;
	double freqtestn_achaz(int,int *,int,double *,double *);
	
	w1 = (double *)calloc(nsam,sizeof(double));
	w2 = (double *)calloc(nsam,sizeof(double));
	freqn = (int *)calloc(nsam,sizeof(double));
	
	for(x=1;x<=floor(nsam/2);x++) {
		if(x == nsam-x) w1[x] = (double)nsam/2.0;
		else w1[x] = (double)nsam/1.0;
		if(x==1) w2[x] = nsam;
		else w2[x] = 0.0;
		if(x == nsam-x) freqn[x] = (int)((freq[x] + freq[nsam-x])/2.0);
		else freqn[x] = freq[x] + freq[nsam-x];
	}
	Test = freqtestn_achaz(nsam,freqn,1,w1,w2);
	
	free(w1);
	free(w2);
	free(freqn);

    ZF_LOGD("END");
	return Test;
}
double Y_achaz(int nsam,int *freq,long int S)
{
    ZF_LOGD("START");
	double Test;
	double *w1,*w2;
	int x;
	double freqtesto_achaz(int,int *,int,double *,double *);
	
	w1 = (double *)calloc(nsam,sizeof(double));
	w2 = (double *)calloc(nsam,sizeof(double));
	
	for(x=1;x<nsam;x++) {
		if(x==1) w1[x] = 0.0;
		else w1[x] = (double)(nsam - x);
		if(x==1) w2[x] = 0.0;
		else w2[x] = 1.0/(double)x;
	}
	Test = freqtesto_achaz(nsam,freq,0,w1,w2);
	
	free(w1);
	free(w2);

    ZF_LOGD("END");
	return Test;
}
double Y2_achaz(int nsam,int *freq,long int S)
{
    ZF_LOGD("START");
	double Test;
	double *w1,*w2;
	int x,*freqn;
	double freqtestn_achaz(int,int *,int,double *,double *);
	
	w1 = (double *)calloc(nsam,sizeof(double));
	w2 = (double *)calloc(nsam,sizeof(double));
	freqn = (int *)calloc(nsam,sizeof(double));
	
	for(x=1;x<=floor(nsam/2);x++) {
		if(x == 1) w1[x] = 0.0;
		else {
			if(x == nsam-x) w1[x] = (double)nsam/2.0;
			else w1[x] = (double)nsam/1.0;
		}
		if(x == 1) w2[x] = 0.0;
		else {
			if(x == nsam-x) w2[x] = (double)nsam/((double)x*(double)(nsam - x)*2.0);
			else w2[x] = (double)nsam/((double)x*(double)(nsam - x)*1.0);
		}
		if(x == nsam-x) freqn[x] = (int)((freq[x] + freq[nsam-x])/2.0);
		else freqn[x] = freq[x] + freq[nsam-x];
	}
	Test = freqtestn_achaz(nsam,freqn,0,w1,w2);
	
	free(w1);
	free(w2);
	free(freqn);

    ZF_LOGD("END");
	return Test;
}

double raggadeness(long int *pwd,long int max_pwd,long int num_comp)
{
	/*
	 Compute the sum of square differences between one position and its neighbour 
	 in a pairwise distribution, using the frequency of each value. Eller et al.
	 MBE 13(8), 1996.
	 */
    ZF_LOGD("START");
	
	int i;
	double r = 0.0;
	
	for (i=1;i<max_pwd;i++)
	{
		r = r + (((double)(*(pwd+i)/(double)num_comp) - ((double)(*(pwd+(i-1)))/(double)num_comp)) * 
		         ((double)(*(pwd+i)/(double)num_comp) - ((double)(*(pwd+(i-1)))/(double)num_comp)));		
	}
	
	/*The last with 0 */
	
	r = r + (((double)0-((double)(*(pwd+(max_pwd-1)))/(double)num_comp)) * 
	         ((double)0-((double)(*(pwd+(max_pwd-1)))/(double)num_comp)));

    ZF_LOGD("END");
	return r;
}
