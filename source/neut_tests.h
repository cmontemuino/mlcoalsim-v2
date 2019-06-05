#pragma once

#ifndef _NEUT_TESTS_H_
#define _NEUT_TESTS_H_

void init_coef(double *p,int sample_size);
double tajima_d(double k_, int S_, double *coef_taj);
double tajima_dvsdmin(double k_, int S_, double *coef_taj,int sample_size);
double fl_d(int sample_size,int fr1,int S, double *coef);
double fl_d2(int sample_size,int fr1w,int S, double *coef);
double fl_f(int sample_size,int fr1, int S, double pi, double *coef);
double fl_f2(int sample_size,int fr1w, int S, double pi, double *coef);
double fay_wu(int sample_size,int *fr,double pi);
double fay_wuvsminH(int sample_size,int *fr,double pi,int S);
double fay_wu_normalized(int n,int *fr,double pi);
double fay_wu_normalized2(int n,double thetaL,double thetaw,double S,double *coef,double pi);
double E_zeng(int n,double thetaL,double thetaw,double S,double *coef);
double Fst(double piwithin, double pibetween,int ntotpop);
double EWtest(int Nsample, int *Freqhap);
double pwh(int n,double *pid);
double testHap(int Nsample, int *Freqhap);
double Fs(int Nsample, double pi, int NumAlelos);
double FunEq23Ewens(int N,int i,double theta, double *qew_);
double estnm(int npop,int nsam,int *config,long int segsites,char **list);
void diff_within_between(int nsam,long int ns,int npop,int *config,char **list,double *pwithin,double *pbetween);
double diff(long int ns,char *gam1,char *gam2);
double R2(long int *unic,double pi,int sample_size,long int S);
double Gxi(int sample_size,int *fr, double thetaw);
double varxi(int i,int sample_size,double theta);
double sigma(int i,int sample_size);
double ai(int i);
double Gximod(int sample_size,int *fr, double thetaw);
double frabs(int sample_size,int *fr, double thetaw);
double freqtesto_achaz(int sample_size,int *fr,int singleton,double *w1,double *w2);
double freqtestn_achaz(int sample_size,int *fr,int singleton,double *w1,double *w2);
double an(int n);
double a2n(int n);
double bn(int n,int i);
double sigmaii(int n,int i);
double sigmaij(int n,int i,int j);
double omegai(int n,int i,double *w1,double *w2);
double psii(int n,int i);
double rhoii(int n,int i);
double rhoij(int n,int i,int j);
double omegain(int n,int i,double *w1,double *w2);
double fl_d_achaz(int nsam,int *freq,long int S);
double fl_d2_achaz(int nsam,int *freq,long int S);
double fl_f_achaz(int nsam,int *freq,long int S);
double fl_f2_achaz(int nsam,int *freq,long int S);
double Y_achaz(int nsam,int *freq,long int S);
double Y2_achaz(int nsam,int *freq,long int S);
double raggadeness(long int *pwd,long int max_pwd,long int num_comp);

#endif
