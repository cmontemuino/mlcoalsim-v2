#ifndef MLCOALSIM_NEUT_TESTS_H
#define MLCOALSIM_NEUT_TESTS_H

double E_zeng(int n,double thetaL,double thetaw,double S,double *coef);
double estnm(int npop,int nsam,int *config,long int segsites,char **list);
double EWtest(int Nsample, int *Freqhap);
void init_coef(double *,int);
double fay_wu(int sample_size,int *fr,double pi);
double fay_wu_normalized2(int n,double thetaL,double thetaw,double S,double *coef,double pi);
double fl_d_achaz(int nsam,int *freq,long int S);
double fl_d2_achaz(int nsam,int *freq,long int S);
double fl_f_achaz(int nsam,int *freq,long int S);
double fl_f2_achaz(int nsam,int *freq,long int S);
double frabs(int sample_size,int *fr, double thetaw);
double Fs(int Nsample, double pi, int NumAlelos);
double Fst(double piwithin, double pibetween,int ntotpop);
double Gxi(int sample_size,int *fr, double thetaw);
double Gximod(int sample_size,int *fr, double thetaw);
double pwh(int n,double *pid);
double R2(long int *unic,double pi,int sample_size,long int S);
double raggadeness(long int *pwd,long int max_pwd,long int num_comp);
double tajima_d(double k_, int S_, double *coef_taj);
double tajima_dvsdmin(double k_, int S_, double *coef_taj,int sample_size);
double testHap(int Nsample, int *Freqhap);
double Y_achaz(int nsam,int *freq,long int S);
double Y2_achaz(int nsam,int *freq,long int S);


#endif //MLCOALSIM_NEUT_TESTS_H
