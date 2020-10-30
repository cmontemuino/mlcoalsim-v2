#ifndef MLCOALSIM_MS2_SM_HELPERS_H
#define MLCOALSIM_MS2_SM_HELPERS_H

#include "mlsp_sm.h"

void biggerlist(int nsam, char **mutations_matrix, long int maxsites);
int compare_(const void *i,const void *j);
int do_heter_gamma_sites(double *categories,double gammashape,double poppar,long int nsites,long int invariable) ;
int do_heter_gamma_sitesrec(double *categories,double gammashape,double poppar,long int nsites);
double fixoutg(int nsam, int Sout, int freq0);
double koutgJC(int nsam, int Sout, int *freq, unsigned long nsites);
long int localize_positiontop(const double *categories,double valuer,long int start,long int end);
void locate(long int n, long int beg, long int *ptr, int mhits, double *weightmut, long int end);
void locate2(long int n, long int beg, long int *ptr, int mhits, int nsam, double *weightmut, long int end, char **mutations_matrix);
void locate_psel(long int n, long int beg, long int *ptr, int mhits, long int sel_nt, double *weightmut, long int end);
void locate2_psel(long int n, long int beg, long int *ptr, int mhits, int nsam, long int sel_nt, double *weightmut, long int end, char **mutations_matrix);
double logPPoisson2(long int Si, double lambda);
void order(long int n, long int *pbuf);
void ordran(long int n, long int *pbuf, int mhits, long int beg, double *weightmut, long int end);
void ranvec(long int n, long int *pbuf, int mhits, long int beg, double *weightmut, long int end);

#endif //MLCOALSIM_MS2_SM_HELPERS_H
