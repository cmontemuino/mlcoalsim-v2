#ifndef MLCOALSIM_MS2_SM_HELPERS_H
#define MLCOALSIM_MS2_SM_HELPERS_H

#include "mlsp_sm.h"

void biggerlist(int nsam, char **mutations_matrix, long int maxsites);
int compare_(const void *i,const void *j);
int do_heter_gamma_sites(double *categories,double gammashape,double poppar,long int nsites,long int invariable) ;
int do_heter_gamma_sitesrec(double *categories,double gammashape,double poppar,long int nsites);
double fixoutg(int nsam, int Sout, int freq0);
int ispolnomhit(long int j,int init, int nsam, int totalsam, char **mutations_matrix, long int *positions);
double koutgJC(int nsam, int Sout, int *freq, unsigned long nsites);
long int localize_positiontop(const double *categories,double valuer,long int start,long int end);
double logPPoisson2(long int Si, double lambda);

#endif //MLCOALSIM_MS2_SM_HELPERS_H
