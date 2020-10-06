#ifndef MLCOALSIM_MS2_SM_HELPERS_H
#define MLCOALSIM_MS2_SM_HELPERS_H

int compare_(const void *i,const void *j);
int do_heter_gamma_sites(double *categories,double gammashape,double poppar,long int nsites,long int invariable) ;
int do_heter_gamma_sitesrec(double *categories,double gammashape,double poppar,long int nsites);
double fixoutg(int nsam, int Sout, int freq0);
double koutgJC(int nsam, int Sout, int *freq, unsigned long nsites);
long int localize_positiontop(const double *categories,double valuer,long int start,long int end);
double logPPoisson2(long int Si, double lambda);

#endif //MLCOALSIM_MS2_SM_HELPERS_H
