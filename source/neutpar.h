#ifndef MLCOALSIM_NEUTPAR_H
#define MLCOALSIM_NEUTPAR_H

#include "mlsp_sm.h"

void calc_neutparSRH(long int segsit, struct var2 **inputp, struct dnapar *ntpar, double valuer, int npopa, int flaghap, char **mutations_matrix, long int *positions);
void calc_neutpar_windowSRH(struct var2 **inputp,struct dnapar *ntpar,long int s0, long int s1, double valuer, int npopa, int flaghap, char **mutations_matrix, long int *positions);
int ispolnomhit(long int j,int init, int nsam, int totalsam, char **mutations_matrix, long int *positions);
int Min_rec(int x, int segsit, int nsam, int inits, int totalsam, char **mutations_matrix, long int *positions);

#endif //MLCOALSIM_NEUTPAR_H
