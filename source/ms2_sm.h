#ifndef MLCOALSIM_MS2_SM_H
#define MLCOALSIM_MS2_SM_H

#include <stdio.h>
int print_matrix_sfixalltheta(struct var **data,FILE *file_thetap,FILE *file_ttotp,struct prob_par **postp);
int print_matrix_rmfix(struct var **data,FILE *file_recp,FILE *file_ttotp,int Sfix,struct prob_par **postp);

#endif //MLCOALSIM_MS2_SM_H
