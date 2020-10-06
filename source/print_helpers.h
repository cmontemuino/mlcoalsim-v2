#ifndef MLCOALSIM_PRINT_HELPERS_H
#define MLCOALSIM_PRINT_HELPERS_H

#include <stdio.h>

int print_matrix_rmfix(struct var **data,FILE *file_recp,FILE *file_ttotp,int Sfix,struct prob_par **postp);
int print_matrix_sfixalltheta(struct var **data,FILE *file_thetap,FILE *file_ttotp,struct prob_par **postp);

#endif //MLCOALSIM_PRINT_HELPERS_H
