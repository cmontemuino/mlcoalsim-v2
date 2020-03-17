#pragma once

#ifndef MLCOALSIM_PRINT_HELPER_H
#define MLCOALSIM_PRINT_HELPER_H

#include <stdio.h>
#include "mlsp_sm.h"

void print_matrix(long int, struct var2 **,FILE *,int,long int,struct var_priors *,int, char **dna_matrix, long int *posit);

void show_progress(double counterp, double *restp, long int *p);

#endif //MLCOALSIM_PRINT_HELPER_H
