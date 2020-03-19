#ifndef MLCOALSIM_NEUTPAR_WINDOW_H
#define MLCOALSIM_NEUTPAR_WINDOW_H

#include "mlsp_sm.h"

void calc_neutpar_window(
        struct var2 **inputp,
        struct dnapar *ntpar,
        long int s0,
        long int s1,
        double valuer,
        int npopa,
        char **dna_matrix,
        long int *posit);

#endif //MLCOALSIM_NEUTPAR_WINDOW_H
