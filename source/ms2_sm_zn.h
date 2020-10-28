#ifndef MLCOALSIM_MS2_SM_ZN_H
#define MLCOALSIM_MS2_SM_ZN_H

double ZnA_(int valuep,long int segsites,struct var2 **inputp,int npopa, char **mutations_matrix, long int *positions);
double ZnA_window(struct var2 **inputp, long int s0, long int s1, int npopa, char **mutations_matrix, long int *positions);
double Zns(int valuep, long int segsites, struct var2 **inputp, int npopa, char **mutations_matrix, long int *positions);
double Zns_window(struct var2 **inputp, long int s0, long int s1, int npopa, char **mutations_matrix, long int *positions);

#endif //MLCOALSIM_MS2_SM_ZN_H
