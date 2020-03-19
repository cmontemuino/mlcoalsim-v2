#ifndef MLCOALSIM_STREEC2_SM_H
#define MLCOALSIM_STREEC2_SM_H

#include "mlsp_sm.h"

struct segl *segtre_mig(long int npop,int nsam,int *inconfig,long int nsites,double r,double f,double track_len,
                        double mig_rate,long int *pnsegs,long int iteration, double *factor,int ifselection, double pop_sel, double sinit,
                        double pop_size,long int sel_nt,int *all_sel,int *selnsam, int *nintn,double **nrec,double **npast, double **tpast,
                        int split_pop, double time_split, double time_scoal, double factor_anc, double *freq, double tlimit,int iflogistic,
                        double *Tts,double factor_chrn,double *weightrec, double **migrate_matrix,int my_rank, int npop_events,
                        struct events *pop_event,int event_forsexratio, double event_sexratio, double sex_ratio, int no_rec_males,
                        double sendt, double sfreqend,double sfreqinit);

#endif //MLCOALSIM_STREEC2_SM_H
