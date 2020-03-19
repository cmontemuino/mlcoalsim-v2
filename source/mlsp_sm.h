/************************************************************************************/
/* Copyright (C) 2006  Sebastian E. Ramos-Onsins                                    */
/*                                                                                  */ 
/* This program is free software; you can redistribute it and/or                    */
/* modify it under the terms of the GNU General Public License                      */
/* as published by the Free Software Foundation; either version 2                   */
/* of the License, or (at your option) any later version.                           */
/*                                                                                  */
/* This program is distributed in the hope that it will be useful,                  */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of                   */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    */
/* GNU General Public License for more details.                                     */
/*                                                                                  */
/* You should have received a copy of the GNU General Public License                */
/* along with this program; if not, write to the Free Software                      */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.  */
/************************************************************************************/

#ifndef MLCOALSIM_MLSP_SM_H
#define MLCOALSIM_MLSP_SM_H

/*20150705*/
#define MLCOALSIM "mlcoalsim version 2.0.0 (20200220)\n"

#define FREQSPECTRUM_TDFSR2	0
/*1 means print S,tajimaD,ZnA,Fs,R2,theta,pi,fwh and freq_spectrum; 2 means freq_spectrum. ONLY WORKS IN OPTION print_neuttest 3, otherwise CRASH*/
#define SHOW_PRIORS	1

/* define the values: SEXRATIOA1 = -1 indicates Ne is in relation 3N. SEXRATIOX1 =-1 indicates Ne is in relation 4N. One MUST be defined */
/* SEXRATIOA1 =  1 && SEXRATIOX1 = -1. 4N and in relation to sexratio=1 */
/* SEXRATIOA1 =  0 && SEXRATIOX1 = -1. 4N and in relation to whatever sexratio ...*/
/* SEXRATIOA1 = -1 && SEXRATIOX1 =  1. 3N and in relation to sexratio=1 */
/* SEXRATIOA1 = -1 && SEXRATIOX1 =  0. 3N and in relation to whatever sexratio ...*/

/* Montse data, Sexratio and data based on X reference*/
/*
#define SEXRATIOA1	-1
#define SEXRATIOX1	0
*/

/**/
#define SEXRATIOA1	0
#define SEXRATIOX1	-1
/**/
#define _CRT_SECURE_NO_DEPRECATE

#ifndef inMPI
	#define inMPI	0
#endif
#ifndef ZNS_ACTIVE
	#define ZNS_ACTIVE	0
#endif

#define MAXREJECT 100000
/*maximum rejection: if accepted values are less than 1/MAXREJECT, jump to next iteration and print 'na'*/
#define SHOWPROGRESS 1
/*set to 0 when the user do not want to see the progress of simulations in stdout */

#define NOBS_STATISTICS 47
#define NEUTVALUES2     47

#define ADDITIONAL_PRINTS 12
/*In case sfixall_thetas == 1 or rmfix == 1*/
/*0: No additional prints. 1: print theta/rec observed. 2: print treelength.*/
/* 12: combinations*/

#define TIMEDURATION 0

#define iEStest 1
#define iHStest 1
#define MAXHAP1 0

#define REFNUMBER 1230000
#define REFNUMBER1 (REFNUMBER+1)
#define PLUSPRIOR 10000

#define NVARSIZESIMPLE 50

struct var
{
    int pr_matrix;
    int mhits;
    double T_out;

    long int  n_iter;
    long int *seed1;
    long int seed2;
    int n_loci;
    double *ratio_sv;
    int linked;
    long int **loci_linked;
    long int despl;
    long int window;

    int *nsam;/*12*/
    double *r;
    long int *nsites;
    double *f;
    double *track_len;
    double *theta_1;
    long int *mutations;

    long int npop;
    int **ssize_pop;
    double mig_rate;
    double **mig_rate_matrix;

    int *ifselection;
    double *pop_sel;/*24*/
    double *sinit;
    double *sendt;
    double *sfreqend;
    double *sfreqinit;
    double *sel_nt;
    double pop_size;

    double *factor_pop;
    int ran_factorpop;
    int same_factorpop;
    int *npop_sampled;
    int max_npop_sampled;

    int neutral_tests;
    int print_neuttest;
    int sfix_allthetas;
    int sfix_pop;
    int mc_jump;

    int range_thetant;/*38*/
    double *thetant_min;
    double *thetant_max;

    int ifgamma;
    double *p_gamma;
    double *alpha_gamma;
    double *correct_gamma;

    int iflogistic;
    double *ts;

    int *nintn;
    double **nrec;
    double **tpast;

    int split_pop;/*50*/
    double time_split;
    double time_scoal;
    double factor_anc;
    double *freq;

    double tlimit;
    int npoprefugia;

    double *factor_chrn;
    double sex_ratio;
    double *event_sexratio;

    int rmfix;
    int rmfix_pop;
    int method_samp;
    int *Rm;
    int *nhapl;

    int range_rnt;/*63*/
    double *recnt_min;
    double *recnt_max;

    int ifgammar;
    double *p_gammar;
    double *alpha_gammar;
    double *correct_gammar;/*69*/

    int no_rec_males;

    int linked_segsites;
    int linked_rm;
    int linked_nhapl;

    int linked_segsites_nregion;
    int linked_rm_nregion;
    int linked_nhapl_nregion;

    double *heter_theta_alphag;
    long int *invariable_mut_sites;
    double *heter_rm_alphag;

    /*observed statistics*/
    int obs_statistics_used[NOBS_STATISTICS];
    double **obs_statistics[NOBS_STATISTICS];

    int likelihood_line;
    double likelihood_error[NOBS_STATISTICS];

    int npriors;

    int npop_events;
    double **pop_event;
    double *fixmax_abstime_event;

    int *ancestral_pol;

    double patcg[5];
    int includeEHHrel;
    long int *ehh_fixnt;

    int pop_outgroup;

    int missing;
};

struct var2
{
    int rmfix;
    int rmfix_pop;
    int Rm;
    int nhapl;

    int range_rnt;
    double recnt_min;
    double recnt_max;

    int ifgammar;
    double p_gammar;
    double alpha_gammar;
    double correct_gammar;

    int range_thetant;
    double thetant_min;
    double thetant_max;

    int ifgamma;
    double p_gamma;
    double alpha_gamma;
    double correct_gamma;

    int Sfix_alltheta;
    int Sfix_pop;
    int mc_jump;

    long int howmany;
    int nsam;
    double r;
    double f;
    double track_len;
    long int nsites;
    double theta;
    long int segsitesin;
    long int npop;
    int *config;
    double migrate;
    double **migrate_matrix;

    int linked;
    long int **loci_linked;
    long int despl;
    long int window;

    int nloci;
    int tloci;
    int pr_matrix;
    int mhits;
    double ratio_sv;
    double T_out;
    int Sout;

    double *factor_pop;
    int ran_factorpop;
    int same_factorpop;
    int npop_sampled;
    int max_npop_sampled;

    int neutral_tests;
    int print_neuttest;

    int ifselection;
    double pop_size;
    double pop_sel;
    double sel_nt;
    double sinit;
    double sendt;
    double sfreqend;
    double sfreqinit;

    int iflogistic;
    double *ts;

    int *nintn;
    double **nrec;
    double **tpast;
    double **tpastS;

    int split_pop;
    double time_split;
    double time_scoal;
    double time_scoalS;
    double factor_anc;
    double *freq;

    double tlimit;

    double factor_chrn;
    double sex_ratio;
    int event_forsexratio;
    double event_sexratio;
    int no_rec_males;

    int linked_segsites;
    int linked_rm;
    int linked_nhapl;

    int linked_segsites_nregion;
    int linked_rm_nregion;
    int linked_nhapl_nregion;

    double heter_theta_alphag;
    long int invariable_mut_sites;
    double heter_rm_alphag;

    int likelihood_line;

    int npriors;
    double **pointtoprior;

    int npop_events;
    struct events *pop_event;
    double *fixmax_abstime_event;

    int type_ancestral;
    int **ancestral_pol;

    double patcg[4];
    int *locspec_prior;
    int includeEHHrel;
    int ifehh_fix;
    long int ehh_margin;
    long int ehh_fixnt;

    int pop_outgroup;

    int missing;
};

struct prob_par {/*posterior probabilities for theta and rec*/
    double thetap;
    double recp;
    double Ttotp;
    double prob;
};

struct var_priors {
    char name_pr[30];/*name of the parameter used*/
    double idprior; /*number of the id prior in data*/
    int kind_var; /*(0)i: long int, (1)d: double*/
    int kind_dist; /*(0)u: uniform, (1)l: log-uniform, (2)g: gamma, (3)b: beta*/
    double dist_par[5]; /*u: 2, l: 2, g: 3, b: 4*/
    long int seed_prior; /*seed for pseudo-random numbers*/
    char prior_file[1000]; /*file with the prior values directly given*/
    double *priordist;
};

struct events { /*events occured, split, mig and factor_pop changes*/
    int number;/*id, sorted*/
    int npop1;/*the new pop has the same name than the name of the smaller id*/
    int npop2;/*if npop1=npop2 means only changes in mig and/or factor are used, no split*/
    double factor_pop;
    double time_event;
    double timeS_event;
    double *mig_fw;/*mig with the rest*/
    double *mig_rv;/*mi of the rest with the new pop*/
    double sex_ratio; /*NOW IS FOR ALL POPS*/
};

/* Els arbres estan fets en nodes (contenen el temps de coalescencia i el nombre del node amb el que conecta per dalt).
 * Tambe tenim els segments-length que contenen l'inici nt(eg 10-500, es 10) del segment, la direccio al nodes (arbre) al qual esta conectat,
 * i el nombre del següent segment (eg seria el segment que conte des del 501). Hi ha un segl per cada fragment produit per recombinacio que te un arbre associat.
 */
struct node {
    int abv;
    int ndes;
    double time;
};
struct segl {
    long int beg;
    struct node *ptree;
    long int next;
};

struct dnapar {
    double k;
    long int S;
    int B1;
    int Q1;
    int *freq;
    double piw;
    double pib;
    long int *unic;
    int maxhapl;
    int maxhapl1;
    int Rm;
    double thetaL;

    double withinw;
    double max_iES;
    double min_uiHS;
    double max_uiHS;

    double *fstallcomp;
    double *piwallcomp;
    double *piaallcomp;

    int nhapl;
    int *fhapl;
    double *fsthapallcomp;
    double *hapwallcomp;
    double *hapaallcomp;

    int *Sanc; /*for 3 pops is: Sx1,Sx2,Sxo,Sf1,Sf2,Sfo,Sx1f2,Sx2f1,Ssh,Sso*/
    int mhsites;

    double pie1;	/*Achaz*/
    double pin1;
    long int Se1;
    long int Sn1;

    double m_sdev;
    double m_skew;
    double m_kurt;
    double ragg;
};

#endif //MLCOALSIM_MLSP_SM_H
