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

#include "mlsp_sm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include "zf_log.h"

#if inMPI==1 
	#include <mpi.h>
#endif

int n_var = 242;
char var_file[242][30] =
{
    {"mhits"},			/* 1 allow multiple hits, 0 not */

    {"n_iterations"},	/* 'howmany' iterations for each simulation */
    {"seed1"},			/* seed for ran1: one value for each locus*/
    {"n_loci"},			/* number of loci */

    {"n_samples"},		/* samples in each locus */
    {"recombination"},	/* recombination in each locus(4Nor)/nt */
    {"n_sites"},		/* number of nt in each locus */
    {"f"},				/* gene conversion  vector*/
    {"track_len"},		/* gene conversion vector*/
    {"thetaw"},			/* theta for each locus/nt, in case speciation, split or subdivision, only for the first group (4Nou) */
    {"mutations"},		/* (segsitesin in case no mhits), this is for the complete tree in each locus */
    {"npop"},			/* the number of pops. The same in each locus */
    {"ssize_pop"},		/* sample size of each population, comma and the next locus */
    {"mig_rate"},		/* rate of migration (en 4Nom). This parameter is only valid now for refugia module. For normal module, use 'mig_rate_matrix'.*/
    {"nlinked_loci"},	/* if linked, indicate the number of linked loci. If 1 means sliding windows */
    {"pos_linked"},		/* the initial position and the end of the fragment for all fragments, then "comma" and next linked loci*/
    {"print_matrixpol"},/* 0 not print, 1 DNA (FASTA), 2 print Hudson modified output(mhits is different),3 dna excluding mhits*/
    {"mhratio_sv"},		/* in case mhits and theta. vector with the ratio vs transitions and transversions */
    {"seed2"},			/* NOW INVALIDATED: Seed2 for generating random numbers */
    
    {"factor_pop"},		/* vector with the Ne proportion of each population */
    {"ran_factorpop"},	/* 0 no random, 1 yes, calculated every iteration */
    {"same_factorpop"},	/* 0 no constant, 1 same Ne for each population */
    {"neutral_tests"},	/* 0 neutral tests not calculated, 1 they are calculated */
    {"print_neuttest"},	/* 0 no print neut tests, 1: print average and variance from total loci, 2: all vales each locus. 3: show all statistic per locus, not avg/var/P-values. 4: show all. Exceptions with linked loci. 3/4 are not possible. 5:Output for DnaSP (special case of type 1)*/
    {"npop_sampled"},  	/* Number of populations sampled, next locus by commas.*/

    { "theta"},			/* REPEATED: theta for each loci/nt, in case speciation, split or subdivision, only for the first group (4Nou) */
    { "sfix_allthetas"},/*Rejection algorithm for Sfix and all thetas (0/1)*/
    { "mc_jump"},		/*separation among chosen values in the markov chain in Sfix_allthetas and rmfix*/
	
						/*Following Braverman et al. 1995:*/
    {"ifselection"},	/*if selection 1, if not 0. Put 0 or 1 for each locus. Default is zero*/
    {"pop_size"},		/*effective population size: N. (Same for all loci)*/
    {"pop_sel"},		/*4Ns for each locus*/
    {"sel_nt"},			/*nt position of the selected mutation. 0 is the left value of the sequence. negative is allowed!!*/
    {"sinit"},			/*time in 4N generations the selection process ended. negative value means sel is not finished..*/
						
						/*limits for theta in case Sfix_all_theta:*/
    {"range_thetant"},	/*Put limits to theta/nt (in case Sfix_all theta) 1 uniform, 2 log10-uniform*/
    {"thetant_min"},	/*theta/nt min. if only one value, all are the same for aeach locus, if not, specify the values for each locus (vector)*/
    {"thetant_max"},	/*theta/nt max. if only one value, all are the same for aeach locus, if not, specify the values for each locus (vector)*/
    
    {"dist_out"},		/* in mhits. Time to the ancestral pop in 4No gen */
    {"displ"},			/*in case sliding windows, displacement*/
    {"window"},			/*in case liding windows, size of the window*/

    { "nintn"},			/* change of the population size. Number of intervals, comma and the next population*/
    { "nrec"},			/* change of the population size. Initial size (en No). First must be 1.0 , comma and the next pop*/
    { "npast"},			/* change of the population size. Final size (en No), comma and the next pop */
    { "tpast"},			/* change of the population size. Duration (en 4No) time is cummulative, starting from the first event with the min value (tpast[0])*/

    {"refugia"},		/*1 if the module refugia 0 if not. Only the parameters specified in this module are then considered (no expansion..etc, no nvents, etc*/
    {"time_split"},		/*time in 4N generations the population split from present to past in n refuges*/
    {"time_scoal"},		/*time (from time_split) the population coalesced in a single pop*/
    {"factor_anc"},		/*factor (relative to N) of the pop size of each refuge */
    {"freq_refugia"},	/*vector of n values (refuges): when time_split, freq of alleles going to each of the refuges*/
    
    {"tlimit"},			/* limit of the trees with recombination*/
	{"npoprefugia"},	/*number of refuges*/

	{"iflogistic"},		/*in case use logistic growth, if not, instantaneous*/
	{"ts"},				/*the first step in a growth process can start (for example) in the middle of the exponential phase*/
	
	{"ifgamma"},		/*in case using a gamma distribution for theta values (1), else 0*/
	{"p_gamma"},		/*p-value for gamma (also named lambda). if only one value, all are the same for aeach locus, if not, specify the values for each locus (vector)*/
	{"alpha_gamma"},	/*alpha value for gamma. if only one value, all are the same for aeach locus, if not, specify the values for each locus (vector)*/
	{"factorn_chr"},	/*value to weight for chrom X, Y, mithocondrial and autosomes..etc.: 1 autosomes, 0.75 is X, 0.25 is Y, mith -0.25 (considering only females reproduce mith) (if one value, all loci are equal)*/

	{"rmfix"},			/*uncertainty in recombination 0/1*/
	
	{"rm"},				/*values of Rm for each locus*/
	{"method_samp"},	/*methodology for sampling trees, in both Sfix_allthetas and rmfix: 1 rejection algorithm*/
	
	{"range_rnt"},		/*in case using a uniform distribution for all loci 1 unform 2 log10-uniform*/
	{"recnt_min"},		/*minimum value of recombination per nt. if only one value, all are the same for aeach locus, if not, specify the values for each locus (vector)*/
	{"recnt_max"},		/*maximum value of recombination per nt. if only one value, all are the same for aeach locus, if not, specify the values for each locus (vector)*/
	
	{"ifgammar"},		/*in case using a gamma distribution*/
	{"alpha_gammar"},	/*alpha value. if only one value, all are the same for aeach locus, if not, specify the values for each locus (vector)*/
	{"p_gammar"},		/*p,also named lambda. if only one value, all are the same for aeach locus, if not, specify the values for each locus (vector)*/
	
	{"no_rec_males"},	/*no recombination in males 1:Drosophila, 0 "rest"*/
	{"nhapl"},			/*number of hapl for each locus*/

	{"linked_fixsegsites"}, /*in case linked loci, indicate the number of segregating sites. They will be fixed in simulations (for a specific region and population)*/
	{"linked_fixrm"},		/*in case linked loci, indicate the number of min. rec. events. They will be fixed in simulations (for a specific region and population)*/
	
	{"heter_theta_alphag"}, /*alpha=p(lambda) parameter from gamma to indicate the heterogeneity in theta values across the locus, (vector) for each independent locus, or 1 value for the entire linked region*/
	{"heter_rec_alphag"},   /*alpha=p(lambda) parameter from gamma to indicate the heterogeneity in recombination values across the locus, (vector) for each independent locus, or 1 value for the entire linked region*/
	
	{"linked_fixnhapl"},	/*in case linked loci, indicate the number of haplotypes. They will be fixed in simulations (for a specific region and population)*/
	{"correct_gamma"},		/*use in case ifgamma = 1. the result of the gamma distribution is multiplied by this value*/
	{"correct_gammar"},		/*use in case ifgammar = 1. the result of the gamma distribution is multiplied by this value*/
	{"invar_mut_sites"},	/*number of positions that are considered invariable. The algorithm will calculate (as an average) the number of invariable positions for each iteration.*/
	
	{"likelihood_line"},	/*In case 1, the likelihood value for each simulation will be shown given the observed data and the interval +-*/
	
	{"sfix_pop"},			/*in case sfix_all_thetas and several populations. indicate the population where the S will be fixed (starting from pop 0 to n-1)*/
	{"rmfix_pop"},			/*in case rmfix and several populations. indicate the population where the S will be fixed (starting from pop 0 to n-1)*/
	
	{"mig_rate_matrix"},	/*in case several pops, indicate the migration rate (4Nm) for each popi vs popj. Indicate the entire matrix, with 0 in the diagonal. Indicate for all pops (not only the pops sampled)*/
	
	{"linked_nregion_fixsegsites"},	/*in case fixing seg sites in liked loci, indicate the region (starting from 0) where the positions will be fixed*/
	{"linked_nregion_fixrm"},		/*in case fixing Rm in liked loci, indicate the region (starting from 0) where the positions will be fixed*/
	{"linked_nhapl_fixnregion"},	/*in case fixing nhapl in liked loci, indicate the region (starting from 0) where the positions will be fixed*/
	
	{"nprior"},				/*indicate the number of the prior (starting with 1). (The same number preceeded by & has to be found in any of the evolutionary parameters included in the file)*/
	{"prior_type_var"},		/*0 integer, 1 is double*/
	{"prior_type_dist"},	/*0 uniform, 1 log10-uniform, 2 gamma, 3 beta*/
	{"prior_parameters"},	/*for 0 and 1 two values (min and max), for 2 three values (alpha, p, and factor_corr), for 3 four values (alpha, beta, min and max)*/
	{"prior_seed"},			/*pseudorandom number for generating the distribution*/
	{"prior_file"},			/*file with prior numbers*/

	{"npop_events"},		/*in order to generate events (changing mig_rates, coalesce pops, changing pop_sizes...) indicate the number of events*/
	{"pop_event"},			/*matrix that contains 'npop_events' columns separated by commas, each with 5 values separated by spaces/tabs + 2 vectors: first the number of the event (starting from 0!), */
							/*second and third, the populations implicated in the event (if only one, type the same number twice, if second is negative, split pop)*/
							/*Next, the the pop_size of the 'new' pop (relative to N). Warning: the 'new' pop must be the number of the smallest pop number when two pops are implicated*/
							/*Next, the time (in 4N gen) the event occurred. Note that this value is relative to the other events occurred before. First, the time is counted from events from 'nint', then, all the npop_events in sucesive order*/
							/*Finally, indicate the migration vector of the new pop(i) i->rest and the vector rest->i in the correct order, with n values (all n population from starting simulation, and writing 0 to the extinct pop and the diagonal)*/
							/*in case split, two events have to be defined: the first event needs pop1 from split and pop2 a negative value (the new pop), factor pop for the current pop and the time of split, and migration vectors for the current pop*/
							/*in case split (cont) second event: pop1=pop2 the new pop, factor pop for the new pop, time (not counted) and the migration vectors for the new pop*/
	
	{"ancestral_pol"},		/*First value: n>1 (n is the number of pops used, the last is the outgroup). For n=3, the values are: Second: number of populations in group 0 ns0: indicate the number of the populations (ns0 values). Third: number of populations in group 1 (ns1): indicate the number of the populations (ns1 values). Fourth: number of populations in outgroup (nso): indicate the number of the populations (nso values).*/
	{"patcg"},				/*if included: proportion of ATCG (from 0 to 1) sum must be 1. Only used in case print_matrixpol is 1*/
	{"include_ehh_rel"},	/*calculate EHH related statistics or not. By default is zero because is very slow.*/
	{"pop_outgroup"},		/*if >= 0 indicates the population defined as outgroup, by defect is -1*/
	{"sehh_fixnt"},			/*The first value indicates if define (1 or not 0). the second value indicates the range to see the closest polymorphism (eg., 100 means +-100 bp from the value indicated) then indicate indicates the nt position ehh (and related statitics) must be used, for each locus. For a windows,it is ony calculated in one vindow*/
	{"pevent_sexratio"},	/*indicate the pop_event number where the factor_chrn changes (e.g., change in male/female ratio) and the values of factor_chrn (two values for each event, no more changes than events)*/
	{"sex_ratio"},			/*ratio homogametic/heterogametic*/
	{"fixmax_abstime_event"},/*in some events it is necessary to fix a maximum absolute time the event occurred*/
	
    {"sendt"}, /*time at the selective event that selection finish (from present to past)*/
    {"sfreqend"}, /*frequence of the selective position at which the selection finish (start from past to present)*/
    
    {"sfreqinit"}, /*frequence of the selective position at which the selection start (finish from past to present)*/
	/*{"ascert_bias"},		TO DO, indicate the number of presamples for each population were analyzed (if 0, not pre-sampled), comma, and the next loci */
	/*{"ascbias_strategy"},	TO DO, 0: no asc bias; 1: asc bias, the entire presample must be polymorphic; 2: asc bias, each presample in each pop must be polymorphic*/
	/*{"mising_values"},	TO DO, 0/1 Consider missing values and not eliminate the column. Only available in the option print_neuttest=3 */

	/*{"include_trees"},	if 1, include the size of lengths of the branches for each population and segment: with option 1, only this file is shown, if 0, branches are not displayed. This is difficult to do...*/
	
    /*
     {"ancestral_pol_active"},
     {"fstn_files_active"},
     {"fsth_files_active"},
     */

    /*observed values: indicate the value for each pop in locus 0, sep by spaces (-10000 means na), comma and the next locus*/
	{"td_obs"},
	{"fs_obs"},
	{"fdn_obs"},
	{"ffn_obs"},
	{"fd_obs"},
	{"ff_obs"},
	{"h_obs"},
	{"b_obs"},
	{"q_obs"},
	{"za_obs"},
	{"fst_obs"},
	{"kw_obs"},
	{"hw_obs"},
	{"r2_obs"},
	{"s_obs"},
	{"pi_w_obs"},
	{"pi_b_obs"},
	{"thetawat_obs"},
	{"thetataj_obs"},
	{"thetafw_obs"},
	{"d_dmin_obs"},
	{"hnorm_obs"},
	{"maxhap_obs"},
	{"maxhap1_obs"},
	{"rm_obs"},
	{"thetafl_obs"},
	{"thetal_obs"},
	{"zenge_obs"},
	{"tew_obs"},
	{"fstw_obs"},
	{"min_uihs_obs"},
	{"max_uihs_obs"},
	{"max_lies_obs"},
	{"fixoutg_obs"},
	{"koutgjc_obs"},
	{"thetawa_obs"},
	{"thetawan_obs"},
	{"thetata_obs"},
	{"thetatan_obs"},
	{"achy_obs"},
	{"achyn_obs"},
	{"msdev_obs"},
	{"mskew_obs"},
	{"mkurt_obs"},
	{"ragg_obs"},
	{"zns_obs"},
	{"zz_obs"},
	/* indicate the interval (+/- the value indicated) to consider the simulated value equal to the observed*/
	{"td_err"},
	{"fs_err"},
	{"fdn_err"},
	{"ffn_err"},
	{"fd_err"},
	{"ff_err"},
	{"h_err"},
	{"b_err"},
	{"q_err"},
	{"za_err"},
	{"fst_err"},
	{"kw_err"},
	{"hw_err"},
	{"r2_err"},
	{"s_err"},
	{"pi_w_err"},
	{"pi_b_err"},
	{"thetawat_err"},
	{"thetataj_err"},
	{"thetafw_err"},
	{"d_dmin_err"},
	{"hnorm_err"},
	{"maxhap_err"},
	{"maxhap1_err"},
	{"rm_err"},
	{"thetafl_err"},
	{"thetal_err"},
	{"zenge_err"},
	{"tew_err"},
	{"fstw_err"},
	{"min_uihs_err"},
	{"max_uihs_err"},
	{"max_lies_err"},
	{"fixoutg_err"},
	{"koutgjc_err"},
	{"thetawa_err"},
	{"thetawan_err"},
	{"thetata_err"},
	{"thetatan_err"},
	{"achy_err"},
	{"achyn_err"},
	{"msdev_err"},
	{"mskew_err"},
	{"mkurt_err"},
	{"ragg_err"},
	{"zns_err"},
	{"zz_err"},
	/*indicate with 1 or 0 if the statistic is considered to calculate the P-values or not*/
	{"td_active"},
	{"fs_active"},
	{"fdn_active"},
	{"ffn_active"},
	{"fd_active"},
	{"ff_active"},
	{"h_active"},
	{"b_active"},
	{"q_active"},
	{"za_active"},
	{"fst_active"},
	{"kw_active"},
	{"hw_active"},
	{"r2_active"},
	{"s_active"},
	{"pi_w_active"},
	{"pi_b_active"},
	{"thetawat_active"},
	{"thetataj_active"},
	{"thetafw_active"},
	{"d_dmin_active"},
	{"hnorm_active"},
	{"maxhap_active"},
	{"maxhap1_active"},
	{"rm_active"},
	{"thetafl_active"},
	{"thetal_active"},
	{"zenge_active"},
	{"tew_active"},
	{"fstw_active"},
	{"min_uihs_active"},
	{"max_uihs_active"},
	{"max_lies_obs_active"},
	{"fixoutg_active"},
	{"koutgjc_active"},
	{"thetawa_active"},
	{"thetawan_active"},
	{"thetata_active"},
	{"thetatan_active"},
	{"achy_active"},
	{"achyn_active"},
	{"msdev_active"},
	{"mskew_active"},
	{"mkurt_active"},
	{"ragg_active"},
	{"zns_active"},
	{"zz_active"},
};

void input_data( FILE *file_input,struct var **data,struct var_priors **priors)
{
    char *f;
	int npopt;
    int c;
    char name_var[30];
    char number[100];
    int v,w,x,y,z,count_pri,ap;
	int totcp=0;
    int nam,numb_1,numb_2/*,val*/;
	int namf=0;
	char name_file[1000];
    double **numb_c;
    /*double end,interval,st;*/
	int maxn=19;
    
    void init_seed1(long);
    double ran1(void);
	
	int totalnloci;
	int defnumber;
	int *defpop;
	int windows;
	
	memset(name_var,0,30*sizeof(char));

    if(!(*data = (struct var *)calloc(1,sizeof(struct var)))) perror("calloc error.0");;
    
    if(!((*data)->nsam = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.1");
    if(!((*data)->r = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.2");
    if(!((*data)->nsites = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.4");    
    if(!((*data)->loci_linked = (long int **)calloc((unsigned)20,sizeof(long int *)))) perror("calloc error.4b");
    for(x=0;x<20;x++) 
        if(!((*data)->loci_linked[x] = (long int *)calloc((unsigned)20,sizeof(long int))))
            perror("calloc error.4c");
    if(!((*data)->f = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.5");
    if(!((*data)->track_len = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.7");
    if(!((*data)->theta_1 = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.9");    

    if(!((*data)->mutations = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.10");
    if(!((*data)->ssize_pop = (int **)calloc((unsigned)20,sizeof(int *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->ssize_pop[x] = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.11b");
    if(!((*data)->factor_pop = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->ratio_sv = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
	if(!((*data)->npop_sampled = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.34");
    
    if(!((*data)->ifselection = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.34");
    if(!((*data)->pop_sel = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sel_nt = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sinit = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sendt = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sfreqend = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sfreqinit = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    /*defaults*/
    (*data)->ifselection[0] = 1;
    (*data)->pop_sel[0] = 1;
    (*data)->sel_nt[0] = 1;
    (*data)->sinit[0] = 1;
    (*data)->sendt[0] = 1;
    (*data)->sfreqend[0] = 1;
    (*data)->sfreqinit[0] = 1;
    for(x=1;x<20;x++) {
        (*data)->ifselection[x] = 0;
        (*data)->pop_sel[x] = 1E6;
        (*data)->sel_nt[x] = -10000;
        (*data)->sinit[x] = 0;
        (*data)->sendt[x] = 1E6;
        (*data)->sfreqend[x] = 0.0;
        (*data)->sfreqinit[x] = 1.0;
    }
	if(!((*data)->ts = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->nintn = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.88");	
	if(!((*data)->nrec = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->nrec[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.11b");
    if(!((*data)->tpast = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->tpast[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.11b");

    if(!((*data)->freq = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.16");
    
	if(!((*data)->p_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
    if(!((*data)->alpha_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
    if(!((*data)->correct_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
    if(!((*data)->factor_chrn = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 58");
    
	if(!((*data)->p_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
	if(!((*data)->alpha_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
    if(!((*data)->correct_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
	if(!((*data)->Rm = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.56");
	if(!((*data)->nhapl = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.56");

	if(!((*data)->heter_theta_alphag = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.91");
	if(!((*data)->heter_rm_alphag = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->invariable_mut_sites = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.9");

	if(!((*data)->thetant_min = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->thetant_max = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");

	if(!((*data)->recnt_min = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->recnt_max = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");

	if(!((*data)->mig_rate_matrix = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->mig_rate_matrix[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.11b");

	if(!((*data)->pop_event = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->pop_event[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.11b");
	
	if(!(*priors = (struct var_priors *)calloc((unsigned)1,sizeof(struct var_priors)))) perror("calloc error.110");
	count_pri = 0;
	
 	if(!((*data)->ancestral_pol = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.561");
 	if(!((*data)->ehh_fixnt = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.562");
 	if(!((*data)->event_sexratio = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.564");

 	if(!((*data)->seed1 = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.563");

	if(!((*data)->fixmax_abstime_event = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.881");
	
	(*data)->pr_matrix = 1;
    (*data)->mhits = 0;
    (*data)->n_iter = 0;
    (*data)->n_loci = 1;
    (*data)->linked = 0;
    (*data)->despl = 0;
    (*data)->window = 0;
    (*data)->npop = 1;
    (*data)->mig_rate = 0.0;
    
    (*data)->seed1[0] = 1;  (*data)->seed1[1] = 12345678;
    (*data)->seed2 = 546373;

    (*data)->ran_factorpop = 0;
    (*data)->same_factorpop = 0;
    (*data)->neutral_tests = 0;
    (*data)->print_neuttest = 0;

    (*data)->sfix_allthetas = 0;
    (*data)->sfix_pop = -1;
    (*data)->mc_jump = 0;
    (*data)->pop_size = 1E+06;
    
    (*data)->range_thetant = 0;
    
    (*data)->T_out = 0.0;
       
    (*data)->split_pop = 0;
    (*data)->time_split = 0.;
    (*data)->time_scoal = (double)1E07;
	
	(*data)->pop_event[0][0] = 1;
	(*data)->pop_event[0][1] = 0;
    (*data)->factor_anc = 1.;

 	(*data)->npop_events = 0;

	(*data)->tlimit = 1000.;
    
	(*data)->iflogistic = 0;

    (*data)->ifgamma = 0;

    (*data)->rmfix = 0;
	(*data)->rmfix_pop = -1;
    (*data)->range_rnt = 0;
    (*data)->ifgammar = 0;
    (*data)->no_rec_males = 0;

	(*data)->linked_segsites = -1;
	(*data)->linked_rm = -1;
	(*data)->linked_nhapl = 0;
	(*data)->linked_segsites_nregion = -1;
	(*data)->linked_rm_nregion = -1;
	(*data)->linked_nhapl_nregion = -1;
	
	(*data)->includeEHHrel = 0;	
	(*data)->pop_outgroup = -1;
	(*data)->sex_ratio = 1.0;

	/*observed values*/
	for(x=0;x<NOBS_STATISTICS;x++) {
		if(!((*data)->obs_statistics[x] = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.113");
		for(y=0;y<20;y++) {
			if(!((*data)->obs_statistics[x][y] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.113b");
			(*data)->obs_statistics[x][y][0] = (double)0;/*redundant*/
		}
		(*data)->likelihood_error[x] = (double)1e-6;
		(*data)->obs_statistics_used[x] = 0;
	}

	(*data)->likelihood_line = 0;
    (*data)->npriors = 0;

    if(!(f = (char *)malloc(BUFSIZ))) perror("malloc error.1");
    setbuf(file_input,f);
    
    x = 0;
    c = 1;
    
    if(!(numb_c = (double **)malloc(20*sizeof(double *)))) perror("malloc error.2");
    for(x=0;x<20;x++) if(!(numb_c[x] = (double *)malloc(20*sizeof(double)))) perror("realloc error.1");

    while(c > 0) {
        
        /* look for the name of the variable */
        nam = 0;
        while(!nam) {
            if(c == 34) { /* 34 is '"' */
                c = fgetc(file_input);
                while((c > 0) && (c != 34))
                    c = fgetc(file_input);
            }
            if(c <= 0) break;
            if(c != 34) {/*Modificat el 8.5.03 */
                if((c < 91 && c > 64) || (c < 123 && c > 96) || c == 95) {/* if letters or '_' */
                    x=0;
                    while(c && ((c < 91 && c > 64) || (c < 123 && c > 96) || c == 95 || (c > 47 && c < 58))) {
                        name_var[x] = c;	/* assign c to name_var */
                        x++;
						if(x==29) break;
                        c = fgetc(file_input);
                        if(c <= 0) break;
						if(c == 58) break;
                    }
                    name_var[x] = '\0';
                    for(y=0;y<x;y++) if(name_var[y] < 91 && name_var[y] > 64) name_var[y] += 32;/* do lowercase */
                    nam = 1;		/* we have new name_var */
                    if(c <= 0) break;
                }
                else {
                    /*if(c < 58 && c > 47) {
                        puts("error in input. Numbers with no variable"); 
                        exit(1);
                    }*/
                    c = fgetc(file_input);
                }
                if(c <= 0) break;
                
                if(name_var[0] == 101 && x == 1) {	/* in case we have E+- */
                    puts("error in input. Numbers with no variable"); 
                    exit(1);
                }
            }
            else c = fgetc(file_input);/*Modificat el 8.5.03 */
		}
        /*Assign the variable*/
        for(w=0;w<n_var;w++)
            if(strcmp(name_var,var_file[w]) == 0)
                break;
		
		if(w < n_var) {
			/* look for the content of the variable */                
			numb_1 = x = 0;
			numb_2 = 1;
			/*if(!(numb_c[0] = realloc(numb_c[0],2*sizeof(double)))) perror("realloc error.1");*/
			numb_c[numb_1][0] = 0;
		
			while(nam && w!= 87) {
                if(c == 34) { /* 34 is '"' */
                    c = fgetc(file_input);
                    while((c > 0) && (c != 34))
                        c = fgetc(file_input);
                }
				while((c < 58 && c > 42) || c == 69 || c == 101 || c == 38) {
					if((c < 58 && c > 47) || c == 69 || c == 101 || c == 43 || c == 46) { 	/* numbers, E, dec, +, also & for priors*/
						number[x] = c;
						x++;
					}
					if(c == '&'/*38*/) {		/* for prior, value always separated by spaces */
						 if(x && number[x-1] != '\t' && number[x-1] != 32 && number[x-1] != ','/*44*/) {
							 puts("Error in input. '&[number_prior]' value must be separated by spaces, tabs or comma"); 
							 exit(1);
						 }
						 else {
							number[x] = c;
							x++;
						 }
					}
					if(c == '-') {		/* for a minus, but also for first-end */
                        number[x] = c;
                        x++;
					}
					if(c == ',') {/* do a new vector */
						if(x) {
							if(number[0] == '&') {
								number[x] = '\0';
								/*new code for enabling priors*/
								/*id number for the prior*/
								if(!(number[1] < 58 && number[1] > 42)) {
									printf("Error: prior mark '&' needs an identifier number linked after &.\n");
									exit(1);
								}
								if(count_pri < atof(number+1)) {
									if(!(*priors = (struct var_priors *)realloc(*priors,(unsigned)((int)(double)atof(number+1))*sizeof(struct var_priors))))
										perror("realloc error.9");
									for(y=count_pri;y<(int)(double)atof(number+1);y++) {
										priors[0][y].idprior = 0;
										priors[0][y].prior_file[0] = '\0';
									}
									count_pri = (int)(double)atof(number+1);
									(*data)->npriors = count_pri;
								}
								if(priors[0][count_pri-1].idprior) {
									printf("Error: prior %d is already defined in another parameter.\n",(int)(double)atof(number+1));
									exit(1);
								}
								memcpy(priors[0][(int)(double)atof(number+1)-1].name_pr,var_file[w],30);
								priors[0][(int)(double)atof(number+1)-1].idprior = atof(number+1);
								numb_c[numb_1][numb_2] = REFNUMBER + (double)priors[0][(int)(double)atof(number+1)-1].idprior;
								numb_c[numb_1][0] ++;
							}
							else {
								number[x] = '\0';
								numb_c[numb_1][numb_2] = atof(number);
								numb_c[numb_1][0] ++;
							}
						}
						numb_2 = 1;
						numb_1++;
						if(numb_1 >= 20) {
							if(!(numb_c = realloc(numb_c,(numb_1+1)*sizeof(double *)))) perror("realloc error.5");
							if(!(numb_c[numb_1] = (double *) malloc(20 * sizeof(double)))) perror("realloc error.6");
						}
						numb_c[numb_1][0] = 0;
						x = 0;
					}
					c = fgetc(file_input);
					if(c == 34) {
						c = fgetc(file_input);
						while((c > 0) && (c != 34))
							c = fgetc(file_input);
					}
					if(!((c < 58 && c > 42) || c == 69 || c == 101 || c == 38)) {
						if((c < 91 && c > 64) || (c < 123 && c > 96) || c == 95) {	/* if letters, stop */
							nam = 1;
							break;
						}
						else {
							if(x) {
								number[x] = '\0';
								if(number[0] == '&') {
									/*new code for enabling priors*/
									/*id number for the prior*/
									if(!(number[1] < 58 && number[1] > 42)) {
										printf("Error: prior mark '&' needs an identifier number linked after &.\n");
										exit(1);
									}
									/*count_pri++;*/
									if(count_pri < atof(number+1)) {
										if(!(*priors = (struct var_priors *)realloc(*priors,(unsigned)((int)(double)atof(number+1))*sizeof(struct var_priors))))
											perror("realloc error.9");
										for(y=count_pri;y<(int)(double)atof(number+1);y++) {
											priors[0][y].idprior = 0;
											priors[0][y].prior_file[0] = '\0';
										}
										count_pri = (int)(double)atof(number+1);
										(*data)->npriors = count_pri;
									}
									if(priors[0][count_pri-1].idprior) {
										printf("Error: prior %d is already defined in another parameter.\n",(int)(double)atof(number+1));
										exit(1);
									}
									memcpy(priors[0][(int)(double)atof(number+1)-1].name_pr,var_file[w],30);
									priors[0][(int)(double)atof(number+1)-1].idprior = atof(number+1);
									numb_c[numb_1][numb_2] = REFNUMBER + (double)priors[0][(int)(double)atof(number+1)-1].idprior;
								}
								else /**/
									numb_c[numb_1][numb_2] = atof(number);
								
								numb_c[numb_1][0]++;
								numb_2++;
								if(numb_2 >= 20) 
									if(!(numb_c[numb_1] = realloc(numb_c[numb_1],(numb_2+1)*sizeof(double)))) perror("realloc error.7");
								x = 0;
							}
						}
					}
				}
				if((c < 91 && c > 64) || (c < 123 && c > 96) || c == 95) {	/* if letters, stop */
					nam = 1;
					break;
				}
				if(c <= 0) break;
				c = fgetc(file_input);
			}
			if(maxn < numb_1) maxn = numb_1;
			
			while(nam && w==87) { /*collect the file name for each prior, if defined*/
				namf = 0;
				while(!namf) {
					if(c == 34) { /* 34 is '"' 9 is tab*/
						c = fgetc(file_input);
						while((c > 0) && (c != 34))
							c = fgetc(file_input);
					}
					if(c <= 0 || c == 10) break;
					if(c != 34) {
						while(c== 9) c = fgetc(file_input);
						x=0;
						while(c) {
							name_file[x] = c;	/* assign c to name_var */
							x++;
							if(x==999) break;
							c = fgetc(file_input);
							if(c <= 0 || c==10 || c<32) break;
						}
						name_file[x] = '\0';
						namf = x;		/* we have new name_file */
						if(c <= 0 || c==10 || c<32) break;
					}
					else c = fgetc(file_input);
				}
				break;
			}
			/* do assignments to the struct variables */
			switch(w) {
				case 0:
					(*data)->mhits = (int)numb_c[0][1];
					break;
				case 1:
					(*data)->n_iter = (long int)numb_c[0][1];
					break;
				case 2:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->seed1 = (long int *)realloc((*data)->seed1,(unsigned)(numb_c[0][0]+1)*sizeof(long int))))
							perror("realloc error.9o");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->seed1[z] = (long int)numb_c[0][z];
					break;
				case 3:
					(*data)->n_loci = (int)numb_c[0][1];
					break;
				case 4:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->nsam = (int *)realloc((*data)->nsam,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
							perror("realloc error.9");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->nsam[z] = (int)numb_c[0][z];
					break;
				case 5:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->r = (double *)realloc((*data)->r,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.9");
					for(z=0;z<numb_c[0][0]+1;z++)
						(*data)->r[z] = numb_c[0][z];
					break;
				case 6:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->nsites = (long int *)realloc((*data)->nsites,
											   (unsigned)(numb_c[0][0]+1)*sizeof(long int))))
						perror("realloc error.13");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->nsites[z] = (long int)numb_c[0][z];
					break;
				case 7:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->f = (double *)realloc((*data)->f,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.9");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->f[z] = numb_c[0][z];
					break;
				case 8:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->track_len = (double *)realloc((*data)->track_len,
											(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.9");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->track_len[z] = numb_c[0][z];
					break;
				case 9:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->theta_1 = (double *)realloc((*data)->theta_1,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.33");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->theta_1[z] = numb_c[0][z];
					break;
				case 10:
					if(numb_c[0][0] >= 20) if(!((*data)->mutations = (long int *)realloc((*data)->mutations,(unsigned)(numb_c[0][0]+1)*sizeof(long int)))) 
						perror("realloc error.21");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->mutations[z] = (long int)numb_c[0][z];
					break;
				case 11:
					(*data)->npop = (long int)numb_c[0][1];
					break;
				case 12:
					if(numb_1+1 != (*data)->n_loci && numb_1 != 0) {/*numb_1 is the number of dimensions. i.e. loci*/
						perror("Error: nloci must be first defined or ssize_pop is different from nloci. ");
						exit(1);
					}
					if(numb_1 >= 20) {
						if(!((*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,(unsigned)(numb_1+1)*sizeof(int *))))
							perror("realloc error.22");
					}
					for(z=0;z<numb_1+1;z++) {
						if(z>=20) {
							if(!((*data)->ssize_pop[z] = (int *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(int))))
								perror("realloc error.11");
						}
						else {
							if(numb_c[z][0] >= 20) 
								if(!((*data)->ssize_pop[z]=(int *)realloc((*data)->ssize_pop[z],(unsigned)(numb_c[z][0]+1)*sizeof(int))))
									perror("realloc error.22b");
						}
						for(y=0;y<numb_c[z][0]+1;y++) (*data)->ssize_pop[z][y] = (int)numb_c[z][y];
					}
					break;
				case 13:
					(*data)->mig_rate = numb_c[0][1];
					break;
				case 14:
					(*data)->linked = (int)numb_c[0][1];
					break;
				case 15:
					if(numb_1 >= 20) {
						if(!((*data)->loci_linked = (long int **)realloc((*data)->loci_linked,(unsigned)(numb_1+1)*sizeof(long int *))))
							perror("realloc error.38");
					}
					for(z=0;z<numb_1+1;z++) {
						if(z>=20) {
							if(!((*data)->loci_linked[z] = (long int *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(long int))))
								perror("realloc error.11");
						}
						else {
							if(numb_c[z][0] >= 20) 
							if(!((*data)->loci_linked[z] = (long int *)realloc((*data)->loci_linked[z],(unsigned)(numb_c[z][0]+1)* sizeof(long int))))
								perror("realloc error.39");
						}
						for(y=0;y<numb_c[z][0]+1;y++) (*data)->loci_linked[z][y] = (long int)numb_c[z][y];
					}
					break;
				case 16:
					(*data)->pr_matrix = (int)numb_c[0][1];
					break;
				case 17:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,
						(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.33");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->ratio_sv[z] = numb_c[0][z];
					break;
				case 18:
					(*data)->seed2 = (long int)numb_c[0][1];
					break;
				case 19:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->factor_pop = (double *)realloc((*data)->factor_pop,
						(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.58");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->factor_pop[z] = (double)numb_c[0][z];                
					break;
				case 20:
					(*data)->ran_factorpop = (int)numb_c[0][1];
					break;
				case 21:
					(*data)->same_factorpop = (int)numb_c[0][1];
					break;
				case 22:
					(*data)->neutral_tests = (int)numb_c[0][1];
					break;
				case 23:
					(*data)->print_neuttest = (int)numb_c[0][1];
					break;            
				case 24:
					if(numb_c[0][0] >= 20) {
						if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,(unsigned)(numb_c[0][0]+1)*sizeof(int)))) 
							perror("realloc error.58");
					}
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->npop_sampled[z] = (int)numb_c[0][z];                
					break;
				case 25: /*Exactly equal than case 9*/
					if(numb_c[0][0] >= 20) 
						if(!((*data)->theta_1 = (double *)realloc((*data)->theta_1,
						(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.33");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->theta_1[z] = numb_c[0][z];
					break;
				case 26:
					(*data)->sfix_allthetas = (int)numb_c[0][1];
					break;
				case 27:
					(*data)->mc_jump = (int)numb_c[0][1];
					break;
				case 28:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->ifselection = (int *)realloc((*data)->ifselection,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
							perror("realloc error.9");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->ifselection[z] = (int)numb_c[0][z];
					break;
				case 29:
					(*data)->pop_size = (double)numb_c[0][1];
					break;
				case 30:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->pop_sel = (double *)realloc((*data)->pop_sel,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.9");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->pop_sel[z] = (double)numb_c[0][z];
					break;
				case 31:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->sel_nt = (double *)realloc((*data)->sel_nt,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.9");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->sel_nt[z] = (double)numb_c[0][z];
					break;
				case 32:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->sinit = (double *)realloc((*data)->sinit,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.9");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->sinit[z] = (double)numb_c[0][z];
					break;
				case 33:
					(*data)->range_thetant = (int)numb_c[0][1];
					break;
				case 34:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->thetant_min = (double *)realloc((*data)->thetant_min,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.69");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->thetant_min[z] = (double)numb_c[0][z];
					break;
				case 35:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->thetant_max = (double *)realloc((*data)->thetant_max,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.69");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->thetant_max[z] = (double)numb_c[0][z];
					break;
				case 36:
					(*data)->T_out = numb_c[0][1];
					break;
				case 37:
					(*data)->despl = (long int)numb_c[0][1];
					break;
				case 38:
					(*data)->window = (long int)numb_c[0][1];
					break;
				case 39:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->nintn = (int *)realloc((*data)->nintn,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
							perror("realloc error.69");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->nintn[z] = (int)numb_c[0][z];
					break;
				case 40:
					if(numb_1 >= 20) {
						if(!((*data)->nrec = (double **)realloc((*data)->nrec,(unsigned)(numb_1+1)*sizeof(double *))))
							perror("realloc error.38");
					}
					for(z=0;z<numb_1+1;z++) {
						if(z>=20) {
							if(!((*data)->nrec[z] = (double *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(double))))
								perror("realloc error.11");
						}
						else {
							if(numb_c[z][0] >= 20) 
								if(!((*data)->nrec[z] = (double *)realloc((*data)->nrec[z],(unsigned)(numb_c[z][0]+1)* sizeof(double))))
									perror("realloc error.39");
						}
						for(y=0;y<numb_c[z][0]+1;y++) (*data)->nrec[z][y] = (double)numb_c[z][y];
					}
					break;
				case 41:
					printf("\n Warning: 'npast' is not enabled in this version.\n");
					break;
				case 42:
					if(numb_1 >= 20) {
						if(!((*data)->tpast = (double **)realloc((*data)->tpast,(unsigned)(numb_1+1)*sizeof(double *))))
							perror("realloc error.38");
					}
					for(z=0;z<numb_1+1;z++) {
						if(z>=20) {
							if(!((*data)->tpast[z] = (double *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(double))))
								perror("realloc error.11");
						}
						else {
							if(numb_c[z][0] >= 20) 
								if(!((*data)->tpast[z] = (double *)realloc((*data)->tpast[z],(unsigned)(numb_c[z][0]+1)* sizeof(double))))
									perror("realloc error.39");
						}
						for(y=0;y<numb_c[z][0]+1;y++) (*data)->tpast[z][y] = (double)numb_c[z][y];
					}
					break;
				case 43:
					(*data)->split_pop = (int)numb_c[0][1];
					break;
				case 44:
					(*data)->time_split = (double)numb_c[0][1];
					break;
				case 45:
					(*data)->time_scoal = (double)numb_c[0][1];
					break;
				case 46:
					(*data)->factor_anc = (double)numb_c[0][1];
					break;
				case 47:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->freq = (double *)realloc((*data)->freq,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.26");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->freq[z] = numb_c[0][z];
					break;
				case 48:
					(*data)->tlimit = numb_c[0][1];
					break;
				case 49:
					(*data)->npoprefugia = (int)numb_c[0][1];
					break;
				case 50:
					(*data)->iflogistic = (int)numb_c[0][1];
					break;
				case 51:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->ts = (double *)realloc((*data)->ts,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.269");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->ts[z] = numb_c[0][z];
					break;
				case 52:
					(*data)->ifgamma = (int)numb_c[0][1];
					break;
				case 53:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->p_gamma = (double *)realloc((*data)->p_gamma,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.53");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->p_gamma[z] = (double)numb_c[0][z];
					break;
				case 54:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->alpha_gamma = (double *)realloc((*data)->alpha_gamma,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.54");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->alpha_gamma[z] = (double)numb_c[0][z];
					break;
				case 55:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->factor_chrn = (double *)realloc((*data)->factor_chrn,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.55");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->factor_chrn[z] = (double)numb_c[0][z];
					break;
				case 56:
					(*data)->rmfix = (int)numb_c[0][1];
					break;
				case 57:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->Rm = (int *)realloc((*data)->Rm,(unsigned)(numb_c[0][0]+1)*sizeof(int)))) 
							perror("realloc error.57");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->Rm[z] = (int)numb_c[0][z];
					break;
				case 58:
					(*data)->method_samp = (int)numb_c[0][1];
					break;
				case 59:
					(*data)->range_rnt = (int)numb_c[0][1];
					break;
				case 60:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->recnt_min = (double *)realloc((*data)->recnt_min,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.69");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->recnt_min[z] = (double)numb_c[0][z];
					break;
				case 61:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->recnt_max = (double *)realloc((*data)->recnt_max,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.69");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->recnt_max[z] = (double)numb_c[0][z];
					break;
				case 62:
					(*data)->ifgammar = (int)numb_c[0][1];
					break;
				case 63:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->alpha_gammar = (double *)realloc((*data)->alpha_gammar,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.62");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->alpha_gammar[z] = (double)numb_c[0][z];
					break;
				case 64:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->p_gammar = (double *)realloc((*data)->p_gammar,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
							perror("realloc error.63");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->p_gammar[z] = (double)numb_c[0][z];
					break;
				case 65:
					(*data)->no_rec_males = (int)numb_c[0][1];
					break;
				case 66:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->nhapl = (int *)realloc((*data)->nhapl,(unsigned)(numb_c[0][0]+1)*sizeof(int)))) 
							perror("realloc error.57");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->nhapl[z] = (int)numb_c[0][z];
					break;
				case 67:
					(*data)->linked_segsites = (int)numb_c[0][1];
					break;
				case 68:
					(*data)->linked_rm = (int)numb_c[0][1];
					break;
				case 69:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->heter_theta_alphag = (double *)realloc((*data)->heter_theta_alphag,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.93");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->heter_theta_alphag[z] = (double)numb_c[0][z];
					break;
				case 70:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->heter_rm_alphag = (double *)realloc((*data)->heter_rm_alphag,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.93");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->heter_rm_alphag[z] = (double)numb_c[0][z];
					break;
				case 71:
					(*data)->linked_nhapl = (int)numb_c[0][1];
					break;
				case 72:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->correct_gamma = (double *)realloc((*data)->correct_gamma,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.97");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->correct_gamma[z] = (double)numb_c[0][z];
					break;
				case 73:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->correct_gammar = (double *)realloc((*data)->correct_gammar,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.98");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->correct_gammar[z] = (double)numb_c[0][z];
					break;
				case 74:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->invariable_mut_sites = (long int *)realloc((*data)->invariable_mut_sites,(unsigned)(numb_c[0][0]+1)*sizeof(long int))))
							perror("realloc error.98");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->invariable_mut_sites[z] = (long int)numb_c[0][z];
					break;
				case 75:
					(*data)->likelihood_line = (int)numb_c[0][1];
					break;
				case 76:
					(*data)->sfix_pop = (int)numb_c[0][1];
					break;
				case 77:
					(*data)->rmfix_pop = (int)numb_c[0][1];
					break;
				case 78:
					if(numb_1 >= 20) {
						if(!((*data)->mig_rate_matrix = (double **)realloc((*data)->mig_rate_matrix,(unsigned)(numb_1+1)*sizeof(double *))))
							perror("realloc error.22");
					}
					for(z=0;z<numb_1+1;z++) {
						if(z>=20) {
							if(!((*data)->mig_rate_matrix[z] = (double *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(double))))
								perror("realloc error.11");
						}
						else {
							if(numb_c[z][0] >= 20) 
								if(!((*data)->mig_rate_matrix[z]=(double *)realloc((*data)->mig_rate_matrix[z],(unsigned)(numb_c[z][0]+1)*sizeof(double))))
									perror("realloc error.22b");
						}
						for(y=0;y<numb_c[z][0]+1;y++) (*data)->mig_rate_matrix[z][y] = (double)numb_c[z][y];
					}
					break;
				case 79:
					(*data)->linked_segsites_nregion = (int)numb_c[0][1];
					break;
				case 80:
					(*data)->linked_rm_nregion = (int)numb_c[0][1];
					break;
				case 81:
					(*data)->linked_nhapl_nregion = (int)numb_c[0][1];
					break;
				case 82:
					for(z=0;z<count_pri;z++) {
						if(priors[0][z].idprior == (int)numb_c[0][1]) {
							totcp = z;
							break;
						}
					}
					if(totcp == count_pri) {
						printf("Error: the number of the nprior parameter must be defined previously\n ");
						exit(1);
					}
					break;
				case 83:
					priors[0][totcp].kind_var = (int)numb_c[0][1];
					break;
				case 84:
					priors[0][totcp].kind_dist = (int)numb_c[0][1];
					break;
				case 85:
					if(numb_c[0][0] > 4) {
						printf("Error: Excessive number of parameters in 'prior_parameters'\n");
						exit(1);
					}
					for(z=0;z<numb_c[0][0]+1;z++) priors[0][totcp].dist_par[z] = (double)numb_c[0][z];
					break;
				case 86:
					priors[0][totcp].seed_prior = (long int)numb_c[0][1];
					break;
				case 87:/*read the name of the file for priors*/
					memcpy(priors[0][totcp].prior_file,name_file,namf);
					break;
				case 88:
					(*data)->npop_events = (int)numb_c[0][1];
					break;
				case 89:
					if(numb_1 >= 20) {
						if(!((*data)->pop_event = (double **)realloc((*data)->pop_event,(unsigned)(numb_1+1)*sizeof(double *))))
							perror("realloc error.22");
					}
					for(z=0;z<numb_1+1;z++) {
						if(z>=20) {
							if(!((*data)->pop_event[z] = (double *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(double))))
								perror("realloc error.11");
						}
						else {
							if(numb_c[z][0] >= 20) 
								if(!((*data)->pop_event[z]=(double *)realloc((*data)->pop_event[z],(unsigned)(numb_c[z][0]+1)*sizeof(double))))
									perror("realloc error.22b");
						}
						for(y=0;y<numb_c[z][0]+1;y++) (*data)->pop_event[z][y] = (double)numb_c[z][y];
					}
					break;
				case 90:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->ancestral_pol = (int *)realloc((*data)->ancestral_pol,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
							perror("realloc error.98");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->ancestral_pol[z] = (int)numb_c[0][z];
					break;
				case 91:
					if(numb_c[0][0] != 4) {
						printf("\nError: Once patcg parameter defined it must contain the proportion of the four nt separated by spaces (only four values)\n");
						exit(1);
					}
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->patcg[z] = (double)numb_c[0][z];
				case 92:
					(*data)->includeEHHrel = (int)numb_c[0][1];
					break;
				case 93:
					(*data)->pop_outgroup = (int)numb_c[0][1];
					break;
				case 94:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->ehh_fixnt = (long int *)realloc((*data)->ehh_fixnt,(unsigned)(numb_c[0][0]+1)*sizeof(long int))))
							perror("realloc error.98");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->ehh_fixnt[z] = (long int)numb_c[0][z];
					break;
				case 95:
					if(numb_c[0][0] >= 20) 
						if(!((*data)->event_sexratio = (double *)realloc((*data)->event_sexratio,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
							perror("realloc error.99");
					for(z=0;z<numb_c[0][0]+1;z++) (*data)->event_sexratio[z] = (double)numb_c[0][z];
					break;					
				case 96:
					(*data)->sex_ratio = (double)numb_c[0][1];
					break;
                case 97:
                    if(numb_c[0][0] >= 20)
                        if(!((*data)->fixmax_abstime_event = (double *)realloc((*data)->fixmax_abstime_event,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                            perror("realloc error.982");
                    for(z=0;z<numb_c[0][0]+1;z++) (*data)->fixmax_abstime_event[z] = (double)numb_c[0][z];
                    break;
                case 98:
                    if(numb_c[0][0] >= 20)
                        if(!((*data)->sendt = (double *)realloc((*data)->sendt,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                            perror("realloc error.9842");
                    for(z=0;z<numb_c[0][0]+1;z++) (*data)->sendt[z] = (double)numb_c[0][z];
                    break;
                case 99:
                    if(numb_c[0][0] >= 20)
                        if(!((*data)->sfreqend = (double *)realloc((*data)->sfreqend,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                            perror("realloc error.9843");
                    for(z=0;z<numb_c[0][0]+1;z++) (*data)->sfreqend[z] = (double)numb_c[0][z];
                    break;
                case 100:
                    if(numb_c[0][0] >= 20)
                        if(!((*data)->sfreqinit = (double *)realloc((*data)->sfreqinit,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                            perror("realloc error.9843");
                    for(z=0;z<numb_c[0][0]+1;z++) (*data)->sfreqinit[z] = (double)numb_c[0][z];
                    break;
				default:
					/*observed values, and errors*/
					defnumber = 101;
					if(w >= defnumber && w < defnumber+NOBS_STATISTICS) {
						if(numb_1 >= 20) {
							if(!((*data)->obs_statistics[w-defnumber] = (double **)realloc((*data)->obs_statistics[w-defnumber],(unsigned)(numb_1+1)*sizeof(double *))))
								perror("realloc error.38");
						}
						for(z=0;z<numb_1+1;z++) {
							if(z>=20) {
								if(!((*data)->obs_statistics[w-defnumber][z] = (double *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(double))))
									perror("realloc error.11");
							}
							else {
								if(numb_c[z][0] >= 20) 
									if(!((*data)->obs_statistics[w-defnumber][z] = (double *)realloc((*data)->obs_statistics[w-defnumber][z],(unsigned)(numb_c[z][0]+1)* sizeof(double))))
										perror("realloc error.39");
							}
							for(y=0;y<numb_c[z][0]+1;y++) (*data)->obs_statistics[w-defnumber][z][y] = (double)numb_c[z][y];
						}
						break;
					}
					else if(w >= defnumber+NOBS_STATISTICS && w < defnumber+NOBS_STATISTICS*2) {
						(*data)->likelihood_error[w-(defnumber+NOBS_STATISTICS)] = (double)numb_c[0][1];
						break;
					}
					else if(w >= defnumber+2*NOBS_STATISTICS && w < defnumber+NOBS_STATISTICS*3) {
						(*data)->obs_statistics_used[w-(defnumber+2*NOBS_STATISTICS)] = (int)numb_c[0][1];
						break;
					}
					break;
			}
		}
		x = 0;
    }
    free(f);
    for(x=0;x<=maxn;x++) free(numb_c[x]);
	free(numb_c);
    
    /* INIT SEED FOR RANDOM VALUES */
    if((*data)->seed1[0] > 0) {
		if((*data)->seed1[0] != (*data)->n_loci) {
			printf("Error: not the same loci than defined in seed1 (it must define nloci seeds)\n");
			exit(1);
		}
		for(x=1;x<(*data)->n_loci+1;x++) {
			if((*data)->seed1[x] < 0) {
				printf("Error: seed1 must be between 0 and 2147483562");
				exit(1);
			}
		}
	}
	
	init_seed1((*data)->seed1[1]);
    /*srand((*data)->seed1);*/	/*INACTIVAT. emprem rands() per outgroup i mhits en especiacio*/
		
    /* FILTERS, AJ! VERY BAD !*/
    if((*data)->print_neuttest > 0) (*data)->neutral_tests = 1;

    if((*data)->mhits > 1) {
        printf("Error: mhits must be 0 or 1\n");
        exit(1);
    }
    if((*data)->T_out < 0.0) {
        printf("Error: dist_out must be >= 0\n");
        exit(1);
    }
    if((*data)->seed2 < 0) {
        printf("Error: seed2 must be between 0 and 2147483398\n");
        exit(1);
    }
    if((*data)->n_loci <  1) {
        printf("Error: n_loci must be at least 1\n");
        exit(1);
    }
    if((*data)->n_loci >  1 && (*data)->linked > 0) {
        printf("Error: When loci are linked, is not allowed more than 1 independent loci. \n");
        exit(1);
    }
    if((*data)->n_loci >  1 && (*data)->despl > 0. && (*data)->window > 0.) {
        printf("Error: When more than one loci are linked is not allowed sliding windows.\n ");
        exit(1);
    }
    if((*data)->linked !=  1 && (*data)->despl > 0. && (*data)->window > 0.) {
        printf("Error: Sliding windows is only allowed with linked value equal to 1. \n");
        exit(1);
    }
    if((*data)->includeEHHrel <  0 || (*data)->includeEHHrel > 1) {
        printf("Error: include_EHH_rel must be 0 or 1. \n");
        exit(1);
    }
    if((*data)->pop_outgroup < -1 || (*data)->pop_outgroup > (*data)->npop-1) {
        printf("Error: pop_outgroup must be -1 (no outgroup) or the number of the population defined as outgroup (from 0 to npop) \n");
        exit(1);
    }
    if((*data)->pop_outgroup >= -1 || (*data)->pop_outgroup < (*data)->npop) {
        (*data)->T_out = 0.0;
    }
	if((*data)->linked > 0) {
		if((*data)->loci_linked[(*data)->linked-1][2] >= (*data)->nsites[1]) {
			printf("Error: size of the linked loci exceeds the total size of the region (0-%ld)\n",(*data)->nsites[1]-1);
			exit(1);
		}
		for(x=0;x<(*data)->linked;x++) {
			if((*data)->loci_linked[x][2] < (*data)->loci_linked[x][1]) {
				printf("Error: size of the linked loci %d is less than zero\n",x);
				exit(1);
			}
			if(x) {
				if((*data)->loci_linked[x][1] < (*data)->loci_linked[x-1][2]) {
					printf("Error: linked locus not sorted\n");
					exit(1);
				}
				if((*data)->loci_linked[x][1] <= (*data)->loci_linked[x-1][1]) {
					printf("Error: linked locus have to be defined at different positions\n");
					exit(1);
				}
			}
		}
		if((*data)->linked > 1) {
			if((*data)->linked_segsites > 0) {
				if((*data)->linked_segsites_nregion > (*data)->linked || (*data)->linked_segsites_nregion < 0) {
						printf("Error: indicate the region to fix segsites in linked_nregion_segsites\n");
						exit(1);
				}
				if((*data)->mutations[0] == (*data)->n_loci && (*data)->mutations[1] != -1) {
						printf("Error: Only allowed linked_segsites or mutations (for the entire region). Not together\n");
						exit(1);
				}
			}
			if(((*data)->linked_segsites != (double)-1) &&
			   ((*data)->sfix_allthetas == 0)) {
				printf("Error: Not allowed to fix segsites in linked fragments wihthout the option 'Sfix_allthetas'.\n");
				exit(1);
			}
			
			if((*data)->linked_rm > -1) {
				if((*data)->linked_rm_nregion >= (*data)->linked) {
						printf("Error: linked_nregion_rm/_nhapl; define the region to fix rm values\n");
						exit(1);
				}
				if((*data)->Rm[0] == (*data)->n_loci && (*data)->Rm[0] != -1) {
						printf("Error: Only allowed linked_rm or rm (for the entire region). Not together\n");
						exit(1);
				}
                if((*data)->linked_nhapl > 0) {
                    if((*data)->linked_rm_nregion != (*data)->linked_nhapl_nregion) {
                        printf("Error: linked__nregion_rm and linked_nregion_nhapl must be in the same region.\n");
                        exit(1);
                    }
                    if((*data)->linked_nhapl_nregion == -1) {
                        (*data)->linked_nhapl_nregion = (*data)->linked_rm_nregion;
                    }
                }
			}
			if((*data)->linked_nhapl > 0) {
				if((*data)->linked_nhapl_nregion >= (*data)->linked) {
						printf("Error: linked_nregion_rm/_nhapl; define the region to fix nhapl values\n\n");
						exit(1);
				}
			}

			if((*data)->linked_rm != (double)-1 && (*data)->rmfix == 0) {
				printf("Error: Not allowed to fix Rm in linked fragments wihthout the option 'rmfix'.\n");
				exit(1);
			}
			if((*data)->linked_nhapl != (double)0 && (*data)->rmfix == 0) {
				printf("Error: Not allowed to fix the number of haplotypes in linked fragments wihthout the option 'rmfix'.\n");
				exit(1);
			}
		}
    }		
	
	windows = (*data)->n_loci;
	if((*data)->linked == 1 && (*data)->window && (*data)->despl) /*sliding windows*/
		windows = (int)ceil(((double)(*data)->nsites[1] - (double)(*data)->window) / ((double)(*data)->despl)) + 1;
	else if((*data)->linked > 1)/*separated but linked regions*/
		windows = (*data)->linked;
	

    if((*data)->nsam[0] != (*data)->n_loci) {
        printf("Error: n_samples needs n_loci inputs\n");
        exit(1);
    }
    if((*data)->mhits) {
        if(((*data)->theta_1[0] > 0 && (*data)->theta_1[1] > 0.) || (*data)->sfix_allthetas || (*data)->range_thetant || (*data)->ifgamma) {
            if((*data)->ratio_sv[0] < (*data)->n_loci) {
                if(!((*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
                    perror("realloc error.33");
                (*data)->ratio_sv[0] = (*data)->n_loci;
                if((*data)->ratio_sv[1] == 0. || (*data)->ratio_sv[1] == 0.5) (*data)->ratio_sv[1] = 0.5;
                for(x=2;x<(*data)->n_loci+1;x++) (*data)->ratio_sv[x] = (*data)->ratio_sv[1];
            }
            for(x=1;x<(*data)->n_loci+1;x++) {
                if((*data)->ratio_sv[x] < 0.) {
                    printf("Error: mhratio_sv must be higher than 0.\n");
                    exit(1);
                }
            }
        }
    }
    if((*data)->r[0] > 0 && !((*data)->r[0] == 1 && (*data)->r[1] <= 0.)) {
		if((*data)->r[0] == 0 || ((*data)->r[0] != (*data)->n_loci && (*data)->r[0] != 1)) {
			printf("Error: recombination needs n_loci inputs or a single input (all loci equal)\n");
			exit(1);
		}
    }
    if((*data)->f[0] > 0) {
		if((*data)->f[0] != (*data)->n_loci) {
			printf("Error: not the same loci than defined in f\n");
			exit(1);
		}
	}
	else {
		if(!((*data)->f = (double *)realloc((*data)->f,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.93");
		(*data)->f[0] = (double)(*data)->n_loci;
		for(x=1;x<(*data)->n_loci+1;x++) (*data)->f[x] = (double)0;
	}
	
    if((*data)->track_len[0] > 0) {
		if((*data)->track_len[0] != (*data)->n_loci) {
			printf("Error: not the same loci than defined in track_len\n");
			exit(1);
		}
	}
	else {
		if(!((*data)->track_len = (double *)realloc((*data)->track_len,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.93");
		(*data)->track_len[0] = (double)(*data)->n_loci;
		for(x=1;x<(*data)->n_loci+1;x++) (*data)->track_len[x] = (double)0;
	}
	
    
	if((*data)->heter_theta_alphag[0] > 0) {
		if((*data)->heter_theta_alphag[0] != (*data)->n_loci) {
				printf("Error: not the same loci than defined in heter_theta_alphag\n");
				exit(1);
		}
	}
	else {
		if(!((*data)->heter_theta_alphag = (double *)realloc((*data)->heter_theta_alphag,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.93");
		(*data)->heter_theta_alphag[0] = (double)(*data)->n_loci;
		for(x=1;x<(*data)->n_loci+1;x++) (*data)->heter_theta_alphag[x] = (double)-1;
	}
		
    if((*data)->invariable_mut_sites[0] > 0) {
		if((*data)->invariable_mut_sites[0] != (*data)->n_loci) {
				printf("Error: not the same loci than defined in invariable_mut_sites\n");
				exit(1);
		}
		else {
			for(x=1;x<(*data)->n_loci+1;x++) {
				if((*data)->invariable_mut_sites[x] >= (*data)->nsites[x]) {
					printf("Error: invariable_mut_sites must be smaller than n_sites.\n");
					exit(1);
				}
			}
		}
	}
	else {
		if(!((*data)->invariable_mut_sites = (long int *)realloc((*data)->invariable_mut_sites,(long int)((*data)->n_loci+1)*sizeof(long int)))) 
			perror("realloc error.99");
		(*data)->invariable_mut_sites[0] = (long int)(*data)->n_loci;
		for(x=1;x<(*data)->n_loci+1;x++) (*data)->invariable_mut_sites[x] = (long int)-1;
	}

    if((*data)->heter_rm_alphag[0] > 0) {
		if((*data)->heter_rm_alphag[0] != (*data)->n_loci) {
				printf("Error: not the same loci than defined in heter_rm_alphag\n");
				exit(1);
		}
	}
	else {
		if(!((*data)->heter_rm_alphag = (double *)realloc((*data)->heter_rm_alphag,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.93");
		(*data)->heter_rm_alphag[0] = (double)(*data)->n_loci;
		for(x=1;x<(*data)->n_loci+1;x++) (*data)->heter_rm_alphag[x] = (double)-1;
	}
    if((int)(*data)->nsites[0] != (*data)->n_loci) {
        printf("Error: nsites needs n_loci inputs\n");
        exit(1);
    }
    if((*data)->f[0] >0 && (*data)->f[1] >0) {
            if((*data)->f[0] == 0) {
                printf("Error: recombination input error\n");
                exit(1);
            }
    }
    if((*data)->track_len[0] && (*data)->track_len[1] > 0) {
            if((*data)->track_len[0] == 0) {
                printf("Error: recombination input error\n");
                exit(1);
            }
    }
    if((*data)->theta_1[0] > 0 && !((*data)->theta_1[0] == 1 && (*data)->theta_1[1] <= 0.)) {
            if((*data)->theta_1[0] == 0 || ((*data)->theta_1[0] != (*data)->n_loci && (*data)->theta_1[0] != 1)) {
                printf("Error: theta_1 needs n_loci inputs or a single input (all loci equal)\n");
                exit(1);
            }
    }
    if((*data)->mutations[0] > 0 && !((*data)->mutations[0] == 1 && (*data)->mutations[1] == -1)) {
        if((int)(*data)->mutations[0] != (*data)->n_loci) {
            printf("Error: mutations needs n_loci inputs\n");
            exit(1);
        }
    }
    if((*data)->mutations[0] == 0 || ((*data)->mutations[0] == 1 && (*data)->mutations[1] == -1)) {
       if(!((*data)->mutations = (long int *)realloc((*data)->mutations,((*data)->n_loci+1)*sizeof(long int)))) 
            perror("realloc error.38bdb");
		(*data)->mutations[0] = (*data)->n_loci;
		for(z=1;z<(*data)->n_loci+1;z++) (*data)->mutations[z] = -1;
    }
    if((*data)->split_pop == 1) {
        if((*data)->n_loci >= 20)
            if(!((*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,((*data)->n_loci+1)*sizeof(int *)))) 
                perror("realloc error.38bb");
        if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,((*data)->n_loci+1)*sizeof(int)))) 
            perror("realloc error.38bb");/*afegit21.4.2003*/
        (*data)->npop_sampled[0] = (*data)->n_loci;/*afegit21.4.2003*/
        for(z=0;z<(*data)->n_loci;z++) {
            if(z<20) {
                if(!((*data)->ssize_pop[z]=(int *)realloc((*data)->ssize_pop[z],2*sizeof(int)))) 
                    perror("realloc error.22bb");
            }
            else {
                if(!((*data)->ssize_pop[z]=(int *)malloc(2*sizeof(int)))) 
                    perror("realloc error.22bb");
            }
        }
        for(x=0;x<(*data)->n_loci;x++) {
            (*data)->ssize_pop[x][0] = 1;
            (*data)->ssize_pop[x][1] = (*data)->nsam[x+1];
            (*data)->npop_sampled[x+1] = 1;/*afegit 21.4.2003*/
        }
    }
    
    if((*data)->npop == 1 && (*data)->split_pop == 0) {
        if(windows >= 20)
            if(!((*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,(windows+1)*sizeof(int *)))) 
                perror("realloc error.38bb");
        if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,(windows+1)*sizeof(int)))) 
            perror("realloc error.38bb");/*afegit21.4.2003*/
        (*data)->npop_sampled[0] = (*data)->n_loci;/*afegit21.4.2003*/
        for(z=0;z<windows;z++) {
            if(z<20) {
                if(!((*data)->ssize_pop[z]=(int *)realloc((*data)->ssize_pop[z],2*sizeof(int)))) 
                    perror("realloc error.22bb");
            }
            else {
                if(!((*data)->ssize_pop[z]=(int *)malloc(2*sizeof(int)))) 
                    perror("realloc error.22bb");
            }
        }
        for(x=0;x<windows;x++) {
            (*data)->ssize_pop[x][0] = 1;
            (*data)->ssize_pop[x][1] = (*data)->nsam[x+1];
            (*data)->npop_sampled[x+1] = 1;/*afegit 21.4.2003*/
        }
        (*data)->ran_factorpop = 0;
        (*data)->same_factorpop = 1;
        (*data)->factor_pop[0] = 1;
        (*data)->factor_pop[1] = 1;
    }
    if((*data)->npop > 0  && (*data)->split_pop == 0) {
        if((*data)->npop > 1) {
            if((*data)->npop_sampled[0] == 1 && windows > 1) { 
                /*Si només hi ha 1 valor és que no s´han definit els altres loci*/
                if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,(windows+1)*sizeof(int))))
                    perror("realloc error.38bb");
                (*data)->npop_sampled[0] = windows; /*definim tots els loci*/
                for(z=1;z<windows;z++) (*data)->npop_sampled[z+1] = (*data)->npop_sampled[1];
                /*tots els loci tenen la mateixa mostra que el primer loci*/
            }
			if((*data)->ssize_pop[0][0] == 0) {
				/*if ssize is not included in input file*/
				(*data)->ssize_pop[0][0] = (*data)->npop_sampled[0];
				if(!((*data)->ssize_pop[0] = (int *)realloc((*data)->ssize_pop[0],(unsigned)((*data)->npop_sampled[0]+1)*sizeof(int)))) perror("calloc error.11b");
				(*data)->ssize_pop[0][1] = (*data)->nsam[0+1];
				for(x=2;x<=(*data)->ssize_pop[0][0];x++) {
					(*data)->ssize_pop[0][x] = 0;
				}
			}
			if(((*data)->ssize_pop[1][0] == 0 && windows > 1) || ((*data)->ssize_pop[0][0] < (*data)->npop_sampled[1])) {
				/*if only the first loci is defined or if ssize_pop is not entirey defined*/
				if(!((*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,(unsigned)windows*sizeof(int *)))) perror("calloc error.11");
				/*define size for ssize_pop per loci*/
				for(x=0;x<windows;x++) {
					if(windows < 20) {
						if(!((*data)->ssize_pop[x] = (int *)realloc((*data)->ssize_pop[x],(unsigned)((*data)->npop_sampled[x+1]+1)*sizeof(int)))) perror("calloc error.11b");
					}
					else {
						if(!((*data)->ssize_pop[x] = (int *)calloc((unsigned)(*data)->npop_sampled[x+1]+1,sizeof(int)))) perror("calloc error.11b");
					}	
					(*data)->ssize_pop[x][0] = (*data)->npop_sampled[x+1];
					(*data)->ssize_pop[x][1] = (*data)->nsam[x+1];
					for(z=2;z<=(*data)->ssize_pop[x][0];z++) {
						(*data)->ssize_pop[x][z] = 0;
					}
				}
			}			
        }
        else {
            if((*data)->npop == 1) {
                if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,(windows+1)*sizeof(int))))
                    perror("realloc error.38bb");
                for(x=0;x<windows;x++) { /*definir valors*/
                    if(!((*data)->ssize_pop[x]=(int *)realloc((*data)->ssize_pop[x],2*sizeof(int)))) 
                        perror("realloc error.22bb");
                    (*data)->ssize_pop[x][0] = 1;
                    (*data)->ssize_pop[x][1] = (*data)->nsam[x+1];
                    (*data)->npop_sampled[x+1] = 1;
                }
            }
        }
    }

	if((*data)->npop <= 0)  {
		printf("Error: npop must be 1 or higher\n");
		exit(1);
	}
	if((*data)->ts[0] > (*data)->npop) {
		printf("Error: ts can not have more than 'npop' values\n");
		exit(1);
	}
	
	if((*data)->npoprefugia > (*data)->npop) 
		npopt = (*data)->npoprefugia;
	else {
		if((*data)->ifselection[0]) npopt = (int)((*data)->npop > 2 ? (*data)->npop:2);
		else npopt = (int)(*data)->npop;
	}
	
	if((*data)->nintn[0] == 0) {
		if(npopt >= 20) {
			if(!((*data)->nintn=(int *)realloc((*data)->nintn,(npopt+1)*sizeof(int)))) 
				perror("realloc error.225b");
			if(!((*data)->ts=(double *)realloc((*data)->ts,(npopt+1)*sizeof(double)))) 
				perror("realloc error.22t5b");
			if(!((*data)->nrec=(double **)realloc((*data)->nrec,(npopt)*sizeof(double *)))) 
				perror("realloc error.226b");
			if(!((*data)->tpast=(double **)realloc((*data)->tpast,(npopt)*sizeof(double *)))) 
				perror("realloc error.227b");
			for(x=0;x<npopt;x++) {
				(*data)->nintn[x] = 0;
				(*data)->ts[x] = 0.;
				if(!((*data)->nrec[x]=(double *)calloc(20,sizeof(double)))) 
					perror("realloc error.226b");
				if(!((*data)->tpast[x]=(double *)calloc(20,sizeof(double)))) 
					perror("realloc error.227b");
			}
		}
	}
	(*data)->nintn[0] = npopt;
	(*data)->ts[0] = npopt;
	
	for(x=0;x<npopt;x++) {
		if((*data)->nintn[x+1] > 0 && (*data)->nintn[0] > x) {
			if((*data)->nrec[x][0] != (*data)->nintn[x+1]) {
				printf("Error: nrec has not the same intervals indicated in nintn\n");
				exit(1);
			}
			if((*data)->nrec[x][1] != 1.0) {
				printf("Error: nrec must be 1.0 in the first value of each pop. Change initial size for populations with 'factor_pop' \n (except pop 0 that must be 1)\n");
				exit(1);
			}
			if((*data)->tpast[x][0] != (*data)->nintn[x+1]) {
				printf("Error: tpast has not the same intervals indicated in nintn\n");
				exit(1);
			}
		}
	}
    if((*data)->ran_factorpop < 0 || (*data)->ran_factorpop > 1) {
        printf("Error ran_factorpop: it should be  0 or 1\n");
        exit(1);
    }
    if((*data)->same_factorpop < 0 || (*data)->same_factorpop > 1) {
        printf("Error ran_factorpop: it should be 0 or 1\n");
        exit(1);
    }
    if((*data)->ran_factorpop != 0 && (*data)->same_factorpop != 0) {
        printf("Error ran_factorpop and same_factorpop can not be activated at the same time\n");
        exit(1);
    }
    if((*data)->neutral_tests < 0 || (*data)->neutral_tests > 1) {
        printf("Error neutral_tests: it should be 0 or 1\n");
        exit(1);
    }
	if((*data)->sfix_allthetas == (*data)->rmfix) {
		if((*data)->sfix_pop != (*data)->rmfix_pop) {
            printf("Error: sfix_pop and rmfix_pop must use the same population\n");
            exit(1);
		}
	}
	if((*data)->sfix_allthetas) {
		if((*data)->sfix_pop == -1) {
            printf("Error: sfix_pop must be defined when sifx_allthetas is defined.\n");
            exit(1);
		}
	}
	if((*data)->rmfix && (*data)->npop > 1) {
		if((*data)->rmfix_pop == -1) {
            printf("Error: rmfix_pop must be defined when npop > 1 and rmfix is defined.\n");
            exit(1);
		}
	}
    
	for(x=0;x<windows;x++) {
        if(((*data)->npop_sampled[x+1] <= (*data)->sfix_pop && (*data)->sfix_pop != (*data)->npop) || ((*data)->npop_sampled[x+1] <= (*data)->rmfix_pop  && (*data)->rmfix_pop != (*data)->npop)) {
            printf("Error: sfix_pop/rmfix_pop must be a population (npop_sampled) included in the analysis, or the whole 'megapop' (named as the number of npops) \n");
            exit(1);
        }
    }
    for(x=0;x<(*data)->n_loci;x++) {
        y = 0;
        for(z=1;z<(*data)->ssize_pop[x][0]+1;z++) y += (*data)->ssize_pop[x][z];
        if((*data)->nsam[x+1] != y) {
            printf("Error: ssize_pop is not coincident with the total number of samples\n");
            exit(1);
        }
    }
    for(x=0;x<windows;x++) {
        if((*data)->npop_sampled[x+1] > 0 && (*data)->npop_sampled[x+1] != (*data)->ssize_pop[x][0]) {
            printf("Error: npop_sampled is not coincident with the #pops in ssize_pop\n");
            exit(1);
        }
    }
    /*vull tenir un nombre de poblacions mostrejades per cada locus i que tingui un nombre de mostres per població*/
    /*tots els loci han de tenir el mateix nombre de poblacions, pero el nombre de mostres poden ser diferents*/
    for(x=0;x<windows;x++)/*aquí es posen els valors de npop_sampled en cas no estigui definit...*/
        if((*data)->npop_sampled[x+1] == 0)
            (*data)->npop_sampled[x+1] = (*data)->ssize_pop[x][0];
    (*data)->npop_sampled[0] = windows;
    /*De fet, el que vull es que en les mostres, en alguns loci es fan servir poblacion diferent i hi han zeros pel mig */
    /*Això no ho se per què ho he definit, miro més endavant... ÉS IMPORTANT!! */
    if((*data)->ssize_pop[0][0] < (*data)->npop) {   
        for(x=0;x<windows;x++) {
            if(!((*data)->ssize_pop[x] = (int *)realloc((*data)->ssize_pop[x],((*data)->npop+1)*sizeof(int))))
                perror("realloc error.78");
            for(z=(*data)->ssize_pop[x][0]+1;z<(*data)->npop+1;z++) (*data)->ssize_pop[x][z] = 0;
            (*data)->ssize_pop[x][0] = (int)(*data)->npop;
        }
    }
    /**/
    if((*data)->factor_pop[0] > 0 && (*data)->factor_pop[1] > 0.0) {
        for(x=1;x<(*data)->factor_pop[0]+1;x++) {
            if((*data)->factor_pop[x] <= 0.0 ) {
                printf("Error factor_pop values: they should be higher than 0\n");
                exit(1);
            }
        }
        if((*data)->factor_pop[1] != 1.0 && (*data)->split_pop == 0) {
            printf("Error factor_pop vector: First value must be 1.0\n");
            exit(1);
        }
        /* popsizes: 1 or between 1 and 10 relative to the first pop if all are not indicated*/
        /**/
        if((*data)->factor_pop[0] < (*data)->npop) {
            if(!((*data)->factor_pop = (double *)realloc((*data)->factor_pop,((*data)->npop+1)*sizeof(double))))
                perror("realloc error.1");
            if((*data)->ran_factorpop == 0 && (*data)->same_factorpop == 0) z = (int)(*data)->factor_pop[0]+1;
            else z = 1;
            (*data)->factor_pop[0] = (*data)->npop;
            for(x=z;x<(*data)->npop+1;x++) {
                if((*data)->same_factorpop == 1 || ((*data)->ran_factorpop == 0 && (*data)->same_factorpop == 0)) (*data)->factor_pop[x] = 1.0;
                else {
                    if(x==1) (*data)->factor_pop[x] = 1.0;
                    else {
                        (*data)->factor_pop[x] = (double)ran1()*9.0 + 1.0; 
                        if((double)ran1() < 0.5)  (*data)->factor_pop[x] = 1./(*data)->factor_pop[x];
                    }
                }
            }
        }
    }
	
    if((*data)->theta_1[0] == 0 || ((*data)->theta_1[0] == 1 && (*data)->theta_1[1] <= (double)0)) {
		if((*data)->mutations[0] == 0 || ((*data)->mutations[0] == 1 && (*data)->mutations[1] == -1 && (int)(*data)->mutations[0] != (*data)->n_loci)) {
			if((*data)->range_thetant == 0 && (*data)->ifgamma == 0) {
                printf("Error: thetaw, mutations, range_thetant and/or ifgamma must be defined\n");
                exit(1);
			}
		}
	}
	
    if((*data)->r[0] == 0 || ((*data)->r[0] == 1 && (*data)->r[1] <= 0.  && (*data)->r[0] != (*data)->n_loci)) {
		if((*data)->range_rnt == 0 && (*data)->ifgammar == 0) {
			printf("Error: recombination, range_rnt or ifgammar must be defined\n");
			exit(1);
		}
	}
	if((*data)->pr_matrix == 1) {
		if((*data)->patcg[0] == 4) {
			if((*data)->patcg[1]+(*data)->patcg[2]+(*data)->patcg[3]+(*data)->patcg[4] != 1.0) {
				printf("\nError: patcg values must sum 1.0\n");
				exit(1);
			}	
		}
	}
    if((*data)->neutral_tests == 1) {
        if((*data)->pr_matrix != 0 && (*data)->print_neuttest > 2) {
            printf("Error: print_neuttest > 2 is not possible when pr_matrix is active\n");
            exit(1);
        }
        if((*data)->print_neuttest < 0 || (*data)->print_neuttest > 5) {
            printf("Error print_neuttest: it should be between 0 and 4\n");
            exit(1);
        }
    }
	if((*data)->neutral_tests == 0 && (*data)->likelihood_line == 1) {
		printf("Error: neutral_test must be active when likelihood_line is active.\n");
		exit(1);
	}
    if((*data)->sfix_allthetas > 1) {
		printf("Error: sfix_allthetas can only be 0 or 1.\n");
		exit(1);
    }
    if((*data)->method_samp > 2) {
		printf("Error: method_samp can only be 0, 1 or 2.\n");
		exit(1);
    }

    if((*data)->ifgamma == 1 && ((*data)->theta_1[0] > 0 && !((*data)->theta_1[0] > 0 && (*data)->theta_1[1] <= 0))) {
        printf("Error: 'theta' and 'ifgamma' can not be defined at the same time.\n");
        exit(1);
    }
    if(((*data)->range_thetant) && ((*data)->theta_1[0] > 0 && ((*data)->theta_1[0] != 1) && (*data)->theta_1[1] <= 0)) {
        printf("Error: 'theta' and 'range_thetant' can not be defined at the same time.\n");
        exit(1);
    }
    if((*data)->ifgammar == 1 && ((*data)->r[0] > 0 && !(((*data)->r[0] > 0) && (*data)->r[1] == 0))) {
        printf("Error: 'recombination' and 'ifgammar' can not be defined at the same time.\n");
        exit(1);
    }
    if(((*data)->range_rnt) && ((*data)->r[0] > 0 && ((*data)->r[0] != 1) && (*data)->r[1] == 0)) {
        printf("Error: 'recombination' and 'range_rnt' can not be defined at the same time.\n");
        exit(1);
    }

    if(((*data)->sfix_allthetas == 1 || (*data)->rmfix == 1) && (*data)->method_samp == 2) {
		if((*data)->mc_jump < 1){
            printf("Error: mc_jump must be a positive integer\n");
            exit(1);
        }
    }
    if((*data)->sfix_allthetas) {
        if((*data)->npop > 1 && ((*data)->sfix_pop > (*data)->npop || (*data)->sfix_pop < 0)) {
            printf("Error: sfix_pop must be in the range [0,'npop'] .\n");
            exit(1);
		}
		if((*data)->range_thetant != 0 && (*data)->range_thetant != 1 && (*data)->range_thetant != 2) {
            printf("Error: range_recnt must be 0, 1 or 2\n");
            exit(1);
        }
        else {
            if((*data)->range_thetant == 1 || (*data)->range_thetant == 2) {
                if((*data)->thetant_min[0] == (*data)->thetant_max[0] && (*data)->thetant_min[0] > 0) {
					for(x=1;x<=(*data)->thetant_min[0];x++) {
						if((*data)->thetant_min[x] < 0) {
							printf("Error: thetant_min must be positive or a zero value\n");
							exit(1);
						}
						if((*data)->thetant_max[x] < (*data)->thetant_min[x]) {
							printf("Error: thetant_max must be positive or a zero value, and lower or equal than thetant_min\n\n");
							exit(1);
						}
					}
				}
				else {
					printf("Error: thetant_min and thetant_max must be defined\n");
					exit(1);
				}
            }
			if((*data)->range_thetant == 0) {
				if((*data)->ifgamma == 0) {
					if((*data)->theta_1[0]) {
						printf("Error in input: sfix_allthetas need theta, range_thetant or ifgamma be defined. \n");
						exit(1);
					}
				}
			}
        }
    }
    if((*data)->rmfix > 0) {
        if((*data)->rmfix_pop > (*data)->npop || (*data)->rmfix_pop < 0) {
            printf("Error: rmfix_pop must be in the range [0,'npop') .\n");
            exit(1);
		}
	}
    if((*data)->ifselection[0] < (*data)->n_loci) { 
        /*Si només hi ha 1 valor és que no s´han definit els altres loci*/
        if(!((*data)->ifselection = (int *)realloc((*data)->ifselection,((*data)->n_loci+1)*sizeof(int))))
            perror("realloc error.38bb");
        for(z=(*data)->ifselection[0];z<(*data)->n_loci;z++) (*data)->ifselection[z+1] = 0;/*no selection*/
        (*data)->ifselection[0] = (*data)->n_loci; /*definim tots els loci*/
    }
    /*Si només hi ha 1 valor és que no s´han definit els altres loci*/
    if((*data)->split_pop == 1 && (*data)->freq[0] != (*data)->npoprefugia) {
        printf("Error: freq_refugia must be defined for all refugia\n");
        exit(1);
    }

    if((*data)->split_pop == 1 && (*data)->factor_pop[0] != (*data)->npoprefugia) {
        printf("Error: factor_pop must be defined for all refugia\n");
        exit(1);
    }

    if((*data)->split_pop == 0 && (*data)->npop > 1 && (*data)->factor_pop[0] != (*data)->npop) {
        printf("Error: factor_pop must be defined for all populations\n");
        exit(1);
    }
	if((*data)->ifgamma == 1) {
		if((*data)->p_gamma[0] != (double)(*data)->n_loci) {
			printf("Error: p gamma parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
		if((*data)->alpha_gamma[0] != (double)(*data)->n_loci) {
			printf("Error: alpha gamma parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
		if((*data)->correct_gamma[0] == (double)0) {
			if(!((*data)->correct_gamma = (double *)realloc((*data)->correct_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
				perror("realloc error.155");
			(*data)->correct_gamma[0] = (double)(*data)->n_loci;
			for(x=1;x<=(*data)->n_loci;x++) {
				(*data)->correct_gamma[x] = (double)1.0;
			}
		}
		if((*data)->correct_gamma[0] != (double)(*data)->n_loci) {
			printf("Error: correction gamma parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
	}
	for(x=1;x<=(*data)->n_loci;x++) {
		if(((*data)->ifgamma == 1 && (*data)->sfix_allthetas > 0) && ((*data)->p_gamma[x] <= (double)0.0 || (*data)->alpha_gamma[x] <= (double)0. || (*data)->correct_gamma[x] <= (double)0.)) {
			printf("Error: gamma parameters must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
	}
	
	if((*data)->factor_chrn[0] == (double)0) {
		if(!((*data)->factor_chrn = (double *)realloc((*data)->factor_chrn,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->factor_chrn[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->factor_chrn[x] = (double)1.0;
		}
	}
	if((*data)->factor_chrn[0] > 1 && (*data)->factor_chrn[0] < (*data)->n_loci) {
		if(!((*data)->factor_chrn = (double *)realloc((*data)->factor_chrn,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->factor_chrn[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->factor_chrn[x] = (double)1.0;
		}
		(*data)->factor_chrn[0] = (double)(*data)->n_loci;
	}
	
#if SEXRATIOA1 > -1
	for(x=1;x<=(int)(*data)->n_loci;x++) {
		if((*data)->factor_chrn[x] != (double)1.0 && (*data)->factor_chrn[x] != (double)0.75 && (*data)->factor_chrn[x] != (double)0.25 && (*data)->factor_chrn[x] != (double)-0.25) {
			printf("Error: factor_chr can only be 1.00 (Autosomes), 0.75 (X,Z), 0.25 (Y,W) or -0.25 (mithocondrial).\n");
			exit(1);
		}
	}
#else 
	for(x=1;x<=(int)(*data)->n_loci;x++) {
		if((*data)->factor_chrn[x] != (double)1.0 && (*data)->factor_chrn[x] != (double)1.33 && (*data)->factor_chrn[x] != (double)0.33 && (*data)->factor_chrn[x] != (double)-0.33) {
			printf("Error: factor_chr can only be 1.00 (Autosomes), 1.33 (X,Z), 0.33 (Y,W) or -0.33 (mithocondrial).\n");
			exit(1);
		}
	}
#endif

	if((*data)->ratio_sv[0] == (double)0) {
		if(!((*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->ratio_sv[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->ratio_sv[x] = (double)0.50;
		}
	}
	if((*data)->ratio_sv[0] < (*data)->n_loci) {
		if(!((*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->ratio_sv[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->ratio_sv[x] = (double)0.50;
		}
		(*data)->ratio_sv[0] = (double)(*data)->n_loci;
	}

	if((*data)->p_gamma[0] < (*data)->n_loci) {
		if(!((*data)->p_gamma = (double *)realloc((*data)->p_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->p_gamma[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->p_gamma[x] = (double)1.0;
		}
	}
	if((*data)->p_gammar[0] < (*data)->n_loci) {
		if(!((*data)->p_gammar = (double *)realloc((*data)->p_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->p_gammar[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->p_gammar[x] = (double)1.0;
		}
		(*data)->p_gammar[0] = (double)(*data)->n_loci;
	}

	if((*data)->alpha_gamma[0] < (*data)->n_loci) {
		if(!((*data)->alpha_gamma = (double *)realloc((*data)->alpha_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->alpha_gamma[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->alpha_gamma[x] = (double)1.0;
		}
	}
	if((*data)->alpha_gammar[0] < (*data)->n_loci) {
		if(!((*data)->alpha_gammar = (double *)realloc((*data)->alpha_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->alpha_gammar[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->alpha_gammar[x] = (double)1.0;
		}
		(*data)->alpha_gammar[0] = (double)(*data)->n_loci;
	}

 	if((*data)->correct_gamma[0] < (*data)->n_loci) {
		if(!((*data)->correct_gamma = (double *)realloc((*data)->correct_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->correct_gamma[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->correct_gamma[x] = (double)1.0;
		}
	}
	if((*data)->correct_gammar[0] < (*data)->n_loci) {
		if(!((*data)->correct_gammar = (double *)realloc((*data)->correct_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->correct_gammar[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->correct_gammar[x] = (double)1.0;
		}
		(*data)->correct_gammar[0] = (double)(*data)->n_loci;
	}
	if((*data)->mutations[0] > 0 && !((*data)->mutations[1] == -1)) {
		for(x=1;x<(int)(*data)->n_loci+1;x++) if((*data)->mutations[x] != 0) break;
		if(x < (*data)->n_loci) {
			if((*data)->sfix_allthetas == 0 && (*data)->mhits != 0) {
				printf("Error. mhits option is not allowed when theta is not defined (define thetaw or a distribution of theta)).\n");
				exit(1);
			}
		}
	}

    if((*data)->rmfix > 1) {
		printf("Error: rmfix can only be 0 or 1.\n");
		exit(1);
    }
    if((*data)->rmfix && (*data)->linked < 2) {
		if((*data)->Rm[0] != (*data)->n_loci) { 
			printf("Error: Once 'rmfix' is defined, 'Rm' must be defined for each locus.\n");
			exit(1);
		}
	}
	if((*data)->nhapl[0] == (int)0 && (*data)->rmfix) {
		if(!((*data)->nhapl = (int *)realloc((*data)->nhapl,(unsigned)((*data)->n_loci+1)*sizeof(int)))) 
			perror("realloc error.1558");
		(*data)->nhapl[0] = (int)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->nhapl[x] = (int)0;
		}
	}
    if((*data)->rmfix && (*data)->linked < 2) {
		if((int)(*data)->nhapl[0] != (int)(*data)->n_loci) { 
			printf("Error: Once 'rmfix' and 'nhapl' is defined, 'nhapl' must be defined for each locus.\n");
			exit(1);
		}
	}
    if((*data)->rmfix) {
        if((*data)->range_rnt != 0 && (*data)->range_rnt != 1 && (*data)->range_rnt != 2) {
            printf("Error: range_recnt must be 0, 1 or 2\n");
            exit(1);
        }
        else {
            if((*data)->range_rnt == 1 || (*data)->range_rnt == 2) {
                if((*data)->recnt_min[0] == (*data)->recnt_max[0] && (*data)->recnt_min[0] > 0) {
					for(x=1;x<=(*data)->recnt_min[0];x++) {
						if((*data)->recnt_min[x] < 0) {
							printf("Error: recnt_min must be positive or a zero value\n");
							exit(1);
						}
						if((*data)->recnt_max[x] < (*data)->recnt_min[x]) {
							printf("Error: recnt_max must be positive or a zero value, and lower or equal than recnt_min\n\n");
							exit(1);
						}
					}
				}
				else {
					printf("Error: recnt_max and recnt_max must be defined\n");
					exit(1);
				}
			}
			if((*data)->range_rnt == 0) {
				if((*data)->ifgammar == 0) {
					if((*data)->r[0] == 0) {
						printf("Error in input: options rmfix, range_rant and ifgammar. \n");
						exit(1);
					}
				}
			}
        }
    }

	if((*data)->ifgammar == 1) {
		if((*data)->p_gammar[0] != (double)(*data)->n_loci) {
			printf("Error: p gammar parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
		if((*data)->alpha_gammar[0] != (double)(*data)->n_loci) {
			printf("Error: alpha gammar parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
		if((*data)->correct_gammar[0] == (double)0) {
			if(!((*data)->correct_gammar = (double *)realloc((*data)->correct_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
				perror("realloc error.155");
			(*data)->correct_gammar[0] = (double)(*data)->n_loci;
			for(x=1;x<=(*data)->n_loci;x++) {
				(*data)->correct_gammar[x] = (double)1.0;
			}
		}
		if((*data)->correct_gammar[0] != (double)(*data)->n_loci) {
			printf("Error: correction gammar parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
	}
	for(x=1;x<=(*data)->n_loci;x++) {
		if(((*data)->ifgammar == 1 && (*data)->rmfix > 0) && ((*data)->p_gammar[x] <= (double)0.0 || (*data)->alpha_gammar[x] <= (double)0. || (*data)->correct_gammar[x] <= (double)0.)) {
			printf("Error: gammar parameters must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
	}

	if((*data)->p_gammar[0] == (double)0) {
		if(!((*data)->p_gammar = (double *)realloc((*data)->p_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.158");
		(*data)->p_gammar[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->p_gammar[x] = (double)1.0;
		}
	}
	if((*data)->p_gammar[0] < (*data)->n_loci) {
		if(!((*data)->p_gammar = (double *)realloc((*data)->p_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.159");
		for(x=(int)(*data)->p_gammar[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->p_gammar[x] = (double)1.0;
		}
		(*data)->p_gammar[0] = (double)(*data)->n_loci;
	}

	if((*data)->alpha_gammar[0] == (double)0) {
		if(!((*data)->alpha_gammar = (double *)realloc((*data)->alpha_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.159n");
		(*data)->alpha_gammar[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->alpha_gammar[x] = (double)1.0;
		}
	}
	if((*data)->alpha_gammar[0] < (*data)->n_loci) {
		if(!((*data)->alpha_gammar = (double *)realloc((*data)->alpha_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.159p");
		for(x=(int)(*data)->alpha_gammar[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->alpha_gammar[x] = (double)1.0;
		}
		(*data)->alpha_gammar[0] = (double)(*data)->n_loci;
	}
	
	if(((*data)->sfix_allthetas == 1 || (*data)->rmfix == 1) && (*data)->method_samp == 0) {
		printf("Error: when sfix_allthetas or rmfix is defined, method_samp must be 1 or 2.\n");
		exit(1);
	}

	if(!((*data)->no_rec_males == 1 || (*data)->no_rec_males == 0)) {
		printf("Error: no_rec_males must be 1 or 0.\n");
		exit(1);
	}

	/*observed data*/
    if((*data)->despl > 0 && (*data)->window > 0) {
        totalnloci = (int)ceil(((double)(*data)->nsites[1] - (double)(*data)->window) / ((double)(*data)->despl)) + (int)1;
	}	
    else {
        if((*data)->linked > 1) totalnloci = (*data)->linked;
        else totalnloci = (*data)->n_loci;
	}

	if((*data)->likelihood_line) {
		y = 0;
		for(x=0;x<NOBS_STATISTICS;x++) {
			if((*data)->obs_statistics_used[x] == (double)1) {
				y = 1;
				break;
			}
		}
		if(y==0) {
			printf("Error: At least one observed value must be defined when likelihood_line is active.\n");
			exit(1);
		}
	}
	
	if((*data)->npop_events) {
		for(x=0;x<(*data)->npop_events;x++) {
			if((int)(*data)->pop_event[x][0] != 5 + 2*(*data)->npop) {
				printf("Error: pop_event must include for each event 5 values (number_of_event pop1 pop2 factor_pop time_event) plus two vectors of migration values of npop values each (the row, pop->others, and the column, others->pop) starting with the 0 value and all separated by spaces\n");
				exit(1);
			}
			if((*data)->pop_event[x][1] >= (*data)->npop_events) {
				printf("Error: pop_event: the number_of_event can not be bigger than npop_events (starting from 0)\n");
				exit(1);
			}
			if((*data)->pop_event[x][2] >= (*data)->npop && (*data)->pop_event[x][2] < 0.) {
				printf("Error: pop_event pop1 value: the value must be between 0 and npop-1\n");
				exit(1);
			}
			if((*data)->pop_event[x][3] >= (*data)->npop && (*data)->pop_event[x][3] < 0.) {
				printf("Error: pop_event pop2 value: the value must be between 0 and npop-1\n");
				exit(1);
			}
			if((*data)->pop_event[x][4] <= 0.) {
				printf("Error: pop_event factor_pop value: the value must be higher than 0.\n");
				exit(1);
			}
			if((*data)->pop_event[x][5] < 0.) {
				printf("Error: pop_event time_event value: the value must be higher than or equal to 0.\n");
				exit(1);
			}
			if((*data)->pop_event[x][6+(int)(*data)->pop_event[x][2]] != 0.) {
				printf("Error: pop_event migration value: the migration for the population to itself must be zero\n");
				exit(1);
			}
			if((*data)->pop_event[x][6+(*data)->npop+(int)(*data)->pop_event[x][2]] != 0.) {
				printf("Error: pop_event migration value: the migration for the population to itself must be zero\n");
				exit(1);
			}
			for(y=0;y<2*(*data)->npop;y++) {
				if((*data)->pop_event[x][y+6] < 0.) {
					printf("Error: pop_event migration values must be positive or zero values\n");
					exit(1);
				}
			}
		}
		if((*data)->fixmax_abstime_event[0]) {
			for(y=1;y<=(*data)->fixmax_abstime_event[0];y++) {
				if((double)y/2.0 != (int)y/2) {
					if((*data)->fixmax_abstime_event[y] >= (*data)->npop_events) {
						printf("Error: fixmax_abstime_event values must be defined as: [number_event] [max abs time] ... [number_event] [max abs time]\n");
						exit(1);
					}
				}
			}
		}
	}
	
	/*check for migration matrix, npop x npop values*/
	if((*data)->npop > 1) {
		if((*data)->npop >= 20) {
			(*data)->mig_rate_matrix = (double **)realloc((*data)->mig_rate_matrix,((*data)->npop+1)*sizeof(double *));
			for(x=0;x<(*data)->npop;x++) {
				if(x<20)
					(*data)->mig_rate_matrix[x] = (double *)realloc((*data)->mig_rate_matrix[x],((*data)->npop+1)*sizeof(double));
				else
					(*data)->mig_rate_matrix[x] = (double *)calloc(((*data)->npop+1),sizeof(double));
			}
		}
		if((*data)->mig_rate_matrix[0][0] == 1) {
			/*all values in matrix are equal*/
            (*data)->mig_rate = (*data)->mig_rate_matrix[0][1];
            (*data)->mig_rate_matrix[0][0] = (*data)->mig_rate_matrix[0][1] = 0.;
		}
		else {
			if((*data)->mig_rate_matrix[0][0] > 1) {
				for(x=0;x<(*data)->npop;x++) {
					if((int)(*data)->mig_rate_matrix[x][0] != (*data)->npop) {
						printf("Error: mig_rate_matrix must have npop values for each row, comma, and again for npop times (npop x npop). Diagonal values are not considered (type 0)\n");
						exit(1);
					}
					if((*data)->mig_rate_matrix[x][x+1] != 0.) {
						printf("Error: mig_rate_matrix Diagonal values must be 0\n");
						exit(1);
					}
					for(y=1;y<=(*data)->npop;y++) {
						if((*data)->mig_rate_matrix[x][y] < 0.) {
							printf("Error: mig_rate_matrix must have positive or zero values\n");
							exit(1);
						}
					}
				}
			}
		}
	}

	/*filters for ancestral_pol*/
	/*first value is n (number of groups here defined)
	second value is the number of pops in group 1
	next values are the numbers of the pops for group 1
	next is the number of pops in the second group 
	... etc
	the last group indicates the group outgroup
	*/
	/*ancestral population is defined*//**/	
	if((*data)->ancestral_pol[0] > 1 && (*data)->ancestral_pol[1] > 1) {
		if((*data)->ancestral_pol[1] > (*data)->npop) {
			printf("Error: ancestral_pol: first value can not be higher than npop\n");
			exit(1);
		}
		if((*data)->ancestral_pol[0] > ((*data)->ancestral_pol[1] + (*data)->npop + 1)) {
			printf("Error: ancestral_pol: too many values defined\n");
			exit(1);
		}
		defpop = (int *)calloc(((*data)->ancestral_pol[0]-(*data)->ancestral_pol[1]-1),sizeof(int));
		
		/*look for repeated values*/
		x = 2;
		w = 0;
		v = 0;
		for(y=0;y<(*data)->ancestral_pol[1];y++) {
			x++;
			w += (*data)->ancestral_pol[x-1];
			ap = (*data)->ancestral_pol[x-1];
			for(z=0;z<ap;z++,x++) {
				defpop[v] = (*data)->ancestral_pol[x];
				v++;
			}
		}
		if(w > (*data)->npop) {
			printf("Error: ancestral_pol: Too much number of populations per group defined\n");
			exit(1);
		}
		for(x=0;x<((*data)->ancestral_pol[0]-(*data)->ancestral_pol[1]-2);x++) {
			for(y=x+1;y<((*data)->ancestral_pol[0]-(*data)->ancestral_pol[1]-1);y++) {
				if(defpop[x] == defpop[y]) {
					printf("Error: ancestral_pol: populations are repeated in the defined groups or wrong defined\n");
					exit(1);
				}
			}
		}
		free(defpop);
	}

	if((*data)->ehh_fixnt[0] > 0 || ((*data)->ehh_fixnt[0] > 2 && (*data)->ehh_fixnt[1] > 1)) {
	   printf("Error: ehh_fixnt: the first value must be 0 (inactive) or 1 (active)\n");
	   exit(1);
	}
	if((*data)->ehh_fixnt[0] > 0 || ((*data)->ehh_fixnt[0] > 1 && (*data)->ehh_fixnt[1] == 1)) {
		if((*data)->ehh_fixnt[2] < 0) {
			printf("Error: ehh_fixnt: the second value (the interval to find a polymorphism in the region) must be more or equal than 0\n");
			exit(1);
		}
		for(x=0;x<(*data)->n_loci;x++) {
			if((*data)->ehh_fixnt[x+3] < 0 || ((*data)->ehh_fixnt[x+3] >= (*data)->nsites[x+1])) {
				printf("Error: ehh_fixnt: the following valuse (third and rest) must be more or equal than 0 and less than the length of the loci\n");
				exit(1);
			}
		}
	}
	
    if((*data)->sex_ratio <= 0.) {
		printf("Error: sex_ratio must be higher than 0.\n");
		exit(1);
    }
	
	/*check observed values*/
	(*data)->max_npop_sampled = 1;
	for(x=0;x<windows;x++) {
		if((*data)->npop_sampled[x+1] > (*data)->max_npop_sampled) 
			(*data)->max_npop_sampled = (*data)->npop_sampled[x+1];
	}
	for(x=0;x<NOBS_STATISTICS;x++) {
		if((*data)->obs_statistics_used[x] == (double)0) {
			if(windows >= 20) {
				if(!((*data)->obs_statistics[x] = (double **)realloc((*data)->obs_statistics[x],(unsigned)(windows+1)*sizeof(double *)))) 
					perror("realloc error.obs_statistics");
			}
		}
		for(z=0;z<windows;z++) {
			if((*data)->obs_statistics_used[x] == (double)0) {
				if(z >= 20) {
					if(!((*data)->obs_statistics[x][z] = (double *)calloc((unsigned)((*data)->max_npop_sampled+1),sizeof(double)))) 
						perror("realloc error.obs_statistics");
				}
				else {
					if(!((*data)->obs_statistics[x][z] = (double *)realloc((*data)->obs_statistics[x][z],(unsigned)((*data)->max_npop_sampled+1)*sizeof(double)))) 
						perror("realloc error.obs_statistics");
				}
				(*data)->obs_statistics[x][z][0] = (double)(*data)->max_npop_sampled;
				for(y=1;y<=(*data)->max_npop_sampled;y++) (*data)->obs_statistics[x][z][y] = (double)-10000;
			}
			else {
				if((*data)->obs_statistics_used[x] == (double)1 && (*data)->obs_statistics[x][z][0] != (double)(*data)->max_npop_sampled) {
					switch(x) {
						case 0:
							printf("Error: TD observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 1:
							printf("Error: Fs observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 2:
							printf("Error: FDn observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 3:
							printf("Error: FFn observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 4:
							printf("Error: FD observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 5:
							printf("Error: FF observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 6:
							printf("Error: H observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 7:
							printf("Error: B observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 8:
							printf("Error: Q observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 9:
							printf("Error: ZA observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 10:
							printf("Error: Fst observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 11:
							printf("Error: Kw observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 12:
							printf("Error: Hw observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 13:
							printf("Error: R2 observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 14:
							printf("Error: Ssites observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 15:
							printf("Error: piw observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 16:
							printf("Error: pib observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 17:
							printf("Error: thetaWatt observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 18:
							printf("Error: thetaTaj observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 19:
							printf("Error: thetaFW observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 20:
							printf("Error: D_Dmin observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 21:
							printf("Error: H_norm observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 22:
							printf("Error: maxhap observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 23:
							printf("Error: maxhap1 observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 24:
							printf("Error: Rm observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 25:
							printf("Error: thetafl observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 26:
							printf("Error: thetal observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 27:
							printf("Error: zengE observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 28:
							printf("Error: EW observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 29:
							printf("Error: Fstw observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 30:
							printf("Error: min_uiHS observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 31:
							printf("Error: max_uiHS observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						case 32:
							printf("Error: max_liES observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
						default:
							printf("Error: observed values (if defined) must have the maximum number of samped populations values for all loci. \n Please include the values for the MAXIMUM number of sampled populations defined (if not exist, include the value -10000), comma, and the next loci. \n");
							exit(1);
							break;
					}
				}
				if((*data)->obs_statistics_used[x] != (double)0 && (*data)->obs_statistics_used[x] != (double)1) {
					switch(x) {
						case 0:
							printf("Error: TD_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 1:
							printf("Error: Fs _obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 2:
							printf("Error: FDn_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 3:
							printf("Error: FFn_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 4:
							printf("Error: FD_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 5:
							printf("Error: FF_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 6:
							printf("Error: H_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 7:
							printf("Error: B_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 8:
							printf("Error: Q_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 9:
							printf("Error: ZA_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 10:
							printf("Error: Fst_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 11:
							printf("Error: Kw_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 12:
							printf("Error: Hw_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 13:
							printf("Error: R2_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 14:
							printf("Error: Ssites_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 15:
							printf("Error: piw_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 16:
							printf("Error: pib_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 17:
							printf("Error: thetaWatt_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 18:
							printf("Error: thetaTaj_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 19:
							printf("Error: thetaFW_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 20:
							printf("Error: D_Dmin_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 21:
							printf("Error: H_norm_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 22:
							printf("Error: maxhap_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 23:
							printf("Error: maxhap1 observed values must have nloci values when is define\nd");
							exit(1);
							break;
						case 24:
							printf("Error: Rm_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 25:
							printf("Error: thetafl_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 26:
							printf("Error: thetal_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 27:
							printf("Error: zengE_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 28:
							printf("Error: EW_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 29:
							printf("Error: Fstw_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 30:
							printf("Error: min_uiHS_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 31:
							printf("Error: max_uiHS_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						case 32:
							printf("Error: max_liES_obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
						default:
							printf("Error: obs_active value must be 0 when is not defined and 1 when defined\n");
							exit(1);
							break;
					}
				}
				if((*data)->likelihood_error[x] <= (double)0 && (*data)->obs_statistics_used[x] == (double)1 && (*data)->likelihood_line) {
					switch(x) {
						case 0:
							printf("Error: TD_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 1:
							printf("Error: Fs_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 2:
							printf("Error: FDn_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 3:
							printf("Error: FFn_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 4:
							printf("Error: FD_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 5:
							printf("Error: FF_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 6:
							printf("Error: H_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 7:
							printf("Error: B_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 8:
							printf("Error: Q_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 9:
							printf("Error: ZA_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 10:
							printf("Error: Fst_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 11:
							printf("Error: Kw_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 12:
							printf("Error: Hw_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 13:
							printf("Error: R2_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 14:
							printf("Error: Ssites_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 15:
							printf("Error: piw_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 16:
							printf("Error: pib_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 17:
							printf("Error: thetaWatt_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 18:
							printf("Error: thetaTaj_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 19:
							printf("Error: thetaFW_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 20:
							printf("Error: D_Dmin_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 21:
							printf("Error: H_norm_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 22:
							printf("Error: maxhap_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 23:
							printf("Error: maxhap1 observed values must have nloci values when is define\nd");
							exit(1);
							break;
						case 24:
							printf("Error: Rm_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 25:
							printf("Error: thetafl_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 26:
							printf("Error: thetal_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 27:
							printf("Error: zengE_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 28:
							printf("Error: EW_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 29:
							printf("Error: Fstw_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 30:
							printf("Error: min_uiHS_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 31:
							printf("Error: max_uiHS_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						case 32:
							printf("Error: max_liES_err must be defined if likelihod_line is active\n");
							exit(1);
							break;
						default:
							printf("Error: err must be defined if likelihod_line is active\n");
							exit(1);
							break;
					}
				}
			}
		}
	}
}

void output_data(FILE *out, char *in, char *file_out, struct var **data,struct var_priors *priors)
{
    void print_var(int,FILE *);

    int x,y,z,zz,flag,npopt;
    double j;
    time_t now;
    struct tm *date;
    char s[80];
	int defnumber;
	int totalnloci;
	long int lx;
	FILE *outputind;
	char dfiles[420];
	char *el,nchar[6];
	div_t npr;
    
    time(&now);
    date = localtime(&now);
    strftime(s,80,"%c",date);
    
	/*observed data*/
    if((*data)->despl > 0 && (*data)->window > 0) {
        totalnloci = (int)ceil(((double)(*data)->nsites[1] - (double)(*data)->window) / ((double)(*data)->despl)) + (int)1;
	}	
    else {
        if((*data)->linked > 1) totalnloci = (*data)->linked;
        else totalnloci = (*data)->n_loci;
	}

    fprintf(out,"\n\"OUTPUT FILE: date %s\" \n\n\"Input data from the file: %s\"\n",s,in);
    print_var(23,out);
    fprintf(out," %d",(*data)->print_neuttest);
	print_var(16,out);
    fprintf(out," %d",(*data)->pr_matrix);
	if((*data)->pr_matrix == 1 && (*data)->patcg[0] == 5) {
		print_var(91,out);
		for(x=1;x<5;x++) fprintf(out," %G",(*data)->patcg[x]);
	}
    if((*data)->mhits && ((*data)->theta_1[1] > 0. || (*data)->sfix_allthetas || (*data)->range_thetant || (*data)->ifgamma)) {
        print_var(0,out);
        fprintf(out," %d",(*data)->mhits);
        if((*data)->T_out != 0) {
            print_var(36,out);
            if((*data)->T_out>=REFNUMBER && (*data)->T_out < REFNUMBER + PLUSPRIOR) fprintf(out," &%d",(int)((double)(*data)->T_out-REFNUMBER));
            else fprintf(out," %f",(*data)->T_out);

            print_var(17,out);
            j = (*data)->ratio_sv[1];
            for(y=2;y<(*data)->n_loci+1;y++) if((*data)->ratio_sv[y] != j) break;
            if(y < (*data)->n_loci + 1) {
                for(z=1;z<(*data)->n_loci+1;z++) {
                    if((*data)->ratio_sv[z]>=REFNUMBER && (*data)->ratio_sv[z] < REFNUMBER + PLUSPRIOR) 
                        fprintf(out," &%d",(int)((double)(*data)->ratio_sv[z]-REFNUMBER));
                    else fprintf(out," %.3G",(*data)->ratio_sv[z]);
                }
            }
            else {
                if(j>=REFNUMBER && j < REFNUMBER + PLUSPRIOR) fprintf(out," &%d",(int)((double)j-REFNUMBER));
                else fprintf(out," %.3G",j);
            }
        }
    }
    print_var(1,out);
    fprintf(out," %ld",(*data)->n_iter);                
    print_var(2,out);        
    for(y=1;y<(*data)->seed1[0]+1;y++) fprintf(out," %ld",(*data)->seed1[y]);
    
	fputs("\n",out);
    
	print_var(3,out);
    fprintf(out," %d",(*data)->n_loci);
    print_var(65,out);
	fprintf(out," %d",(*data)->no_rec_males);
	if((*data)->invariable_mut_sites[1] > (double)0) {
		print_var(74,out);        
		for(y=1;y<(*data)->invariable_mut_sites[0]+1;y++) 
			if((*data)->invariable_mut_sites[y]>=REFNUMBER && (*data)->invariable_mut_sites[y] < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->invariable_mut_sites[y]-REFNUMBER));
			else fprintf(out," %ld",(*data)->invariable_mut_sites[y]);
	}	
    print_var(55,out);
	for(y=1;y<(*data)->factor_chrn[0]+1;y++) fprintf(out," %f",(*data)->factor_chrn[y]);
    print_var(95,out);
	for(y=1;y<(*data)->event_sexratio[0]+1;y++) fprintf(out," %f",(*data)->event_sexratio[y]);
	print_var(97,out);
	for(y=1;y<(*data)->fixmax_abstime_event[0]+1;y++) fprintf(out," %f",(*data)->fixmax_abstime_event[y]);
	print_var(96,out);
	if((*data)->sex_ratio>=REFNUMBER && (*data)->sex_ratio < REFNUMBER + PLUSPRIOR) 
		fprintf(out," &%d",(int)((double)(*data)->sex_ratio-REFNUMBER));
	else fprintf(out," %f",(*data)->sex_ratio);
	print_var(6,out);        
    for(y=1;y<(int)(*data)->nsites[0]+1;y++) fprintf(out," %ld",(*data)->nsites[y]);
    print_var(4,out);        
    for(y=1;y<(*data)->nsam[0]+1;y++) fprintf(out," %d",(*data)->nsam[y]);
	print_var(11,out);
	fprintf(out," %ld",(*data)->npop);
	
    fputs("\n",out);
	
	if((*data)->method_samp) {
		print_var(58,out);
		fprintf(out," %d",(*data)->method_samp);
		if((*data)->method_samp > 1) {
			print_var(27,out);
			fprintf(out," %d",(*data)->mc_jump);
		}
	}
	
    if((*data)->linked) {

		fputs("\n",out);
		
        print_var(14,out);
        fprintf(out," %d",(*data)->linked);
        if((*data)->linked > 1) {
            print_var(15,out);
            for(x=0;x<(*data)->linked;x++) {
                if(x) fputs(",",out);
                for(z=1;z<3;z++) fprintf(out," %ld",(*data)->loci_linked[x][z]); 
            }
			/**/
			if((*data)->linked_segsites != -1) {
				print_var(67,out);
				fprintf(out," %d",(*data)->linked_segsites);
				print_var(79,out);
				fprintf(out," %d",(*data)->linked_segsites_nregion);
			}
			if((*data)->linked_rm != -1) {
				print_var(68,out);
				fprintf(out," %d",(*data)->linked_rm);
				print_var(80,out);
				fprintf(out," %d",(*data)->linked_rm_nregion);
			}
			if((*data)->linked_nhapl != 0) {
				print_var(71,out);
				fprintf(out," %d",(*data)->linked_nhapl);
				print_var(81,out);
				fprintf(out," %d",(*data)->linked_nhapl_nregion);
			}
			/**/
        }
        else {
            print_var(37,out);
            fprintf(out," %ld",(*data)->despl);
            print_var(38,out);
            fprintf(out," %ld",(*data)->window);
        }
    }    

    fputs("\n",out);

	if((*data)->range_thetant) {
		print_var(33,out);
		fprintf(out," %d",(*data)->range_thetant);
		if((*data)->ifgamma == 0) {
			print_var(34,out);
			for(y=1;y<(*data)->thetant_min[0]+1;y++) fprintf(out," %f",(*data)->thetant_min[y]);
			print_var(35,out);
			for(y=1;y<(*data)->thetant_max[0]+1;y++) fprintf(out," %f",(*data)->thetant_max[y]);
		}
	}
	if((*data)->ifgamma) {
		print_var(52,out);
		fprintf(out," %d",(*data)->ifgamma);
		print_var(54,out);        
		for(y=1;y<(*data)->alpha_gamma[0]+1;y++) fprintf(out," %G",(*data)->alpha_gamma[y]);
		print_var(53,out);        
		for(y=1;y<(*data)->p_gamma[0]+1;y++) fprintf(out," %g",(*data)->p_gamma[y]);
		print_var(72,out);        
		for(y=1;y<(*data)->correct_gamma[0]+1;y++) fprintf(out," %g",(*data)->correct_gamma[y]);
	}
	if(!((*data)->mutations[0] == 0 || ((*data)->mutations[0] == 1 && (*data)->mutations[1] == -1))) {
        for(y=1;y<(int)(*data)->mutations[0]+1;y++) if((*data)->mutations[y] > -1) break;
		if(y<(int)(*data)->mutations[0]+1) {
			print_var(10,out);
			for(y=1;y<(int)(*data)->mutations[0]+1;y++) fprintf(out," %ld",(*data)->mutations[y]);
		}
    }	
	if((*data)->sfix_allthetas && (*data)->method_samp) {
		print_var(26,out);
		fprintf(out," %d",(*data)->sfix_allthetas);
		print_var(76,out);
		fprintf(out," %d",(*data)->sfix_pop);
	}
	if(!((*data)->theta_1[0] == (double)0 || ((*data)->theta_1[0] == (double)1 && (*data)->theta_1[1] == (double)0))) {
        print_var(9,out);
        for(y=1;y<(*data)->theta_1[0]+1;y++) 
			if((*data)->theta_1[y]>=REFNUMBER && (*data)->theta_1[y] < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->theta_1[y]-REFNUMBER));
			else fprintf(out," %G",(*data)->theta_1[y]);
    }
	if((*data)->heter_theta_alphag[1] > (double)0) {
		print_var(69,out);
		for(y=1;y<(*data)->n_loci+1;y++) 
			if((*data)->heter_theta_alphag[y]>=REFNUMBER && (*data)->heter_theta_alphag[y] < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->heter_theta_alphag[y]-REFNUMBER));
			else fprintf(out," %f",(*data)->heter_theta_alphag[y]);
	}
    fputs("\n",out);

	if((*data)->range_rnt == 1 || (*data)->range_rnt == 2) {
		print_var(59,out);
		fprintf(out," %d",(*data)->range_rnt);
		print_var(60,out);
		for(y=1;y<(*data)->recnt_min[0]+1;y++) fprintf(out," %f",(*data)->recnt_min[y]);
		print_var(61,out);
		for(y=1;y<(*data)->recnt_max[0]+1;y++) fprintf(out," %f",(*data)->recnt_max[y]);
	}
	if((*data)->ifgammar == 1) {
		print_var(62,out);
		fprintf(out," %d",(*data)->ifgammar);
		print_var(63,out);        
		for(y=1;y<(*data)->alpha_gammar[0]+1;y++) fprintf(out," %G",(*data)->alpha_gammar[y]);
		print_var(64,out);        
		for(y=1;y<(*data)->p_gammar[0]+1;y++) fprintf(out," %g",(*data)->p_gammar[y]);
		print_var(73,out);        
		for(y=1;y<(*data)->correct_gammar[0]+1;y++) fprintf(out," %g",(*data)->correct_gammar[y]);
	}
	if((*data)->ifgammar == 0 && (*data)->range_rnt == 0) {
		print_var(5,out);        
		for(y=1;y<(*data)->r[0]+1;y++) 
			if((*data)->r[y]>=REFNUMBER && (*data)->r[y] < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->r[y]-REFNUMBER));
			else fprintf(out," %G",(*data)->r[y]); 
    }
	if((*data)->rmfix && (*data)->linked == 0 && (*data)->method_samp) {
		print_var(57,out);
		for(y=1;y<(int)(*data)->Rm[0]+1;y++) fprintf(out," %d",(*data)->Rm[y]);
		for(y=1;y<(int)(*data)->nhapl[0]+1;y++) {
			if((*data)->nhapl[y] > 0) {
				print_var(66,out);
				break;
			}
		}
		for(y=1;y<(int)(*data)->nhapl[0]+1;y++) {
			if((*data)->nhapl[y] > 0)
				fprintf(out," %d",(*data)->nhapl[y]);
		}
	}	
	if((*data)->rmfix && (*data)->method_samp) {
		print_var(56,out);
		fprintf(out," %d",(*data)->rmfix);
		print_var(77,out);
		fprintf(out," %d",(*data)->rmfix_pop);
	}
	if((*data)->heter_rm_alphag[1] > (double)0) {
		print_var(70,out);
		for(y=1;y<(*data)->n_loci+1;y++) 
			if((*data)->heter_rm_alphag[y]>=REFNUMBER && (*data)->heter_rm_alphag[y] < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->heter_rm_alphag[y]-REFNUMBER));
			else fprintf(out," %f",(*data)->heter_rm_alphag[y]);
	}
	
    z = 0;
    for(y=1;y<=(*data)->f[0];y++) {/*only print f and track_len if defined*/
        if((*data)->f[y]>0.) {
            if((*data)->ifselection == 0) z=1;
            break;
        }
    }
    if(z==1) {
        fputs("\n",out);
        if(!((*data)->f[0] == 0 || (*data)->f[1] == 0.)) {
            print_var(7,out);
            for(y=1;y<=(*data)->n_loci;y++) { 
                if((*data)->f[y]>=REFNUMBER && (*data)->f[y] < REFNUMBER + PLUSPRIOR) 
					fprintf(out," &%d",(int)((double)(*data)->f[y]-REFNUMBER));
				else fprintf(out," %G",(*data)->f[y]); 
            }
        }
        if(!((*data)->track_len[0] == 0 || (*data)->track_len[1] == 0.)) {
            print_var(8,out);        
            for(y=1;y<=(*data)->n_loci;y++) { 
                if((*data)->track_len[y]>=REFNUMBER && (*data)->track_len[y] < REFNUMBER + PLUSPRIOR) 
					fprintf(out," &%d",(int)((double)(*data)->track_len[y]-REFNUMBER));
				else fprintf(out," %G",(*data)->track_len[y]); 
            }
        }
    }

    z = 0;
	if((*data)->ifselection[0]>0) {
		for(y=1;y<=(*data)->ifselection[0];y++) {
            if((*data)->ifselection[y]==1) {
				z=1;
				break;
			}
        }
    }
    if(z==0) {
        if((*data)->npop > 1 && (*data)->split_pop == 0) {
            fputs("\n",out);
			if((*data)->mig_rate_matrix[0][0]) {
				print_var(78,out);
				for(x=0;x<(*data)->npop;x++) {
					if(x!=0) fprintf(out,",");
					for(y=0;y<(*data)->npop;y++) {
						if((*data)->mig_rate_matrix[x][y+1]>=REFNUMBER && (*data)->mig_rate_matrix[x][y+1] < REFNUMBER + PLUSPRIOR) 
							fprintf(out," &%d",(int)((double)(*data)->mig_rate_matrix[x][y+1]-REFNUMBER));
						else fprintf(out," %G",(*data)->mig_rate_matrix[x][y+1]);
					}
				}
			}
			else {
				if((*data)->mig_rate) {
					print_var(13,out);
					if((*data)->mig_rate >=REFNUMBER && (*data)->mig_rate < REFNUMBER + PLUSPRIOR) 
						fprintf(out," &%d",(int)((double)(*data)->mig_rate-REFNUMBER));
					else fprintf(out," %G",(*data)->mig_rate);
				}
			}
            print_var(24,out);
            for(y=1;y<=(*data)->npop_sampled[0];y++)
				fprintf(out," %d",(*data)->npop_sampled[y]);
            
            if((*data)->ran_factorpop == 0 && (*data)->same_factorpop == 0) {
                zz= 1;
                if((*data)->split_pop == 0) {
                    for(y=2;y<(*data)->npop_sampled[0]+1;y++) 
                        if((*data)->npop_sampled[y] > (*data)->npop_sampled[zz]) zz = y;
                    if(!((*data)->factor_pop[0] == 0 || ((*data)->factor_pop[0] == 1 && (*data)->factor_pop[1] == 0.0))){
                        print_var(19,out);
                        for(y=1;y<(*data)->factor_pop[0]+1;y++) {
							if((*data)->factor_pop[y]>=REFNUMBER && (*data)->factor_pop[y] < REFNUMBER + PLUSPRIOR) 
								fprintf(out," &%d",(int)((double)(*data)->factor_pop[y]-REFNUMBER));
							else fprintf(out," %.3G",(*data)->factor_pop[y]); 
						}
                    }
                }
                else {
					print_var(19,out);
					for(y=1;y<(*data)->npop_sampled[zz]+1;y++) {
						if((*data)->factor_pop[y]>=REFNUMBER && (*data)->factor_pop[y] < REFNUMBER + PLUSPRIOR) 
							fprintf(out," &%d",(int)((double)(*data)->factor_pop[y]-REFNUMBER));
						else fprintf(out," %.3G",(*data)->factor_pop[y]);
					}
				}
            }
            print_var(20,out);
            fprintf(out," %d",(*data)->ran_factorpop);
            print_var(21,out);
            fprintf(out," %d",(*data)->same_factorpop);
            if(!((*data)->ssize_pop[0][0] == 0 || ((*data)->ssize_pop[0][0] == 1 && (*data)->ssize_pop[0][1] == 0))){
                print_var(12,out);
                for(y=0;y<(*data)->n_loci;y++) { 
                    /*però sembla ho tinc restringit a nomes la mateixa mostra per tots locus...*/
                    if(y) fputs(",",out);
                    for(zz=1;zz<(*data)->npop_sampled[y+1]+1;zz++)
                        fprintf(out," %d",(*data)->ssize_pop[y][zz]); 
                }
            }            fputs("\n",out);
        }
		
		flag = 0;

		if((*data)->npoprefugia > (*data)->npop) 
			npopt = (*data)->npoprefugia;
		else {
			if((*data)->ifselection[0]) npopt = (int)((*data)->npop > 2 ? (*data)->npop:2);
			else npopt = (int)(*data)->npop;
		}
		
        if(!(*data)->npop_events) {
            for(x=0;x<npopt;x++) {
                if((*data)->nintn[x] != 0) {
                    print_var(50,out);
                    fprintf(out," %d",(*data)->iflogistic);
                    flag = 1;
                    break;
                }
            }
            if(flag == 1) {
                print_var(39,out);
                for(x=0;x<(*data)->npop;x++) {
                    fprintf(out," %d",(*data)->nintn[x+1]);
                }
                if((*data)->iflogistic) {
                    print_var(51,out);
                    if((*data)->ts[x+1]>=REFNUMBER && (*data)->ts[x+1] < REFNUMBER + PLUSPRIOR) 
                        fprintf(out," &%d",(int)((double)(*data)->ts[x+1]-REFNUMBER));
                    else fprintf(out," %G",(*data)->ts[x+1]);
                }
                print_var(40,out);
                for(x=0;x<(*data)->npop;x++) {
                    if(x>0) fprintf(out,", ");
                    for(y=1;y<=(*data)->nintn[x+1];y++)
                        if((*data)->nrec[x][y]>=REFNUMBER && (*data)->nrec[x][y] < REFNUMBER + PLUSPRIOR) 
                            fprintf(out," &%d",(int)((double)(*data)->nrec[x][y]-REFNUMBER));
                        else  fprintf(out," %G",(*data)->nrec[x][y]);
                }
                print_var(42,out);        
                for(x=0;x<(*data)->npop;x++) {
                    if(x>0) fprintf(out,", ");
                    for(y=1;y<=(*data)->nintn[x+1];y++) 
                        if((*data)->tpast[x][y]>=REFNUMBER && (*data)->tpast[x][y] < REFNUMBER + PLUSPRIOR) 
                            fprintf(out," &%d",(int)((double)(*data)->tpast[x][y]-REFNUMBER));
                        else  fprintf(out," %G",(*data)->tpast[x][y]);
                }
                
                fputs("\n",out);
            }
        }

		if((*data)->split_pop == 1) {
            print_var(43,out);
            fprintf(out," %d",(*data)->split_pop);
            print_var(49,out);
            fprintf(out," %d",(*data)->npoprefugia);
            print_var(13,out);
            if((*data)->mig_rate>=REFNUMBER && (*data)->mig_rate < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->mig_rate-REFNUMBER));
			else fprintf(out," %G",(*data)->mig_rate);
            print_var(44,out);
            if((*data)->time_split>=REFNUMBER && (*data)->time_split < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->time_split-REFNUMBER));
			else fprintf(out," %f",(*data)->time_split);
			print_var(45,out);
			if((*data)->time_scoal>=REFNUMBER && (*data)->time_scoal < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->time_scoal-REFNUMBER));
			else fprintf(out," %f",(*data)->time_scoal);
            print_var(47,out);        
            for(y=1;y<(*data)->npoprefugia+1;y++) 
				if((*data)->freq[y]>=REFNUMBER && (*data)->freq[y] < REFNUMBER + PLUSPRIOR) 
					fprintf(out," &%d",(int)((double)(*data)->freq[y]-REFNUMBER));
				else fprintf(out," %G",(*data)->freq[y]);
			print_var(46,out);
			if((*data)->factor_anc>=REFNUMBER && (*data)->factor_anc < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->factor_anc-REFNUMBER));
			else fprintf(out," %f",(*data)->factor_anc);
            if((*data)->ran_factorpop == 0 && (*data)->same_factorpop == 0) {
				print_var(19,out);
                for(y=1;y<(*data)->npoprefugia+1;y++) 
					if((*data)->factor_pop[y]>=REFNUMBER && (*data)->factor_pop[y] < REFNUMBER + PLUSPRIOR) 
						fprintf(out," &%d",(int)((double)(*data)->factor_pop[y]-REFNUMBER));
					else fprintf(out," %.3G",(*data)->factor_pop[y]);
            }
            print_var(20,out);
            fprintf(out," %d",(*data)->ran_factorpop);
            print_var(21,out);
            fprintf(out," %d",(*data)->same_factorpop);
        }
    }
    else{ /*in case selection*/
        fputs("\n",out);
		print_var(29,out);
        if((*data)->pop_size>=REFNUMBER && (*data)->pop_size < REFNUMBER + PLUSPRIOR) 
			fprintf(out," &%d",(int)((double)(*data)->pop_size-REFNUMBER));
		else fprintf(out," %.2E",(*data)->pop_size);
        print_var(28,out);
        for(y=1;y<(*data)->ifselection[0]+1;y++) fprintf(out," %d",(*data)->ifselection[y]);
        print_var(30,out);
        for(y=1;y<(*data)->pop_sel[0]+1;y++) 
			if((*data)->pop_sel[y]>=REFNUMBER && (*data)->pop_sel[y] < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->pop_sel[y]-REFNUMBER));
			else fprintf(out," %G",(*data)->pop_sel[y]);
        print_var(31,out);
        for(y=1;y<(*data)->sel_nt[0]+1;y++) 
			if((*data)->sel_nt[y]>=REFNUMBER && (*data)->sel_nt[y] < REFNUMBER + PLUSPRIOR) 
				fprintf(out," &%d",(int)((double)(*data)->sel_nt[y]-REFNUMBER));
			else fprintf(out," %ld",(long int)(*data)->sel_nt[y]);
        print_var(32,out);
        for(y=1;y<(*data)->sinit[0]+1;y++)
            if((*data)->sinit[y]>=REFNUMBER && (*data)->sinit[y] < REFNUMBER + PLUSPRIOR)
                fprintf(out," &%d",(int)((double)(*data)->sinit[y]-REFNUMBER));
            else fprintf(out," %G",(*data)->sinit[y]);
        print_var(98,out);
        for(y=1;y<(*data)->sendt[0]+1;y++)
            if((*data)->sendt[y]>=REFNUMBER && (*data)->sendt[y] < REFNUMBER + PLUSPRIOR)
                fprintf(out," &%d",(int)((double)(*data)->sendt[y]-REFNUMBER));
            else fprintf(out," %G",(*data)->sendt[y]);
        print_var(99,out);
        for(y=1;y<(*data)->sfreqend[0]+1;y++)
            if((*data)->sfreqend[y]>=REFNUMBER && (*data)->sfreqend[y] < REFNUMBER + PLUSPRIOR)
                fprintf(out," &%d",(int)((double)(*data)->sfreqend[y]-REFNUMBER));
            else fprintf(out," %G",(*data)->sfreqend[y]);
        print_var(100,out);
        for(y=1;y<(*data)->sfreqinit[0]+1;y++)
            if((*data)->sfreqinit[y]>=REFNUMBER && (*data)->sfreqinit[y] < REFNUMBER + PLUSPRIOR)
                fprintf(out," &%d",(int)((double)(*data)->sfreqinit[y]-REFNUMBER));
            else fprintf(out," %G",(*data)->sfreqinit[y]);
        fputs("\n\n\"Recombination parameter value for the studied region and for all the region up the selective position.\"\n",out);
    }
	
	/*print pop_events*/
	fputs("\n",out);
	if((*data)->npop_events) {
		print_var(88,out);
		fprintf(out," %d",(int)(*data)->npop_events);
		print_var(89,out);
		for(x=0;x<(*data)->npop_events;x++) {
			if(x!=0) fprintf(out,",");
			for(y=0;y<(*data)->pop_event[x][0];y++) {
				if((*data)->pop_event[x][y+1]>=REFNUMBER && (*data)->pop_event[x][y+1] < REFNUMBER + PLUSPRIOR) {
					fprintf(out," &%d",(int)((double)(*data)->pop_event[x][y+1]-REFNUMBER));
				}
				else {
					if(y<3) fprintf(out," %d",(int)(*data)->pop_event[x][y+1]);
					else fprintf(out," %G",(*data)->pop_event[x][y+1]);
				}
			}
		}
	}
	
	/*ancestral_pol*/
	fputs("\n",out);
	if((*data)->npop > 1 && (*data)->ancestral_pol[1] == 1) {
		print_var(90,out);
		for(y=1;y<(*data)->ancestral_pol[0]+1;y++) 
			fprintf(out," %d",(*data)->ancestral_pol[y]);
	}
	
	if((*data)->pop_outgroup) {
		print_var(91,out);
		fprintf(out," %d",(*data)->pop_outgroup);
	}
	
	/*printing priors*/
	fputs("\n",out);
	for(x=0;x<(*data)->npriors;x++) {
		fprintf(out,"\n\"Prior referred to parameter %s\"",priors[x].name_pr);
		print_var(82,out);
		fprintf(out," %d",(int)priors[x].idprior);
		if(priors[x].prior_file[0] != '\0') {
			print_var(87,out);
			fprintf(out,"%s",priors[x].prior_file);
		}
		else {
			print_var(83,out);
			fprintf(out," %d",priors[x].kind_var);
			if(priors[x].kind_var == 0) fprintf(out," \"long int\"");
			if(priors[x].kind_var == 1) fprintf(out," \"double\"");
			print_var(84,out);
			fprintf(out," %d",priors[x].kind_dist);
			if(priors[x].kind_dist == 0) fprintf(out," \"uniform\"");
			if(priors[x].kind_dist == 1) fprintf(out," \"log10-uniform\"");
			if(priors[x].kind_dist == 2) fprintf(out," \"gamma\"");
			if(priors[x].kind_dist == 3) fprintf(out," \"beta\"");
			print_var(85,out);
			for(y=1;y<=priors[x].dist_par[0];y++) fprintf(out," %G",priors[x].dist_par[y]);
			if(priors[x].kind_dist == 0) fprintf(out," \"min max\"");
			if(priors[x].kind_dist == 1) fprintf(out," \"min max\"");
			if(priors[x].kind_dist == 2) fprintf(out," \"alpha p factor_corr\"");
			if(priors[x].kind_dist == 3) fprintf(out," \"alpha beta min max\"");
			print_var(86,out);
			fprintf(out," %ld",priors[x].seed_prior);		
			fputs("\n",out);
		}
		
		dfiles[0] = '\0';
		strcat(dfiles,file_out);
		el = strrchr(dfiles,'.');
		*el = '\0';
		strcat(dfiles,"_prior");
		/*there is no function??*/
		zz = (int)(priors[x].idprior);
		z = 10000;
		for(y=0;y<5;y++) {
			npr = div(zz,z);
			nchar[y] = 48 + npr.quot;
			zz = npr.rem;
			z /= 10;
		}
		nchar[5] = '\0';
		strcat(dfiles,nchar);
		strcat(dfiles,".out");
		if(!(outputind = fopen(dfiles,"w"))) {
			printf(" Error creating prior values files.\n");
			exit(1);
		}
		fprintf(outputind,"\n\"Prior referred to parameter %s\"",priors[x].name_pr);
		print_var(82,outputind);
		fprintf(outputind," %d",(int)priors[x].idprior);
		if(priors[x].prior_file[0] != '\0') {
			print_var(87,outputind);
			fprintf(outputind,"%s",priors[x].prior_file);
		}
		else {
			print_var(83,outputind);
			fprintf(outputind," %d",priors[x].kind_var);
			if(priors[x].kind_var == 0) fprintf(outputind," \"long int\"");
			if(priors[x].kind_var == 1) fprintf(outputind," \"double\"");
			print_var(84,outputind);
			fprintf(outputind," %d",priors[x].kind_dist);
			if(priors[x].kind_dist == 0) fprintf(outputind," \"uniform\"");
			if(priors[x].kind_dist == 1) fprintf(outputind," \"log10-uniform\"");
			if(priors[x].kind_dist == 2) fprintf(outputind," \"gamma\"");
			if(priors[x].kind_dist == 3) fprintf(outputind," \"beta\"");
			print_var(85,outputind);
			for(y=1;y<=priors[x].dist_par[0];y++) fprintf(outputind," %G",priors[x].dist_par[y]);
			if(priors[x].kind_dist == 0) fprintf(outputind," \"min max\"");
			if(priors[x].kind_dist == 1) fprintf(outputind," \"min max\"");
			if(priors[x].kind_dist == 2) fprintf(outputind," \"alpha p factor_corr\"");
			if(priors[x].kind_dist == 3) fprintf(outputind," \"alpha beta min max\"");
			print_var(86,outputind);
			fprintf(outputind," %ld",priors[x].seed_prior);	
		}
		fputs("\n\nPrior distribution values:\n",outputind);
		
		for(lx=0;lx<(*data)->n_iter;lx++) {
			fprintf(outputind,"%G\n",priors[x].priordist[lx]);
		}
		fputs("\n",outputind);
		fclose(outputind);
	}

	if((*data)->includeEHHrel) {
		print_var(92,out);
		fprintf(out," %d",(*data)->includeEHHrel);
		if((*data)->ehh_fixnt[1] == 1) {
			print_var(94,out);
			for(y=1;y<(*data)->ehh_fixnt[0]+1;y++) 
				fprintf(out," %ld",(*data)->ehh_fixnt[y]);
		}
	}
	
	defnumber = 101;
	for(x=0;x<NOBS_STATISTICS;x++) {
		if((*data)->obs_statistics_used[x]) {
			print_var(defnumber+x,out);
			for(y=0;y<totalnloci/*(*data)->n_loci*/;y++) {
				if(y>0) fputs(", ",out);
				for(z=1;z<=(*data)->max_npop_sampled;z++) {
					fprintf(out," %G",(*data)->obs_statistics[x][y][z]);
				}
			}
		}
	}
	fputs("\n",out);
	for(x=0;x<NOBS_STATISTICS;x++) {
		if((*data)->obs_statistics_used[x]) {
			print_var(defnumber+2*NOBS_STATISTICS+x,out);
			fprintf(out," %d",(*data)->obs_statistics_used[x]);
		}
	}
	fputs("\n",out);
	if((*data)->likelihood_line) {
		for(x=0;x<NOBS_STATISTICS;x++) {
			if((*data)->obs_statistics_used[x]) {
				print_var(defnumber+NOBS_STATISTICS+x,out);
				fprintf(out," %G",(*data)->likelihood_error[x]);
			}
		}
	}
		
	if((*data)->sfix_allthetas == 1 && (*data)->method_samp == 1) 
        fputs("\n\"REJECTION ALGORITHM. Fix S and screening thetas given a prior.\"",out);
    if((*data)->rmfix == 1 && (*data)->method_samp == 1) 
        fputs("\n\"REJECTION ALGORITHM. Fix Rm (and optionally nhapl) and screening R given a prior.\"",out);
    fputs("\n\n",out);
    
    if(fflush(out) != 0) {
        puts("\nError. Buffer print error\n");
        exit(1);
    }

}

void print_var(int x,FILE *out)
{
    int y,z;
    char til20[30];

    til20[0] = var_file[x][0] - 32;
    for(y=1;y<29;y++) {
        til20[y] = var_file[x][y];
        if(til20[y] == '\0') break;
    }
    til20[y] = ':';
    for(z=y+1;z<29;z++) til20[z] = 32;
    til20[29] = '\0';
    fprintf(out,"\n%s   ",til20);
    return;
}

/*************** FREE *********************/
void free_inputdata(struct var **data,struct var_priors *priors,int my_rank)
{
    int x,y;
    int loc;

	for(x=0;x<(*data)->npriors;x++) {
		/*printf("\nx=%d,npriors=%d",x,(*data)->npriors);*/
		free(priors[x].priordist);
	}
	
	free((*data)->nsam);
    free((*data)->r);
    free((*data)->nsites);
    if((*data)->linked > 20) loc = (*data)->linked;
    else loc = 20;   
    for(x=0;x<loc;x++) free((*data)->loci_linked[x]);
    free((*data)->loci_linked);
    free((*data)->f);
    free((*data)->track_len);
    free((*data)->theta_1);    
    free((*data)->mutations);
    if((*data)->n_loci > 20) loc = (*data)->n_loci;
    else loc = 20;   
    for(y=0;y<loc;y++) free((*data)->ssize_pop[y]);
    free((*data)->ssize_pop);
    free((*data)->factor_pop);
    free((*data)->ratio_sv);
    free((*data)->npop_sampled);
    free((*data)->ifselection);
    free((*data)->pop_sel);
    free((*data)->sel_nt);
    free((*data)->sinit);
    free((*data)->sendt);
    free((*data)->sfreqend);
    free((*data)->sfreqinit);
    free((*data)->ts);

	if((*data)->npop > 20) loc = (int)(*data)->npop;
    else loc = 20;   
    for(x=0;x<loc;x++) {
		free((*data)->nrec[x]);
		free((*data)->tpast[x]);
		free((*data)->mig_rate_matrix[x]);
		free((*data)->pop_event[x]);
	}
    free((*data)->nrec);
    free((*data)->tpast);
	free((*data)->mig_rate_matrix);
	free((*data)->pop_event);
    
	free((*data)->freq);
	
    free((*data)->p_gamma);
    free((*data)->alpha_gamma);
    free((*data)->correct_gamma);
    free((*data)->factor_chrn);

    free((*data)->p_gammar);
    free((*data)->alpha_gammar);
    free((*data)->correct_gammar);
    free((*data)->Rm);
    free((*data)->nhapl);
	
	free((*data)->heter_theta_alphag);
    free((*data)->invariable_mut_sites);
    free((*data)->heter_rm_alphag);

    free((*data)->thetant_min);
    free((*data)->thetant_max);

    free((*data)->recnt_min);
    free((*data)->recnt_max);

    free((*data)->ancestral_pol);
    free((*data)->ehh_fixnt);
    free((*data)->event_sexratio);
	free((*data)->fixmax_abstime_event);

    if(my_rank == 0) {
		if((*data)->n_loci > 20) loc = (*data)->n_loci;
		else loc = 20;   
		for(x=0;x<NOBS_STATISTICS;x++) {
			for(y=0;y<loc;y++)
				free((*data)->obs_statistics[x][y]);
			free((*data)->obs_statistics[x]);
		}
	}
    return;
}

#if inMPI == 1
void init_data_priors(struct var **data,struct var_priors **priors) 
{
    int x,y;

	/*initialize memory*/
	if(!(*data = (struct var *)calloc(1,sizeof(struct var)))) perror("calloc error.0");;
    
    if(!((*data)->nsam = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.1");
    if(!((*data)->r = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.2");
    if(!((*data)->nsites = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.4");    
    if(!((*data)->loci_linked = (long int **)calloc((unsigned)20,sizeof(long int *)))) perror("calloc error.4b");
    for(x=0;x<20;x++) 
        if(!((*data)->loci_linked[x] = (long int *)calloc((unsigned)20,sizeof(long int))))
            perror("calloc error.4c");
    if(!((*data)->f = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.5");
    if(!((*data)->track_len = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.7");
    if(!((*data)->theta_1 = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.9");    
	
    if(!((*data)->mutations = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.10");
    if(!((*data)->ssize_pop = (int **)calloc((unsigned)20,sizeof(int *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->ssize_pop[x] = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.11b");
    if(!((*data)->factor_pop = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->ratio_sv = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->npop_sampled = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.34");
    if(!((*data)->ifselection = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.34");
    if(!((*data)->pop_sel = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sel_nt = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sinit = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sendt = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sfreqend = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sfreqinit = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
	
	if(!((*data)->ts = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->nintn = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.88");	
	if(!((*data)->nrec = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->nrec[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.11b");
    if(!((*data)->tpast = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->tpast[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.11b");
	
    if(!((*data)->freq = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.16");
    
	if(!((*data)->p_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
    if(!((*data)->alpha_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
    if(!((*data)->correct_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
    if(!((*data)->factor_chrn = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 58");
    
	if(!((*data)->p_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
	if(!((*data)->alpha_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
    if(!((*data)->correct_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
	if(!((*data)->Rm = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.56");
	if(!((*data)->nhapl = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.56");
	
	if(!((*data)->heter_theta_alphag = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.91");
	if(!((*data)->heter_rm_alphag = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->invariable_mut_sites = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.9");
	
	if(!((*data)->thetant_min = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->thetant_max = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	
	if(!((*data)->recnt_min = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->recnt_max = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	
	if(!((*data)->mig_rate_matrix = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->mig_rate_matrix[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.11b");
	
	if(!((*data)->pop_event = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->pop_event[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.11b");
		
	if(!(*priors = (struct var_priors *)calloc((unsigned)1,sizeof(struct var_priors)))) perror("calloc error.110");
	
 	if(!((*data)->ancestral_pol = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.561");
 	if(!((*data)->ehh_fixnt = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.561");
 	if(!((*data)->seed1 = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.561o");
 	if(!((*data)->event_sexratio = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.561o");
 	if(!((*data)->fixmax_abstime_event = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.561o");
	
    (*data)->pr_matrix = 1;
    (*data)->mhits = 0;
    (*data)->n_iter = 0;
    (*data)->n_loci = 1;
    (*data)->linked = 0;
    (*data)->despl = 0;
    (*data)->window = 0;
    (*data)->npop = 1;
    (*data)->mig_rate = 0.0;
    
    (*data)->seed1[0] = 1;  (*data)->seed1[1] = 12345678;
    (*data)->seed2 = 546373;
	
    (*data)->ran_factorpop = 0;
    (*data)->same_factorpop = 0;
    (*data)->neutral_tests = 0;
    (*data)->print_neuttest = 0;
	
    (*data)->sfix_allthetas = 0;
    (*data)->sfix_pop = -1;
    (*data)->mc_jump = 0;
    (*data)->pop_size = 1E+06;
    
    (*data)->range_thetant = 0;
    
    (*data)->T_out = 0.0;
	
    (*data)->split_pop = 0;
    (*data)->time_split = 0.;
    (*data)->time_scoal = (double)1E07;
    (*data)->factor_anc = 1.;

	(*data)->pop_event[0][0] = 1;
	(*data)->pop_event[0][1] = 0;
	
    (*data)->tlimit = 1000.;
    
	(*data)->iflogistic = 0;
	
    (*data)->ifgamma = 0;
	
    (*data)->rmfix = 0;
	(*data)->rmfix_pop = -1;
    (*data)->range_rnt = 0;
    (*data)->ifgammar = 0;
    (*data)->no_rec_males = 0;
	
	(*data)->linked_segsites = -1;
	(*data)->linked_rm = -1;
	(*data)->linked_nhapl = 0;
	(*data)->linked_segsites_nregion = -1;
	(*data)->linked_rm_nregion = -1;
	(*data)->linked_nhapl_nregion = -1;
	
	(*data)->includeEHHrel = 0;
	(*data)->pop_outgroup = -1;
	(*data)->sex_ratio = 1;

	/*observed values*/
	for(x=0;x<NOBS_STATISTICS;x++) {
		if(!((*data)->obs_statistics[x] = (double **)calloc((unsigned)20,sizeof(double *)))) perror("calloc error.113");
		for(y=0;y<20;y++) {
			if(!((*data)->obs_statistics[x][y] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.113b");
			(*data)->obs_statistics[x][y][0] = (double)0;/*redundant*/
		}
		(*data)->likelihood_error[x] = (double)1e-6;
		(*data)->obs_statistics_used[x] = 0;
	}
	
	(*data)->likelihood_line = 0;
    (*data)->npriors = 0;
	
	return;
}


void BuildBcast_derived_data(struct var **data,struct var_priors **priors,int my_rank,int *npmpi,int **nppr,long int **niterpr)
{
	MPI_Datatype mesg_mpi_data;
	MPI_Datatype mesg_mpi_prior;
	
	int block_lengths[NVARSIZESIMPLE];
	MPI_Aint displacements[NVARSIZESIMPLE];
	MPI_Datatype typelist[NVARSIZESIMPLE];
	MPI_Aint start_address;
	MPI_Aint address;
	/*
	MPI_Datatype column_mpi_t;
	MPI_Aint *displacements2;
	MPI_Datatype *typelist2;
	int *block_lengths2;
	*/
	div_t npl;
	long int npi;
	int restl,resti;

	int i,j/*,k*/;
	void init_data_priors(struct var **,struct var_priors **);
	void rankprint(struct var **,struct var_priors **,int);
	
	if(my_rank != 0) 
		init_data_priors(data,priors); 
		
	/*communicate simple variables from struct var*/
	for(i=0;i<NVARSIZESIMPLE;i++) block_lengths[i] = 1;
	block_lengths[41] = NOBS_STATISTICS;
	block_lengths[43] = NOBS_STATISTICS;
	block_lengths[46] = 5;
	
	i=0;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_LONG;
	typelist[i++] = MPI_LONG;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_LONG;
	typelist[i++] = MPI_LONG;
	typelist[i++] = MPI_LONG;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_DOUBLE;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_INT;
	typelist[i++] = MPI_DOUBLE;
	
	i=0;
	displacements[i++] = 0;
	MPI_Get_address(&((*data)->pr_matrix),&start_address);
	MPI_Get_address(&((*data)->mhits),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->T_out),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->n_iter),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->seed2),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->n_loci),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->linked),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->despl),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->window),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->npop),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->mig_rate),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->pop_size),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->ran_factorpop),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->same_factorpop),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->max_npop_sampled),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->neutral_tests),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->print_neuttest),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->sfix_allthetas),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->sfix_pop),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->mc_jump),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->range_thetant),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->ifgamma),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->iflogistic),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->split_pop),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->time_split),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->time_scoal),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->factor_anc),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->tlimit),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->npoprefugia),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->rmfix),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->rmfix_pop),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->method_samp),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->range_rnt),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->ifgammar),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->no_rec_males),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->linked_segsites),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->linked_rm),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->linked_nhapl),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->linked_segsites_nregion),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->linked_rm_nregion),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->linked_nhapl_nregion),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address((*data)->obs_statistics_used,&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->likelihood_line),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address((*data)->likelihood_error,&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->npriors),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->npop_events),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->patcg),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->includeEHHrel),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->pop_outgroup),&address);
	displacements[i++] = address - start_address;
	MPI_Get_address(&((*data)->sex_ratio),&address);
	displacements[i++] = address - start_address;
		
	MPI_Type_create_struct(NVARSIZESIMPLE,block_lengths,displacements,typelist,&mesg_mpi_data);
	MPI_Type_commit(&mesg_mpi_data);
	MPI_Bcast(*data,1,mesg_mpi_data,0,MPI_COMM_WORLD);
	
	/*ABORT IF NUMBER OF PROCESSES are higher than nloci or niter*/
	if(*npmpi > (*data)->n_loci && (*data)->n_loci > 1) {
		if(my_rank == 0) printf("\nSorry, it is not allowed more processes than nloci (for nloci > 1)\n");
		MPI_Abort(MPI_COMM_WORLD,1);
		exit(0); 
	}
	if(((long int)*npmpi > ((*data)->n_iter + (*data)->mc_jump)) && (*data)->n_loci == 1) {
		if(my_rank == 0) printf("\nSorry, it is not allowed more processes than niter (for nloci == 1)\n");
		MPI_Abort(MPI_COMM_WORLD,1);
		exit(0); 
	}
	/*DEFINE THE NUMBER OF LOCI AND ITERATIONS PER PROCESS*/
	if(!(*nppr = (int *)realloc(*nppr,(*npmpi+1)*sizeof(int)))) {
		perror("calloc error 097");
		exit(1);
	}
	if(!(*niterpr = (long int *)realloc(*niterpr,(*npmpi+1)*sizeof(long int)))) {
		perror("calloc error 098");
		exit(1);
	}

	/*we do cummulative, then we set the value 0 as initial (0)*/
	nppr[0][0] = 0;
	niterpr[0][0] = 0;
	
	if((*data)->n_loci > 1) {
		if(*npmpi <= (*data)->n_loci) {
			/*for more loci than processes, divide the nloci*/
			npl = div((*data)->n_loci,*npmpi);
			restl = npl.rem;
			for(i=1;i<*npmpi+1;i++) {
				nppr[0][i] = npl.quot;
				if(restl) {
					nppr[0][i] += 1;
					restl -= 1;
				}
				niterpr[0][i] = (*data)->n_iter + (*data)->mc_jump;
			}
		}
		else {
			for(i=1;i<(*data)->n_loci+1;i++) {
				nppr[0][i] = 1;
				niterpr[0][i] = (*data)->n_iter + (*data)->mc_jump;
			}
			for(i=(*data)->n_loci;i<*npmpi+1;i++) {
				nppr[0][i] = 0;
				niterpr[0][i] = 0;
			}
		}
	}
	else {
		/*for a single locus, divide the iterations...*/
		npi = (long int)floor(((double)(*data)->n_iter + (*data)->mc_jump)/(double)(*npmpi));
		resti = ((long int)(*data)->n_iter + (*data)->mc_jump) - (long int)npi*(long int)(*npmpi);
		for(i=1;i<*npmpi+1;i++) {
			nppr[0][i] = 0;
			niterpr[0][i] = npi;
			if(resti) {
				niterpr[0][i] += 1;
				resti -= 1;
			}
		}
	}
	/*Cummulative*/
	for(i=1;i<*npmpi+1;i++) {
		nppr[0][i] += nppr[0][i-1];
		niterpr[0][i] += niterpr[0][i-1];
	}
	/*
	 printf("\nmy_rank: %d, npmpi: %d ",my_rank,*npmpi);
	 for(i=1;i<*npmpi+1;i++) {
		printf("nppr[0][%d]: %d ",i,nppr[0][i]);
		printf("niterpr[0][%d]: %ld ",i,niterpr[0][i]);
	 }
	 fflush(stdout);
	 */
	/* comunicate pointers from struct var*/
	if(my_rank != 0) {
		/*reallocate pointers*/
		(*data)->nsam = (int *)realloc((*data)->nsam,(unsigned)((*data)->n_loci+1)*sizeof(int));
		(*data)->nsites = (long int *)realloc((*data)->nsites,(unsigned)((*data)->n_loci+1)*sizeof(long int));
	}
	MPI_Bcast((*data)->nsam,(*data)->n_loci+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->nsites,(*data)->n_loci+1,MPI_LONG,0,MPI_COMM_WORLD);
	
	/*take/send the first value of the vector. This is the size of the vector*/
	MPI_Bcast((*data)->seed1,1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->npop_sampled,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ratio_sv,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->r,1,MPI_DOUBLE,0,MPI_COMM_WORLD);	
	MPI_Bcast((*data)->f,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->track_len,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->theta_1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->mutations,1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->factor_pop,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ifselection,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->pop_sel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->sel_nt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((*data)->sinit,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((*data)->sendt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((*data)->sfreqend,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((*data)->sfreqinit,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ts,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->nintn,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->freq,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->p_gamma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->alpha_gamma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->correct_gamma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);	
	MPI_Bcast((*data)->factor_chrn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->p_gammar,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->alpha_gammar,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->correct_gammar,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->Rm,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->nhapl,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->heter_theta_alphag,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->heter_rm_alphag,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->invariable_mut_sites,1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->thetant_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->thetant_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->recnt_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->recnt_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ancestral_pol,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ehh_fixnt,1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->event_sexratio,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->fixmax_abstime_event,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if(my_rank != 0) {
		/*reallocate pointers*/
		(*data)->seed1 = (long int *)realloc((*data)->seed1,(unsigned)((*data)->seed1[0]+1)*sizeof(long int));
		(*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,(unsigned)((*data)->npop_sampled[0]+1)*sizeof(int));
		(*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,(unsigned)((*data)->ratio_sv[0]+1)*sizeof(double));
		(*data)->r = (double *)realloc((*data)->r,(unsigned)((*data)->r[0]+1)*sizeof(double));
		(*data)->f = (double *)realloc((*data)->f,(unsigned)((*data)->f[0]+1)*sizeof(double));
		(*data)->track_len = (double *)realloc((*data)->track_len,(unsigned)((*data)->track_len[0]+1)*sizeof(double));
		(*data)->theta_1 = (double *)realloc((*data)->theta_1,(unsigned)((*data)->theta_1[0]+1)*sizeof(double));
		(*data)->mutations = (long int *)realloc((*data)->mutations,(unsigned)((*data)->mutations[0]+1)*sizeof(long int));
		(*data)->factor_pop = (double *)realloc((*data)->factor_pop,(unsigned)((*data)->factor_pop[0]+1)*sizeof(double));
		(*data)->ifselection = (int *)realloc((*data)->ifselection,(unsigned)((*data)->ifselection[0]+1)*sizeof(int));
		(*data)->pop_sel = (double *)realloc((*data)->pop_sel,(unsigned)((*data)->pop_sel[0]+1)*sizeof(double));
		(*data)->sel_nt = (double *)realloc((*data)->sel_nt,(unsigned)((*data)->sel_nt[0]+1)*sizeof(double));
        (*data)->sinit = (double *)realloc((*data)->sinit,(unsigned)((*data)->sinit[0]+1)*sizeof(double));
        (*data)->sendt = (double *)realloc((*data)->sendt,(unsigned)((*data)->sendt[0]+1)*sizeof(double));
        (*data)->sfreqend = (double *)realloc((*data)->sfreqend,(unsigned)((*data)->sfreqend[0]+1)*sizeof(double));
        (*data)->sfreqinit = (double *)realloc((*data)->sfreqinit,(unsigned)((*data)->sfreqinit[0]+1)*sizeof(double));
		(*data)->ts = (double *)realloc((*data)->ts,(unsigned)((*data)->ts[0]+1)*sizeof(double));
		(*data)->nintn = (int *)realloc((*data)->nintn,(unsigned)((*data)->nintn[0]+1)*sizeof(int));
		(*data)->freq = (double *)realloc((*data)->freq,(unsigned)((*data)->freq[0]+1)*sizeof(double));
		(*data)->p_gamma = (double *)realloc((*data)->p_gamma,(unsigned)((*data)->p_gamma[0]+1)*sizeof(double));
		(*data)->alpha_gamma = (double *)realloc((*data)->alpha_gamma,(unsigned)((*data)->alpha_gamma[0]+1)*sizeof(double));
		(*data)->correct_gamma = (double *)realloc((*data)->correct_gamma,(unsigned)((*data)->correct_gamma[0]+1)*sizeof(double));
		(*data)->factor_chrn = (double *)realloc((*data)->factor_chrn,(unsigned)((*data)->factor_chrn[0]+1)*sizeof(double));
		(*data)->p_gammar = (double *)realloc((*data)->p_gammar,(unsigned)((*data)->p_gammar[0]+1)*sizeof(double));
		(*data)->alpha_gammar = (double *)realloc((*data)->alpha_gammar,(unsigned)((*data)->alpha_gammar[0]+1)*sizeof(double));
		(*data)->correct_gammar = (double *)realloc((*data)->correct_gammar,(unsigned)((*data)->correct_gammar[0]+1)*sizeof(double));
		(*data)->Rm = (int *)realloc((*data)->Rm,(unsigned)((*data)->Rm[0]+1)*sizeof(int));
		(*data)->nhapl = (int *)realloc((*data)->nhapl,(unsigned)((*data)->nhapl[0]+1)*sizeof(int));
		(*data)->heter_theta_alphag = (double *)realloc((*data)->heter_theta_alphag,(unsigned)((*data)->heter_theta_alphag[0]+1)*sizeof(double));
		(*data)->heter_rm_alphag = (double *)realloc((*data)->heter_rm_alphag,(unsigned)((*data)->heter_rm_alphag[0]+1)*sizeof(double));
		(*data)->invariable_mut_sites = (long int *)realloc((*data)->invariable_mut_sites,(unsigned)((*data)->invariable_mut_sites[0]+1)*sizeof(long int));
		(*data)->thetant_min = (double *)realloc((*data)->thetant_min,(unsigned)((*data)->thetant_min[0]+1)*sizeof(double));
		(*data)->thetant_max = (double *)realloc((*data)->thetant_max,(unsigned)((*data)->thetant_max[0]+1)*sizeof(double));
		(*data)->recnt_min = (double *)realloc((*data)->recnt_min,(unsigned)((*data)->recnt_min[0]+1)*sizeof(double));
		(*data)->recnt_max = (double *)realloc((*data)->recnt_max,(unsigned)((*data)->recnt_max[0]+1)*sizeof(double));
		(*data)->ancestral_pol = (int *)realloc((*data)->ancestral_pol,(unsigned)((*data)->ancestral_pol[0]+1)*sizeof(int));
		(*data)->ehh_fixnt = (long int *)realloc((*data)->ehh_fixnt,(unsigned)((*data)->ehh_fixnt[0]+1)*sizeof(long int));
		(*data)->event_sexratio = (double *)realloc((*data)->event_sexratio,(unsigned)((*data)->event_sexratio[0]+1)*sizeof(double));
		(*data)->fixmax_abstime_event = (double *)realloc((*data)->fixmax_abstime_event,(unsigned)((*data)->fixmax_abstime_event[0]+1)*sizeof(double));
	}
	MPI_Bcast((*data)->seed1,(*data)->seed1[0]+1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->npop_sampled,(*data)->npop_sampled[0]+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ratio_sv,(*data)->ratio_sv[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->r,(*data)->r[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);	
	MPI_Bcast((*data)->f,(*data)->f[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->track_len,(*data)->track_len[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->theta_1,(*data)->theta_1[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->mutations,(*data)->mutations[0]+1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->factor_pop,(*data)->factor_pop[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ifselection,(*data)->ifselection[0]+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->pop_sel,(*data)->pop_sel[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->sel_nt,(*data)->sel_nt[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((*data)->sinit,(*data)->sinit[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((*data)->sendt,(*data)->sendt[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((*data)->sfreqend,(*data)->sfreqend[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((*data)->sfreqinit,(*data)->sfreqinit[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ts,(*data)->ts[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->nintn,(*data)->nintn[0]+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->freq,(*data)->freq[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->p_gamma,(*data)->p_gamma[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->alpha_gamma,(*data)->alpha_gamma[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->correct_gamma,(*data)->correct_gamma[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);	
	MPI_Bcast((*data)->factor_chrn,(*data)->factor_chrn[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->p_gammar,(*data)->p_gammar[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->alpha_gammar,(*data)->alpha_gammar[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->correct_gammar,(*data)->correct_gammar[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->Rm,(*data)->Rm[0]+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->nhapl,(*data)->nhapl[0]+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->heter_theta_alphag,(*data)->heter_theta_alphag[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->heter_rm_alphag,(*data)->heter_rm_alphag[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->invariable_mut_sites,(*data)->invariable_mut_sites[0]+1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->thetant_min,(*data)->thetant_min[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->thetant_max,(*data)->thetant_max[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->recnt_min,(*data)->recnt_min[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->recnt_max,(*data)->recnt_max[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ancestral_pol,(*data)->ancestral_pol[0]+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->ehh_fixnt,(*data)->ehh_fixnt[0]+1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->event_sexratio,(*data)->event_sexratio[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((*data)->fixmax_abstime_event,(*data)->fixmax_abstime_event[0]+1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	/* double pointers */
	
	/*linked*/
	if((*data)->linked > 1) {		
		if(my_rank != 0) {
			if((*data)->linked>=20) (*data)->loci_linked = (long int **)realloc((*data)->loci_linked,(*data)->linked*sizeof(long int *));	
			for(j=0;j<(*data)->linked;j++) {
				if(j<20) 
					(*data)->loci_linked[j] = (long int *)realloc((*data)->loci_linked[j],3*sizeof(long int));
				else 
					(*data)->loci_linked[j] = (long int *)malloc(3*sizeof(long int));
			}
		}
		for(j=0;j<(*data)->linked;j++) {
			MPI_Bcast((*data)->loci_linked[j],3,MPI_LONG,0,MPI_COMM_WORLD);
		}
	}
	
	/*ssize_pop*/
	if(my_rank != 0) {
		if((*data)->n_loci>=20) (*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,(*data)->n_loci*sizeof(int *));
		for(j=0;j<(*data)->n_loci;j++) {
			if(j<20)
				(*data)->ssize_pop[j] = (int *)realloc((*data)->ssize_pop[j],((*data)->npop+1)*sizeof(int));
			else
				(*data)->ssize_pop[j] = (int *)malloc(((*data)->npop+1)*sizeof(int));
		}
	}	
	for(j=0;j<(*data)->n_loci;j++) {
		MPI_Bcast((*data)->ssize_pop[j],(*data)->npop+1,MPI_INT,0,MPI_COMM_WORLD);
	}
	
	/*mig_rate_matrix*/
	if(my_rank != 0) {
		if((*data)->npop>=20) (*data)->mig_rate_matrix = (double **)realloc((*data)->mig_rate_matrix,(*data)->npop*sizeof(double *));
		for(j=0;j<(*data)->npop;j++) {
			if(j<20)
				(*data)->mig_rate_matrix[j] = (double *)realloc((*data)->mig_rate_matrix[j],((*data)->npop+1)*sizeof(double));
			else
				(*data)->mig_rate_matrix[j] = (double *)malloc(((*data)->npop+1)*sizeof(double));
		}
	}	
	for(j=0;j<(*data)->npop;j++) {
		MPI_Bcast((*data)->mig_rate_matrix[j],(*data)->npop+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	
	/*pop_event*/
	if(my_rank != 0) {
		if((*data)->npop_events>=20) (*data)->pop_event = (double **)realloc((*data)->pop_event,(*data)->npop_events*sizeof(double *));
		for(j=0;j<(*data)->npop_events;j++) {
			if(j<20)
				(*data)->pop_event[j] = (double *)realloc((*data)->pop_event[j],(6+2*(*data)->npop)*sizeof(double));
			else
				(*data)->pop_event[j] = (double *)malloc((6+2*(*data)->npop)*sizeof(double));
		}
	}
	for(j=0;j<(*data)->npop_events;j++) {
		MPI_Bcast((*data)->pop_event[j],6+2*(*data)->npop,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	
	/*nrec*/
	if(my_rank != 0) {
		if((*data)->npop>=20) (*data)->nrec = (double **)realloc((*data)->nrec,(*data)->npop*sizeof(double *));
		for(j=0;j<(*data)->npop;j++) {
			if(j<20)
				(*data)->nrec[j] = (double *)realloc((*data)->nrec[j],((*data)->nintn[j+1]+2)*sizeof(double));
			else
				(*data)->nrec[j] = (double *)malloc(((*data)->nintn[j+1]+2)*sizeof(double));
		}
	}
	for(j=0;j<(*data)->npop;j++) {
		MPI_Bcast((*data)->nrec[j],(*data)->nintn[j+1]+2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	
	/*tpast*/
	if(my_rank != 0) {
		if((*data)->npop>=20) (*data)->tpast = (double **)realloc((*data)->tpast,(*data)->npop*sizeof(double *));
		for(j=0;j<(*data)->npop;j++) {
			if(j<20)
				(*data)->tpast[j] = (double *)realloc((*data)->tpast[j],((*data)->nintn[j+1]+2)*sizeof(double));
			else
				(*data)->tpast[j] = (double *)malloc(((*data)->nintn[j+1]+2)*sizeof(double));
		}
	}
	for(j=0;j<(*data)->npop;j++) {
		MPI_Bcast((*data)->tpast[j],(*data)->nintn[j+1]+2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	
	/*obs_statistics*//*
	if(my_rank != 0) {
		for(i=0;i<NOBS_STATISTICS;i++) {
			if((*data)->n_loci>=20) (*data)->obs_statistics[i] = (double **)realloc((*data)->obs_statistics[i],(*data)->n_loci*sizeof(double *));
			for(j=0;j<(*data)->n_loci;j++) {
				if(j<20)
					(*data)->obs_statistics[i][j] = (double *)realloc((*data)->obs_statistics[i][j],((*data)->max_npop_sampled+1)*sizeof(double));
				else
					(*data)->obs_statistics[i][j] = (double *)malloc(((*data)->max_npop_sampled+1)*sizeof(double));
			}
		}
	}
	for(i=0;i<NOBS_STATISTICS;i++)
		for(j=0;j<(*data)->n_loci;j++)
			MPI_Bcast((*data)->obs_statistics[i][j],(*data)->max_npop_sampled+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	*/
	/* communicate simple variables from struct var_priors*/
	if((*data)->npriors) {
		if(my_rank != 0) {
			/*reallocate and initialize memory*/
			if(!(priors[0] = (struct var_priors *)realloc(priors[0],(unsigned)((*data)->npriors)*sizeof(struct var_priors))))
				perror("realloc error.49");
			for(j=0;j<(*data)->npriors;j++) {
				if(!(priors[0][j].priordist = (double *)calloc((long unsigned)(*data)->n_iter,sizeof(double)))) 
					perror("realloc error.99");
			}
		}
			
		for(j=0;j<(*data)->npriors;j++) {
			i=0;
			block_lengths[i++] = 30;
			block_lengths[i++] = 1;
			block_lengths[i++] = 1;		
			block_lengths[i++] = 1;
			block_lengths[i++] = 5;
			block_lengths[i++] = 1;
			block_lengths[i++] = 1000;
			
			i=0;
			typelist[i++] = MPI_CHAR;
			typelist[i++] = MPI_DOUBLE;
			typelist[i++] = MPI_INT;
			typelist[i++] = MPI_INT;
			typelist[i++] = MPI_DOUBLE;
			typelist[i++] = MPI_LONG;
			typelist[i++] = MPI_CHAR;

			i=0;
			displacements[i++] = 0;
			MPI_Get_address(&(priors[0][j].name_pr),&start_address);
			MPI_Get_address(&(priors[0][j].idprior),&address);
			displacements[i++] = address - start_address;
			MPI_Get_address(&(priors[0][j].kind_var),&address);
			displacements[i++] = address - start_address;
			MPI_Get_address(&(priors[0][j].kind_dist),&address);
			displacements[i++] = address - start_address;
			MPI_Get_address(&(priors[0][j].dist_par),&address);
			displacements[i++] = address - start_address;
			MPI_Get_address(&(priors[0][j].seed_prior),&address);
			displacements[i++] = address - start_address;
			MPI_Get_address(&(priors[0][j].prior_file),&address);
			displacements[i++] = address - start_address;
			
			MPI_Type_create_struct(7,block_lengths,displacements,typelist,&mesg_mpi_prior);
			MPI_Type_commit(&mesg_mpi_prior);
			MPI_Bcast(&(priors[0][j]),1,mesg_mpi_prior,0,MPI_COMM_WORLD);
		}
		
		/* communicate pointers from struct var_priors*/	
		for(j=0;j<(*data)->npriors;j++) {
			MPI_Bcast(priors[0][j].priordist,(*data)->n_iter,MPI_DOUBLE,0,MPI_COMM_WORLD);
		}
	}
	return;
}

void rankprint(struct var **data,struct var_priors **priors,int my_rank)
{
	int i,j/*,k*/;
	long int x;
	/**/	
	printf("rank: %d, pr_matrix: %d\n",my_rank,(*data)->pr_matrix);
	printf("rank: %d, mhits: %d\n",my_rank,(*data)->mhits);
	printf("rank: %d, T_out: %f\n",my_rank,(*data)->T_out);
	printf("rank: %d, n_iter: %ld\n",my_rank,(*data)->n_iter);
	printf("rank: %d, seed2: %ld\n",my_rank,(*data)->seed2);
	printf("rank: %d, n_loci: %d\n",my_rank,(*data)->n_loci);
	printf("rank: %d, linked: %d\n",my_rank,(*data)->linked);
	printf("rank: %d, despl: %ld\n",my_rank,(*data)->despl);
	printf("rank: %d, window: %ld\n",my_rank,(*data)->window);
	printf("rank: %d, npop: %ld\n",my_rank,(*data)->npop);
	printf("rank: %d, mig_rate: %f\n",my_rank,(*data)->mig_rate);
	printf("rank: %d, pop_size: %f\n",my_rank,(*data)->pop_size);
	printf("rank: %d, ran_factorpop: %d\n",my_rank,(*data)->ran_factorpop);
	printf("rank: %d, same_factorpop: %d\n",my_rank,(*data)->same_factorpop);
	printf("rank: %d, max_npop_sampled: %d\n",my_rank,(*data)->max_npop_sampled);
	printf("rank: %d, neutral_tests: %d\n",my_rank,(*data)->neutral_tests);
	printf("rank: %d, print_neuttests: %d\n",my_rank,(*data)->print_neuttest);
	printf("rank: %d, sfix_allthetas: %d\n",my_rank,(*data)->sfix_allthetas);
	printf("rank: %d, sfix_pop: %d\n",my_rank,(*data)->sfix_pop);
	printf("rank: %d, mc_jump: %d\n",my_rank,(*data)->mc_jump);
	printf("rank: %d, range_thetant: %d\n",my_rank,(*data)->range_thetant);
	printf("rank: %d, ifgamma: %d\n",my_rank,(*data)->ifgamma);
	printf("rank: %d, iflogistic: %d\n",my_rank,(*data)->iflogistic);
	printf("rank: %d, split_pop: %d\n",my_rank,(*data)->split_pop);
	printf("rank: %d, time_split: %f\n",my_rank,(*data)->time_split);
	printf("rank: %d, time_scoal: %f\n",my_rank,(*data)->time_scoal);
	printf("rank: %d, factor_anc: %f\n",my_rank,(*data)->factor_anc);
	printf("rank: %d, tlimit: %f\n",my_rank,(*data)->tlimit);
	printf("rank: %d, npoprefugia: %d\n",my_rank,(*data)->npoprefugia);
	printf("rank: %d, rmfix: %d\n",my_rank,(*data)->rmfix);
	printf("rank: %d, rmfix_pop: %d\n",my_rank,(*data)->rmfix_pop);
	printf("rank: %d, method_samp: %d\n",my_rank,(*data)->method_samp);
	printf("rank: %d, range_rnt: %d\n",my_rank,(*data)->range_rnt);
	printf("rank: %d, ifgammar: %d\n",my_rank,(*data)->ifgammar);
	printf("rank: %d, no_rec_males: %d\n",my_rank,(*data)->no_rec_males);
	printf("rank: %d, linked_segsites: %d\n",my_rank,(*data)->linked_segsites);
	printf("rank: %d, linked_rm: %d\n",my_rank,(*data)->linked_rm);
	printf("rank: %d, linked_nhapl: %d\n",my_rank,(*data)->linked_nhapl);
	printf("rank: %d, linked_segsites_nregion: %d\n",my_rank,(*data)->linked_segsites_nregion);
	printf("rank: %d, linked_rm_nregion: %d\n",my_rank,(*data)->linked_rm_nregion);
	printf("rank: %d, linked_nhapl_nregion: %d\n",my_rank,(*data)->linked_nhapl_nregion);
	printf("rank: %d, likelihood_line: %d\n",my_rank,(*data)->likelihood_line);
	printf("rank: %d, npriors: %d\n",my_rank,(*data)->npriors);
	printf("rank: %d, npriors: %d\n",my_rank,(*data)->includeEHHrel);
	printf("rank: %d, npriors: %d\n",my_rank,(*data)->pop_outgroup);
	printf("rank: %d, sex_ratio: %f\n",my_rank,(*data)->sex_ratio);

	for(i=0;i<=NOBS_STATISTICS;i++) printf("rank: %d, obs_statistics_used[%d]: %d\n",my_rank,i,(*data)->obs_statistics_used[i]);
	for(i=0;i<=NOBS_STATISTICS;i++) printf("rank: %d, likelihood_error[%d]: %f\n",my_rank,i,(*data)->likelihood_error[i]);
	for(i=0;i<=5;i++) printf("rank: %d, patcg[%d]: %f\n",my_rank,i,(*data)->patcg[i]);

	for(i=0;i<(*data)->nsam[0]+1;i++) printf("rank: %d, nsam[%d]: %d\n",my_rank,i,(*data)->nsam[i]);
	for(i=0;i<(*data)->nsites[0]+1;i++) printf("rank: %d, nsites[%d]: %ld\n",my_rank,i,(*data)->nsites[i]);
	for(i=0;i<(*data)->seed1[0]+1;i++) printf("rank: %d, seed1[%d]: %ld\n",my_rank,i,(*data)->seed1[i]);
	for(i=0;i<(*data)->npop_sampled[0]+1;i++) printf("rank: %d, npop_sampled[%d]: %d\n",my_rank,i,(*data)->npop_sampled[i]);
	for(i=0;i<(*data)->ratio_sv[0]+1;i++) printf("rank: %d, ratio_sv[%d]: %f\n",my_rank,i,(*data)->ratio_sv[i]);
	for(i=0;i<(*data)->r[0]+1;i++) printf("rank: %d, r[%d]: %f\n",my_rank,i,(*data)->r[i]);
	for(i=0;i<(*data)->f[0]+1;i++) printf("rank: %d, f[%d]: %f\n",my_rank,i,(*data)->f[i]);
	for(i=0;i<(*data)->track_len[0]+1;i++) printf("rank: %d, track_len[%d]: %f\n",my_rank,i,(*data)->track_len[i]);
	for(i=0;i<(*data)->theta_1[0]+1;i++) printf("rank: %d, theta_1[%d]: %f\n",my_rank,i,(*data)->theta_1[i]);
	for(i=0;i<(*data)->mutations[0]+1;i++) printf("rank: %d, mutations[%d]: %ld\n",my_rank,i,(*data)->mutations[i]);
	for(i=0;i<(*data)->factor_pop[0]+1;i++) printf("rank: %d, factor_pop[%d]: %f\n",my_rank,i,(*data)->factor_pop[i]);
	for(i=0;i<(*data)->ifselection[0]+1;i++) printf("rank: %d, ifselection[%d]: %d\n",my_rank,i,(*data)->ifselection[i]);
	for(i=0;i<(*data)->pop_sel[0]+1;i++) printf("rank: %d, pop_sel[%d]: %f\n",my_rank,i,(*data)->pop_sel[i]);
	for(i=0;i<(*data)->sel_nt[0]+1;i++) printf("rank: %d, sel_nt[%d]: %f\n",my_rank,i,(*data)->sel_nt[i]);
    for(i=0;i<(*data)->sinit[0]+1;i++) printf("rank: %d, sinit[%d]: %f\n",my_rank,i,(*data)->sinit[i]);
    for(i=0;i<(*data)->sendt[0]+1;i++) printf("rank: %d, sendt[%d]: %f\n",my_rank,i,(*data)->sendt[i]);
    for(i=0;i<(*data)->sfreqend[0]+1;i++) printf("rank: %d, sfreqend[%d]: %f\n",my_rank,i,(*data)->sfreqend[i]);
    for(i=0;i<(*data)->sfreqinit[0]+1;i++) printf("rank: %d, sfreqinit[%d]: %f\n",my_rank,i,(*data)->sfreqinit[i]);
	for(i=0;i<(*data)->ts[0]+1;i++) printf("rank: %d, ts[%d]: %f\n",my_rank,i,(*data)->ts[i]);
	for(i=0;i<(*data)->nintn[0]+1;i++) printf("rank: %d, nintn[%d]: %d\n",my_rank,i,(*data)->nintn[i]);
	for(i=0;i<(*data)->freq[0]+1;i++) printf("rank: %d, freq[%d]: %f\n",my_rank,i,(*data)->freq[i]);
	for(i=0;i<(*data)->p_gamma[0]+1;i++) printf("rank: %d, p_gamma[%d]: %f\n",my_rank,i,(*data)->p_gamma[i]);
	for(i=0;i<(*data)->alpha_gamma[0]+1;i++) printf("rank: %d, alpha_gamma[%d]: %f\n",my_rank,i,(*data)->alpha_gamma[i]);
	for(i=0;i<(*data)->correct_gamma[0]+1;i++) printf("rank: %d, correct_gamma[%d]: %f\n",my_rank,i,(*data)->correct_gamma[i]);
	for(i=0;i<(*data)->factor_chrn[0]+1;i++) printf("rank: %d, factor_chrn[%d]: %f\n",my_rank,i,(*data)->factor_chrn[i]);
	for(i=0;i<(*data)->p_gammar[0]+1;i++) printf("rank: %d, p_gammar[%d]: %f\n",my_rank,i,(*data)->p_gammar[i]);
	for(i=0;i<(*data)->alpha_gammar[0]+1;i++) printf("rank: %d, alpha_gammar[%d]: %f\n",my_rank,i,(*data)->alpha_gammar[i]);
	for(i=0;i<(*data)->correct_gammar[0]+1;i++) printf("rank: %d, correct_gammar[%d]: %f\n",my_rank,i,(*data)->correct_gammar[i]);
	for(i=0;i<(*data)->Rm[0]+1;i++) printf("rank: %d, Rm[%d]: %d\n",my_rank,i,(*data)->Rm[i]);
	for(i=0;i<(*data)->nhapl[0]+1;i++) printf("rank: %d, nhapl[%d]: %d\n",my_rank,i,(*data)->nhapl[i]);
	for(i=0;i<(*data)->heter_theta_alphag[0]+1;i++) printf("rank: %d, heter_theta_alphag[%d]: %f\n",my_rank,i,(*data)->heter_theta_alphag[i]);
	for(i=0;i<(*data)->heter_rm_alphag[0]+1;i++) printf("rank: %d, heter_rm_alphag[%d]: %f\n",my_rank,i,(*data)->heter_rm_alphag[i]);
	for(i=0;i<(*data)->invariable_mut_sites[0]+1;i++) printf("rank: %d, invariable_mut_sites[%d]: %ld\n",my_rank,i,(*data)->invariable_mut_sites[i]);
	for(i=0;i<(*data)->thetant_min[0]+1;i++) printf("rank: %d, thetant_min[%d]: %f\n",my_rank,i,(*data)->thetant_min[i]);
	for(i=0;i<(*data)->thetant_max[0]+1;i++) printf("rank: %d, thetant_max[%d]: %f\n",my_rank,i,(*data)->thetant_max[i]);
	for(i=0;i<(*data)->recnt_min[0]+1;i++) printf("rank: %d, recnt_min[%d]: %f\n",my_rank,i,(*data)->recnt_min[i]);
	for(i=0;i<(*data)->recnt_max[0]+1;i++) printf("rank: %d, recnt_max[%d]: %f\n",my_rank,i,(*data)->recnt_max[i]);
	for(i=0;i<(*data)->ancestral_pol[0]+1;i++) printf("rank: %d, recnt_max[%d]: %d\n",my_rank,i,(*data)->ancestral_pol[i]);
	for(i=0;i<(*data)->ehh_fixnt[0]+1;i++) printf("rank: %d, recnt_max[%d]: %ld\n",my_rank,i,(*data)->ehh_fixnt[i]);
	for(i=0;i<(*data)->event_sexratio[0]+1;i++) printf("rank: %d, event_sexratio[%d]: %f\n",my_rank,i,(*data)->event_sexratio[i]);
	for(i=0;i<(*data)->fixmax_abstime_event[0]+1;i++) printf("rank: %d, fixmax_abstime_event[%d]: %f\n",my_rank,i,(*data)->fixmax_abstime_event[i]);
	
	for(i=0;i<(*data)->linked;i++) for(j=0;j<3;j++) printf("rank: %d, loci_linked[%d][%d]: %ld\n",my_rank,i,j,(*data)->loci_linked[i][j]);
	for(i=0;i<(*data)->n_loci;i++) for(j=0;j<(*data)->npop+1;j++) printf("rank: %d, ssize_pop[%d][%d]: %d\n",my_rank,i,j,(*data)->ssize_pop[i][j]);
	for(i=0;i<(*data)->npop;i++) for(j=0;j<(*data)->npop+1;j++) printf("rank: %d, mig_rate_matrix[%d][%d]: %f\n",my_rank,i,j,(*data)->mig_rate_matrix[i][j]);
	for(i=0;i<(*data)->npop_events;i++) for(j=0;j<6+2*(*data)->npop;j++) printf("rank: %d, pop_event[%d][%d]: %f\n",my_rank,i,j,(*data)->pop_event[i][j]);
	for(i=0;i<(*data)->npop;i++) for(j=0;j<(*data)->nintn[i+1]+2;j++) printf("rank: %d, nrec[%d][%d]: %f\n",my_rank,i,j,(*data)->nrec[i][j]);
	for(i=0;i<(*data)->npop;i++) for(j=0;j<(*data)->nintn[i+1]+2;j++) printf("rank: %d, tpast[%d][%d]: %f\n",my_rank,i,j,(*data)->tpast[i][j]);
	/**//*for(i=0;i<NOBS_STATISTICS;i++) for(j=0;j<(*data)->n_loci;j++) for(k=0;k<(*data)->max_npop_sampled+1;k++) printf("rank: %d, obs_statistics[%d][%d][%d]: %f\n",my_rank,i,j,k,(*data)->obs_statistics[i][j][k]);*/
	
	for(i=0;i<(*data)->npriors;i++) {
		printf("rank: %d, name_pr[%d]: %s\n",my_rank,i,priors[0][i].name_pr);
		printf("rank: %d, idprior[%d]: %f\n",my_rank,i,priors[0][i].idprior);
		printf("rank: %d, kind_var[%d]: %d\n",my_rank,i,priors[0][i].kind_var);
		printf("rank: %d, kind_dist[%d]: %d\n",my_rank,i,priors[0][i].kind_dist);
		for(j=0;j<5;j++) printf("rank: %d, dist_par[%d][%d]: %f\n",my_rank,i,j,priors[0][i].dist_par[j]);
		printf("rank: %d, seed_prior[%d]: %ld\n",my_rank,i,priors[0][i].seed_prior);
		for(x=0;x<(*data)->n_iter;x++) printf("rank: %d, priordist[%d][%ld]: %f\n",my_rank,i,x,priors[0][i].priordist[x]);
		printf("rank: %d, prior_file[%d]: %s\n",my_rank,i,priors[0][i].prior_file);
	}
	
	return;
}
#endif

void rankprint_inputp(struct var2 **inputp,int my_rank,int x)
{
	int i,j;
	
	printf("rank: %d, locus: %d, rmfix: %d\n",my_rank,x,(*inputp)->rmfix);
	printf("rank: %d, locus: %d, rmfix_pop: %d\n",my_rank,x,(*inputp)->rmfix_pop);
	printf("rank: %d, locus: %d, Rm: %d\n",my_rank,x,(*inputp)->Rm);
	printf("rank: %d, locus: %d, nhapl: %d\n",my_rank,x,(*inputp)->nhapl);
	printf("rank: %d, locus: %d, range_rnt: %d\n",my_rank,x,(*inputp)->range_rnt);
	printf("rank: %d, locus: %d, recnt_min: %f\n",my_rank,x,(*inputp)->recnt_min);
	printf("rank: %d, locus: %d, recnt_max: %f\n",my_rank,x,(*inputp)->recnt_max);
	printf("rank: %d, locus: %d, ifgammar: %d\n",my_rank,x,(*inputp)->ifgammar);
	printf("rank: %d, locus: %d, p_gammar: %f\n",my_rank,x,(*inputp)->p_gammar);
	printf("rank: %d, locus: %d, alpha_gammar: %f\n",my_rank,x,(*inputp)->alpha_gammar);
	printf("rank: %d, locus: %d, correct_gammar: %f\n",my_rank,x,(*inputp)->correct_gammar);
	printf("rank: %d, locus: %d, range_thetant: %d\n",my_rank,x,(*inputp)->range_thetant);
	printf("rank: %d, locus: %d, thetant_min: %f\n",my_rank,x,(*inputp)->thetant_min);
	printf("rank: %d, locus: %d, thetant_max: %f\n",my_rank,x,(*inputp)->thetant_max);
	printf("rank: %d, locus: %d, ifgamma: %d\n",my_rank,x,(*inputp)->ifgamma);
	printf("rank: %d, locus: %d, p_gamma: %f\n",my_rank,x,(*inputp)->p_gamma);
	printf("rank: %d, locus: %d, alpha_gamma: %f\n",my_rank,x,(*inputp)->alpha_gamma);
	printf("rank: %d, locus: %d, correct_gamma: %f\n",my_rank,x,(*inputp)->correct_gamma);
	printf("rank: %d, locus: %d, Sfix_alltheta: %d\n",my_rank,x,(*inputp)->Sfix_alltheta);
	printf("rank: %d, locus: %d, Sfix_pop: %d\n",my_rank,x,(*inputp)->Sfix_pop);
	printf("rank: %d, locus: %d, mc_jump: %d\n",my_rank,x,(*inputp)->mc_jump);
	printf("rank: %d, locus: %d, howmany: %ld\n",my_rank,x,(*inputp)->howmany);
	printf("rank: %d, locus: %d, npop: %ld\n",my_rank,x,(*inputp)->npop);
	printf("rank: %d, locus: %d, migrate: %f\n",my_rank,x,(*inputp)->migrate);
	printf("rank: %d, locus: %d, nsam: %d\n",my_rank,x,(*inputp)->nsam);
	printf("rank: %d, locus: %d, nsites: %ld\n",my_rank,x,(*inputp)->nsites);
	printf("rank: %d, locus: %d, r: %f\n",my_rank,x,(*inputp)->r);
	printf("rank: %d, locus: %d, f: %f\n",my_rank,x,(*inputp)->f);
	printf("rank: %d, locus: %d, track_len: %f\n",my_rank,x,(*inputp)->track_len);
	printf("rank: %d, locus: %d, theta: %f\n",my_rank,x,(*inputp)->theta);
	printf("rank: %d, locus: %d, segsitesin: %ld\n",my_rank,x,(*inputp)->segsitesin);
	if((*inputp)->npop>1)
		for(i=0;i<(*inputp)->npop;i++) 
			for(j=0;j<(*inputp)->npop;j++) 
				printf("rank: %d, locus: %d, migrate_matrix[%d][%d]: %f\n",my_rank,x,i,j,(*inputp)->migrate_matrix[i][j]);
	printf("rank: %d, locus: %d, linked: %d\n",my_rank,x,(*inputp)->linked);
	printf("rank: %d, locus: %d, despl: %ld\n",my_rank,x,(*inputp)->despl);
	printf("rank: %d, locus: %d, window: %ld\n",my_rank,x,(*inputp)->window);
	for(i=0;i<(*inputp)->linked;i++) 
		for(j=0;j<2;j++) 
			printf("rank: %d, locus: %d, loci_linked[%d][%d]: %ld\n",my_rank,x,i,j,(*inputp)->loci_linked[i][j]);
	printf("rank: %d, locus: %d, tloci: %d\n",my_rank,x,(*inputp)->tloci);	
	printf("rank: %d, locus: %d, pr_matrix: %d\n",my_rank,x,(*inputp)->pr_matrix);
	printf("rank: %d, locus: %d, mhits: %d\n",my_rank,x,(*inputp)->mhits);
	printf("rank: %d, locus: %d, T_out: %f\n",my_rank,x,(*inputp)->T_out);
	printf("rank: %d, locus: %d, ratio_sv: %f\n",my_rank,x,(*inputp)->ratio_sv);
	if((*inputp)->npop>1)
		for(i=0;i<(*inputp)->npop;i++) 
			printf("rank: %d, locus: %d, factor_pop[%d]: %f\n",my_rank,x,i,(*inputp)->factor_pop[i]);
	printf("rank: %d, locus: %d, ran_factorpop: %d\n",my_rank,x,(*inputp)->ran_factorpop);
	printf("rank: %d, locus: %d, same_factorpop: %d\n",my_rank,x,(*inputp)->same_factorpop);
	printf("rank: %d, locus: %d, npop_sampled: %d\n",my_rank,x,(*inputp)->npop_sampled);
	printf("rank: %d, locus: %d, max_npop_sampled: %d\n",my_rank,x,(*inputp)->max_npop_sampled);
	printf("rank: %d, locus: %d, neutral_tests: %d\n",my_rank,x,(*inputp)->neutral_tests);
	printf("rank: %d, locus: %d, print_neuttests: %d\n",my_rank,x,(*inputp)->print_neuttest);
	printf("rank: %d, locus: %d, iflogistic: %d\n",my_rank,x,(*inputp)->iflogistic);
	for(i=0;i<(*inputp)->npop;i++) 
		printf("rank: %d, locus: %d, ts[%d]: %f\n",my_rank,x,i,(*inputp)->ts[i]);
	for(i=0;i<(*inputp)->npop;i++) 
		printf("rank: %d, locus: %d, nintn[%d]: %d\n",my_rank,x,i,(*inputp)->nintn[i]);
	for(i=0;i<(*inputp)->npop;i++) 
		for(j=0;j<(*inputp)->nintn[i]+1;j++) 
			printf("rank: %d, locus: %d, nrec[%d][%d]: %f\n",my_rank,x,i,j,(*inputp)->nrec[i][j]);
	for(i=0;i<(*inputp)->npop;i++) 
		for(j=0;j<(*inputp)->nintn[i]+1;j++) 
			printf("rank: %d, locus: %d, tpast[%d][%d]: %f\n",my_rank,x,i,j,(*inputp)->tpast[i][j]);
	printf("rank: %d, locus: %d, split_pop: %d\n",my_rank,x,(*inputp)->split_pop);
	printf("rank: %d, locus: %d, time_split: %f\n",my_rank,x,(*inputp)->time_split);
	printf("rank: %d, locus: %d, time_scoal: %f\n",my_rank,x,(*inputp)->time_scoal);
	if((*inputp)->npop_events) {
		for(i=0;i<(*inputp)->npop_events;i++) {
			printf("rank: %d, event: %d, pop_event[%d].number: %d\n",my_rank,x,i,(int)(*inputp)->pop_event[i].number);
			printf("rank: %d, event: %d, pop_event[%d].npop1: %d\n",my_rank,x,i,(int)(*inputp)->pop_event[i].npop1);
			printf("rank: %d, event: %d, pop_event[%d].npop2: %d\n",my_rank,x,i,(int)(*inputp)->pop_event[i].npop2);
			printf("rank: %d, event: %d, pop_event[%d].factor_pop: %f\n",my_rank,x,i,(*inputp)->pop_event[i].factor_pop);
			printf("rank: %d, event: %d, pop_event[%d].time_event: %f\n",my_rank,x,i,(*inputp)->pop_event[i].time_event);
			printf("rank: %d, event: %d, pop_event[%d].timeS_event: %f\n",my_rank,x,i,(*inputp)->pop_event[i].timeS_event);
			for(j=0;j<(*inputp)->npop;j++) 
				printf("rank: %d, event: %d, pop_event[%d].mig_fw[%d]: %f\n",my_rank,x,i,j,(*inputp)->pop_event[i].mig_fw[j]);
			for(j=0;j<(*inputp)->npop;j++) 
				printf("rank: %d, event: %d, pop_event[%d].mig_rv[%d]: %f\n",my_rank,x,i,j,(*inputp)->pop_event[i].mig_rv[j]);
		}
		for(i=0;i<(*inputp)->fixmax_abstime_event[0]+1;i++) 
			printf("rank: %d, event: %d, fixmax_abstime_event[%d]: %f\n",my_rank,x,i,(*inputp)->fixmax_abstime_event[i]);
	}
	printf("rank: %d, locus: %d, factor_anc: %f\n",my_rank,x,(*inputp)->factor_anc);
	if((*inputp)->split_pop)
		for(i=0;i<(*inputp)->npop;i++) 
			printf("rank: %d, locus: %d, freq[%d]: %f\n",my_rank,x,i,(*inputp)->freq[i]);
	printf("rank: %d, locus: %d, tlimit: %f\n",my_rank,x,(*inputp)->tlimit);
	printf("rank: %d, locus: %d, factor_chrn: %f\n",my_rank,x,(*inputp)->factor_chrn);
	printf("rank: %d, locus: %d, sex_ratio: %f\n",my_rank,x,(*inputp)->sex_ratio);
	printf("rank: %d, locus: %d, no_rec_males: %d\n",my_rank,x,(*inputp)->no_rec_males);
	printf("rank: %d, locus: %d, linked_segsites: %d\n",my_rank,x,(*inputp)->linked_segsites);
	printf("rank: %d, locus: %d, linked_rm: %d\n",my_rank,x,(*inputp)->linked_rm);
	printf("rank: %d, locus: %d, linked_nhapl: %d\n",my_rank,x,(*inputp)->linked_nhapl);
	printf("rank: %d, locus: %d, linked_segsites_nregion: %d\n",my_rank,x,(*inputp)->linked_segsites_nregion);
	printf("rank: %d, locus: %d, linked_rm_nregion: %d\n",my_rank,x,(*inputp)->linked_rm_nregion);
	printf("rank: %d, locus: %d, linked_nhapl_nregion: %d\n",my_rank,x,(*inputp)->linked_nhapl_nregion);
	printf("rank: %d, locus: %d, heter_theta_alphag: %f\n",my_rank,x,(*inputp)->heter_theta_alphag);
	printf("rank: %d, locus: %d, heter_rm_alphag: %f\n",my_rank,x,(*inputp)->heter_rm_alphag);
	printf("rank: %d, locus: %d, invariable_mut_sites: %ld\n",my_rank,x,(*inputp)->invariable_mut_sites);
	printf("rank: %d, locus: %d, likelihood_line: %d\n",my_rank,x,(*inputp)->likelihood_line);
	printf("rank: %d, locus: %d, npriors: %d\n",my_rank,x,(*inputp)->npriors);
	printf("rank: %d, locus: %d, ifselection: %d\n",my_rank,x,(*inputp)->ifselection);
	printf("rank: %d, locus: %d, pop_size: %f\n",my_rank,x,(*inputp)->pop_size);
	printf("rank: %d, locus: %d, pop_sel: %f\n",my_rank,x,(*inputp)->pop_sel);
	printf("rank: %d, locus: %d, sel_nt: %f\n",my_rank,x,(*inputp)->sel_nt);
    printf("rank: %d, locus: %d, sinit: %f\n",my_rank,x,(*inputp)->sinit);
    printf("rank: %d, locus: %d, sendt: %f\n",my_rank,x,(*inputp)->sendt);
    printf("rank: %d, locus: %d, sfreqend: %f\n",my_rank,x,(*inputp)->sfreqend);
    printf("rank: %d, locus: %d, sfreqend: %f\n",my_rank,x,(*inputp)->sfreqinit);
	for(i=0;i<(*inputp)->npop;i++)
		printf("rank: %d, locus: %d, config[%d]: %d\n",my_rank,x,i,(*inputp)->config[i]);
	for(i=0;i<4;i++) 
		printf("rank: %d, locus: %d, patcg[%d]: %g\n",my_rank,x,i,(*inputp)->patcg[i]);
	printf("rank: %d, locus: %d, includeEHHrel: %d\n",my_rank,x,(*inputp)->includeEHHrel);
	printf("rank: %d, locus: %d, pop_outgroup: %d\n",my_rank,x,(*inputp)->pop_outgroup);

	return;
}

