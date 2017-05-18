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
#include <stdlib.h>
#include <stdio.h>

/* fer una nova estructura que contingui totes les variables en proces */
void getpars_fix(struct var **data, struct var2 **inputp)
{    
    int x,y,z,ap;
	/*double maxtpast;*/
    double ran1(void);
	double tmax;

	(*inputp)->includeEHHrel = (*data)->includeEHHrel;
	(*inputp)->max_npop_sampled = (*data)->max_npop_sampled;
	
	(*inputp)->Sfix_alltheta = (*data)->sfix_allthetas;
	(*inputp)->Sfix_pop = (*data)->sfix_pop ;
	(*inputp)->mc_jump = (*data)->mc_jump;
	
	(*inputp)->rmfix = (*data)->rmfix;
	(*inputp)->rmfix_pop = (*data)->rmfix_pop ;
	
    (*inputp)->howmany = (*data)->n_iter;			
    (*inputp)->pr_matrix = (*data)->pr_matrix;						
    (*inputp)->mhits = (*data)->mhits;					
	
    (*inputp)->tloci = (*data)->n_loci;
	
    if((*data)->split_pop == 1) (*inputp)->npop = (*data)->npoprefugia;
	else (*inputp)->npop = (*data)->npop;
    (*inputp)->config = (int *)calloc((*inputp)->npop,sizeof(int));

    (*inputp)->neutral_tests = (*data)->neutral_tests;
    (*inputp)->print_neuttest = (*data)->print_neuttest;
	
    (*inputp)->tlimit = (*data)->tlimit;    
    (*inputp)->no_rec_males = (*data)->no_rec_males;
	
    (*inputp)->iflogistic = (*data)->iflogistic;
	
	(*inputp)->likelihood_line = (*data)->likelihood_line;
	
    (*inputp)->range_thetant = (*data)->range_thetant;
    (*inputp)->range_rnt = (*data)->range_rnt;
	(*inputp)->ifgamma = (*data)->ifgamma;
	(*inputp)->ifgammar = (*data)->ifgammar;
	
	if((*inputp)->tloci == 1 && (*data)->window == 0 && (*data)->linked < 2) (*inputp)->linked = 1; /*** Here I use a single locus like a linked loci for MPI ***/
	else (*inputp)->linked = (*data)->linked;

    if((*inputp)->linked > 1) {
        (*inputp)->loci_linked = (long int **)malloc((*inputp)->linked*sizeof(long int *));		
        for(x=0;x<(*inputp)->linked;x++) {
            (*inputp)->loci_linked[x] = (long int *)malloc(2*sizeof(long int));
            (*inputp)->loci_linked[x][0] = (*data)->loci_linked[x][1];
            (*inputp)->loci_linked[x][1] = (*data)->loci_linked[x][2];
        }
		(*inputp)->linked_segsites = (*data)->linked_segsites;
		(*inputp)->linked_rm = (*data)->linked_rm;
		(*inputp)->linked_nhapl = (*data)->linked_nhapl;
		(*inputp)->linked_segsites_nregion = (*data)->linked_segsites_nregion;
		(*inputp)->linked_rm_nregion = (*data)->linked_rm_nregion;
		(*inputp)->linked_nhapl_nregion = (*data)->linked_nhapl_nregion;
    }
	
	if((*inputp)->linked == 1) {
		if((*data)->window) {
			(*inputp)->despl = (*data)->despl;
			(*inputp)->window = (*data)->window;
		}
		else {
			(*inputp)->despl = (*data)->nsites[1];
			(*inputp)->window = (*data)->nsites[1]/*-1*/;
			
			(*inputp)->linked_segsites_nregion = -1;
			(*inputp)->linked_rm_nregion = -1;
			(*inputp)->linked_nhapl_nregion = -1;

			if((*data)->theta_1[0] == (double)1 && (*data)->theta_1[1] == (double)0) {		
				if((*data)->linked_segsites != -1) (*inputp)->linked_segsites = (*data)->linked_segsites;
				else (*inputp)->linked_segsites = (int)(*data)->mutations[1];
				(*inputp)->linked_segsites_nregion = 0;
			}
			if((*inputp)->Sfix_alltheta) {
				if((*data)->linked_segsites != -1) (*inputp)->linked_segsites = (*data)->linked_segsites;
				else (*inputp)->linked_segsites = (int)(*data)->mutations[1];
				(*inputp)->linked_segsites_nregion = 0;
			}		
			if((*inputp)->rmfix) {
				if((*data)->linked_rm != -1) (*inputp)->linked_rm = (*data)->linked_rm;
				else (*inputp)->linked_rm = (*data)->Rm[1];
				
				if((*data)->linked_nhapl != 0) (*inputp)->linked_nhapl = (*data)->linked_nhapl;
				else (*inputp)->linked_nhapl = (*data)->nhapl[1];
				
				(*inputp)->linked_rm_nregion = 0;
				(*inputp)->linked_nhapl_nregion = 0;		
			}
		}
	}
	else {
		(*inputp)->despl = 0;
		(*inputp)->window = 0;
	}
    
	(*inputp)->npriors = (*data)->npriors;
	if((*inputp)->npriors) {
		if(!((*inputp)->pointtoprior = (double **)calloc((*inputp)->npriors,sizeof(double *)))) {
			printf("Error allocating memory. get_mod 33.\n");
			exit(1);
		}
	}
	
	if((*data)->npop > 1 || (*data)->split_pop == 1 ) {
        (*inputp)->factor_pop = (double *) malloc(((*inputp)->npop)*sizeof(double));
        for(x=0;x<(*data)->factor_pop[0];x++) {
			(*inputp)->factor_pop[x] = (*data)->factor_pop[x+1];
			if((*data)->factor_pop[x+1]>=REFNUMBER && (*data)->factor_pop[x+1] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->factor_pop[x+1] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = (*inputp)->factor_pop+x;
			}
        }
		(*inputp)->ran_factorpop = (*data)->ran_factorpop;
        (*inputp)->same_factorpop = (*data)->same_factorpop;
		
        if((*inputp)->ran_factorpop == 1) {
            (*inputp)->factor_pop[0] = (double)1;
            for(x=1;x<(*inputp)->npop;x++) {
                (*inputp)->factor_pop[x] = (double)ran1()*(double)9 + (double)1;
                if((double)ran1() < 0.5) (*inputp)->factor_pop[x] = (double)1/(*inputp)->factor_pop[x];
            }
        }
        if((*inputp)->same_factorpop == 1) {
            for(x=0;x<(*inputp)->npop;x++) {
                (*inputp)->factor_pop[x] = (double)1;
            }
        }        
    }
    else {
        (*inputp)->factor_pop = (double *)malloc(1*sizeof(double));
        (*inputp)->factor_pop[0] = 1.;
        (*inputp)->ran_factorpop = 0;
        (*inputp)->same_factorpop = 1;
    }
	
	(*inputp)->ts = (double *)malloc(((*inputp)->npop)*sizeof(double));
	for(x=0;x<(*inputp)->npop;x++) {
		(*inputp)->ts[x] = (*data)->ts[x+1];
		if((*data)->ts[x+1]>=REFNUMBER && (*data)->ts[x+1] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->ts[x+1] - REFNUMBER1);
			(*inputp)->pointtoprior[z] = (*inputp)->ts+x;
		}
	}
    (*inputp)->nintn = (int *)malloc(((*inputp)->npop)*sizeof(int));
	for(x=0;x<(*inputp)->npop;x++) (*inputp)->nintn[x] = (*data)->nintn[x+1];
	
    (*inputp)->nrec  = (double **)malloc((*inputp)->npop*sizeof(double *));
    (*inputp)->tpast = (double **)malloc((*inputp)->npop*sizeof(double *));
    (*inputp)->tpastS = (double **)malloc((*inputp)->npop*sizeof(double *));
	
	/*maxtpast = 0.;*/
	for(x=0;x<(*inputp)->npop;x++) {
		(*inputp)->nrec[x] = (double *)malloc((2+(*inputp)->nintn[x])*sizeof(double));
		(*inputp)->tpast[x] = (double *)malloc((2+(*inputp)->nintn[x])*sizeof(double));
		(*inputp)->tpastS[x] = (double *)malloc((2+(*inputp)->nintn[x])*sizeof(double));
		
		(*inputp)->nrec[x][0] = (double)1;
		(*inputp)->nrec[x][1] = (double)1;
		(*inputp)->tpast[x][0] = (double)0;
		(*inputp)->tpast[x][1] = (double)1;
		(*inputp)->tpastS[x][0] = (double)0;
		(*inputp)->tpastS[x][1] = (double)1;
		
		for(y=1;y<=(*inputp)->nintn[x];y++) {
			(*inputp)->nrec[x][y] = (*data)->nrec[x][y];    
			if((*data)->nrec[x][y]>=REFNUMBER && (*data)->nrec[x][y] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->nrec[x][y] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->nrec[x][y]);
			}
			(*inputp)->tpastS[x][y] = (*data)->tpast[x][y];
			if((*data)->tpast[x][y]>=REFNUMBER && (*data)->tpast[x][y] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->tpast[x][y] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->tpastS[x][y]);
			}
			(*inputp)->tpast[x][y] = (*inputp)->tpastS[x][y];
		}
	}
	
    (*inputp)->pop_size = (*data)->pop_size;/*selection*/
	if((*data)->pop_size>=REFNUMBER && (*data)->pop_size < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->pop_size- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->pop_size);
	}
	
	(*inputp)->T_out = (*data)->T_out;					
	if((*data)->T_out>=REFNUMBER && (*data)->T_out < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->T_out- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->T_out);
	}
	
    (*inputp)->split_pop = (*data)->split_pop;
	
    (*inputp)->time_split = (*data)->time_split;
	if((*data)->time_split>=REFNUMBER && (*data)->time_split < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->time_split- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->time_split);
	}
	
	(*inputp)->time_scoalS = (*data)->time_scoal;
	if((*data)->time_scoal>=REFNUMBER && (*data)->time_scoal < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->time_scoal - REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->time_scoalS);
	}
	(*inputp)->time_scoal = (*inputp)->time_scoalS;
	
	(*inputp)->npop_events = (*data)->npop_events;
	if((*inputp)->npop_events) {
		(*inputp)->pop_event = (struct events *)malloc((*inputp)->npop_events*sizeof(struct events));
		for(x=0;x<(*inputp)->npop_events;x++) {
			(*inputp)->pop_event[x].sex_ratio = -1.;
		}
		if((*data)->event_sexratio[0]<=2*(*inputp)->npop_events) {
			for(x=1;x<(*data)->event_sexratio[0]+1;x+=2) {
				(*inputp)->pop_event[(int)(*data)->event_sexratio[x]].sex_ratio = (*data)->event_sexratio[x+1];
				if((*data)->event_sexratio[x+1]>=REFNUMBER && (*data)->event_sexratio[x+1] < REFNUMBER + PLUSPRIOR) {
					z = (int)((*data)->event_sexratio[x+1]- REFNUMBER1);
					(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[(int)(*data)->event_sexratio[x]].sex_ratio);
				}
			}
		}
		for(x=0;x<(*inputp)->npop_events;x++) {
			(*inputp)->pop_event[x].number = (int)(*data)->pop_event[x][1];
			(*inputp)->pop_event[x].npop1 = (int)(*data)->pop_event[x][2];
			(*inputp)->pop_event[x].npop2 = (int)(*data)->pop_event[x][3];
			(*inputp)->pop_event[x].factor_pop = (*data)->pop_event[x][4];
			if((*data)->pop_event[x][4]>=REFNUMBER && (*data)->pop_event[x][4] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->pop_event[x][4]- REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[x].factor_pop);
			}
			(*inputp)->pop_event[x].timeS_event = (*data)->pop_event[x][5];
			if((*data)->pop_event[x][5]>=REFNUMBER && (*data)->pop_event[x][5] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->pop_event[x][5] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[x].timeS_event);
			}
			(*inputp)->pop_event[x].time_event = (*inputp)->pop_event[x].timeS_event;

			(*inputp)->pop_event[x].mig_fw = (double *)malloc(((*inputp)->npop)*sizeof(double));
			(*inputp)->pop_event[x].mig_rv = (double *)malloc(((*inputp)->npop)*sizeof(double));
			for(y=0;y<(*inputp)->npop;y++) {
				(*inputp)->pop_event[x].mig_fw[y] = (*data)->pop_event[x][6+y];
				if((*data)->pop_event[x][6+y]>=REFNUMBER && (*data)->pop_event[x][6+y] < REFNUMBER + PLUSPRIOR) {
					z = (int)((*data)->pop_event[x][6+y]- REFNUMBER1);
					(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[x].mig_fw[y]);
				}
				(*inputp)->pop_event[x].mig_rv[y] = (*data)->pop_event[x][6+(*inputp)->npop+y];
				if((*data)->pop_event[x][6+(*inputp)->npop+y]>=REFNUMBER && (*data)->pop_event[x][6+(*inputp)->npop+y] < REFNUMBER + PLUSPRIOR) {
					z = (int)((*data)->pop_event[x][6+(*inputp)->npop+y] - REFNUMBER1);
					(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[x].mig_rv[y]);
				}
			}
		}
        (*inputp)->fixmax_abstime_event = (double *)malloc((*inputp)->npop_events*sizeof(double));		
        for(x=0;x<(*inputp)->npop_events;x+=1) {
			(*inputp)->fixmax_abstime_event[x] = 1000.;
			for(y=1;y<(*data)->fixmax_abstime_event[0]+1;y+=2) {
				if((*data)->fixmax_abstime_event[y] == x) {
					(*inputp)->fixmax_abstime_event[x] = (*data)->fixmax_abstime_event[y+1];
					break;
				}				
			}
		}
		tmax = (*inputp)->fixmax_abstime_event[(*inputp)->npop_events-1];
		for(x=(*inputp)->npop_events-2;x>=0;x-=1) {
			if((*inputp)->fixmax_abstime_event[x] > tmax) {
				(*inputp)->fixmax_abstime_event[x] = tmax;
			}
			else {
				tmax = (*inputp)->fixmax_abstime_event[x];
			}
		}
	}
	else {
		(*inputp)->pop_event = (struct events *)malloc(1*sizeof(struct events));
		(*inputp)->pop_event[0].number = 0;
		(*inputp)->pop_event[0].npop1 = 0;
		(*inputp)->pop_event[0].npop2 = 0;
		(*inputp)->pop_event[0].factor_pop = 1;
		(*inputp)->pop_event[0].time_event = 1E07;
		(*inputp)->pop_event[0].timeS_event = 1E07;
	}
    (*inputp)->factor_anc = (*data)->factor_anc;
	if((*data)->factor_anc>=REFNUMBER && (*data)->factor_anc < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->factor_anc- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->factor_anc);
	}
	
	if((*data)->split_pop) {
        (*inputp)->freq = (double *)malloc((*inputp)->npop*sizeof(double));
        for(x=0;x<(*inputp)->npop;x++) {
            (*inputp)->freq[x] = (*data)->freq[x+1];
			if((*data)->freq[x+1]>=REFNUMBER && (*data)->freq[x+1] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->freq[x+1] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = (*inputp)->freq+x;
			}
        }
    }
    (*inputp)->migrate = (*data)->mig_rate;    
	if((*data)->mig_rate>=REFNUMBER && (*data)->mig_rate < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->mig_rate- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->migrate);
	}
	
	(*inputp)->migrate_matrix = (double **)malloc((*inputp)->npop*sizeof(double *));
    for(x=0;x<(*inputp)->npop;x++) {
		(*inputp)->migrate_matrix[x] = (double *)malloc(((*inputp)->npop)*sizeof(double));
		for(y=0;y<(*inputp)->npop;y++) {
			(*inputp)->migrate_matrix[x][y] = (*data)->mig_rate_matrix[x][y+1];
			if((*data)->mig_rate_matrix[x][y+1]>=REFNUMBER && (*data)->mig_rate_matrix[x][y+1] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->mig_rate_matrix[x][y+1] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->migrate_matrix[x][y]);
			}
		}
	}	
	
	(*inputp)->type_ancestral = (*data)->ancestral_pol[1];
	if((*inputp)->type_ancestral) {
		(*inputp)->ancestral_pol = (int **)malloc((*inputp)->type_ancestral*sizeof(int *));
		x = 2;
		ap = (*data)->ancestral_pol[x];
		for(y=0;y<(*inputp)->type_ancestral;y++) {
			(*inputp)->ancestral_pol[y] = (int *)malloc(((*inputp)->npop+1)*sizeof(int));
			for(z=0;z<ap+1;z++,x++) {
				(*inputp)->ancestral_pol[y][z] = (*data)->ancestral_pol[x];
			}
		}
	}
    (*inputp)->pop_outgroup = (*data)->pop_outgroup;
	
	if((*data)->patcg[0] == 0) {
		for(y=0;y<4;y++) 
			(*inputp)->patcg[y] = -1.;
	}
	else {
		for(y=0;y<(*data)->patcg[0];y++) 
			(*inputp)->patcg[y] = (*data)->patcg[y+1];
	}

    (*inputp)->sex_ratio = (*data)->sex_ratio;
	if((*data)->sex_ratio>=REFNUMBER && (*data)->sex_ratio < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->sex_ratio- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->sex_ratio);
	}
	
	/*in order to detect locus-specific priors*/
	(*inputp)->ifehh_fix  = (int)(*data)->ehh_fixnt[1];
	(*inputp)->ehh_margin = (*data)->ehh_fixnt[2];
	
	if((*inputp)->npriors) {
		(*inputp)->locspec_prior = (int *)calloc((*inputp)->npriors,sizeof(int));
		for(x=0;x<(*inputp)->npriors;x++)
			if((*inputp)->pointtoprior[x] == (double *)0) /*locus-specific prior*/
				(*inputp)->locspec_prior[x] = -1;
	}
	
	return;
}    

void getpars_mod(struct var **data, struct var2 **inputp,int p1)
{
    int x,y,z;
		
    /*WE HAVE TO TAKE INTO ACCOUNT THE PRIOR DISTRIBUTIONS*/
	
	/*check for the locus-specific priors*/
	for(z=0;z<(*inputp)->npriors;z++) 
		if((*inputp)->locspec_prior[z] == -1) 
			(*inputp)->pointtoprior[z] = (double *)0;
	
	/*defining all parameters*/
    if((*data)->npop > 1 || (*data)->split_pop == 1 ) {
        for(x=0;x<(*data)->factor_pop[0];x++) {
			(*inputp)->factor_pop[x] = (*data)->factor_pop[x+1];
			if((*data)->factor_pop[x+1]>=REFNUMBER && (*data)->factor_pop[x+1] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->factor_pop[x+1] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = (*inputp)->factor_pop+x;
			}
        }
    }
	
	for(x=0;x<(*inputp)->npop;x++) {
		(*inputp)->ts[x] = (*data)->ts[x+1];
		if((*data)->ts[x+1]>=REFNUMBER && (*data)->ts[x+1] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->ts[x+1] - REFNUMBER1);
			(*inputp)->pointtoprior[z] = (*inputp)->ts+x;
		}
	}
 	for(x=0;x<(*inputp)->npop;x++) (*inputp)->nintn[x] = (*data)->nintn[x+1];
	
	/*maxtpast = 0.;*/
	for(x=0;x<(*inputp)->npop;x++) {
		(*inputp)->nrec[x][0] = (double)1;
		(*inputp)->nrec[x][1] = (double)1;
		(*inputp)->tpast[x][0] = (double)0;
		(*inputp)->tpast[x][1] = (double)1;
		(*inputp)->tpastS[x][0] = (double)0;
		(*inputp)->tpastS[x][1] = (double)1;
		
		for(y=1;y<=(*inputp)->nintn[x];y++) {
			(*inputp)->nrec[x][y] = (*data)->nrec[x][y];    
			if((*data)->nrec[x][y]>=REFNUMBER && (*data)->nrec[x][y] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->nrec[x][y] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->nrec[x][y]);
			}
			(*inputp)->tpastS[x][y] = (*data)->tpast[x][y];
			if((*data)->tpast[x][y]>=REFNUMBER && (*data)->tpast[x][y] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->tpast[x][y] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->tpastS[x][y]);
			}
			(*inputp)->tpast[x][y] = (*inputp)->tpastS[x][y];
		}
	}
	
    (*inputp)->pop_size = (*data)->pop_size;/*selection*/
	if((*data)->pop_size>=REFNUMBER && (*data)->pop_size < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->pop_size- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->pop_size);
	}
	
	(*inputp)->T_out = (*data)->T_out;					
	if((*data)->T_out>=REFNUMBER && (*data)->T_out < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->T_out- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->T_out);
	}
	
    (*inputp)->time_split = (*data)->time_split;
	if((*data)->time_split>=REFNUMBER && (*data)->time_split < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->time_split- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->time_split);
	}
	
	(*inputp)->time_scoalS = (*data)->time_scoal;
	if((*data)->time_scoal>=REFNUMBER && (*data)->time_scoal < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->time_scoal- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->time_scoalS);
	}
	(*inputp)->time_scoal = (*inputp)->time_scoalS;
	
    for(x=0;x<(*inputp)->npop_events;x++) {
		(*inputp)->pop_event[x].factor_pop = (*data)->pop_event[x][4];
		if((*data)->pop_event[x][4]>=REFNUMBER && (*data)->pop_event[x][4] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->pop_event[x][4]- REFNUMBER1);
			(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[x].factor_pop);
		}
		(*inputp)->pop_event[x].timeS_event = (*data)->pop_event[x][5];
		if((*data)->pop_event[x][5]>=REFNUMBER && (*data)->pop_event[x][5] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->pop_event[x][5] - REFNUMBER1);
			(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[x].timeS_event);
		}
		(*inputp)->pop_event[x].time_event = (*inputp)->pop_event[x].timeS_event;
		
		for(y=0;y<(*inputp)->npop;y++) {
			(*inputp)->pop_event[x].mig_fw[y] = (*data)->pop_event[x][6+y];
			if((*data)->pop_event[x][6+y]>=REFNUMBER && (*data)->pop_event[x][6+y] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->pop_event[x][6+y]- REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[x].mig_fw[y]);
			}
			(*inputp)->pop_event[x].mig_rv[y] = (*data)->pop_event[x][6+(*inputp)->npop+y];
			if((*data)->pop_event[x][6+(*inputp)->npop+y]>=REFNUMBER && (*data)->pop_event[x][6+(*inputp)->npop+y] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->pop_event[x][6+(*inputp)->npop+y] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->pop_event[x].mig_rv[y]);
			}
		}
	}	
	
    (*inputp)->factor_anc = (*data)->factor_anc;
	if((*data)->factor_anc>=REFNUMBER && (*data)->factor_anc < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->factor_anc- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->factor_anc);
		/*printf("\nlocus: %d, z: %d, (*inputp)->pointtoprior[z]: %ld",p1,z,(*inputp)->pointtoprior[z]);*/
	}

    if((*data)->split_pop) {
        for(x=0;x<(*inputp)->npop;x++) {
            (*inputp)->freq[x] = (*data)->freq[x+1];
			if((*data)->freq[x+1]>=REFNUMBER && (*data)->freq[x+1] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->freq[x+1] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = (*inputp)->freq+x;
			}
        }
    }
    (*inputp)->migrate = (*data)->mig_rate;    
	if((*data)->mig_rate>=REFNUMBER && (*data)->mig_rate < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->mig_rate- REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->migrate);
	}
	
    for(x=0;x<(*inputp)->npop;x++) {
		for(y=0;y<(*inputp)->npop;y++) {
			(*inputp)->migrate_matrix[x][y] = (*data)->mig_rate_matrix[x][y+1];
			if((*data)->mig_rate_matrix[x][y+1]>=REFNUMBER && (*data)->mig_rate_matrix[x][y+1] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->mig_rate_matrix[x][y+1] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->migrate_matrix[x][y]);
			}
		}
	}	
	
	(*inputp)->nloci = p1;
    (*inputp)->nsam = (*data)->nsam[p1+1];
    (*inputp)->nsites = (*data)->nsites[p1+1];    
	    
	/*selection*/
    (*inputp)->ifselection = (*data)->ifselection[p1+1];
    if((*data)->ifselection[p1+1] == 1) {
        (*inputp)->npop = 2;
        (*inputp)->config = (int *)realloc((*inputp)->config,2*sizeof(int));
        (*inputp)->pop_sel = (*data)->pop_sel[p1+1];
		if((*data)->pop_sel[p1+1]>=REFNUMBER && (*data)->pop_sel[p1+1] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->pop_sel[p1+1] - REFNUMBER1);
			(*inputp)->pointtoprior[z] = &((*inputp)->pop_sel);
		}
        (*inputp)->sel_nt = (*data)->sel_nt[p1+1];
		if((*data)->sel_nt[p1+1]>=REFNUMBER && (*data)->sel_nt[p1+1] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->sel_nt[p1+1] - REFNUMBER1);
			(*inputp)->pointtoprior[z] = &((*inputp)->sel_nt);
		}
        (*inputp)->sinit = (*data)->sinit[p1+1];
        if((*data)->sinit[p1+1]>=REFNUMBER && (*data)->sinit[p1+1] < REFNUMBER + PLUSPRIOR) {
            z = (int)((*data)->sinit[p1+1] - REFNUMBER1);
            (*inputp)->pointtoprior[z] = &((*inputp)->sinit);
        }
        (*inputp)->sendt = (*data)->sendt[p1+1];
        if((*data)->sendt[p1+1]>=REFNUMBER && (*data)->sendt[p1+1] < REFNUMBER + PLUSPRIOR) {
            z = (int)((*data)->sendt[p1+1] - REFNUMBER1);
            (*inputp)->pointtoprior[z] = &((*inputp)->sendt);
        }
        (*inputp)->sfreqend = (*data)->sfreqend[p1+1];
        if((*data)->sfreqend[p1+1]>=REFNUMBER && (*data)->sfreqend[p1+1] < REFNUMBER + PLUSPRIOR) {
            z = (int)((*data)->sfreqend[p1+1] - REFNUMBER1);
            (*inputp)->pointtoprior[z] = &((*inputp)->sfreqend);
        }
        (*inputp)->sfreqinit = (*data)->sfreqinit[p1+1];
        if((*data)->sfreqinit[p1+1]>=REFNUMBER && (*data)->sfreqinit[p1+1] < REFNUMBER + PLUSPRIOR) {
            z = (int)((*data)->sfreqinit[p1+1] - REFNUMBER1);
            (*inputp)->pointtoprior[z] = &((*inputp)->sfreqinit);
        }

        /***************/
        (*inputp)->nintn = (int *)realloc((*inputp)->nintn,((*inputp)->npop)*sizeof(int));
        for(x=1;x<(*inputp)->npop;x++) (*inputp)->nintn[x] = 0;
        /***************/
    }
    if((*data)->ifselection[p1+1] == 0 && (p1 > 0 && (*data)->ifselection[p1] == 1)) {
        (*inputp)->npop = 1;
	}

    for(x=0;x<(*data)->ssize_pop[p1][0];x++)
        (*inputp)->config[x] = (*data)->ssize_pop[p1][x+1];
    for(x=(*data)->ssize_pop[p1][0];x<(*inputp)->npop;x++)
        (*inputp)->config[x] = 0;
            
	(*inputp)->ratio_sv = (*data)->ratio_sv[p1+1];
	if((*data)->ratio_sv[p1+1]>=REFNUMBER && (*data)->ratio_sv[p1+1] < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->ratio_sv[p1+1] - REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->ratio_sv);
	}
	
	/*theta vs S*/
	if((*data)->theta_1[0]) {
		if(!((*data)->theta_1[0]==1 && (*data)->theta_1[1] == (double)0)) {
			if((*data)->theta_1[0] == 1) {
				(*inputp)->theta = (*data)->theta_1[0+1];
				if((*data)->theta_1[0+1]>=REFNUMBER && (*data)->theta_1[0+1] < REFNUMBER + PLUSPRIOR) {
					z = (int)((*data)->theta_1[0+1] - REFNUMBER1);
					(*inputp)->pointtoprior[z] = &((*inputp)->theta);
				}
				(*inputp)->segsitesin = -1;
			}
			else {
				(*inputp)->theta = (*data)->theta_1[p1+1];
				if((*data)->theta_1[p1+1]>=REFNUMBER && (*data)->theta_1[p1+1] < REFNUMBER + PLUSPRIOR) {
					z = (int)((*data)->theta_1[p1+1] - REFNUMBER1);
					(*inputp)->pointtoprior[z] = &((*inputp)->theta);
				}
				(*inputp)->segsitesin = -1;
			}
		}
		else {
			(*inputp)->segsitesin = (*data)->mutations[p1+1];
			(*inputp)->theta = (double)0;
		}
	}
    else {
		(*inputp)->segsitesin = (*data)->mutations[p1+1];
		(*inputp)->theta = (double)0;
	}
	
	if((*inputp)->Sfix_alltheta) {
		(*inputp)->segsitesin = (*data)->mutations[p1+1];
		(*inputp)->theta = (double)0;
	}
	
	if((*inputp)->rmfix) {
		(*inputp)->Rm = (*data)->Rm[p1+1];
		(*inputp)->nhapl = (*data)->nhapl[p1+1];
		(*inputp)->r = (double)-1;
	}
	else {
		(*inputp)->Rm = -1;
		(*inputp)->nhapl = -1;
		if((*data)->ifgammar==0 && (*data)->range_rnt==0) {
			if((*data)->r[0] == 1) {
				(*inputp)->r = (*data)->r[0+1];
				if((*data)->r[0+1]>=REFNUMBER && (*data)->r[0+1] < REFNUMBER + PLUSPRIOR) {
					z = (int)((*data)->r[0+1] - REFNUMBER1);
					(*inputp)->pointtoprior[z] = &((*inputp)->r);
				}
				/*(*inputp)->Rm = -1;*/
			}
			else {
				(*inputp)->r = (*data)->r[p1+1];
				if((*data)->r[p1+1]>=REFNUMBER && (*data)->r[p1+1] < REFNUMBER + PLUSPRIOR) {
					z = (int)((*data)->r[p1+1] - REFNUMBER1);
					(*inputp)->pointtoprior[z] = &((*inputp)->r);
				}
			}
		}
		/*
		else 
			(*inputp)->Rm = (double)-1;
		*/
	}

    (*inputp)->npop_sampled = (*data)->npop_sampled[p1+1];    
    
    if((*data)->thetant_min[0] >= p1+1) x = p1+1;
	else x = 1;
	(*inputp)->thetant_min   = (*data)->thetant_min[x] * (double)(*data)->nsites[p1+1];
    (*inputp)->thetant_max   = (*data)->thetant_max[x] * (double)(*data)->nsites[p1+1];
	
    if((*data)->recnt_min[0] >= p1+1) x = p1+1;
	else x = 1;
    (*inputp)->recnt_min   = (*data)->recnt_min[x] * (double)(*data)->nsites[p1+1];
    (*inputp)->recnt_max   = (*data)->recnt_max[x] * (double)(*data)->nsites[p1+1];
	
	/*if((*inputp)->ifgamma) {*/
		(*inputp)->p_gamma = (*data)->p_gamma[p1+1];
		(*inputp)->alpha_gamma = (*data)->alpha_gamma[p1+1];
		(*inputp)->correct_gamma = (*data)->correct_gamma[p1+1];
	/*}*/
	/*if((*inputp)->ifgammar) {*/
		(*inputp)->p_gammar = (*data)->p_gammar[p1+1];
		(*inputp)->alpha_gammar = (*data)->alpha_gammar[p1+1];
		(*inputp)->correct_gammar = (*data)->correct_gammar[p1+1];
	/*}*/
	
	/*(*inputp)->factor_chrn = (double)1/(double)(*data)->factor_chrn[p1+1];*/
    if((*data)->factor_chrn[0] == 1) {
		(*inputp)->factor_chrn = (*data)->factor_chrn[1];
		/*
		if((*data)->factor_chrn[1] >=REFNUMBER && (*data)->factor_chrn[1] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->factor_chrn[1] - REFNUMBER1);
			(*inputp)->pointtoprior[z] = &((*inputp)->factor_chrn);
		}
		*/
	}
	else {
		(*inputp)->factor_chrn = (*data)->factor_chrn[p1+1];
		/*
		if((*data)->factor_chrn[p1+1] >=REFNUMBER && (*data)->factor_chrn[p1+1] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->factor_chrn[p1+1] - REFNUMBER1);
			(*inputp)->pointtoprior[z] = &((*inputp)->factor_chrn);
		}
		*/
	}
/*
	if((*data)->event_sexratio[0] == 2) {
		(*inputp)->event_forsexratio = (int)(*data)->event_sexratio[1];
		(*inputp)->event_sexratio = (*data)->event_sexratio[2];
		if((*data)->event_sexratio[2] >=REFNUMBER && (*data)->event_sexratio[2] < REFNUMBER + PLUSPRIOR) {
			z = (int)((*data)->event_sexratio[2] - REFNUMBER1);
			(*inputp)->pointtoprior[z] = &((*inputp)->event_sexratio);
		}
	}
	else {
		if((*data)->event_sexratio[0] < 2) {
			(*inputp)->event_forsexratio = -1;
			(*inputp)->event_sexratio = 1.0;
		}
		else {
			(*inputp)->event_forsexratio = (int)(*data)->event_sexratio[1];
			(*inputp)->event_sexratio = (*data)->event_sexratio[p1+2];
			if((*data)->event_sexratio[p1+2] >=REFNUMBER && (*data)->event_sexratio[p1+2] < REFNUMBER + PLUSPRIOR) {
				z = (int)((*data)->event_sexratio[p1+2] - REFNUMBER1);
				(*inputp)->pointtoprior[z] = &((*inputp)->event_sexratio);
			}
		}
	}	
*/	
    (*inputp)->heter_theta_alphag = (*data)->heter_theta_alphag[p1+1];
    (*inputp)->invariable_mut_sites = (*data)->invariable_mut_sites[p1+1];
    (*inputp)->heter_rm_alphag = (*data)->heter_rm_alphag[p1+1];

    (*inputp)->f = (*data)->f[p1+1];
	if((*data)->f[p1+1] >=REFNUMBER && (*data)->f[p1+1] < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->f[p1+1] - REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->f);
	}
    (*inputp)->track_len = (*data)->track_len[p1+1];
	if((*data)->track_len[p1+1]>=REFNUMBER && (*data)->track_len[p1+1] < REFNUMBER + PLUSPRIOR) {
		z = (int)((*data)->track_len[p1+1] - REFNUMBER1);
		(*inputp)->pointtoprior[z] = &((*inputp)->track_len);
	}
	
	if((*data)->ehh_fixnt[0] == 0)
		(*inputp)->ehh_fixnt = 0;
	else
		(*inputp)->ehh_fixnt = (*data)->ehh_fixnt[p1+3];

    return;
}
/*************************************** FREE **************************************/
void free_getpars_fix(struct var **data, struct var2 **inputp) 
{
    int x;
    
    free((*inputp)->config);
	free((*inputp)->ts);
	free((*inputp)->nintn);
	for(x=0;x<(*data)->npop;x++) {
		free((*inputp)->nrec[x]);
		free((*inputp)->tpast[x]);	
		free((*inputp)->tpastS[x]);
		free((*inputp)->migrate_matrix[x]);
	}
	free((*inputp)->nrec);
	free((*inputp)->tpast);	
	free((*inputp)->tpastS);
	free((*inputp)->migrate_matrix);
	
	for(x=0;x<(*inputp)->npop_events;x++) {
		free((*inputp)->pop_event[x].mig_fw);
		free((*inputp)->pop_event[x].mig_rv);
	}
	free((*inputp)->pop_event);
	free((*inputp)->fixmax_abstime_event);
    
    if((*data)->linked > 1) {
        for(x=0;x<(*data)->linked;x++) 
            free((*inputp)->loci_linked[x]);
    }
    if((*inputp)->linked > 1) {
        free((*inputp)->loci_linked);
	}
    free((*inputp)->factor_pop);
    
    if((*inputp)->split_pop)
        free((*inputp)->freq);
	
	if((*inputp)->npriors) free((*inputp)->pointtoprior);
	if((*inputp)->npriors) free((*inputp)->locspec_prior);
	
	if((*inputp)->type_ancestral) {
		for(x=0;x<(*inputp)->type_ancestral;x++) {
			free((*inputp)->ancestral_pol[x]);
		}
		free((*inputp)->ancestral_pol);
	}

    return;
}
