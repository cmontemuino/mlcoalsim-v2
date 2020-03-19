/************************************************************************************/
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

/*HUDSON ms.c file modified*/
/*ALSO INCLUDED A FUNCTION FOR RM FROM J. Wall.*/

#include "mlsp_sm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "print_helper.h"
#include "neut_tests.h"
#include "neutpar_common.h"
#include "neutpar_window.h"

#define FACTOR (double)1.5
/*please, NIT_PRIOR no more than 32000*/
#define SITESINC 100
#define PRINTTHETAS 0
#define DEBUGSOUT 0

long int maxsites = 1000;	/* la llargada de la regió */

double **coef = NULL;
struct dnapar *neutpar = NULL;    
long int *posit;	/*important! externes perque al reallocar memoria no localitza la nova posicio. */
                        /*Tambe es fa posant la posicio de memoria*/
char **list;

char dfiles[NEUTVALUES2+1][420] = {{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},
	{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},
	{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},
	{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},
	{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},
	{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"}
};
char ndfiles[NEUTVALUES2+1][20] = {{"_TD.out"},{"_Fs.out"},{"_FDn.out"},{"_FFn.out"},{"_FD.out"},{"_FF.out"},{"_H.out"},
	{"_B.out"},{"_Q.out"},{"_ZnA.out"},{"_Fst.out"},{"_Kw.out"},{"_Hw.out"},{"_R2.out"},
	{"_S.out"},{"_piw.out"},{"_pib.out"},{"_ThetaW.out"},{"_ThetaT.out"},
	{"_ThetaFW.out"},{"_D_Dmin.out"},{"_H_norm.out"},{"_maxhap.out"},{"_maxhap1.out"},{"_Rm.out"},
	{"_ThetaFL.out"},{"_ThetaL.out"},{"_ZengE.out"},{"_EW.out"},{"_Fstw.out"},{"_minuiHS.out"},{"_maxuiHS.out"},{"_maxliES.out"},{"_fixoutg.out"},{"_KoutgJC.out"},
	{"_ThetaWA.out"},{"_ThetaWAn.out"},{"_ThetaTA.out"},{"_ThetaTAn.out"},{"_AchY.out"},{"_AchYn.out"},
	{"_msdev.out"},{"_mskew.out"},{"_mkurt.out"},{"_ragg.out"},{"_ZnS.out"},{"_ZZ.out"},
	{"_PPercentiles.out"}
};
char names[NEUTVALUES2][20]  = {{"TD"},{"Fs"},{"FDn"},{"FFn"},{"FD"},{"FF"},{"H"},
	{"B"},{"Q"},{"ZnA"},{"Fst"},{"Kw"},{"Hw"},{"R2"},
	{"S"},{"piw"},{"pib"},{"ThetaW"},{"ThetaT"},
	{"ThetaFW"},{"D_Dmin"},{"H_norm"},{"maxhap"},{"maxhap1"},{"Rm"},
	{"ThetaFL"},{"ThetaL"},{"ZengE"},{"EW"},{"Fstw"},{"minuiHS"},{"maxuiHS"},{"maxliES"},
	{"fixoutg"},{"KoutgJC"},{"ThetaWA"},{"ThetaWAn"},{"ThetaTA"},{"ThetaTAn"},{"AchY"},{"AchYn"},
	{"msdev"},{"mskew"},{"mkurt"},{"ragg"},{"ZnS"},{"ZZ"}
};

int ms(struct var2 **inputp,char *file_out,double **matrix_test,struct prob_par **postp,long int count0,long int countmax,long int *listnumbers, long int *jcount2,long int *mcount2,struct var_priors *priors,int my_rank,long int seed1)
{
	long int segsites=0;
    int i;

    char **cmatrix(int,long int);
    long int gensam(long int,int,int *,long int,double,long int,double,double,double,double,int,
        long int,double *,double *,int,double,double,double,long int,double,int *,int *,double **,double **,double **,
        int, double, double, double, double *, double,int,double *,double,double,double *,double *,double **,int,int,struct events *,
		int,long int **,int,double,double,int,double,double,double);
    void mod_mhits(long int, struct var2 **,double *);
    void mod_outgroup(long int, struct var2 **,int);
    double logp(double);
    double ran1(void);
	void init_seed1(long int);

    void calc_neutpar(int,long int,struct var2 **,struct dnapar *,double,int);
    void calc_neutparSRH(int,long int,struct var2 **,struct dnapar *,double,int,int);
    double Zns(int,long int,struct var2 **,int);
    double Zns_window(struct var2 **,long int,long int,int);
    double ZnA_(int,long int,struct var2 **,int);
    double ZnA_window(struct var2 **,long int,long int,int);
    double Fstw(double *, int *,double,int);
    double logPPoisson2(long int, double);
    double koutg(int,int,int *);
	double koutgJC(int,int,int *,unsigned long);
	double fixoutg(int,int,int);

    double thetae=0.;
	double thetaemin=0.;
	double thetaemax= 0.;
	double lengtht=0.;
    int j,k,xw,neutv,sep=0;
    long int jcount=0;
	long int kcount=0;
	long int mcount=0;
	double u,logv;
	long int ui;
    double logPoissonkk=0.;
    double div=0.; double divt,divr;
    int burn_in=0;
    int compare_(const void *,const void *);
    int nwindow,aa;
    long int s0,s1,kk,ll;
	long int inputS;

	int Rmi;
	int nhi;
	double rece=0.;
	double recemin=0.;
	double recemax=0.;
	double recombinationv=0.;
    double correction_recabs(double,double,int);
    void calc_neutpar_windowSRH(struct var2 **,struct dnapar *,long int,long int,double,int,int);

	double gammadist(double);


	/*counting simulations in stdout*/
    static double counterp;
	static double restp;
    static long int p;
	static int npm = 0;
	int x;
	
	/*heterogeneity in mutation and recombination*/
	double *weightmut=0;
	double *weightrec=0;
	int do_heter_gamma_sites(double *,double,double,long int,long int);
	int do_heter_gamma_sitesrec(double *,double,double,long int);
	double thetaSv;
    	
	long int nrej; /*avoid too much rejections*/
	long int rejflag;
	int nv;
	
	/*count populations separately*/
	int npops,npopa;
	
	/******priors******/
	double maxtpast;
	int y;
			
	/**************** FILES   **************/
	FILE *output=0;
	FILE **file_outputm=0;
	
	FILE *outfstallcomp=0;
	FILE *outpiwallcomp=0;
	FILE *outpiaallcomp=0;

	FILE **outfstallcompm = 0;
	FILE **outpiwallcompm = 0;
	FILE **outpiaallcompm = 0;
	
	FILE *outfsthapallcomp=0;
	FILE *outhapwallcomp=0;
	FILE *outhapaallcomp=0;
	
	FILE **outfsthapallcompm = 0;
	FILE **outhapwallcompm = 0;
	FILE **outhapaallcompm = 0;
	
	FILE *outancestral=0;
	FILE **outancestralm = 0;
	
	FILE *file_outputpost=0;
	
	char namefile_[420];
	char file_out_[420];
	char *el;

	#if DEBUGSOUT
    FILE *file_debug;
    if (!(file_debug = fopen("debug.out","w"))) perror("Error in input/output");
	#endif
    #if PRINTTHETAS
    FILE *file_debug_sel;
    if (!(file_debug_sel = fopen("debug_sel.out","w"))) perror("Error in input/output");
    fputs("piTaj\tpiWatt\tpiFW\n",file_debug_sel);
    #endif

	/***** count for MPI *****/
	long int countmpi = count0;

	/******************************* PRINTING OUTPUT FILES ************************************/	
	file_out_[0] = '\0';		
	strncat(file_out_,file_out,420);
	el = strrchr(file_out_,'.');
	*el = '\0';
	
	if((*inputp)->likelihood_line == 0 && (*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
		if(((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt)) {
			x = (*inputp)->nloci;
			if(count0 == 0) {
				for(j=0;j<420;j++) namefile_[j] = '\0';
				sprintf(namefile_,"%s_locus_%05d_postp.out",file_out_,x);
				if(!(file_outputpost = fopen(namefile_,"w"))) perror("Error in input/output");
				fputs(MLCOALSIM,file_outputpost);
				fprintf(file_outputpost,"n_iter\t");
				if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
					fputs("theta_post\t",file_outputpost);
				if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
					fputs("rec_post\t",file_outputpost);
				if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt)
					fputs("lengtht\tprob_accept\n",file_outputpost);
			}
			else {
				for(j=0;j<420;j++) namefile_[j] = '\0';
				sprintf(namefile_,"%s_locus_%05d_postp.out",file_out_,x);
				if(!(file_outputpost = fopen(namefile_,"a+"))) perror("Error in input/output");
			}
		}
	}

	if((*inputp)->likelihood_line == 0 && (*inputp)->pr_matrix && (*inputp)->linked) {		
		namefile_[0] = '\0';		
		strncat(namefile_,file_out,420);
		el = strrchr(namefile_,'.');
		*el = '\0';
		strncat(namefile_,"_linkedlocus_\0",420);
		sprintf(namefile_,"%s_rank%03d.out",namefile_,my_rank);
		if(count0 == countmpi) {
			if(!(output = fopen(namefile_,"w"))) perror("Error in input/output");
			if(count0 == 0) fputs(MLCOALSIM,output);
		}
		else {
			if(!(output = fopen(namefile_,"a+"))) perror("Error in input/output");
		}
	}
	if(((*inputp)->likelihood_line == 0 &&  (*inputp)->pr_matrix          && (*inputp)->linked == 0) ||
	   ((*inputp)->likelihood_line == 0 && ((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5 && (*inputp)->linked == 0))) {
		x = (*inputp)->nloci;
		if(count0 == 0) {
			for(j=0;j<420;j++) namefile_[j] = '\0';	
			sprintf(namefile_,"%s_locus_%05d.out",file_out_,x);
			if(!(output = fopen(namefile_,"w"))) perror("Error in input/output");
			fputs(MLCOALSIM,output);
		}
		else {
			for(j=0;j<420;j++) namefile_[j] = '\0';	
			sprintf(namefile_,"%s_locus_%05d.out",file_out_,x);
			if(!(output = fopen(namefile_,"a+"))) perror("Error in input/output");
		}
	}
	if(((*inputp)->likelihood_line == 0 &&  (*inputp)->pr_matrix          && (*inputp)->linked == 0) ||
	   ((*inputp)->likelihood_line == 0 && ((*inputp)->print_neuttest >/*=*/ 2 && (*inputp)->print_neuttest < 5 && (*inputp)->linked == 0))) {
		if(((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) && (*inputp)->print_neuttest >/*=*/ 2 && (*inputp)->print_neuttest < 5) {
			x = (*inputp)->nloci;
			if(count0 == 0) {
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_Fst_allcomp.out",file_out_,x);
				if(!(outfstallcomp = fopen(namefile_,"w"))) perror("Error in input/output");
				fputs(MLCOALSIM,outfstallcomp);
				fflush(outfstallcomp);
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_PiW_allcomp.out",file_out_,x);
				if(!(outpiwallcomp = fopen(namefile_,"w"))) perror("Error in input/output");
				fputs(MLCOALSIM,outpiwallcomp);
				fflush(outpiwallcomp);
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_PiA_allcomp.out",file_out_,x);
				if(!(outpiaallcomp = fopen(namefile_,"w"))) perror("Error in input/output");
				fputs(MLCOALSIM,outpiaallcomp);
				fflush(outpiaallcomp);
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_Fsthap_allcomp.out",file_out_,x);
				if(!(outfsthapallcomp = fopen(namefile_,"w"))) perror("Error in input/output");
				fputs(MLCOALSIM,outfsthapallcomp);
				fflush(outfsthapallcomp);
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_HapW_allcomp.out",file_out_,x);
				if(!(outhapwallcomp = fopen(namefile_,"w"))) perror("Error in input/output");
				fputs(MLCOALSIM,outhapwallcomp);
				fflush(outhapwallcomp);
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_HapA_allcomp.out",file_out_,x);
				if(!(outhapaallcomp = fopen(namefile_,"w"))) perror("Error in input/output");
				fputs(MLCOALSIM,outhapaallcomp);
				fflush(outhapaallcomp);
				
				if((*inputp)->type_ancestral) {
					for(j=0;j<420;j++) namefile_[j] = '\0';	
					sprintf(namefile_,"%s_locus_%05d_Sancestral.out",file_out_,x);
					if(!(outancestral = fopen(namefile_,"w"))) perror("Error in input/output");
					fputs(MLCOALSIM,outancestral);
					fflush(outancestral);
				}
			}
			else {
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_Fst_allcomp.out",file_out_,x);
				if(!(outfstallcomp = fopen(namefile_,"a+"))) perror("Error in input/output");
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_PiW_allcomp.out",file_out_,x);
				if(!(outpiwallcomp = fopen(namefile_,"a+"))) perror("Error in input/output");
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_PiA_allcomp.out",file_out_,x);
				if(!(outpiaallcomp = fopen(namefile_,"a+"))) perror("Error in input/output");
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_Fsthap_allcomp.out",file_out_,x);
				if(!(outfsthapallcomp = fopen(namefile_,"a+"))) perror("Error in input/output");

				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_HapW_allcomp.out",file_out_,x);
				if(!(outhapwallcomp = fopen(namefile_,"a+"))) perror("Error in input/output");
				
				for(j=0;j<420;j++) namefile_[j] = '\0';	
				sprintf(namefile_,"%s_locus_%05d_HapA_allcomp.out",file_out_,x);
				if(!(outhapaallcomp = fopen(namefile_,"a+"))) perror("Error in input/output");
				
				if((*inputp)->type_ancestral) {
					for(j=0;j<420;j++) namefile_[j] = '\0';	
					sprintf(namefile_,"%s_locus_%05d_Sancestral.out",file_out_,x);
					if(!(outancestral = fopen(namefile_,"a+"))) perror("Error in input/output");
				}
			}
		}

	}
	if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0 && (*inputp)->linked == 0) {
		if(count0 == 0) {
			fprintf(output,"n_iter\t");
		#if SHOW_PRIORS
			for(j=1;j<=(*inputp)->npriors;j++) fprintf(output,"prior[%03d]\t",j);
		#endif
			for(j=0;j<(*inputp)->max_npop_sampled;j++) {
		#if FREQSPECTRUM_TDFSR2 == 0
				for(x=0;x<NEUTVALUES2;x++) fprintf(output,"%s[pop=%d]\t",names[x],j);
		#endif
		#if FREQSPECTRUM_TDFSR2 == 1
				fprintf(output,"%s[pop=%d]\t",names[0],j);
				fprintf(output,"%s[pop=%d]\t",names[1],j);
				fprintf(output,"%s[pop=%d]\t",names[7],j);
				fprintf(output,"%s[pop=%d]\t",names[8],j);
				fprintf(output,"%s[pop=%d]\t",names[9],j);
				fprintf(output,"%s[pop=%d]\t",names[13],j);
				fprintf(output,"%s[pop=%d]\t",names[14],j);
				fprintf(output,"%s[pop=%d]\t",names[17],j);
				fprintf(output,"%s[pop=%d]\t",names[18],j);
				fprintf(output,"%s[pop=%d]\t",names[21],j);
				fprintf(output,"%s[pop=%d]\t",names[33],j);
				fprintf(output,"%s[pop=%d]\t",names[34],j);
				fprintf(output,"%s[pop=%d]\t",names[41],j);
				fprintf(output,"%s[pop=%d]\t",names[42],j);
				fprintf(output,"%s[pop=%d]\t",names[43],j);
		#endif
				fprintf(output,"nmhits[pop=%d]\t",j);
				for(x=1;x<(*inputp)->config[j];x++) fprintf(output,"fr[pop=%d,%d]\t",j,x);
			}
			fputs("\n",output); 
			fflush(output);
		}
	}
	if((*inputp)->print_neuttest >/*=*/ 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0 && (*inputp)->linked == 0) {
		if(((*inputp)->max_npop_sampled > 2 && (*inputp)->print_neuttest < 5) || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
			if(count0 == 0) {
				fprintf(outfstallcomp,"n_iter\t");
				for(j=0;j<(*inputp)->max_npop_sampled;j++) for(x=j+1;x<(*inputp)->max_npop_sampled;x++) fprintf(outfstallcomp,"Fst[pop=%d][pop=%d]\t",j,x);
				fputs("\n",outfstallcomp); 
				fflush(outfstallcomp);
				
				fprintf(outpiwallcomp,"n_iter\t");
				for(j=0;j<(*inputp)->max_npop_sampled;j++) fprintf(outpiwallcomp,"PiW[pop=%d]\t",j);
				fputs("\n",outpiwallcomp); 
				fflush(outpiwallcomp);
				
				fprintf(outpiaallcomp,"n_iter\t");
				for(j=0;j<(*inputp)->max_npop_sampled;j++) for(x=j+1;x<(*inputp)->max_npop_sampled;x++) fprintf(outpiaallcomp,"PiA[pop=%d][pop=%d]\t",j,x);
				fputs("\n",outpiaallcomp); 
				fflush(outpiaallcomp);
				
				fprintf(outfsthapallcomp,"n_iter\t");
				for(j=0;j<(*inputp)->max_npop_sampled;j++) for(x=j+1;x<(*inputp)->max_npop_sampled;x++) fprintf(outfsthapallcomp,"FstH[pop=%d][pop=%d]\t",j,x);
				fputs("\n",outfsthapallcomp); 
				fflush(outfsthapallcomp);
				
				fprintf(outhapwallcomp,"n_iter\t");
				for(j=0;j<(*inputp)->max_npop_sampled;j++) fprintf(outhapwallcomp,"HapW[pop=%d]\t",j);
				fputs("\n",outhapwallcomp); 
				fflush(outhapwallcomp);
				
				fprintf(outhapaallcomp,"n_iter\t");
				for(j=0;j<(*inputp)->max_npop_sampled;j++) for(x=j+1;x<(*inputp)->max_npop_sampled;x++) fprintf(outhapaallcomp,"HapA[pop=%d][pop=%d]\t",j,x);
				fputs("\n",outhapaallcomp); 
				fflush(outhapaallcomp);
				
				if((*inputp)->type_ancestral) {
					fprintf(outancestral,"n_iter\t");
					if((*inputp)->type_ancestral > 3) {
						for(j=0;j<(*inputp)->type_ancestral-1;j++) {
							if(j>0) fprintf(outancestral,"\t");
							fprintf(outancestral,"Sx[%d]\tSf[%d]\tSxf[%d,rest]\tSs[%d,rest]",j,j,j,j);
						}
						fprintf(outancestral,"\tSx[out]\tSf[out]\tSxf[out,rest]\tSs[out,rest]");
					}
					if((*inputp)->type_ancestral == 3) fprintf(outancestral,"Sx1\tSx2\tSxo\tSf1\tSf2\tSfo\tSx1f2\tSx2f1\tSsh\tSso");
					if((*inputp)->type_ancestral == 2) fprintf(outancestral,"Sx1\tSx2\tSf\tSsh");
					if((*inputp)->type_ancestral == 1) fprintf(outancestral,"Sx1");
					fputs("\n",outancestral); 
					fflush(outancestral);
				}
			}
		}
	}
	/*PRINT LINKED REGIONS*/
	if((*inputp)->print_neuttest >/*=*/ 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0 && (*inputp)->linked) {
		if((*inputp)->linked == 1) {
			if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5)
				if(!(file_outputm = (FILE **)calloc(1,sizeof(FILE *)))) perror("Error in input/output");
			if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
				if(!(outfstallcompm = (FILE **)calloc(1,sizeof(FILE *)))) perror("Error in input/output");
				if(!(outpiwallcompm = (FILE **)calloc(1,sizeof(FILE *)))) perror("Error in input/output");
				if(!(outpiaallcompm = (FILE **)calloc(1,sizeof(FILE *)))) perror("Error in input/output");
				if(!(outancestralm = (FILE **)calloc(1,sizeof(FILE *)))) perror("Error in input/output");
			}
			x = 0;
			for(kk=0,ll=(*inputp)->window;kk<(*inputp)->nsites;kk += (long int)(*inputp)->despl) {
				if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
					if(!(file_outputm = (FILE **)realloc(file_outputm,(x+1)*sizeof(FILE *)))) perror("Error in input/output");
					for(j=0;j<420;j++) namefile_[j] = '\0';
					sprintf(namefile_,"%s_locus_%05d_rank%03d.out",file_out_,x,my_rank);
					if(count0 == countmpi) {
						if(!(file_outputm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,file_outputm[x]);
							fprintf(file_outputm[x],"n_iter\t");
					#if SHOW_PRIORS
							for(j=1;j<=(*inputp)->npriors;j++) fprintf(file_outputm[x],"prior[%03d]\t",j);
					#endif
							for(j=0;j<(*inputp)->max_npop_sampled;j++) {
					#if FREQSPECTRUM_TDFSR2 == 0
								for(k=0;k<NEUTVALUES2;k++) fprintf(file_outputm[x],"%s[pop=%d]\t",names[k],j);
					#endif
					#if FREQSPECTRUM_TDFSR2 == 1
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[0],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[1],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[7],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[8],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[9],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[13],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[14],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[17],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[18],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[21],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[33],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[34],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[41],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[42],j);
								fprintf(file_outputm[x],"%s[pop=%d]\t",names[43],j);
					#endif
								fprintf(file_outputm[x],"nmhits[pop=%d]\t",j);
								for(xw=1;xw<(*inputp)->config[j];xw++) fprintf(file_outputm[x],"fr[pop=%d,%d]\t",j,xw);
							}
							fputs("\n",file_outputm[x]); 
							fflush(file_outputm[x]);
						}
					}
				}
				if(count0 == countmpi) {
					if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
						if(!(outfstallcompm = (FILE **)realloc(outfstallcompm,(x+1)*sizeof(FILE *)))) perror("Error in input/output");
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_Fst_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outfstallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outfstallcompm[x]);
							fprintf(outfstallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) for(xw=j+1;xw<(*inputp)->max_npop_sampled;xw++) fprintf(outfstallcompm[x],"Fst[pop=%d][pop=%d]\t",j,xw);
							fputs("\n",outfstallcompm[x]); 
							fflush(outfstallcompm[x]);
						}
						
						if(!(outpiwallcompm = (FILE **)realloc(outpiwallcompm,(x+1)*sizeof(FILE *)))) perror("Error in input/output");
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_PiW_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outpiwallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outpiwallcompm[x]);
							fprintf(outpiwallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) fprintf(outpiwallcompm[x],"PiW[pop=%d]\t",j);
							fputs("\n",outpiwallcompm[x]); 
							fflush(outpiwallcompm[x]);
						}
						
						if(!(outpiaallcompm = (FILE **)realloc(outpiaallcompm,(x+1)*sizeof(FILE *)))) perror("Error in input/output");
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_PiA_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outpiaallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outpiaallcompm[x]);
							fprintf(outpiaallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) for(xw=j+1;xw<(*inputp)->max_npop_sampled;xw++) fprintf(outpiaallcompm[x],"PiA[pop=%d][pop=%d]\t",j,xw);
							fputs("\n",outpiaallcompm[x]); 
							fflush(outpiaallcompm[x]);
						}
						
						if(!(outfsthapallcompm = (FILE **)realloc(outfsthapallcompm,(x+1)*sizeof(FILE *)))) perror("Error in input/output");
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_Fsthap_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outfsthapallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outfsthapallcompm[x]);
							fprintf(outfsthapallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) for(xw=j+1;xw<(*inputp)->max_npop_sampled;xw++) fprintf(outfsthapallcompm[x],"FstH[pop=%d][pop=%d]\t",j,xw);
							fputs("\n",outfsthapallcompm[x]); 
							fflush(outfsthapallcompm[x]);
						}
						
						if(!(outhapwallcompm = (FILE **)realloc(outhapwallcompm,(x+1)*sizeof(FILE *)))) perror("Error in input/output");
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_HapW_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outhapwallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outhapwallcompm[x]);
							fprintf(outhapwallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) fprintf(outhapwallcompm[x],"HapW[pop=%d]\t",j);
							fputs("\n",outhapwallcompm[x]); 
							fflush(outhapwallcompm[x]);
						}
						
						if(!(outhapaallcompm = (FILE **)realloc(outhapaallcompm,(x+1)*sizeof(FILE *)))) perror("Error in input/output");
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_HapA_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outhapaallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outhapaallcompm[x]);
							fprintf(outhapaallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) for(xw=j+1;xw<(*inputp)->max_npop_sampled;xw++) fprintf(outhapaallcompm[x],"HapA[pop=%d][pop=%d]\t",j,xw);
							fputs("\n",outhapaallcompm[x]); 
							fflush(outhapaallcompm[x]);
						}
						
						if((*inputp)->type_ancestral) {
							if(!(outancestralm = (FILE **)realloc(outancestralm,(x+1)*sizeof(FILE *)))) perror("Error in input/output");
							for(j=0;j<420;j++) namefile_[j] = '\0';	
							sprintf(namefile_,"%s_locus_%05d_Sancestral_rank%03d.out",file_out_,x,my_rank);
							if(!(outancestralm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
							if(count0 == 0) {
								fputs(MLCOALSIM,outancestralm[x]);
								fprintf(outancestralm[x],"n_iter\t");
								if((*inputp)->type_ancestral > 3) {
									for(j=0;j<(*inputp)->type_ancestral-1;j++) {
										if(j>0) fprintf(outancestralm[x],"\t");
										fprintf(outancestralm[x],"Sx[%d]\tSf[%d]\tSxf[%d,rest]\tSs[%d,rest]",j,j,j,j);
									}
									fprintf(outancestralm[x],"\tSx[out]\tSf[out]\tSxf[out,rest]\tSs[out,rest]");
								}
								if((*inputp)->type_ancestral == 3) fprintf(outancestralm[x],"Sx1\tSx2\tSxo\tSf1\tSf2\tSfo\tSx1f2\tSx2f1\tSsh\tSso");
								if((*inputp)->type_ancestral == 2) fprintf(outancestralm[x],"Sx1\tSx2\tSf\tSsh");
								if((*inputp)->type_ancestral == 1) fprintf(outancestralm[x],"Sx1");
								fputs("\n",outancestralm[x]); 
								fflush(outancestralm[x]);
							}
						}
					}
				}
				else {
					if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) 
						if(!(file_outputm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
					
					if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_Fst_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outfstallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_PiW_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outpiwallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_PiA_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outpiaallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_FstHap_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outfsthapallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_HapW_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outhapwallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_HapA_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outhapaallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						if((*inputp)->type_ancestral) {
							for(j=0;j<420;j++) namefile_[j] = '\0';	
							sprintf(namefile_,"%s_locus_%05d_Sancestral_rank%03d.out",file_out_,x,my_rank);
							if(!(outancestralm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						}
					}
				}
				if(ll == (*inputp)->nsites) break;
				else ll += (*inputp)->despl;
				if(ll > (*inputp)->nsites) ll = (*inputp)->nsites;
				x++;
			}
		}
		else {
			if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5)
				if(!(file_outputm = (FILE **)calloc((*inputp)->linked,sizeof(FILE *)))) perror("Error in input/output");
			if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
				if(!(outfstallcompm = (FILE **)calloc((*inputp)->linked,sizeof(FILE *)))) perror("Error in input/output");
				if(!(outpiwallcompm = (FILE **)calloc((*inputp)->linked,sizeof(FILE *)))) perror("Error in input/output");
				if(!(outpiaallcompm = (FILE **)calloc((*inputp)->linked,sizeof(FILE *)))) perror("Error in input/output");
				if(!(outancestralm = (FILE **)calloc((*inputp)->linked,sizeof(FILE *)))) perror("Error in input/output");
			}
			for(x=0;x<(*inputp)->linked;x++) {
				if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
					for(j=0;j<420;j++) namefile_[j] = '\0';
					sprintf(namefile_,"%s_locus_%05d_rank%03d.out",file_out_,x,my_rank);
					if(count0 == countmpi) {
						if(!(file_outputm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,file_outputm[x]);
							fprintf(file_outputm[x],"n_iter\t");
						#if SHOW_PRIORS
							for(j=1;j<=(*inputp)->npriors;j++) fprintf(file_outputm[x],"prior[%03d]\t",j);
						#endif
							for(j=0;j<(*inputp)->max_npop_sampled;j++) {
								for(k=0;k<NEUTVALUES2;k++) fprintf(file_outputm[x],"%s[pop=%d]\t",names[k],j);
								fprintf(file_outputm[x],"nmhits[pop=%d]\t",j);
								for(xw=1;xw<(*inputp)->config[j];xw++) fprintf(file_outputm[x],"fr[pop=%d,%d]\t",j,xw);
							}
							fputs("\n",file_outputm[x]); 
							fflush(file_outputm[x]);
						}
					}
				}
				if(count0 == countmpi) {
					if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_Fst_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outfstallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outfstallcompm[x]);
							fprintf(outfstallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) for(xw=j+1;xw<(*inputp)->max_npop_sampled;xw++) fprintf(outfstallcompm[x],"Fst[pop=%d][pop=%d]\t",j,xw);
							fputs("\n",outfstallcompm[x]); 
							fflush(outfstallcompm[x]);
						}
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_PiW_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outpiwallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outpiwallcompm[x]);
							fprintf(outpiwallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) fprintf(outpiwallcompm[x],"PiW[pop=%d]\t",j);
							fputs("\n",outpiwallcompm[x]);
							fflush(outpiwallcompm[x]);
						}
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_PiA_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outpiaallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outpiaallcompm[x]);
							fprintf(outpiaallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) for(xw=j+1;xw<(*inputp)->max_npop_sampled;xw++) fprintf(outpiaallcompm[x],"PiA[pop=%d][pop=%d]\t",j,xw);
							fputs("\n",outpiaallcompm[x]);
							fflush(outpiaallcompm[x]);
						}

						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_FstHAp_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outfsthapallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outfsthapallcompm[x]);
							fprintf(outfsthapallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) for(xw=j+1;xw<(*inputp)->max_npop_sampled;xw++) fprintf(outfsthapallcompm[x],"FstH[pop=%d][pop=%d]\t",j,xw);
							fputs("\n",outfsthapallcompm[x]); 
							fflush(outfsthapallcompm[x]);
						}
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_HapW_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outhapwallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outhapwallcompm[x]);
							fprintf(outhapwallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) fprintf(outhapwallcompm[x],"HapW[pop=%d]\t",j);
							fputs("\n",outhapwallcompm[x]);
							fflush(outhapwallcompm[x]);
						}
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_HApA_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outhapaallcompm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
						if(count0 == 0) {
							fputs(MLCOALSIM,outhapaallcompm[x]);
							fprintf(outhapaallcompm[x],"n_iter\t");
							for(j=0;j<(*inputp)->max_npop_sampled;j++) for(xw=j+1;xw<(*inputp)->max_npop_sampled;xw++) fprintf(outhapaallcompm[x],"HApA[pop=%d][pop=%d]\t",j,xw);
							fputs("\n",outhapaallcompm[x]);
							fflush(outhapaallcompm[x]);
						}
						
						if((*inputp)->type_ancestral) {
							for(j=0;j<420;j++) namefile_[j] = '\0';	
							sprintf(namefile_,"%s_locus_%05d_Sancestral_rank%03d.out",file_out_,x,my_rank);
							if(!(outancestralm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
							if(count0 == 0) {
								fputs(MLCOALSIM,outancestralm[x]);
								fprintf(outancestralm[x],"n_iter\t");
								if((*inputp)->type_ancestral > 3) {
									for(j=0;j<(*inputp)->type_ancestral-1;j++) {
										if(j>0) fprintf(outancestralm[x],"\t");
										fprintf(outancestralm[x],"Sx[%d]\tSf[%d]\tSxf[%d,rest]\tSs[%d,rest]",j,j,j,j);
									}
									fprintf(outancestralm[x],"\tSx[out]\tSf[out]\tSxf[out,rest]\tSs[out,rest]");
								}
								if((*inputp)->type_ancestral == 3) fprintf(outancestralm[x],"Sx1\tSx2\tSxo\tSf1\tSf2\tSfo\tSx1f2\tSx2f1\tSsh\tSso");
								if((*inputp)->type_ancestral == 2) fprintf(outancestralm[x],"Sx1\tSx2\tSf\tSsh");
								if((*inputp)->type_ancestral == 1) fprintf(outancestralm[x],"Sx1");
								fputs("\n",outancestralm[x]); 
								fflush(outancestralm[x]); 
							}
						}
					}
				}
				else {
					if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5)
						if(!(file_outputm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
					
					if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_Fst_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outfstallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_PiW_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outpiwallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_PiA_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outpiaallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_FstHap_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outfsthapallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_HapW_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outhapwallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						for(j=0;j<420;j++) namefile_[j] = '\0';	
						sprintf(namefile_,"%s_locus_%05d_HapA_allcomp_rank%03d.out",file_out_,x,my_rank);
						if(!(outhapaallcompm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						
						if((*inputp)->type_ancestral) {
							for(j=0;j<420;j++) namefile_[j] = '\0';	
							sprintf(namefile_,"%s_locus_%05d_Sancestral_rank%03d.out",file_out_,x,my_rank);
							if(!(outancestralm[x] = fopen(namefile_,"a+"))) perror("Error in input/output");
						}
					}
				}
			}
		}
	}
	
	/*********************** heterogeneity in mutation and recombination *****************/
	/*init*/
	if(!(weightmut = (double *)calloc((*inputp)->nsites+1,sizeof(double)))) {
		perror("calloc error ms.89");
		exit(1);
	}
	if(!(weightrec = (double *)calloc((*inputp)->nsites,sizeof(double)))) {
		perror("calloc error ms.89");
		exit(1);
	}

	/*counting simulations in stdout*/
	if((*inputp)->Sfix_alltheta == 2 || (*inputp)->rmfix == 2) {
		if(npm == 0) {
			counterp = (double)(((*inputp)->howmany+(*inputp)->mc_jump) * (*inputp)->tloci)/50.;
			p = 1;
			restp = 0.;
			npm = 1;
		}
	}
	else {
		if(npm == 0) {
			counterp = (double)((*inputp)->howmany * (*inputp)->tloci)/50.;
			p = 1;
			restp = 0.;
			npm = 1;
		}
	}
    /*two important parameters for Sfixallthetas = 2 !! ELIMINATED*/
    if((*inputp)->Sfix_alltheta == 2 || (*inputp)->rmfix == 2) {
        sep = (*inputp)->mc_jump; 	/*separation of the values we pick in the markov chain: 1 fastests...*/
        if((*inputp)->Sfix_alltheta == 2 || (*inputp)->rmfix == 2)
            div = 0.01;
        else div = 0.01;
        /*Fixed to 1 works well. Range of the chosen value in the uniform distribution: (0,1]*/
    }
	/*Define the matrix for posterior probabilities for theta and other values*/
    
    if((*inputp)->theta > 0.0 || (*inputp)->ifgamma == 1 || (*inputp)->range_thetant) {        
        list = cmatrix((*inputp)->nsam,maxsites+1);
        if(!(posit = (long int *)malloc((long int)(maxsites*sizeof(long int)))))
            perror("ms error.1");
    }
    else {
        list = cmatrix((*inputp)->nsam,(long int)(*inputp)->segsitesin+1);
        if(!(posit = (long int *)malloc((long int)(((long int)(*inputp)->segsitesin+1)*sizeof(long int)))))
            perror("ms error.2");
    }
    /* matriu double test taj,fs,fd,ff,h,B,Q,ZnS,Fst,#hap,divhap... per tots loci */
    if(coef == NULL) {
        if((neutpar = (struct dnapar *)calloc(1,sizeof(struct dnapar))) == NULL) {
            perror("calloc error ms.1d1");
            exit(1);
        }
        if((coef = (double **)calloc(1,sizeof(double *))) == NULL) {
            perror("calloc error ms.1d1");
            exit(1);
        }
        for(i=0;i<1;i++) {
            if((coef[i] = (double *)calloc(18,sizeof(double))) == NULL) {
                perror("calloc error ms.1d1");
                exit(1);
            }
        }
    }

    if((neutpar[0].freq = (int *)calloc((*inputp)->nsam,sizeof(int))) == NULL) {
        perror("calloc error ms.1d3");
        exit(1);
    }
    if((neutpar[0].fhapl = (int *)calloc((*inputp)->nsam,sizeof(int))) == NULL) {
        perror("calloc error ms.1d3");
        exit(1);
    }
    if((neutpar[0].unic = (long int *)calloc((*inputp)->nsam,sizeof(long int))) == NULL) {
        perror("calloc error ms.1d3");
        exit(1);
    }
    if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
        if((neutpar[0].fstallcomp = (double *)calloc((*inputp)->max_npop_sampled,sizeof(double))) == NULL) {
            perror("calloc error ms.1d30");
            exit(1);
        }
        if((neutpar[0].piwallcomp = (double *)calloc((*inputp)->max_npop_sampled,sizeof(double))) == NULL) {
            perror("calloc error ms.1d31");
            exit(1);
        }
        if((neutpar[0].piaallcomp = (double *)calloc((*inputp)->max_npop_sampled,sizeof(double))) == NULL) {
            perror("calloc error ms.1d32");
            exit(1);
        }
        if((neutpar[0].fsthapallcomp = (double *)calloc((*inputp)->max_npop_sampled,sizeof(double))) == NULL) {
            perror("calloc error ms.1d30");
            exit(1);
        }
        if((neutpar[0].hapwallcomp = (double *)calloc((*inputp)->max_npop_sampled,sizeof(double))) == NULL) {
            perror("calloc error ms.1d31");
            exit(1);
        }
        if((neutpar[0].hapaallcomp = (double *)calloc((*inputp)->max_npop_sampled,sizeof(double))) == NULL) {
            perror("calloc error ms.1d32");
            exit(1);
        }
    }
	if((*inputp)->type_ancestral > 3) {
		if((neutpar[0].Sanc = (int *)calloc(4*(*inputp)->type_ancestral,sizeof(int))) == NULL) {
			perror("calloc error ms.1d321");
			exit(1);
		}
	}
	else {
		if((neutpar[0].Sanc = (int *)calloc(10,sizeof(int))) == NULL) {
			perror("calloc error ms.1d321");
			exit(1);
		}
	}
	
	if((*inputp)->range_thetant) {
		thetaemax = (*inputp)->thetant_max;
		thetaemin = (*inputp)->thetant_min;
	}
	if((*inputp)->range_rnt) {
		recemax = (*inputp)->recnt_max;
		recemin = (*inputp)->recnt_min;
	}
	
	/*other parameters*/        
    if(((*inputp)->Sfix_alltheta || (*inputp)->rmfix) && (*inputp)->linked < 2) { /*Calculate approx. min and max for theta. Do NIT_PRIOR iterations and estimate theta...*/
		if(!((*inputp)->Sfix_alltheta && (*inputp)->rmfix == 0)) {
			thetaemax = (*inputp)->thetant_max;
			thetaemin = (*inputp)->thetant_min;
		}
		if((*inputp)->rmfix && (*inputp)->linked < 2) {
			recemax = (*inputp)->recnt_max;
			recemin = (*inputp)->recnt_min;
		}
		
		if((*inputp)->Sfix_alltheta==2 || (*inputp)->rmfix==2) {
			divt = div*(thetaemax-thetaemin); /*range of the chosen value in the uniform distribution for theta not used?*/
			divr = div*(recemax-recemin); /*range of the chosen value in the uniform distribution for rec*/
			if((*inputp)->range_thetant == 2) divt = div*((double)log(thetaemax)-(double)log(thetaemin)); /*range of the chosen value in the uniform distribution for theta*/
			if((*inputp)->range_rnt == 2) divr = div*((double)log(recemax)-(double)log(recemin)); /*range of the chosen value in the uniform distribution for rec*/
			burn_in = 1000;
		}
        if((*inputp)->Sfix_alltheta) {
			if((*inputp)->segsitesin > 0) logPoissonkk = (double)logPPoisson2((long int)(*inputp)->segsitesin,(double)(*inputp)->segsitesin);
			else logPoissonkk = (double)logPPoisson2((long int)(*inputp)->segsitesin,(double)0.1);
		}
	}
	
	if((*inputp)->Sfix_alltheta || (*inputp)->rmfix) {
        j = 0; 
        k = 0; 
        jcount = *jcount2; 
        mcount = *mcount2; 
        kcount = 0; 
		neutv = NEUTVALUES2;
    }
	if((*inputp)->Sfix_alltheta == 0 && (*inputp)->rmfix == 0  && (*inputp)->linked < 2) {
		if((*inputp)->range_thetant) {
			thetaemax = (*inputp)->thetant_max;
			thetaemin = (*inputp)->thetant_min;
		}
		if((*inputp)->range_rnt) {
			recemax = (*inputp)->recnt_max;
			recemin = (*inputp)->recnt_min;
		}
	}
	/*assign to npopa the pop to do rejection or mcmc*/
	npopa = 0;
	if((*inputp)->Sfix_alltheta) npopa = (*inputp)->Sfix_pop;
	else if((*inputp)->rmfix) npopa = (*inputp)->rmfix_pop;

	rejflag = 0;
	
	/*change the seed for each process in case linked*/
	if((*inputp)->linked) init_seed1((long int)-(2*count0+seed1));
	
	/*vectors for heterogeneity: the SAME for ALL iterations (but each locus is different). In each iteration they are just weighted by the theta or rec (now over 1)*/
	do_heter_gamma_sites   (weightmut,(*inputp)->heter_theta_alphag,1.0,(*inputp)->nsites+1,(*inputp)->invariable_mut_sites);
	do_heter_gamma_sitesrec(weightrec,(*inputp)->heter_rm_alphag,   1.0,(*inputp)->nsites);
    
	/************************************************  ITERATIONS: Routine from Hudson and modified ****************************************************/
    while(countmax - count0++) {    
		/*take into account the prior distributions*/
		for(x=0;x<(*inputp)->npriors;x++) {
			if((*inputp)->pointtoprior[x]) 
				*((*inputp)->pointtoprior[x]) = priors[x].priordist[count0-1];
		}
		
		maxtpast = 0.;
		for(x=0;x<(*inputp)->npop/*-((*inputp)->ifselection)*/;x++) { /*NOT ALLOWED MORE THAN ONE POPULATION NOW*/
			for(y=1;y<=(*inputp)->nintn[x];y++) {
				(*inputp)->tpast[x][y] = (*inputp)->tpastS[x][y] + (*inputp)->tpastS[x][y-1]; 
				if(maxtpast < (*inputp)->tpast[x][y]) maxtpast = (*inputp)->tpast[x][y];
			}
		}
		if((*inputp)->split_pop) {
			/*time_split is the first event in split_pop*/
			(*inputp)->time_scoal = (*inputp)->time_scoalS + (*inputp)->time_split;
		}
		else {
			for(x=0;x<(*inputp)->npop_events;x++) {/*events*/
				if((*inputp)->fixmax_abstime_event[x] > (*inputp)->pop_event[x].timeS_event + maxtpast)
					(*inputp)->pop_event[x].time_event = (*inputp)->pop_event[x].timeS_event + maxtpast; 
				else
					(*inputp)->pop_event[x].time_event = (*inputp)->fixmax_abstime_event[x];
				maxtpast = (*inputp)->pop_event[x].time_event;
			}
		}		
		/* ranfactor: between 0.1 and 1, or between 1 and 10. Equally divided */
        if((*inputp)->ran_factorpop == 2) {
            for(i=1;i<(*inputp)->npop;i++) {
                (*inputp)->factor_pop[i] = (double)ran1()* (double)9 + (double)1;
                if((double)ran1() < 0.5) (*inputp)->factor_pop[i] = 1./(*inputp)->factor_pop[i];
            }
        }
        if((*inputp)->Sfix_alltheta == 0 && (*inputp)->rmfix == 0) {/*"normal" simulations: Fix theta (or fix S, for the wole populations together) and fix R*/
			if((*inputp)->ifgamma == 1) thetae = gammadist((*inputp)->alpha_gamma)/(*inputp)->p_gamma * (*inputp)->correct_gamma;
			else {
				if((*inputp)->range_thetant == 1) thetae = thetaemin + (thetaemax - thetaemin) * ran1();
				else {
					if((*inputp)->range_thetant == 2) thetae = (double)exp((double)log((double)thetaemin) + ((double)log((double)thetaemax) - (double)log((double)thetaemin)) * ran1());
					else thetae = (double)(*inputp)->theta * (double)(*inputp)->nsites;
				}
			}
			if((*inputp)->ifgammar == 1) rece = gammadist((*inputp)->alpha_gammar)/(*inputp)->p_gammar * (*inputp)->correct_gammar;
			else {
				if((*inputp)->range_rnt == 1) rece = recemin + (recemax - recemin) * ran1();
				else {
					if((*inputp)->range_rnt == 2) rece = (double)exp((double)log((double)recemin) + ((double)log((double)recemax) - (double)log((double)recemin)) * ran1());
					else rece = (double)(*inputp)->r * (double)(*inputp)->nsites;
				}
			}			
            recombinationv = rece * correction_recabs((*inputp)->factor_chrn,(*inputp)->sex_ratio,(*inputp)->no_rec_males);/*recombinationv is only useful for statistics calculations*/
			if((*inputp)->segsitesin == -1 || (*inputp)->range_thetant || (*inputp)->ifgamma) {
				thetaSv = thetae /** correction_theta((*inputp)->factor_chrn,(*inputp)->sex_ratio)*/;
				inputS = -1;
			}
			else thetaSv = (double)1.0; /*indicated 1.0 (and not 0) because it has to be taken into account the heterogeneity in relation to thetaSv*/

			segsites = gensam((*inputp)->npop,(*inputp)->nsam, (*inputp)->config, (*inputp)->nsites,
            thetaSv, (*inputp)->segsitesin, rece, (*inputp)->f, (*inputp)->track_len,
            (*inputp)->migrate,(*inputp)->mhits,count0,(*inputp)->factor_pop, 
            &lengtht,(*inputp)->ifselection, 
            (*inputp)->pop_sel,(*inputp)->sinit,(*inputp)->pop_size,(long int)(*inputp)->sel_nt,(*inputp)->T_out,&(*inputp)->Sout,
            (*inputp)->nintn, (*inputp)->nrec, (*inputp)->nrec, (*inputp)->tpast,
            (*inputp)->split_pop,(*inputp)->time_split,(*inputp)->time_scoal,(*inputp)->factor_anc,(*inputp)->freq,
            (*inputp)->tlimit,(*inputp)->iflogistic,(*inputp)->ts,(*inputp)->factor_chrn,(*inputp)->ratio_sv,
			weightmut,weightrec,(*inputp)->migrate_matrix,my_rank,(*inputp)->npop_events,(*inputp)->pop_event,
			(*inputp)->linked,(*inputp)->loci_linked,(*inputp)->event_forsexratio,(*inputp)->event_sexratio,
			(*inputp)->sex_ratio,(*inputp)->no_rec_males,(*inputp)->sendt,(*inputp)->sfreqend,(*inputp)->sfreqinit);
            if((*inputp)->mhits) mod_mhits(segsites,inputp,weightmut); /******** mhits ******************/
            if((*inputp)->pop_outgroup != -1) mod_outgroup(segsites,inputp,(*inputp)->pop_outgroup); /******** outgroup ******************/
			if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
				postp[0][/*listnumbers[*/count0-1/*]*/].thetap = thetae;
			if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
				postp[0][/*listnumbers[*/count0-1/*]*/].recp = rece;
			if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
				postp[0][/*listnumbers[*/count0-1/*]*/].Ttotp  = lengtht;
				postp[0][/*listnumbers[*/count0-1/*]*/].prob  = 1.0;
			}

			/*printing a point every 2% of the total iterations*/
			#if SHOWPROGRESS == 1
			show_progress(counterp, &restp, &p);
			#endif
        }
        /*else do the rejection algorithm (RA) of Tavaré et al. 1997.*/
        else {
			if((*inputp)->Sfix_alltheta == 1 || (*inputp)->rmfix == 1) { 
				nrej = jcount;
				do {
					rejflag = 0;/*flag for jumping when rejections be > MAXREJECT*/					
					if((*inputp)->ifgamma == 1) thetae = gammadist((*inputp)->alpha_gamma)/(*inputp)->p_gamma * (*inputp)->correct_gamma;
					else {
						if((*inputp)->range_thetant == 1) thetae = thetaemin + (thetaemax - thetaemin) * ran1();
						else {
							if((*inputp)->range_thetant == 2) thetae = (double)exp((double)log((double)thetaemin) + ((double)log((double)thetaemax) - (double)log((double)thetaemin)) * ran1());
							else thetae = (double)(*inputp)->theta * (double)(*inputp)->nsites;
						}
					}
					if((*inputp)->ifgammar == 1) rece = gammadist((*inputp)->alpha_gammar)/(*inputp)->p_gammar * (*inputp)->correct_gammar;
					else {
						if((*inputp)->range_rnt == 1) rece = recemin + (recemax - recemin) * ran1();
						else {
							if((*inputp)->range_rnt == 2) rece = (double)exp((double)log((double)recemin) + ((double)log((double)recemax) - (double)log((double)recemin)) * ran1());
							else rece = (double)(*inputp)->r * (double)(*inputp)->nsites;
						}
					}			
					recombinationv = rece * correction_recabs((*inputp)->factor_chrn,(*inputp)->sex_ratio,(*inputp)->no_rec_males); /*recombinationv is only useful for statistics calculations*/
					
					/*be careful with mhits*/
					inputS = (*inputp)->segsitesin;
					if((*inputp)->mhits == 1) inputS = -1;
					if((*inputp)->segsitesin == -1 || (*inputp)->mhits == 1 || ((*inputp)->npop > 1 && (*inputp)->Sfix_pop < (*inputp)->npop)) {
						thetaSv = thetae /** correction_theta((*inputp)->factor_chrn,(*inputp)->sex_ratio)*/;
						inputS = -1;
					}
					else {
						thetaSv = (double)1.0;
					}
					segsites = gensam((*inputp)->npop,(*inputp)->nsam, (*inputp)->config, (*inputp)->nsites,
					thetaSv, /*(*inputp)->segsitesin*/inputS, rece, (*inputp)->f, (*inputp)->track_len,
					(*inputp)->migrate,(*inputp)->mhits,count0,(*inputp)->factor_pop, 
					&lengtht,(*inputp)->ifselection, 
					(*inputp)->pop_sel,(*inputp)->sinit,(*inputp)->pop_size,(long int)(*inputp)->sel_nt,(*inputp)->T_out,&(*inputp)->Sout,
					(*inputp)->nintn, (*inputp)->nrec, (*inputp)->nrec, (*inputp)->tpast,
					(*inputp)->split_pop,(*inputp)->time_split,(*inputp)->time_scoal,(*inputp)->factor_anc,(*inputp)->freq,
					(*inputp)->tlimit,(*inputp)->iflogistic,(*inputp)->ts,(*inputp)->factor_chrn,(*inputp)->ratio_sv,
					weightmut,weightrec,(*inputp)->migrate_matrix,my_rank,(*inputp)->npop_events,(*inputp)->pop_event,
					(*inputp)->linked,(*inputp)->loci_linked,(*inputp)->event_forsexratio,(*inputp)->event_sexratio,
					(*inputp)->sex_ratio,(*inputp)->no_rec_males,(*inputp)->sendt,(*inputp)->sfreqend,(*inputp)->sfreqinit);
					if((*inputp)->mhits) mod_mhits(segsites,inputp,weightmut); /******** mhits ******************/
					if((*inputp)->pop_outgroup != -1) mod_outgroup(segsites,inputp,(*inputp)->pop_outgroup); /******** outgroup ******************/
					
					if((*inputp)->linked > 1) {/*linked fragments with subset fixed values*/
						s0=s1=0;
						nwindow=0;
						for(aa=0;aa<(*inputp)->linked;aa++) {
							kk=(*inputp)->loci_linked[aa][0];
							ll=(*inputp)->loci_linked[aa][1]+1;
							while(s0 < segsites && posit[s0] < kk) s0++;
							while(s1 < segsites && posit[s1] < ll) s1++;
							/*calc_estadistics*/
							if(((*inputp)->rmfix == 1 && (*inputp)->linked_rm_nregion == aa) || ((*inputp)->Sfix_alltheta == 1 && (*inputp)->linked_segsites_nregion == aa)) {
								calc_neutpar_windowSRH(inputp,neutpar,s0,s1,(double)recombinationv,npopa,(*inputp)->linked_nhapl);
								if((*inputp)->rmfix == 1) {
									if((*inputp)->linked_rm_nregion == aa) {
										Rmi = neutpar[0].Rm;
										nhi = neutpar[0].nhapl;
										if(Rmi != (*inputp)->linked_rm || ((*inputp)->linked_nhapl != 0 && nhi != (*inputp)->linked_nhapl)) {
											jcount += 1;
											if(jcount - nrej >= MAXREJECT) {
												rejflag = 1;
												aa = (*inputp)->linked;
											}
											break; /*reject*/
										}
									}
								}
								if((*inputp)->Sfix_alltheta == 1 && inputS != -1) {
									u = (double)logPPoisson2((long int)(*inputp)->segsitesin,thetae*lengtht) - logPoissonkk;
									logv = (double)log((double)ran1());
									if(u < logv) {
										jcount++;
										if(jcount - nrej >= MAXREJECT) {
											rejflag = 1;
											aa = (*inputp)->linked;
										}
										break; /*reject*/
									}
								}
								else {
									if((*inputp)->Sfix_alltheta == 1) {
										if((*inputp)->linked_segsites_nregion == aa) {
											ui = neutpar[0].S;
											if(ui != (long int)(*inputp)->linked_segsites) {
												jcount += 1;
												if(jcount - nrej >= MAXREJECT) {
													rejflag = 1;
													aa = (*inputp)->linked;
												}
												break; /*reject*/
											}
										}
									}
								}
							}
						}
						if(aa<(*inputp)->linked && rejflag == 0) continue; /*reject*/
						else {
							if(rejflag == 0) {
								/*accept*/ /*We will be in this loop until we have "howmany" success*/
								if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
									postp[0][count0-1].thetap = thetae;
								if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
									postp[0][count0-1].recp = rece;
								if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
									postp[0][count0-1].Ttotp  = lengtht;
									postp[0][count0-1].prob  = 1.0/(double)(1.0 + jcount - nrej);
								}
								mcount++;
								jcount++;
							}
							else {
								if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
									postp[0][count0-1].thetap = -10000.;
								if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
									postp[0][count0-1].recp = -10000.;
								if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
									postp[0][count0-1].Ttotp  = -10000.;
									postp[0][count0-1].prob  = 0.0;
								}
							}
							/*printing a point every 2% of the total iterations*/
							#if SHOWPROGRESS == 1
                            show_progress(counterp, &restp, &p);
							#endif

							break;
						}
					}
					else {/*like not linked fragments*/
						if((*inputp)->rmfix == 1 || (*inputp)->Sfix_alltheta == 1)
							calc_neutparSRH(0,segsites,inputp,neutpar+0,(double)recombinationv,npopa,(*inputp)->nhapl);
						if((*inputp)->rmfix == 1) {
							Rmi = neutpar[0].Rm;
							nhi = neutpar[0].nhapl;
							if(Rmi != (*inputp)->Rm || ((*inputp)->nhapl != 0 && nhi != (*inputp)->nhapl)) {
								jcount += 1;
								if(jcount - nrej >= MAXREJECT) {
									rejflag = 1;
								}
								if(rejflag==0) continue; /*reject*/
							}
						}
						if((*inputp)->Sfix_alltheta == 1 && inputS != -1 && rejflag == 0) {
							u = (double)logPPoisson2((long int)(*inputp)->segsitesin,thetae*lengtht) - logPoissonkk;
							logv = (double)log((double)ran1());
							if(u < logv) {
								jcount++;
								if(jcount - nrej >= MAXREJECT) {
									rejflag = 1;
									aa = (*inputp)->linked;
								}
								if(rejflag==0) continue; /*reject*/
							}
						}
						else {
							if((*inputp)->Sfix_alltheta == 1 && rejflag == 0/*&& (*inputp)->mhits == 1*/) {
								ui = neutpar[0].S;
								if(ui != (long int)(*inputp)->segsitesin) {
									jcount += 1;
									if(jcount - nrej >= MAXREJECT) {
										rejflag = 1;
									}
									if(rejflag==0) continue; /*reject*/
								}
							}
						}
						if(rejflag == 0) {
							/*accept*/ /*We will be in this loop until we have "howmany" success*/
							if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
								postp[0][count0-1].thetap = thetae;
							if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
								postp[0][count0-1].recp = rece;
							if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
								postp[0][count0-1].Ttotp  = lengtht;
								postp[0][count0-1].prob  = 1.0/(double)(jcount - nrej);
							}
							mcount++;
							jcount++;
						}
						else {
							if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
								postp[0][count0-1].thetap = -10000.;
							if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
								postp[0][count0-1].recp = -10000.;
							if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
								postp[0][count0-1].Ttotp  = -10000.;
								postp[0][count0-1].prob  = 1.0/(double)(jcount - nrej);
							}
						}
						/*printing a point every 2% of the total iterations*/
						#if SHOWPROGRESS == 1
                        show_progress(counterp, &restp, &p);
						#endif

						break;
					}
				}while(1);
			}    
		}
		/*MHMCMC routine eliminated*/

		if((*inputp)->pr_matrix) /******** print the whole matrix **********/
			if(rejflag==0) print_matrix(segsites,inputp,output,(*inputp)->pr_matrix,count0-1,priors,(*inputp)->mhits, list, posit);
		/******************************* STATISTICS *******************************/

		if((*inputp)->print_neuttest) {
			for(npops=0;npops<(*inputp)->max_npop_sampled;npops++) {
				if((*inputp)->linked > 0) {
					/*define the limits for each region (min and max values)*/
					if((*inputp)->linked == 1) {
						/*sliding windows*/
						s0=s1=0;
						nwindow=0;
						for(kk=0,ll=(*inputp)->window;kk<(*inputp)->nsites;kk += (long int)(*inputp)->despl) {
							while(s0 < segsites && posit[s0]<kk) s0++;
							while(s1 < segsites && posit[s1]<ll) s1++;
							/*calc_estadistics*/
							if(rejflag==0 && npops < (*inputp)->npop_sampled) {
								init_coef(coef[0],(*inputp)->config[npops]);
								calc_neutpar_window(inputp,neutpar,s0,s1,(double)recombinationv,npops, list, posit);
								/*calc_neut_tests*/
								/*we need to define the number of windows and the current window number*/
							#if FREQSPECTRUM_TDFSR2 < 2		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]] 
									= (double)tajima_d(neutpar[0].k,(int)neutpar[0].S,coef[0]);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]] 
									= (double)Fs((*inputp)->config[npops],neutpar[0].k,neutpar[0].nhapl);
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+2][listnumbers[count0-1]]
									= (double)fl_d2_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+3][listnumbers[count0-1]] 
									= (double)fl_f2_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+4][listnumbers[count0-1]] 
									= (double)fl_d_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+5][listnumbers[count0-1]]
									= (double)fl_f_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+6][listnumbers[count0-1]] 
									= (double)fay_wu((*inputp)->config[npops],neutpar[0].freq,neutpar[0].k);
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2
								if(neutpar[0].S > 1) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]] 
										= (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]] 
										= (double)neutpar[0].Q1/((double)neutpar[0].S);
								}
								else {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]] = -10000;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]] = -10000;
								}
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								if(neutpar[0].S > 1) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]] 
									= (double)ZnA_window(inputp,s0,s1,npops);
								}
								else {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]] = -10000;
								}
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+10][listnumbers[count0-1]] 
										= (double)Fst(neutpar[0].piw,neutpar[0].pib,(int)(*inputp)->npop);
								
								if(neutpar[0].nhapl > -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+11][listnumbers[count0-1]] 
									= (double)neutpar[0].nhapl;
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+11][listnumbers[count0-1]] = -10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+12][listnumbers[count0-1]] 
									= (double)testHap((*inputp)->config[npops],neutpar[0].fhapl);
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								if(neutpar[0].S > 0)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]] 
										= (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp)->config[npops],(int)neutpar[0].S);
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]] = -10000;
							#endif
							#if FREQSPECTRUM_TDFSR2 < 3		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]] 
									= (double)neutpar[0].S;
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+15][listnumbers[count0-1]] 
									= (double)neutpar[0].k;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+16][listnumbers[count0-1]] 
									= (double)neutpar[0].pib;
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								/*S, pi, thetaFW*/
								if(neutpar[0].S > -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]] 
									= (double)neutpar[0].S/(double)coef[0][0];
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]] = -10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]] 
									= (double)neutpar[0].k;
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								if(neutpar[0].k > -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+19][listnumbers[count0-1]] 
									= (double)neutpar[0].k - 
									  (double)(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+6][listnumbers[count0-1]]);
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+19][listnumbers[count0-1]] = -10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+20][listnumbers[count0-1]] 
									= (double)tajima_dvsdmin(neutpar[0].k,(int)neutpar[0].S,coef[0],(*inputp)->config[npops]);
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]] 
									= (double)fay_wu_normalized2((int)(*inputp)->config[npops],(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0],neutpar[0].k);
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								if(neutpar[0].k > -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+22][listnumbers[count0-1]] 
									= (double)neutpar[0].maxhapl/(double)(*inputp)->config[npops];
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+22][listnumbers[count0-1]] = -10000;
								if(neutpar[0].maxhapl1 != -10000) 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+23][listnumbers[count0-1]] 
									= (double)neutpar[0].maxhapl1/(double)(*inputp)->config[npops];
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+23][listnumbers[count0-1]] 
									= (double)neutpar[0].maxhapl1;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+24][listnumbers[count0-1]] 
									= (double)neutpar[0].Rm;
								if(neutpar[0].S != -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+25][listnumbers[count0-1]] 
									= (double)neutpar[0].freq[1];
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+25][listnumbers[count0-1]] = -10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+26][listnumbers[count0-1]] 
									= (double)neutpar[0].thetaL;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+27][listnumbers[count0-1]] 
									= (double)E_zeng((int)(*inputp)->config[npops],(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0]);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+28][listnumbers[count0-1]] 
									= (double)EWtest((int)(*inputp)->config[npops],neutpar[0].fhapl);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+29][listnumbers[count0-1]] 
									= (double)Fst(neutpar[0].withinw,neutpar[0].pib,(int)(*inputp)->npop);
								if(neutpar[0].S > 1) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+30][listnumbers[count0-1]] 
									= neutpar[0].min_uiHS;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+31][listnumbers[count0-1]] 
									= neutpar[0].max_uiHS;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+32][listnumbers[count0-1]]
									= neutpar[0].max_iES;
								}
								else {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+30][listnumbers[count0-1]] = (double)-10000;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+31][listnumbers[count0-1]] = (double)-10000;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+32][listnumbers[count0-1]] = (double)-10000;
								}
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								if((*inputp)->T_out > 0. || (*inputp)->pop_outgroup != -1) {
									if((*inputp)->pop_outgroup != -1) (*inputp)->Sout = 0;
									else neutpar[0].freq[0] = 0;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]] 
									= (double)fixoutg((*inputp)->config[npops],(*inputp)->Sout,neutpar[0].freq[0]);
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]] 
									= (double)koutgJC((*inputp)->config[npops],(*inputp)->Sout,neutpar[0].freq,(*inputp)->nsites);
								}
								else {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]] = -10000;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]] = -10000;
								}
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								if((*inputp)->config[npops] > 2)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+35][listnumbers[count0-1]]
									= (double)neutpar[0].Se1/(coef[0][0]-1.);
                                else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+35][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 3)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+36][listnumbers[count0-1]] 
									= (double)neutpar[0].Sn1/(coef[0][0]-((double)(*inputp)->config[npops]/((double)(*inputp)->config[npops]-1.0)));
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+36][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 2)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+37][listnumbers[count0-1]] 
									= (double)neutpar[0].pie1 * (double)(*inputp)->config[npops]/((double)(*inputp)->config[npops]-2.0);
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+37][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 3)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+38][listnumbers[count0-1]] 
									= (double)neutpar[0].pin1 * ((double)(*inputp)->config[npops] - 1.0)/((double)(*inputp)->config[npops] - 3.0);
                                else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+38][listnumbers[count0-1]] = -10000;
                                
								if((*inputp)->config[npops] > 2)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+39][listnumbers[count0-1]] 
									= (double)Y_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+39][listnumbers[count0-1]] = -10000;
								if((*inputp)->config[npops] > 2)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+40][listnumbers[count0-1]] 
									= (double)(double)Y2_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+40][listnumbers[count0-1]] = -10000;
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
                                if((*inputp)->config[npops] > 2) {
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]]
                                    = (double)neutpar[0].m_sdev;
                                }
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 2) {
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]]
                                    = (double)neutpar[0].m_skew;
                                }
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 3) {
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]]
                                    = (double)neutpar[0].m_kurt;
                                }
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]] = -10000;
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0
								if(neutpar[0].S > 0) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+44][listnumbers[count0-1]] 
									= (double)neutpar[0].ragg;
								}
								#if ZNS_ACTIVE == 1
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]] 
								= (double)Zns_window(inputp,s0,s1,npops);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+46][listnumbers[count0-1]] 
								= matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9] [listnumbers[count0-1]] -
								  matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]];
								#else
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]] 
								= (double)-10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+46][listnumbers[count0-1]] 
								= (double)-10000;
								#endif
							#endif
							}
							else {
								for(nv=0;nv<NEUTVALUES2;nv++) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]] = (double)-10000;
								}
							}
							/*printing values*/
							if((*inputp)->print_neuttest >/*= */2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0) {
								if(npops==0) {
									if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) fprintf(file_outputm[nwindow],"%ld\t",count0-1);
								#if SHOW_PRIORS
									for(nv=0;nv<(*inputp)->npriors;nv++) {
										fprintf(file_outputm[nwindow],"%f\t",priors[nv].priordist[count0-1]);
									}
								#endif
									if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
										fprintf(outfstallcompm[nwindow],"%ld\t",count0-1);
										fprintf(outpiwallcompm[nwindow],"%ld\t",count0-1);
										fprintf(outpiaallcompm[nwindow],"%ld\t",count0-1);
										fflush(outfstallcompm[nwindow]);
										fflush(outpiwallcompm[nwindow]);
										fflush(outpiaallcompm[nwindow]);
										
										fprintf(outfsthapallcompm[nwindow],"%ld\t",count0-1);
										fprintf(outhapwallcompm[nwindow],"%ld\t",count0-1);
										fprintf(outhapaallcompm[nwindow],"%ld\t",count0-1);
										fflush(outfsthapallcompm[nwindow]);
										fflush(outhapwallcompm[nwindow]);
										fflush(outhapaallcompm[nwindow]);
										
										if((*inputp)->type_ancestral) {
											fprintf(outancestralm[nwindow],"%ld\t",count0-1);
											if((*inputp)->type_ancestral > 3) {
												for(nv=0;nv<4*(*inputp)->type_ancestral;nv++) {
													if(neutpar[0].Sanc[nv] == (double)-10000) fprintf(outancestralm[nwindow],"na\t");
													else fprintf(outancestralm[nwindow],"%d\t",neutpar[0].Sanc[nv]);
												}
											}
											if((*inputp)->type_ancestral == 3) {
												for(nv=0;nv<10;nv++) {
													if(neutpar[0].Sanc[nv] == (double)-10000) fprintf(outancestralm[nwindow],"na\t");
													else fprintf(outancestralm[nwindow],"%d\t",neutpar[0].Sanc[nv]);
												}
											}
											if((*inputp)->type_ancestral == 2) {
												if(neutpar[0].Sanc[0] == (double)-10000) fprintf(outancestralm[nwindow],"na\t");
												else fprintf(outancestralm[nwindow],"%d\t",neutpar[0].Sanc[0]);
												if(neutpar[0].Sanc[1] == (double)-10000) fprintf(outancestralm[nwindow],"na\t");
												else fprintf(outancestralm[nwindow],"%d\t",neutpar[0].Sanc[1]);
												if(neutpar[0].Sanc[3] == (double)-10000) fprintf(outancestralm[nwindow],"na\t");
												else fprintf(outancestralm[nwindow],"%d\t",neutpar[0].Sanc[3]);
												if(neutpar[0].Sanc[8] == (double)-10000) fprintf(outancestralm[nwindow],"na\t");
												else fprintf(outancestralm[nwindow],"%d\t",neutpar[0].Sanc[8]);
											}
											if((*inputp)->type_ancestral == 1) {
												if(neutpar[0].Sanc[0] == (double)-10000) fprintf(outancestralm[nwindow],"na\t");
												else fprintf(outancestralm[nwindow],"%d\t",neutpar[0].Sanc[0]);
											}
											if((*inputp)->type_ancestral) fputs("\n",outancestralm[nwindow]);
											fflush(outancestralm[nwindow]);
										}
									}
								}
								if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
							#if FREQSPECTRUM_TDFSR2 == 0		
									for(nv=0;nv<NEUTVALUES2;nv++) {
										if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]] == (double)-10000) 
											fprintf(file_outputm[nwindow],"na\t");
										else 
											fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]]);
									}
							#endif
							#if FREQSPECTRUM_TDFSR2 == 1	
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]]);
									if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[nwindow],"na\t");
									else 
										fprintf(file_outputm[nwindow],"%g\t",matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]]);
								#endif
									fprintf(file_outputm[nwindow],"%d\t",neutpar[0].mhsites);
									for(nv=1;nv<(*inputp)->config[npops];nv++) {
										fprintf(file_outputm[nwindow],"%d\t",neutpar[0].freq[nv]);
									}
									fflush(file_outputm[nwindow]);
								}
								if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
									if(neutpar[0].piwallcomp[npops] == (double)-10000) fprintf(outpiwallcompm[nwindow],"na\t");
									else fprintf(outpiwallcompm[nwindow],"%g\t",neutpar[0].piwallcomp[npops]);
									
									for(nv=npops+1;nv<(*inputp)->max_npop_sampled;nv++) {
										if(neutpar[0].piaallcomp[nv] == (double)-10000) fprintf(outpiaallcompm[nwindow],"na\t");
										else fprintf(outpiaallcompm[nwindow],"%g\t",neutpar[0].piaallcomp[nv]);
										
										if(neutpar[0].fstallcomp[nv] == (double)-10000) fprintf(outfstallcompm[nwindow],"na\t");
										else fprintf(outfstallcompm[nwindow],"%g\t",neutpar[0].fstallcomp[nv]);
									}
									fflush(outfstallcompm[nwindow]);
									fflush(outpiwallcompm[nwindow]);
									fflush(outpiaallcompm[nwindow]);

									if(neutpar[0].hapwallcomp[npops] == (double)-10000) fprintf(outhapwallcompm[nwindow],"na\t");
									else fprintf(outhapwallcompm[nwindow],"%g\t",neutpar[0].hapwallcomp[npops]);
									
									for(nv=npops+1;nv<(*inputp)->max_npop_sampled;nv++) {
										if(neutpar[0].hapaallcomp[nv] == (double)-10000) fprintf(outhapaallcompm[nwindow],"na\t");
										else fprintf(outhapaallcompm[nwindow],"%g\t",neutpar[0].hapaallcomp[nv]);
										
										if(neutpar[0].fsthapallcomp[nv] == (double)-10000) fprintf(outfsthapallcompm[nwindow],"na\t");
										else fprintf(outfsthapallcompm[nwindow],"%g\t",neutpar[0].fsthapallcomp[nv]);
									}
									fflush(outfsthapallcompm[nwindow]);
									fflush(outhapwallcompm[nwindow]);
									fflush(outhapaallcompm[nwindow]);
								}
								
								if(npops == (*inputp)->max_npop_sampled-1) {
									if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
										fputs("\n",file_outputm[nwindow]);
										fflush(file_outputm[nwindow]);
									}
									if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
										fputs("\n",outfstallcompm[nwindow]);
										fputs("\n",outpiwallcompm[nwindow]);
										fputs("\n",outpiaallcompm[nwindow]);

										fflush(outfstallcompm[nwindow]);
										fflush(outpiwallcompm[nwindow]);
										fflush(outpiaallcompm[nwindow]);										

										fputs("\n",outfsthapallcompm[nwindow]);
										fputs("\n",outhapwallcompm[nwindow]);
										fputs("\n",outhapaallcompm[nwindow]);
										
										fflush(outfstallcompm[nwindow]);
										fflush(outpiwallcompm[nwindow]);
										fflush(outpiaallcompm[nwindow]);
									}
								}
							}

							/*end calc neutrality test*/
							nwindow += 1;
							if(ll == (*inputp)->nsites) break;
							else ll += (*inputp)->despl;
							if(ll > (*inputp)->nsites) ll = (*inputp)->nsites;
							s0 = s1;
						}
					}
					else {/*linked regions*/
						s0=s1=0;
						nwindow=0;
						for(aa=0;aa<(*inputp)->linked;aa++) {
							kk=(*inputp)->loci_linked[aa][0];
							ll=(*inputp)->loci_linked[aa][1]+1;
							
							while(s0 < segsites && posit[s0] < kk) {
								s0++;
							}
							
							while(s1 < segsites && posit[s1] < ll) {
								s1++;
							}
							
							if(rejflag==0 && npops < (*inputp)->npop_sampled) {
								init_coef(coef[0],(*inputp)->config[npops]);
								/*calc_estadistics*/
								calc_neutpar_window(inputp,neutpar,s0,s1,(double)recombinationv,npops, list, posit);
								/*calc_neut_tests*/
							#if FREQSPECTRUM_TDFSR2 < 2		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]] 
									= (double)tajima_d(neutpar[0].k,(int)neutpar[0].S,coef[0]);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]] 
									= (double)Fs((*inputp)->config[npops],neutpar[0].k,neutpar[0].nhapl);
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+2][listnumbers[count0-1]] 
									= (double)fl_d2_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+3][listnumbers[count0-1]] 
									= (double)fl_f2_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+4][listnumbers[count0-1]] 
									= (double)fl_d_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+5][listnumbers[count0-1]]
									= (double)fl_f_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+6][listnumbers[count0-1]] 
									= (double)fay_wu((*inputp)->config[npops],neutpar[0].freq,neutpar[0].k);
							#endif							
							#if FREQSPECTRUM_TDFSR2 < 2
								if(neutpar[0].S > 1) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]] 
										= (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]] 
										= (double)neutpar[0].Q1/((double)neutpar[0].S);
								}
								else {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]] = -10000;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]] = -10000;
								}
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								if(neutpar[0].S > 1) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]] 
									= (double)ZnA_window(inputp,s0,s1,npops);
								}
								else {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]] = -10000;
								}
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+10][listnumbers[count0-1]] 
										= (double)Fst(neutpar[0].piw,neutpar[0].pib,(int)(*inputp)->npop);
								
								if(neutpar[0].nhapl > -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+11][listnumbers[count0-1]] 
									= (double)neutpar[0].nhapl;
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+11][listnumbers[count0-1]] = -10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+12][listnumbers[count0-1]] 
									= (double)testHap((*inputp)->config[npops],neutpar[0].fhapl);
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								if(neutpar[0].S > 0)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]] 
										= (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp)->config[npops],(int)neutpar[0].S);
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]] = -10000;
							#endif		
							#if FREQSPECTRUM_TDFSR2 < 3		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]] 
									= (double)neutpar[0].S;
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+15][listnumbers[count0-1]] 
									= (double)neutpar[0].k;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+16][listnumbers[count0-1]] 
									= (double)neutpar[0].pib;
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								/*S, pi, thetaFW*/
								if(neutpar[0].S > -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]] 
									= (double)neutpar[0].S/(double)coef[0][0];
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]] = -10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]] 
									= (double)neutpar[0].k;
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								if(neutpar[0].k > -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+19][listnumbers[count0-1]] 
									= (double)neutpar[0].k - 
									  (double)(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+6][listnumbers[count0-1]]);
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+19][listnumbers[count0-1]] = -10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+20][listnumbers[count0-1]] 
									= (double)tajima_dvsdmin(neutpar[0].k,(int)neutpar[0].S,coef[0],(*inputp)->config[npops]);
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]] 
									= (double)fay_wu_normalized2((int)(*inputp)->config[npops],(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0],neutpar[0].k);
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								if(neutpar[0].k > -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+22][listnumbers[count0-1]] 
									= (double)neutpar[0].maxhapl/(double)(*inputp)->config[npops];
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+22][listnumbers[count0-1]] = -10000;
								if(neutpar[0].maxhapl1 != -10000) 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+23][listnumbers[count0-1]] 
									= (double)neutpar[0].maxhapl1/(double)(*inputp)->config[npops];
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+23][listnumbers[count0-1]] 
									= (double)neutpar[0].maxhapl1;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+24][listnumbers[count0-1]] 
									= (double)neutpar[0].Rm;
								if(neutpar[0].S != -10000)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+25][listnumbers[count0-1]] 
									= (double)neutpar[0].freq[1];
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+25][listnumbers[count0-1]] = -10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+26][listnumbers[count0-1]] 
								= (double)neutpar[0].thetaL;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+27][listnumbers[count0-1]] 
								= (double)E_zeng((int)(*inputp)->config[npops],(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0]);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+28][listnumbers[count0-1]] 
								= (double)EWtest((int)(*inputp)->config[npops],neutpar[0].fhapl);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+29][listnumbers[count0-1]] 
								= (double)Fst(neutpar[0].withinw,neutpar[0].pib,(int)(*inputp)->npop);
								if(neutpar[0].S > 1) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+30][listnumbers[count0-1]] 
									= neutpar[0].min_uiHS;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+31][listnumbers[count0-1]] 
									= neutpar[0].max_uiHS;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+32][listnumbers[count0-1]]
									= neutpar[0].max_iES;
								}
								else {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+30][listnumbers[count0-1]] = (double)-10000;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+31][listnumbers[count0-1]] = (double)-10000;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+32][listnumbers[count0-1]] = (double)-10000;
								}
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
								if((*inputp)->T_out > 0. || (*inputp)->pop_outgroup != -1) {
									if((*inputp)->pop_outgroup != -1) (*inputp)->Sout = 0;
									else neutpar[0].freq[0] = 0;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]] 
									= (double)fixoutg((*inputp)->config[npops],(*inputp)->Sout,neutpar[0].freq[0]);
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]] 
									= (double)koutgJC((*inputp)->config[npops],(*inputp)->Sout,neutpar[0].freq,(*inputp)->nsites);
								}
								else {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]] = -10000;
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]] = -10000;
								}
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0		
								if((*inputp)->config[npops] > 2)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+35][listnumbers[count0-1]] 
									= (double)neutpar[0].Se1/(coef[0][0]-1.);
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+35][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 3)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+36][listnumbers[count0-1]] 
									= (double)neutpar[0].Sn1/(coef[0][0]-((double)(*inputp)->config[npops]/((double)(*inputp)->config[npops]-1.0)));
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+36][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 2)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+37][listnumbers[count0-1]] 
									= (double)neutpar[0].pie1 * (double)(*inputp)->config[npops]/((double)(*inputp)->config[npops]-2.0);
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+37][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 3)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+38][listnumbers[count0-1]] 
									= (double)neutpar[0].pin1 * ((double)(*inputp)->config[npops] - 1.0)/((double)(*inputp)->config[npops] - 3.0);
								else
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+38][listnumbers[count0-1]] = -10000;
								
								if((*inputp)->config[npops] > 2)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+39][listnumbers[count0-1]] 
									= (double)Y_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+39][listnumbers[count0-1]] = -10000;
								if((*inputp)->config[npops] > 2)
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+40][listnumbers[count0-1]] 
									= (double)Y_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
								else 
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+40][listnumbers[count0-1]] = -10000;
							#endif
							#if FREQSPECTRUM_TDFSR2 < 2		
                                if((*inputp)->config[npops] > 2) {
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]]
                                    = (double)neutpar[0].m_sdev;
                                }
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 2) {
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]]
                                    = (double)neutpar[0].m_skew;
                                }
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]] = -10000;
                                if((*inputp)->config[npops] > 3) {
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]]
                                    = (double)neutpar[0].m_kurt;
                                }
                                else
                                    matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]] = -10000;
							#endif
							#if FREQSPECTRUM_TDFSR2 == 0
								if(neutpar[0].S > 0) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+44][listnumbers[count0-1]] 
									= (double)neutpar[0].ragg;
								}
								#if ZNS_ACTIVE == 1
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]] 
								= (double)Zns_window(inputp,s0,s1,npops);
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+46][listnumbers[count0-1]] 
								= matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9] [listnumbers[count0-1]] -
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]];
								#else
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]] 
								= (double)-10000;
								matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+46][listnumbers[count0-1]] 
								= (double)-10000;
								#endif	
							#endif	
								/*end calc neutrality test*/
							}
							else {
								for(nv=0;nv<NEUTVALUES2;nv++) {
									matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]] = (double)-10000;
								}
							}
							/*printing values*/
							if((*inputp)->print_neuttest >/*=*/ 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0) {
								if(npops==0) {
									if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5)
										fprintf(file_outputm[aa],"%ld\t",/*listnumbers[*/count0-1/*]*/);
								#if SHOW_PRIORS
									for(nv=0;nv<(*inputp)->npriors;nv++) {
										fprintf(file_outputm[aa],"%f\t",priors[nv].priordist[count0-1]);
									}
								#endif
									if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
										fprintf(outfstallcompm[aa],"%ld\t",/*listnumbers[*/count0-1/*]*/);
										fprintf(outpiwallcompm[aa],"%ld\t",/*listnumbers[*/count0-1/*]*/);
										fprintf(outpiaallcompm[aa],"%ld\t",/*listnumbers[*/count0-1/*]*/);
										fflush(outfstallcompm[aa]);
										fflush(outpiwallcompm[aa]);
										fflush(outpiaallcompm[aa]);
										
										fprintf(outfsthapallcompm[aa],"%ld\t",/*listnumbers[*/count0-1/*]*/);
										fprintf(outhapwallcompm[aa],"%ld\t",/*listnumbers[*/count0-1/*]*/);
										fprintf(outhapaallcompm[aa],"%ld\t",/*listnumbers[*/count0-1/*]*/);
										fflush(outfsthapallcompm[aa]);
										fflush(outhapwallcompm[aa]);
										fflush(outhapaallcompm[aa]);
									}
									if((*inputp)->type_ancestral) {
										fprintf(outancestralm[aa],"%ld\t",/*listnumbers[*/count0-1/*]*/);
										if((*inputp)->type_ancestral > 3) {
											for(nv=0;nv<4*(*inputp)->type_ancestral;nv++) {
												if(neutpar[0].Sanc[nv] == (double)-10000) fprintf(outancestralm[aa],"na\t");
												else fprintf(outancestralm[aa],"%d\t",neutpar[0].Sanc[nv]);
											}
										}
										if((*inputp)->type_ancestral == 3) {
											for(nv=0;nv<10;nv++) {
												if(neutpar[0].Sanc[nv] == (double)-10000) fprintf(outancestralm[aa],"na\t");
												else fprintf(outancestralm[aa],"%d\t",neutpar[0].Sanc[nv]);
											}
										}
										if((*inputp)->type_ancestral == 2) {
											if(neutpar[0].Sanc[0] == (double)-10000) fprintf(outancestralm[aa],"na\t");
											else fprintf(outancestralm[aa],"%d\t",neutpar[0].Sanc[0]);
											if(neutpar[0].Sanc[1] == (double)-10000) fprintf(outancestralm[aa],"na\t");
											else fprintf(outancestralm[aa],"%d\t",neutpar[0].Sanc[1]);
											if(neutpar[0].Sanc[3] == (double)-10000) fprintf(outancestralm[aa],"na\t");
											else fprintf(outancestralm[aa],"%d\t",neutpar[0].Sanc[3]);
											if(neutpar[0].Sanc[8] == (double)-10000) fprintf(outancestralm[aa],"na\t");
											else fprintf(outancestralm[aa],"%d\t",neutpar[0].Sanc[8]);
										}
										if((*inputp)->type_ancestral == 1) {
											if(neutpar[0].Sanc[0] == (double)-10000) fprintf(outancestralm[aa],"na\t");
											else fprintf(outancestralm[aa],"%d\t",neutpar[0].Sanc[0]);
										}
										if((*inputp)->type_ancestral) fputs("\n",outancestralm[aa]);
										fflush(outancestralm[aa]);
									}
								}
								if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
							#if FREQSPECTRUM_TDFSR2 == 0		
									for(nv=0;nv<NEUTVALUES2;nv++) {
										if(matrix_test[nwindow*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]] == (double)-10000) 
											fprintf(file_outputm[aa],"na\t");
										else 
											fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]]);
									}
							#endif
							#if FREQSPECTRUM_TDFSR2 == 1		
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]]);
									if(matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]] == (double)-10000) 
										fprintf(file_outputm[aa],"na\t");
									else 
										fprintf(file_outputm[aa],"%g\t",matrix_test[aa*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]]);
#endif
									fprintf(file_outputm[aa],"%d\t",neutpar[0].mhsites);
									for(nv=1;nv<(*inputp)->config[npops];nv++) {
										fprintf(file_outputm[aa],"%d\t",neutpar[0].freq[nv]);
									}
									fflush(file_outputm[aa]);
								}
								if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
									if(neutpar[0].piwallcomp[npops] == (double)-10000) fprintf(outpiwallcompm[aa],"na\t");
									else fprintf(outpiwallcompm[aa],"%g\t",neutpar[0].piwallcomp[npops]);

									for(nv=npops+1;nv<(*inputp)->max_npop_sampled;nv++) {										
										if(neutpar[0].piaallcomp[nv] == (double)-10000) fprintf(outpiaallcompm[aa],"na\t");
										else fprintf(outpiaallcompm[aa],"%g\t",neutpar[0].piaallcomp[nv]);
										
										if(neutpar[0].fstallcomp[nv] == (double)-10000) fprintf(outfstallcompm[aa],"na\t");
										else fprintf(outfstallcompm[aa],"%g\t",neutpar[0].fstallcomp[nv]);
									}
									fflush(outfstallcompm[aa]);
									fflush(outpiwallcompm[aa]);
									fflush(outpiaallcompm[aa]);

									if(neutpar[0].hapwallcomp[npops] == (double)-10000) fprintf(outhapwallcompm[aa],"na\t");
									else fprintf(outhapwallcompm[aa],"%g\t",neutpar[0].hapwallcomp[npops]);
									
									for(nv=npops+1;nv<(*inputp)->max_npop_sampled;nv++) {
										if(neutpar[0].hapaallcomp[nv] == (double)-10000) fprintf(outhapaallcompm[aa],"na\t");
										else fprintf(outhapaallcompm[aa],"%g\t",neutpar[0].hapaallcomp[nv]);
										
										if(neutpar[0].fsthapallcomp[nv] == (double)-10000) fprintf(outfsthapallcompm[aa],"na\t");
										else fprintf(outfsthapallcompm[aa],"%g\t",neutpar[0].fsthapallcomp[nv]);
									}
									fflush(outfsthapallcompm[aa]);
									fflush(outhapwallcompm[aa]);
									fflush(outhapaallcompm[aa]);
								}

								if(npops == (*inputp)->max_npop_sampled-1) {
									if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
										fputs("\n",file_outputm[aa]);
										fflush(file_outputm[aa]);
									}
									if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
										fputs("\n",outfstallcompm[aa]);
										fputs("\n",outpiwallcompm[aa]);
										fputs("\n",outpiaallcompm[aa]);
										
										fflush(outfstallcompm[aa]);
										fflush(outpiwallcompm[aa]);
										fflush(outpiaallcompm[aa]);										

										fputs("\n",outfsthapallcompm[aa]);
										fputs("\n",outhapwallcompm[aa]);
										fputs("\n",outhapaallcompm[aa]);
										
										fflush(outfstallcompm[aa]);
										fflush(outpiwallcompm[aa]);
										fflush(outpiaallcompm[aa]);
									}
								}
							}
							
							nwindow += 1;
							s0 = s1;
						}
					}
				}
				else {
					if(rejflag==0 && npops < (*inputp)->npop_sampled) {
						init_coef(coef[0],(*inputp)->config[npops]);
						calc_neutpar(0,segsites,inputp,neutpar+0,(double)recombinationv,npops);
					#if FREQSPECTRUM_TDFSR2 < 2		
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]]
							= (double)tajima_d(neutpar[0].k,(int)neutpar[0].S,coef[0]);
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]]
							= (double)Fs((*inputp)->config[npops],neutpar[0].k,neutpar[0].nhapl);
					#endif
					#if FREQSPECTRUM_TDFSR2 == 0
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+2][listnumbers[count0-1]]
							= (double)fl_d2_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+3][listnumbers[count0-1]]
							= (double)fl_f2_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+4][listnumbers[count0-1]]
							= (double)fl_d_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+5][listnumbers[count0-1]]
							= (double)fl_f_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+6][listnumbers[count0-1]]
							= (double)fay_wu((*inputp)->config[npops],neutpar[0].freq,neutpar[0].k);
					#endif					
					#if FREQSPECTRUM_TDFSR2 < 2
						if(neutpar[0].S > 1) {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]]
								= (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]]
								= (double)neutpar[0].Q1/((double)neutpar[0].S);
						}
						else {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]] = -10000;
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]] = -10000;
						}
					#endif
					#if FREQSPECTRUM_TDFSR2 < 2		
						if(neutpar[0].S > 1) {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]]
							= (double)ZnA_(0,segsites,inputp,npops);/*INACTIVE*/
						}
						else {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]] = -10000;
						}
					#endif
					#if FREQSPECTRUM_TDFSR2 == 0
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+10][listnumbers[count0-1]]
								= (double)Fst(neutpar[0].piw,neutpar[0].pib,(int)(*inputp)->npop);
						
						if(neutpar[0].nhapl > -10000)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+11][listnumbers[count0-1]]
							= (double)neutpar[0].nhapl;
						else
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+11][listnumbers[count0-1]] =-10000;
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+12][listnumbers[count0-1]]
							= (double)testHap((*inputp)->config[npops],neutpar[0].fhapl);
					#endif
					#if FREQSPECTRUM_TDFSR2 < 2		
						if(neutpar[0].S > 0)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]]
								= (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp)->config[npops],(int)neutpar[0].S);
						else 
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]] = -10000;
					#endif
					#if FREQSPECTRUM_TDFSR2 < 3
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]]
							= (double)neutpar[0].S;/*MODIFIED*/
					#endif
					#if FREQSPECTRUM_TDFSR2 == 0
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+15][listnumbers[count0-1]]
							= (double)neutpar[0].k;
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+16][listnumbers[count0-1]]
							= (double)neutpar[0].pib;
					#endif
					#if FREQSPECTRUM_TDFSR2 < 2
						if(neutpar[0].S > -10000)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]]
							= (double)neutpar[0].S/(double)coef[0][0];
						else
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]] = -10000;
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]]
							= (double)neutpar[0].k;
					#endif
					#if FREQSPECTRUM_TDFSR2 == 0
						if(neutpar[0].k > -10000)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+19][listnumbers[count0-1]]
							= (double)neutpar[0].k - (double)(matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+6][listnumbers[count0-1]]);
						else
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+19][listnumbers[count0-1]] = -10000;
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+20][listnumbers[count0-1]]
							= (double)tajima_dvsdmin(neutpar[0].k,(int)neutpar[0].S,coef[0],(*inputp)->config[npops]);
					#endif
					#if FREQSPECTRUM_TDFSR2 < 2
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]]
							= (double)fay_wu_normalized2((int)(*inputp)->config[npops],(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0],neutpar[0].k);
					#endif
					#if FREQSPECTRUM_TDFSR2 == 0
						if(neutpar[0].k > -10000)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+22][listnumbers[count0-1]]
							= (double)neutpar[0].maxhapl/(double)(*inputp)->config[npops];
						else
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+22][listnumbers[count0-1]] = -10000;
						if(neutpar[0].maxhapl1 != -10000)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+23][listnumbers[count0-1]]
							= (double)neutpar[0].maxhapl1/(double)(*inputp)->config[npops];
						else
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+23][listnumbers[count0-1]]
							= (double)neutpar[0].maxhapl1;
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+24][listnumbers[count0-1]]
							= (double)neutpar[0].Rm;
						if(neutpar[0].S != -10000)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+25][listnumbers[count0-1]]
							= (double)neutpar[0].freq[1];
						else
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+25][listnumbers[count0-1]] = -10000;
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+26][listnumbers[count0-1]]
						= (double)neutpar[0].thetaL;
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+27][listnumbers[count0-1]]
						= (double)E_zeng((int)(*inputp)->config[npops],(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0]);
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+28][listnumbers[count0-1]]
						= (double)EWtest((int)(*inputp)->config[npops],neutpar[0].fhapl);
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+29][listnumbers[count0-1]]
						= (double)Fst(neutpar[0].withinw,neutpar[0].pib,(int)(*inputp)->npop);
						if(neutpar[0].S > 1) {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+30][listnumbers[count0-1]]
							= neutpar[0].min_uiHS;
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+31][listnumbers[count0-1]]
							= neutpar[0].max_uiHS;
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+32][listnumbers[count0-1]]
							= neutpar[0].max_iES;
						}
						else {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+30][listnumbers[count0-1]] = (double)-10000;
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+31][listnumbers[count0-1]] = (double)-10000;
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+32][listnumbers[count0-1]] = (double)-10000;
						}
					#endif
					#if FREQSPECTRUM_TDFSR2 < 2
						if((*inputp)->T_out > 0. || (*inputp)->pop_outgroup != -1) {
							if((*inputp)->pop_outgroup != -1) (*inputp)->Sout = 0;
							else neutpar[0].freq[0] = 0;
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]]
							= (double)fixoutg((*inputp)->config[npops],(*inputp)->Sout,neutpar[0].freq[0]);
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]]
							= (double)koutgJC((*inputp)->config[npops],(*inputp)->Sout,neutpar[0].freq,(*inputp)->nsites);
						}
						else {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]] = -10000;
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]] = -10000;
						}
					#endif
					#if FREQSPECTRUM_TDFSR2 == 0
						if((*inputp)->config[npops] > 2)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+35][listnumbers[count0-1]]
							= (double)neutpar[0].Se1/(coef[0][0]-1.);
                        else
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+35][listnumbers[count0-1]] = -10000;
						if((*inputp)->config[npops] > 3)
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+36][listnumbers[count0-1]]
							= (double)neutpar[0].Sn1/(coef[0][0]-((double)(*inputp)->config[npops]/((double)(*inputp)->config[npops]-1.0)));
                        else
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+36][listnumbers[count0-1]] = -10000;
						if((*inputp)->config[npops] > 2)
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+37][listnumbers[count0-1]]
							= (double)neutpar[0].pie1 * (double)(*inputp)->config[npops]/((double)(*inputp)->config[npops]-2.0);
                        else
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+37][listnumbers[count0-1]] = -10000;
						if((*inputp)->config[npops] > 3)
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+38][listnumbers[count0-1]]
							= (double)neutpar[0].pin1 * ((double)(*inputp)->config[npops] - 1.0)/((double)(*inputp)->config[npops] - 3.0);
						else
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+38][listnumbers[count0-1]] = -10000;
						
						if((*inputp)->config[npops] > 2)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+39][listnumbers[count0-1]]
							= (double)Y_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
						else 
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+39][listnumbers[count0-1]] = -10000;
						if((*inputp)->config[npops] > 2)
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+40][listnumbers[count0-1]]
							= (double)Y2_achaz((*inputp)->config[npops],neutpar[0].freq,neutpar[0].S);
						else 
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+40][listnumbers[count0-1]] = -10000;
					#endif
					#if FREQSPECTRUM_TDFSR2 < 2	
                        if((*inputp)->config[npops] > 2)
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]]
                            = (double)neutpar[0].m_sdev;
                        else
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]] = -10000;
                        if((*inputp)->config[npops] > 2)
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]]
                                = (double)neutpar[0].m_skew;
                        else
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]] = -10000;
                        if((*inputp)->config[npops] > 3)
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]]
                            = (double)neutpar[0].m_kurt;
                        else
                            matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]] = -10000;
					#endif
					#if FREQSPECTRUM_TDFSR2 == 0
						if(neutpar[0].S > 0) {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+44][listnumbers[count0-1]]
							= (double)neutpar[0].ragg;
						}
						#if ZNS_ACTIVE == 1
						matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]] 
						= (double)Zns(0,segsites,inputp,npops);
						matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+46][listnumbers[count0-1]] 
						= matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9] [listnumbers[count0-1]] -
						  matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]];
						#else
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+45][listnumbers[count0-1]]
						= (double)-10000;
						matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+46][listnumbers[count0-1]]
						= (double)-10000;
						#endif
					#endif	
					#if DEBUGSOUT
						fprintf(file_debug,"%d\n",(*inputp)->Sout);
					#endif
					#if PRINTTHETAS
						fprintf(file_debug_sel,"%f\t%f\t%f\n",neutpar[0].k,(double)neutpar[0].S/coef[0][0], 
							(double)neutpar[0].k - (double)(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+6][listnumbers[count0-1]]));
					#endif
					}
					else {
						for(nv=0;nv<NEUTVALUES2;nv++) {
							matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]] = (double)-10000;
						}
					}
					/*printing values*/
					if((*inputp)->print_neuttest >/*=*/ 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0) {
						if(npops==0) {
							if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5)
								fprintf(output,"%ld\t",count0-1);
						#if SHOW_PRIORS
							for(nv=0;nv<(*inputp)->npriors;nv++) {
								fprintf(output,"%f\t",priors[nv].priordist[count0-1]);
							}
						#endif
							if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
								fprintf(outfstallcomp,"%ld\t",count0-1);
								fprintf(outpiwallcomp,"%ld\t",count0-1);
								fprintf(outpiaallcomp,"%ld\t",count0-1);
								fflush(outfstallcomp);
								fflush(outpiwallcomp);
								fflush(outpiaallcomp);
								
								fprintf(outfsthapallcomp,"%ld\t",count0-1);
								fprintf(outhapwallcomp,"%ld\t",count0-1);
								fprintf(outhapaallcomp,"%ld\t",count0-1);
								fflush(outfsthapallcomp);
								fflush(outhapwallcomp);
								fflush(outhapaallcomp);
							}
							if((*inputp)->type_ancestral) {
								fprintf(outancestral,"%ld\t",count0-1);
								if((*inputp)->type_ancestral > 3) {
									for(nv=0;nv<4*(*inputp)->type_ancestral;nv++) {
										if(neutpar[0].Sanc[nv] == (double)-10000) fprintf(outancestral,"na\t");
										else fprintf(outancestral,"%d\t",neutpar[0].Sanc[nv]);
									}
								}
								if((*inputp)->type_ancestral == 3) {
									for(nv=0;nv<10;nv++) {
										if(neutpar[0].Sanc[nv] == (double)-10000) fprintf(outancestral,"na\t");
										else fprintf(outancestral,"%d\t",neutpar[0].Sanc[nv]);
									}
								}
								if((*inputp)->type_ancestral == 2) {
									if(neutpar[0].Sanc[0] == (double)-10000) fprintf(outancestral,"na\t");
									else fprintf(outancestral,"%d\t",neutpar[0].Sanc[0]);
									if(neutpar[0].Sanc[1] == (double)-10000) fprintf(outancestral,"na\t");
									else fprintf(outancestral,"%d\t",neutpar[0].Sanc[1]);
									if(neutpar[0].Sanc[3] == (double)-10000) fprintf(outancestral,"na\t");
									else fprintf(outancestral,"%d\t",neutpar[0].Sanc[3]);
									if(neutpar[0].Sanc[8] == (double)-10000) fprintf(outancestral,"na\t");
									else fprintf(outancestral,"%d\t",neutpar[0].Sanc[8]);
								}
								if((*inputp)->type_ancestral == 1) {
									if(neutpar[0].Sanc[0] == (double)-10000) fprintf(outancestral,"na\t");
									else fprintf(outancestral,"%d\t",neutpar[0].Sanc[0]);
								}
								if((*inputp)->type_ancestral) fputs("\n",outancestral);
								fflush(outancestral);
							}
						}
						if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
					#if FREQSPECTRUM_TDFSR2 == 0
							for(nv=0;nv<NEUTVALUES2;nv++) {
								if(matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]] == (double)-10000)
									fprintf(output,"na\t");
								else 
									fprintf(output,"%g\t",matrix_test[0*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+nv][listnumbers[count0-1]]);
							}
					#endif
					#if FREQSPECTRUM_TDFSR2 == 1 /*Tajima's D, Fu's Fs, theta and pi*/
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+0][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+1][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+7][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+8][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+9][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+13][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+14][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+17][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+18][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+21][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+33][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+34][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+41][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+42][listnumbers[count0-1]]);
							if(matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]] == (double)-10000) 
								fprintf(output,"na\t");
							else 
								fprintf(output,"%g\t",matrix_test[0/*(*inputp)->nloci*/*NEUTVALUES2*(*inputp)->max_npop_sampled+npops*NEUTVALUES2+43][listnumbers[count0-1]]);
						#endif
							fprintf(output,"%d\t",neutpar[0].mhsites);
							for(nv=1;nv<(*inputp)->config[npops];nv++) {
								fprintf(output,"%d\t",neutpar[0].freq[nv]);
							}
							fflush(output);
						}
						if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {							
							if(neutpar[0].piwallcomp[npops] == (double)-10000) fprintf(outpiwallcomp,"na\t");
							else fprintf(outpiwallcomp,"%g\t",neutpar[0].piwallcomp[npops]);
							
							for(nv=npops+1;nv<(*inputp)->max_npop_sampled;nv++) {
								if(neutpar[0].piaallcomp[nv] == (double)-10000) fprintf(outpiaallcomp,"na\t");
								else fprintf(outpiaallcomp,"%g\t",neutpar[0].piaallcomp[nv]);

								if(neutpar[0].fstallcomp[nv] == (double)-10000) fprintf(outfstallcomp,"na\t");
								else fprintf(outfstallcomp,"%g\t",neutpar[0].fstallcomp[nv]);
							}
							fflush(outfstallcomp);
							fflush(outpiwallcomp);
							fflush(outpiaallcomp);
							
							if(neutpar[0].hapwallcomp[npops] == (double)-10000) fprintf(outhapwallcomp,"na\t");
							else fprintf(outhapwallcomp,"%g\t",neutpar[0].hapwallcomp[npops]);
							
							for(nv=npops+1;nv<(*inputp)->max_npop_sampled;nv++) {
								if(neutpar[0].hapaallcomp[nv] == (double)-10000) fprintf(outhapaallcomp,"na\t");
								else fprintf(outhapaallcomp,"%g\t",neutpar[0].hapaallcomp[nv]);
								
								if(neutpar[0].fsthapallcomp[nv] == (double)-10000) fprintf(outfsthapallcomp,"na\t");
								else fprintf(outfsthapallcomp,"%g\t",neutpar[0].fsthapallcomp[nv]);
							}
							fflush(outfsthapallcomp);
							fflush(outhapwallcomp);
							fflush(outhapaallcomp);
						}

						if(npops == (*inputp)->max_npop_sampled-1) {
							if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5) {
								fputs("\n",output);
								fflush(output);
							}
							if(((*inputp)->max_npop_sampled > 2 && (*inputp)->print_neuttest < 5) || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
								fputs("\n",outfstallcomp);
								fputs("\n",outpiwallcomp);
								fputs("\n",outpiaallcomp);
								
								fflush(outfstallcomp);
								fflush(outpiwallcomp);
								fflush(outpiaallcomp);								
								
								fputs("\n",outfsthapallcomp);
								fputs("\n",outhapwallcomp);
								fputs("\n",outhapaallcomp);
								
								fflush(outfstallcomp);
								fflush(outpiwallcomp);
								fflush(outpiaallcomp);
							}
						}
					}
				}
			}
			if((*inputp)->likelihood_line == 0) {
				if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
					fprintf(file_outputpost,"%ld\t",count0-1);
				if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
					fprintf(file_outputpost,"%g\t",thetae);
				if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
					fprintf(file_outputpost,"%g\t",rece);
				if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
					fprintf(file_outputpost,"%g\t",lengtht);
					fprintf(file_outputpost,"%g\n",postp[0][count0-1].prob);
				}
			}
		}
	}
	*jcount2 = jcount;
	*mcount2 = mcount;		
	
    /* alliberar les matrius i vectors */
    free(posit);
    for(i=0;i<(*inputp)->nsam;i++)
        free(list[i]);
    free(list);
    
	free(neutpar[0].freq);
	free(neutpar[0].fhapl);
	free(neutpar[0].unic);
	free(neutpar[0].Sanc);
	if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
		free(neutpar[0].fstallcomp);
		free(neutpar[0].piwallcomp);
		free(neutpar[0].piaallcomp);
		free(neutpar[0].fsthapallcomp);
		free(neutpar[0].hapwallcomp);
		free(neutpar[0].hapaallcomp);
	}
	free(weightmut);
	free(weightrec);	

	/*Close file/s*/
	if(((*inputp)->likelihood_line == 0 && (*inputp)->pr_matrix && (*inputp)->linked)) {
		fclose(output);
	}
	if(((*inputp)->likelihood_line == 0 && (*inputp)->pr_matrix && (*inputp)->linked == 0) ||
	   ((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0 && (*inputp)->linked == 0)) {
		fclose(output);
	}
	if(((*inputp)->likelihood_line == 0 && (*inputp)->pr_matrix && (*inputp)->linked == 0) ||
	   ((*inputp)->print_neuttest >/*=*/ 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0 && (*inputp)->linked == 0)) {
		if(((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) && (*inputp)->print_neuttest >/*=*/ 2 && (*inputp)->print_neuttest < 5) {
			fclose(outfstallcomp);
			fclose(outpiwallcomp);
			fclose(outpiaallcomp);
			fclose(outfsthapallcomp);
			fclose(outhapwallcomp);
			fclose(outhapaallcomp);
			if((*inputp)->type_ancestral) fclose(outancestral);
		}
	}
	if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0 && (*inputp)->linked) {
		if((*inputp)->linked == 1) {
			x = 0;
			for(kk=0,ll=(*inputp)->window;kk<(*inputp)->nsites;kk += (long int)(*inputp)->despl) {
				if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5)
					fclose(file_outputm[x]);
			    if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
			        fclose(outfstallcompm[x]);
			        fclose(outpiwallcompm[x]);
					fclose(outpiaallcompm[x]);
			        fclose(outfsthapallcompm[x]);
			        fclose(outhapwallcompm[x]);
					fclose(outhapaallcompm[x]);
					if((*inputp)->type_ancestral) fclose(outancestralm[x]);
			    }
			    if(ll == (*inputp)->nsites) break;
				else ll += (*inputp)->despl;
				if(ll > (*inputp)->nsites) ll = (*inputp)->nsites;
				x++;
			}
		}
		else {
			for(x=0;x<(*inputp)->linked;x++) {
				if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5)
					fclose(file_outputm[x]);
			    if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
					fclose(outfstallcompm[x]);
					fclose(outpiwallcompm[x]);
					fclose(outpiaallcompm[x]);
					fclose(outfsthapallcompm[x]);
					fclose(outhapwallcompm[x]);
					fclose(outhapaallcompm[x]);
					if((*inputp)->type_ancestral) fclose(outancestralm[x]);
				}
			}
		}
		if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5)
			free(file_outputm);
		if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
			free(outfstallcompm);
			free(outpiwallcompm);
			free(outpiaallcompm);
			free(outfsthapallcompm);
			free(outhapwallcompm);
			free(outhapaallcompm);
			if((*inputp)->type_ancestral) free(outancestralm);
		}
	}
	if((*inputp)->print_neuttest > 2 && (*inputp)->print_neuttest < 5 && (*inputp)->likelihood_line == 0) {
		if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
			fclose(file_outputpost);
		}
	}
	
	#if DEBUGSOUT
    fclose(file_debug);
	#endif
    #if PRINTTHETAS
    fclose(file_debug_sel);
    #endif
    
    return 0;
}

void mod_mhits(long int segsites, struct var2 **inputp,double *weightmut)
{
    long int x,y,z;
    int r,h,i,j,k;
    char a[1];
    int Sout;
    double rr,ratio,r_transv,r_transc;
    double ran1(void);
    int *mhsout;
	double poissondist(double);
	double valuer;
	long int localize_positiontop(double *,double,long int,long int);
    
	ratio = (*inputp)->ratio_sv;
    r_transc = ratio/(ratio + 1.);
    r_transv = (0.5)/(ratio + 1.); /*suma de transc + 1nt a transversio, quan sumen els 2nt es total = 1*/
    Sout = (*inputp)->Sout;
    x = 0;
	*a = '0';
				        
    if(segsites != 0) {
		while(x<(int)segsites-1) {
			k = 0;
			while(posit[x] == posit[x+1]) {	/* buscar els mhits */
				k++;
				x++;
				if(x == (long int)segsites-1) break;
			}
			if(k) {				/* ordenar mhits de mes antic a mes recent (...)*/
				for(y=(x-k);y<x;y++) {
					j = 0;
					for(i=0;i<(*inputp)->nsam;i++) if(list[i][y] != '0') j++;
					for(z=y+1;z<(x+1);z++) {
						h = 0;
						for(i=0;i<(*inputp)->nsam;i++) if(list[i][z] != '0') h++;
						if(j < h) {
							for(i=0;i<(*inputp)->nsam;i++) {
								*a = list[i][y];
								list[i][y] = list[i][z];
								list[i][z] = *a;
							}
							j = h;
						}
					}
				}
				if((*inputp)->segsitesin == -1 || ((*inputp)->segsitesin >= 0 && (*inputp)->mhits == 1 && (*inputp)->Sfix_alltheta)) {/*in case no fixed mutations*/
					for(y=(x-k+1);y<x+1;y++) {	/* afegir les mutacions a la posicio mes antiga */
						r = -1;
						for(i=0;i<(*inputp)->nsam;i++) {
							if(list[i][y] != '0') {
								if(r == -1) {
									rr = (double)ran1();	/* inclou ratio trans/transv */
									if(rr < r_transc) r = 0;	/* transicio, la resta transversions */
									else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 1;
									else r = 2;
									if(list[i][x-k] == '0') {
										if(r==0) list[i][x-k] = *a = '1';
										else if(r==1) list[i][x-k] = *a = '2';
										else if(r==2) list[i][x-k] = *a = '3';
									}
									else
										if(list[i][x-k] == '1') {
											if(r==0) list[i][x-k] = *a = '0';
											else if(r==1) list[i][x-k] = *a = '2';
											else if(r==2) list[i][x-k] = *a = '3';
										}
										else
											if(list[i][x-k] == '2') {
												if(r==0) list[i][x-k] = *a = '3';
												else if(r==1) list[i][x-k] = *a = '0';
												else if(r==2) list[i][x-k] = *a = '1';
											}
											else
												if(list[i][x-k] == '3') {
													if(r==0) list[i][x-k] = *a = '2';
													else if(r==1) list[i][x-k] = *a = '1';
													else if(r==2) list[i][x-k] = *a = '0';
												}
								}
								else list[i][x-k] = *a;
							}
						}
					}
				}
				else {
					for(y=(x+1-k),r=2;y<x+1;y++,r++) {	/* mutacions fix, totes s'han de veure (fins a tres mutacions) */
						for(i=0;i<(*inputp)->nsam;i++) 
							if(list[i][y] != '0') 
								list[i][x-k] = r + '0';
					}
					if(r>4) {
						perror("Error: more than 3 mutations in one position");
						exit(1);
					}
				}
			}
			/*in case no outgroup but mhits: assign '0' randomly to the mhits positions (when no '0' is observed but in putative outgroup)*/
			if((*inputp)->T_out == 0.) {
				r = 0;
				for(i=0;i<(*inputp)->nsam;i++) if(list[i][x-k] == '0') r++;
				if(r==0) {
					r = list[0][x-k];
					for(i=0;i<(*inputp)->nsam;i++) {
						if(list[i][x-k] == r) list[i][x-k] = '0';
					}
				}
			}
			x++;
		}
	}
    /* Treure els mhits de la branca de l'outgroup (no es mostra l'outgroup i per tant es veuen menys). */
    /* Per Sout mutacions en (nsites - segsites) */
	/* incorporacio per mhits en especiacio.Els mhits al outgroup no s'observen, nomes canviem el nt ancestral, pero no estan indicades!*/
	if(Sout) {
        if(!(mhsout = (int *)calloc((long int)((*inputp)->nsites),sizeof(int)))) {
            perror("calloc error in mhits.0b");
            exit(1);
        }
        for(x=0;x<Sout;x++) {
			valuer = ran1()*(double)(weightmut[(*inputp)->nsites-1]);
			y = localize_positiontop(weightmut,valuer,(long int)0,(*inputp)->nsites);
            
			rr = (double)ran1();
            if(rr < r_transc) r = 1;	/* transicio, la resta transversions */
            else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 2;
			else r = 3 ;
			
			if(mhsout[y] == 0) 
				mhsout[y] = r;
            else {
				if(mhsout[y] == 1) {
					if(r==1) mhsout[y] = 0;
					else mhsout[y] = r;
				}
				else {
					if(mhsout[y] == 2) {
						if(r==1) mhsout[y] = 3;
						else mhsout[y] = r-2;
					}
					else {
						if(mhsout[y] == 3) {
							if(r==1) mhsout[y] = 2;
							else mhsout[y] = r-2;
						}
					}
				}
			}
			/*mutation outgroup-ingroup: Els mhits al outgroup no s'observen, nomes canviem el nt ancestral, pero no estan indicades!*/
			for(k=0;k<segsites;k++) if(posit[k] == y) break;
			if(k<segsites) {
				mhsout[y] = 0;
				
				for(i=0;i<(*inputp)->nsam;i++) {
					if(list[i][k] == '0') list[i][k] = (char)r + '0';
					else if (list[i][k] == (char)r + '0') list[i][k] = '0';
				}
			}
        }
        Sout = 0;
        for(x=0;x<(*inputp)->nsites;x++) if(mhsout[x] > 0) Sout++;
        free(mhsout);
    }
	
    (*inputp)->Sout = Sout;
}
void mod_outgroup(long int segsit, struct var2 **inputp,int outgroup) 
{
	int j,h,polc,inito,nt0;
	int oldn;
	
	for(j=0;j<segsit;j++) {
		nt0 = 0;
		for(h=0;h<(*inputp)->pop_outgroup;h++) nt0 += (*inputp)->config[h];	
		inito = list[nt0][j]; /*nt of the outgroup*/
		polc = 0;
		for(h=nt0+1;h<nt0+(*inputp)->config[(*inputp)->pop_outgroup];h++) {
			if(list[h][j] != inito) {
				polc = 1; /*if outgroup polymorphic*/
				break;
			}
		}
		if(polc == 1) {
			if(inito < list[h][j]) oldn = inito;
			else oldn = list[h][j];
		}
		else oldn = inito;
		
		if(oldn != '0') {
			for(h=0;h<(*inputp)->nsam;h++) {
				if(list[h][j] == (char)oldn) list[h][j] = '0';
				else if(list[h][j] == '0') list[h][j] = (char)oldn;
			}
		}
	}
	
	return;
}

char **cmatrix(int nsam,long int len)	/* defineix l'espai per col.locar els polimorfismes */
{
    int i;
    char **m;
    
    if(!(m=(char **)malloc((unsigned)nsam*sizeof(char *)))) perror("alloc error in cmatrix");
    for(i=0;i<nsam;i++)
        if(!(m[i] = (char *)malloc((long int)len*sizeof(char)))) perror("alloc error in cmatrix.2");
    return(m);
}

long int gensam(long int npop,int nsam,int inconfig[],long int nsites,double theta,long int segsites,
	double r,double f,double track_len,double mig_rate,int mhits,long int iteration, double *factor, double *lengtht,
	int ifselection, double pop_sel, double sinit,double pop_size,long int sel_nt,double T_out, int *Sout,int *nintn,double **nrec,
	double **nrec2,double **tpast,int split_pop, double time_split, double time_scoal, double factor_anc, double *freq, 
	double tlimit,int iflogistic,double *ts, double factor_chrn, double rsv,double *weightmut,double *weightrec,double **migrate_matrix,
	int my_rank,int npop_events,struct events *pop_event,int linked,long int **loci_linked,int event_forsexratio,double event_sexratio,double sex_ratio,int no_rec_males,double sendt,double sfreqend,double sfreqinit)
{
	int i,ii;
    long int nsegs,seg,ns,start,end,len,segsit,k; 
    struct segl *seglst;
    double tseg,tt;
    double *pk,tout=0.;double tout2=0.;
    long int *ss;
    long int *len2;
    long int mmax;
	double r_transc,r_transv;
    struct segl *segtre_mig(long int ,int,int *,long int,double,double,double,
        double ,long int *,long int,
        double *,int,double,double,double,long int,int *,int *,int *,double **,double **,double **,
        int, double, double, double, double *,double,int,double *,double,double *,double **,int,
		int,struct events *,int, double,double,int, double, double,double);/* used to be: [MAXSEG]; */
    double ttime(struct node *, int);
    double poissondist(double), ran1(void);
    void biggerlist(int);
    void make_gametes(int,struct node *,double,long int,long int,int,double,double);
    void locate(long int,long int,long int *,int,double *,long int);
    void locate2(long int,long int,long int *,int,int,double *,long int);
    void mnmial2(long int,long int,double *,long int *,long int *);
    void mnmial2_psel(long int,long int,double *,long int *,long int *,long int);
    /*partial selection*/
    int all_sel,*selnsam; /*number of lines under selection and the vector with the lines (the first all_sel lines)*/
    int segsit_sel=0;/*parameter for partial selection*/    
    void locate_psel(long int,long int,long int *,int,long int,double *,long int);
    void locate2_psel(long int,long int,long int *,int,int,long int,double *,long int);
	double wstartm1,ttt;
	long int *len_nozero,nsites_nozero,nz;
	int loopcount;
	int aa,getinto;
	long int kk,ll,mm,bb,cc,dd;
    long int nsites_mod_recinf=nsites;
    
	if(!(selnsam = (int *)malloc((nsam)*sizeof(int)))) perror("malloc error sel. gensam.");
	all_sel = 0;

    seglst = segtre_mig(npop,nsam,inconfig,nsites_mod_recinf,r,f,track_len,mig_rate,&nsegs,
        iteration,factor,ifselection,pop_sel,sinit,pop_size,sel_nt,&all_sel,selnsam,
        nintn, nrec, nrec2, tpast,split_pop,time_split,time_scoal,factor_anc,freq,tlimit,
		iflogistic,ts,factor_chrn,weightrec,migrate_matrix,my_rank,npop_events,pop_event,
		event_forsexratio,event_sexratio,sex_ratio,no_rec_males,sendt,sfreqend,sfreqinit);
    r_transc = rsv/(rsv + (double)1.);
    r_transv = ((double)0.5)/(rsv + (double)1.); /*suma de transc + 1nt a transversio, quan sumen els 2nt es total = 1*/

	/*heterogeneity*/
	if(!(len_nozero = (long int *)malloc((nsegs)*sizeof(long int)))) perror("malloc error len. gensam.");
	nsites_nozero = 0;
	for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {
		end = (k<nsegs-1 ? (long int)seglst[seglst[seg].next].beg -1 : nsites_mod_recinf-1); /*next és l'index, beg és el punt físic */
		start = seglst[seg].beg;
		for(len_nozero[k]=0,nz=start;nz<=end;nz++) {
			if(nz==0) wstartm1 = (double)0;
			else wstartm1 = weightmut[nz-1];
			if(((weightmut[nz]-wstartm1)*theta) != (double)0) {
				len_nozero[k] += 1;
				nsites_nozero += 1;
			}
		}
	}
	
	/*Include mutations*/
    if(segsites == -1) {
        ns = 0;
        *Sout = 0;
		*lengtht = (double)0;
        for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {
			end = (k<nsegs-1 ? (long int)seglst[seglst[seg].next].beg -1 : nsites-1); /*next és l'index, beg és el punt físic */
			start = seglst[seg].beg;
			getinto = 0;
			if(linked > 1) {/*avoid calculation of mutations in not interested regions*/
				for(aa=0;aa<linked;aa++) {
					kk=loci_linked[aa][0];
					ll=loci_linked[aa][1]+1;
					if((kk >= start && kk <  end) ||
					   (ll >= start && ll <  end) || 
                       (kk >= start && ll <  end) || 
                       (kk <  start && ll >= end)) {
						getinto = 1;
						break;
					}
				}
			}
			else {
				getinto = 1;
			}
			if(getinto) {
				/*including heterogeneity*/
				len = end - start + 1;
				if(start==0) wstartm1 = (double)0;
				else wstartm1 = weightmut[start-1];
				tseg = (double)(weightmut[end] - wstartm1)*theta;
				tt = ttime(seglst[seg].ptree,nsam); /* Ttot pel segment, en funció de 4No respecte la pob 0*/
				*lengtht += (double)tt*((double)len/(double)nsites); /*length in 4N generations, no matter the mutation rate*/
				if(mhits) {/*T_out is the time to the ancestor in 4N*//*modified Jan2009*/
					if(T_out == 0.) {
						*Sout = 0;
					}
					else {
						tout = 2.0 * T_out;/*T_out is considered as a fixed value, not a parameter.*/
						tout += -1.0 *(double)log((double)((double)1-ran1()));/*sum the time after divergence, assuming No=1*/
						tout -= (seglst[seg].ptree + 2*nsam-2)->time; /*substract the distance of the sample*/
						if(tout < (double)0) tout = (double)0;	/*Outgroup can not accumulate negative mutations*/
						*Sout += (int)poissondist((double)(tseg*tout));	/* Sout needed to calculate hka and mhits */
					}
				}
				else *Sout = 0;
				
				loopcount = -1;
                if(len==1 && mhits==0)  {
                    if(tseg*tt > ran1()) segsit = 1;
                    else segsit = 0;
                }
                else {
                    segsit = (long int)poissondist((double)(tseg*tt));		/* nombre de mutacions al segment */
                }
                if(segsit == (long int)0 && all_sel > (int)0 && all_sel < nsam && (long int)sel_nt >= (long int)start && (long int)sel_nt <= (long int)end)
                    segsit = 1;/*we force the selective mut*/
                loopcount += 1;
                if(loopcount > 100) {
                    printf("\nSorry, the length of the sequence is too short to include so much mutations: try mhits 1.\n");
                    exit(1);
                }
                if(segsit > len_nozero[k] && mhits == 0)
                    segsit = len_nozero[k];
				if((segsit + ns) >= maxsites) {	/* refem la matriu dels polimorfismes */
					maxsites = segsit + ns + SITESINC;
					posit = (long int *)realloc(posit,(long int)maxsites*sizeof(long int));
					/* canvia mida del vector dels nombres dels polimorfismes */
					if(posit==NULL) perror("realloc error. gensam.1");
					biggerlist(nsam);	/* refem la llista dels polimorfismes */
				}
				/*partial selection*//*not well debugged yet*/ 
				if(all_sel > (int)0 && all_sel < (int)nsam && (long int)sel_nt >= (long int)start && (long int)sel_nt <= (long int)end
                   && (sfreqend == 0.1 || sendt == 1E6)) {
					segsit_sel = 1;
				}
				/*partial selection*/
				if(segsit_sel == 1) {
					locate_psel(segsit,start,posit+ns,mhits,sel_nt,weightmut,end);
					ii = (int)ns;
					while((long int)posit[ii] != (long int)sel_nt) ii++; 
					for(i=0;i<nsam;i++) {
						list[i][ns] = list[i][ii];
						if(all_sel>i) list[i][ii] = '1';
						else list[i][ii] = '0';
					}
					segsit_sel = 0;
				}
				else locate(segsit,start,/*len,*/posit+ns,mhits,weightmut,end);/* posa el nombre de les mutacions a la matriu */
				
				/*begin of substracting mutations out of range (only for the regions included in the study)*/
				/**/
				if(linked > 1) {
					cc = 0;
					dd = 0;
					mm = (long int)0;
					for(aa=0;aa<linked;aa++) {
						kk=loci_linked[aa][0];
						ll=loci_linked[aa][1]+1;
						for(bb=dd+cc;bb<segsit;bb++) {
							if((posit[bb] < kk) && (posit[bb] >= mm)) {
								dd++;
							}
							else {
								if((posit[bb] >= kk) && (posit[bb] < ll)) {
									posit[cc] = posit[bb];
									cc++;
								}
							}
						}
						mm = ll;
					}
					segsit = cc;
				}
				/**/
				/*end of substracting mutations out of range*/
				
				/*changed the order of the functions!! now here down*/
				make_gametes(nsam,seglst[seg].ptree,tt,segsit-segsit_sel,ns+segsit_sel,mhits,r_transc,r_transv);	/* posa les mutacions  a list*/
				free(seglst[seg].ptree);

				ns += segsit;
			}
		/**/	
		}
    }
    else {/*THE PARTIAL SELECTION WITH SEG. SITE ARE NOT WELL DEBUGGED YET IN THE Sfix METHOD*/
        /* en cas de S fix */
		pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
        ss = (long int *)malloc((unsigned)(nsegs*sizeof(long int)));
        len2 = (long int *)malloc((unsigned)(nsegs*sizeof(long int)));
        if((pk==NULL)||(ss==NULL)||(len2 ==NULL)) perror("malloc error. gensam.2");
		/*multiple hits for fixed mutations*/
		if(mhits) mmax = (3 < nsam ? (long int)3 : (long int)(nsam-1));
		else mmax = (long int)1; /*if no multiple hits available*/
        /*set time  and nsites_nozero to zero*/
		tt = ttt = (double)0;
		if(mhits) tout = (double)0;
		/*calcular primer la mida total de tot l'arbre*/
		for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {		
            end = (k<nsegs-1 ? (long int)seglst[seglst[seg].next].beg -1 : nsites_mod_recinf-1);
            start = seglst[seg].beg;
			/*when it is not possible to fix the specified segsites*/
			if(nsites_nozero*mmax < segsites) {
				printf("\n*****WARNING: It is not possible to fix %ld mutations!!. Instead used 0 mutations.*****\n",segsites);
				segsites = 0;
			}
			
			len = end - start + 1;
			if(start==0) wstartm1 = (double)0;
			else wstartm1 = weightmut[start-1];
			tseg = (double)(weightmut[end] - wstartm1)*theta;
			len2[k] = len_nozero[k]*mmax;
			pk[k] = ttime(seglst[seg].ptree,nsam) * tseg;/*time per chromosome section (in function of mutational rate)*/
			tt += pk[k];
			ttt += pk[k]/tseg * (double)len/(double)nsites_mod_recinf;/*time to be counted by lengtht (not in function of mutational rate)*/
			if(mhits && T_out > 0) {		/* incorporacio per mhits */
				tout2 = 2.0*T_out;
				tout2 += -1.0*(double)log((double)((double)1-ran1()));/*time after divergence, assuming equal No*/
				tout2 -= (seglst[seg].ptree + 2*nsam-2)->time; /*substract the distance of the sample*/
				if(tout2 < (double)0) tout2 = (double)0;	/*Outgroup can not accumulate negative mutations*/
				tout += tout2 * tseg;/*time to the outgroup (in function of mutational rate)*/
			}
			else tout = 0.;
        }
        *lengtht = (double)ttt; /*afegit per Sfix_allthetas*/
        for(k=0;k<nsegs;k++) pk[k] /= tt;	/* aleshores dividir el temps proporcionalment per situar les mutacions */
        if(mhits && T_out > 0 && theta) {		/* incorporacio per mhits en especiacio, only if theta is defined*/
            *Sout = (int)poissondist(/**/(double)(theta*tout));
        }
		else *Sout = 0;
		/*partial selection*/
		if(all_sel > (int)0 && all_sel < (int)nsam)
			mnmial2_psel((long int)segsites,nsegs,pk,ss,len2,(long int)sel_nt*mmax);
        else 
			mnmial2((long int)segsites,nsegs,pk,ss,len2);/* afegit, per evitar mes mutacions que posicions, en mhits mes de 3xpos */
        ns = 0;/*mnmial2 distribueix segsites al llarg de la secuencia*/
        for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {
            end = (k<nsegs-1 ? (long int)seglst[seglst[seg].next].beg -1 : nsites_mod_recinf-1);
            start = seglst[seg].beg;
            len = end - start + 1;
			if(start==0) wstartm1 = (double)0;
			else wstartm1 = weightmut[start-1];
			tseg = (double)(weightmut[end] - wstartm1)*theta;
            /*partial selection*/
			segsit_sel = 0;
            if(all_sel > (int)0 && all_sel < (int)nsam && (long int)sel_nt >= (long int)start && sel_nt <= (long int)end && segsites > 0 && (sfreqend == 0.1 || sendt == 1E6)) {
                segsit_sel = 1;
            }
            make_gametes(nsam,seglst[seg].ptree,tt*pk[k]/tseg,ss[k]-segsit_sel,ns+segsit_sel,mhits,r_transc,r_transv);/*posa a la matriu list les mutacions*/
            /*partial selection*/
            if(segsit_sel == 1) {
                locate2_psel(ss[k],start,/*len,*/posit+ns,mhits,nsam,(long int)sel_nt,weightmut,end);
				ii = (int)ns;
				while((long int)posit[ii] != (long int)sel_nt) ii++; 
				for(i=0;i<nsam;i++) {
					list[i][ns] = list[i][ii];
					if(all_sel>i) list[i][ii] = '1';
					else list[i][ii] = '0';
				}
                segsit_sel = 0;
            }
            else locate2(ss[k],start,posit+ns,mhits,nsam,weightmut,end); /* posa el nombre de les mutacions a la matriu */
            /* modificat per mhits i fix muts */
            free(seglst[seg].ptree);
            ns += ss[k];
        }
        free(pk);
        free(ss);
        free(len2);
    }
    free(selnsam);
	free(len_nozero);
    return(ns);
}
void biggerlist(int nsam)	/* fa més gran la matriu dels polimorfismes */
{
    int i;
    for(i=0;i<nsam;i++) {
        list[i] = (char *)realloc(list[i],maxsites*sizeof(char));
        if(list[i] == NULL) perror("realloc error. biggerlist");
    }
}
double ttime(struct node *ptree, int nsam)	/* la Ttot de l'arbre */
{
    double t;
    int i;
    
    t = (ptree + 2*nsam-2)->time;
    for(i=nsam;i<2*nsam-1;i++)
        t += (ptree + i)->time;
    return(t);
}

void make_gametes(int nsam, struct node *ptree, double tt,long int newsites,long int ns, int mhits, double r_transc,double r_transv)
{					/* posa les mutacions a la matriu list */
    long int j;
    int tip,node;
    int pickb(int, struct node *,double);
    int tdesn(struct node *,int,int);
	
	double rr;
	double ran1(void);
	char r;
	
	
    for(j=ns;j<ns+newsites;j++) {
		if(mhits) {
			rr = (double)ran1();
			if(rr < r_transc) r = '1';/* transicio, la resta transversions */
			else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = '2';
			else r = '3';
		}
		else r='1';
        node = pickb(nsam,ptree,tt);	/* busca una branca al'atzar en funcio de la mida de t*/
        for(tip=0;tip<nsam;tip++) {
            if(tdesn(ptree,tip,node))	{ /* posa mutació si la mostra té relació amb la branca */
                list[tip][j] = r; 
			}
            else
                list[tip][j] = '0';
        }
    }
}
int pickb(int nsam, struct node *ptree,double tt) /* agafa la branca a on ha caigut la mutació */
{
    double x,y,z;
    int i;
    double ran1(void);
    
    x = (double)ran1()*tt;
    for(i=0,y=0;i<2*nsam-2;i++) {
        z = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
        if(z < 0.) z = 0.; /* en cas d'extincio, per si hi han problemes de precissio */
        y += z;
        if(y >= x) return(i);
    }
    return(i);
}
int tdesn(struct node *ptree, int tip, int node) /* mira si la mostra que mirem esta relacionada amb la branca */
{
    int k;
    
    for(k=tip;k<node;k = (ptree+k)->abv);
    
    if(k==node) return(1);
    else return(0);
}
void mnmial(long int n,long int nclass,double *p,long int *rv)
{
    
    double x,s;
    long int i,j;
    double ran1(void);
    
    for(i=0;i<(long int)nclass;i++) rv[i]=0;	/* inicialitzar */
    for(i=0;i<(long int)n;i++) {	/* posa les n mutacions */
        x = (double)ran1();	/* agafa un nombre [0,1) */
        j = 0;
        s = p[0];	/* p és el temps total de cada segment relatiu a 1 */
        while((x>s) && (j<((long int)nclass-(long int)1))) s += p[++j]; /* busca el segment fins que ran sigui major que s, posa a j+1 */
        rv[j]++; 	/* afegeix la mutació al segment j */
    }
}
void mnmial2(long int n,long int nclass,double *p,long int *rv,long int *len)
{
    double x,s;
    long int i,j;
    double ran1(void);
    
    for(i=0;i<(long int)nclass;i++) rv[i]=0;	/* inicialitzar */
    i = 0;
    while(i<(long int)n) {	/* posa les n mutacions */
        x = (double)ran1();	/* agafa un nombre [0,1) */
        j = 0;
        s = p[0];	/* p és el temps total de cada segment relatiu a 1 */
        while((x>s) && (j<((long int)nclass-(long int)1))) s += p[++j]; /* busca el segment fins que ran sigui major que s, posa a j+1 */
        if(len[j] > 0) {
            len[j]--;	
            rv[j]++; 	/* afegeix la mutació al segment j */
            i++;
        }
    }
}
void mnmial2_psel(long int n,long int nclass,double *p,long int *rv,long int *len,long int sel_nt)
{
    double x,s;
    long int i,j;
    double ran1(void);
	long int sumlen;
    
    for(i=0;i<(long int)nclass;i++) rv[i]=0;	/* inicialitzar */
    j=0;/*locate the selected position*/
	sumlen = len[j];
	while(sumlen<sel_nt) {
		j++;
		sumlen +=len[j];
	}
	rv[j] += 1;
	len[j] -=1;
	i = 1;	
    while(i<(long int)n) {	/* posa les altres n-1 mutacions */
        x = (double)ran1();	/* agafa un nombre [0,1) */
        j = 0;
        s = p[0];	/* p és el temps total de cada segment relatiu a 1 */
        while((x>s) && (j<((long int)nclass-(long int)1))) s += p[++j]; /* busca el segment fins que ran sigui major que s, posa a j+1 */
        if(len[j] > 0) {
            len[j]--;	
            rv[j]++; 	/* afegeix la mutació al segment j */
            i++;
        }
    }
}
void locate_psel(long int n,long int beg,long int *ptr,int mhits,long int sel_nt,double *weightmut,long int end)
/* localitza les mutacions en un fragment */
{
     /*long int i;*/
     void ordran_psel(long int,long int *,int,long int,long int,double *,long int);
     
     ordran_psel(n,ptr,/*len,*/mhits,sel_nt,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
}

void ordran_psel(long int n,long int *pbuf,int mhits,long int sel_nt,long int beg,double *weightmut,long int end)
{
    void ranvec_psel(long int,long int *,/*long int,*/int,long int,long int,double *,long int);
    void order(long int,long int *);

    ranvec_psel(n,pbuf,/*len,*/mhits,sel_nt,beg,weightmut,end);
    order(n,pbuf);
}
void ranvec_psel(long int n,long int *pbuf,int mhits,long int sel_nt,long int beg,double *weightmut,long int end)
 /* posa un nombre entre [0,len) */
{
    long int i,x;
    double ran1(void);
	double valuer,wstartm1;
	long int localize_positiontop(double *,double,long int,long int);
    
	/*include nsites*/
    pbuf[0] = (long int)sel_nt;
	if(beg==0) wstartm1 = (double)0;
	else wstartm1 = (double)weightmut[beg-1];
    for(i=1;i<(long int)n;i++) {
        valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
		pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
		if(pbuf[i] == pbuf[0]) {
			i--;
			continue;
		}   
        if(!mhits) {	/* bucle per no mhits (mhits=0) */
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

void locate(long int n,long int beg,long int *ptr,int mhits,double *weightmut,long int end)
/* localitza les mutacions en un fragment */
{
     /*long int i;*/
     void ordran(long int,long int *,int,long int,double *,long int);
     
     ordran(n,ptr,/*len,*/mhits,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
}

void ordran(long int n, long int *pbuf,int mhits,long int beg,double *weightmut,long int end)
{
    void ranvec(long int,long int *,int,long int,double *,long int);
    /* posa un nombre entre [0,len) */
    void order(long int,long int *);/* ordena els valors */

    ranvec(n,pbuf,mhits,beg,weightmut,end);
    order(n,pbuf);
}

void ranvec(long int n, /*double */long int *pbuf,int mhits,long int beg,double *weightmut,long int end) /* posa un nombre entre [0,len) */
{
    long int i,x;
    double ran1(void);
	double valuer,wstartm1;
	long int localize_positiontop(double *,double,long int,long int);
    
    
	if(beg==0) wstartm1 = (double)0;
	else wstartm1 = (double)weightmut[beg-1];
    for(i=0;i<(long int)n;i++) {
        valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
		pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
        if(!mhits) {	/* bucle per no mhits (mhits=0) */
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}
void locate2(long int n,long int beg,long int *ptr,int mhits,int nsam,double *weightmut,long int end)/* localitza les mutacions en un fragment */
{
     /*long int i;*/
     void ordran2(long int,long int *,int,int,long int,double *,long int);
     
     ordran2(n,ptr,mhits,nsam,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
}
void ordran2(long int n, long int *pbuf,int mhits,int nsam,long int beg,double *weightmut,long int end)
{
    void ranvec2(long int,long int *,int,int,long int,double *,long int);
    /* posa un nombre entre [0,len) */
    void order(long int,long int *);/* ordena els valors */

    ranvec2(n,pbuf,mhits,nsam,beg,weightmut,end);
    order(n,pbuf);
}
void ranvec2(long int n, long int *pbuf,int mhits,int nsam,long int beg,double *weightmut,long int end)
/* posa un nombre entre [0,len) */
{
    long int i,x;
    int y;
    int z,f;
	
    double dlen;
    double a;
    double ran1(void);
	double valuer,wstartm1;
	long int localize_positiontop(double *,double,long int,long int);
    
	if(beg==0) wstartm1 = (double)0;
	else wstartm1 = (double)weightmut[beg-1];
    if(mhits) {
        dlen = (double)(weightmut[end] - wstartm1);
        for(i=0;i<(long int)n;i++) {
            a = (double)ran1()*dlen + (double)wstartm1;
            pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);   
            x = i-1;
            y = 1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    y++;
                    if(y == 3 || y == nsam) {/*no mes de 4 nt per posicio*/
                        a = (double)ran1()*dlen + (double)wstartm1;
                        pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);   
                        x = i;
                        y = 1;
                    }
					f=0;
					for(z=0;z<nsam;z++) /*all mutations should be observed. equal pattern in different position*/
						if((list[z][i] == '0' && list[z][x] == '0') ||
						   (list[z][i] != '0' && list[z][x] != '0')) 
							f++;
					if(f==nsam) {
						a = (double)ran1()*dlen + (double)wstartm1;
						pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);   
						x = i;
						y = 1;
						break;
					}
                }
                x--;
            }
		}
    }
    else { /* per no mhits */
        for(i=0;i<(long int)n;i++) {
			valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
			pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

void order(long int n,long int *pbuf)/* ordena els valors: es molt lent per un gran nombre de polimofismes! */
{
    long int gap,i,j;
    long int temp;
    
    for(gap= n/2; gap>0;gap /= 2)
        for(i=gap;i<(long int)n;i++)
            for(j=i-gap;j>=0 && pbuf[j]>pbuf[j+gap];j -= gap) {
                temp = pbuf[j];
                pbuf[j] = pbuf[j+gap];
                pbuf[j+gap] = temp;
            }
}

void locate2_psel(long int n,long int beg,long int *ptr,int mhits,int nsam,long int sel_nt,double *weightmut,long int end)/* localitza les mutacions en un fragment */
{
     /*long int i;*/
     void ordran2_psel(long int,long int *,int,int,long int,long int,double *,long int);
     
     ordran2_psel(n,ptr,mhits,nsam,sel_nt,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
}
void ordran2_psel(long int n,long int *pbuf,int mhits,int nsam,long int sel_nt,long int beg,double *weightmut,long int end)
{
    void ranvec2_psel(long int,long int *,int,int,long int,long int,double *,long int);/* posa un nombre entre [0,len) */
    void order(long int,long int *);/* ordena els valors */

    ranvec2_psel(n,pbuf,mhits,nsam,sel_nt,beg,weightmut,end);
    order(n,pbuf);
}
void ranvec2_psel(long int n,long int *pbuf,int mhits,int nsam,long int sel_nt,long int beg,double *weightmut,long int end)
/* posa un nombre entre [0,len) */
{
    long int i,x;
    int y;
    int z,f;
    double dlen;
    double a;
	double ran1(void);
	double valuer,wstartm1;
	long int localize_positiontop(double *,double,long int,long int);
    
	if(beg==0) wstartm1 = (double)0;
	else wstartm1 = (double)weightmut[beg-1];
    if(mhits) {
        dlen = (double)(weightmut[end] - wstartm1);
        for(i=0;i<(long int)n;i++) {
            if(i) {
				a = (double)ran1()*dlen + (double)wstartm1;
				pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);
			}
			else pbuf[0] = (long int)sel_nt;
			
            x = i-1;
            y = 1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    y++;
                    if(y == 3 || y == nsam || x == 0) {/*no mes de 4 nt per posicio*/
                        a = (double)ran1()*dlen + (double)wstartm1;
                        pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);   
                        x = i;
                        y = 1;
                    }
					f=0;
					for(z=0;z<nsam;z++) /*all mutations should be observed. equal pattern in different position*/
						if(list[z][pbuf[i]] == list[z][pbuf[x]]) f++;
					if(f==nsam) {
						a = (double)ran1()*dlen + (double)wstartm1;
						pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);  
						x = i;
						y = 1;
						break;
					}
                }
                x--;
            }
		}
    }
    else { /* per no mhits */
		pbuf[0] = (long int)sel_nt;
        for(i=1;i<(long int)n;i++) {
			valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
			pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

void order2(long int n, double **pbuf)/* ordena els valors */
{
    long int gap,i,j;
    double temp0;
    double temp1;

    for(gap= n/2; gap>0;gap /= 2)
        for(i=gap;i<(long int)n;i++)
            for(j=i-gap;j>=0 && pbuf[j][0]>pbuf[j+gap][0];j -= gap) {
                temp0 = pbuf[j][0];
                temp1 = pbuf[j][1];
                pbuf[j][0] = pbuf[j+gap][0];
                pbuf[j][1] = pbuf[j+gap][1];
                pbuf[j+gap][0] = temp0;
                pbuf[j+gap][1] = temp1;
            }
}

int print_neuttest(struct var **data,FILE *output,char *file_out,char *file_in,double **matrix_test,struct prob_par **postp)
{
    long int x=0;
    long int y;
    long int a,b,c,d;
    double sum,sum2;
    int numloc;
    int totalnloci=1;
    FILE *outputind[NEUTVALUES2];
    char *el;
    
    /*observed values*/
    FILE *filePvalue=0;
    FILE *filePvalueDnaSP=0;
    char filenameDnaSP[420]="\0";
    double **matrix_avg=0;
    double **matrix_var=0;
    double **matrix_Pvalues=0;
    double **matrix_Pvaluesequal=0;
    double ***percentiles=0;
    int n=0;
    int observed=1;
    int pv=0;
    int compare_(const void *,const void *);
    long int count,countequal,totalobs;
    long int total=0;
    double fcount=0.;
    double *average=0; double *variance=0;
    long int **validiter=0;
    FILE *duplik;
    char file_out1[420];
    char DnaSPtxt[20]="_DnaSP.txt\0";
    double thetaptot,recptot;

    double calc_quantile(double *, long int, double, long int);
                         
    if((*data)->print_neuttest == 3) return 0;
    
    c = NEUTVALUES2;
    if((*data)->linked == 1 && (*data)->despl > 0 && (*data)->window > 0) {
        totalnloci = (int)ceil(((double)(*data)->nsites[1] - (double)(*data)->window) / ((double)(*data)->despl)) + (int)1;
    }
    else {
        if((*data)->linked > 1) totalnloci = (*data)->linked;
        else totalnloci = (*data)->n_loci;
    }
    
    /*observed values*/
    if(observed == 1) {
        if(!(average = (double *)malloc((unsigned)(c*(*data)->max_npop_sampled*sizeof(double))))) {
            perror("calloc error ms.matrix_avg*");
            exit(1);
        }
        if(!(variance = (double *)malloc((unsigned)(c*(*data)->max_npop_sampled*sizeof(double))))) {
            perror("calloc error ms.matrix_var*");
            exit(1);
        }
        if(!(validiter = (long int **)malloc((unsigned)(c*(*data)->max_npop_sampled*sizeof(long int *))))) {
            perror("calloc error ms.matrix_avg*");
            exit(1);
        }
        if(!(matrix_avg = (double **)malloc((unsigned)(c*(*data)->max_npop_sampled*sizeof(double *))))) {
            perror("calloc error ms.matrix_avg*");
            exit(1);
        }
        if(!(matrix_var = (double **)malloc((unsigned)(c*(*data)->max_npop_sampled*sizeof(double *))))) {
            perror("calloc error ms.matrix_var*");
            exit(1);
        }
        if(!(matrix_Pvalues = (double **)malloc((unsigned)(c*(*data)->max_npop_sampled*sizeof(double *))))) {
            perror("calloc error ms.matrix_var*");
            exit(1);
        }
        if(!(matrix_Pvaluesequal = (double **)malloc((unsigned)(c*(*data)->max_npop_sampled*sizeof(double *))))) {
            perror("calloc error ms.matrix_var*");
            exit(1);
        }
        if(!(percentiles = (double ***)malloc((unsigned)(c*(*data)->max_npop_sampled*sizeof(double **))))) {
            perror("calloc error ms.matrix_var*");
            exit(1);
        }
        for(a=0;a<c*(*data)->max_npop_sampled;a++) {
            if(!(validiter[a] = (long int *)malloc((unsigned)((totalnloci+2)*sizeof(long int))))) {
                perror("calloc error ms.matrix_avg*");
                exit(1);
            }
            if(!(matrix_avg[a] = (double *)malloc((unsigned)((*data)->n_iter*sizeof(double))))) {
                perror("calloc error ms.matrix_avg*");
                exit(1);
            }
            if(!(matrix_var[a] = (double *)malloc((unsigned)((*data)->n_iter*sizeof(double))))) {
                perror("calloc error ms.matrix_var*");
                exit(1);
            }
            if(!(matrix_Pvalues[a] = (double *)malloc((unsigned)((totalnloci+2)*sizeof(double))))) {
                perror("calloc error ms.matrix_var*");
                exit(1);
            }
            if(!(matrix_Pvaluesequal[a] = (double *)malloc((unsigned)((totalnloci+3)*sizeof(double))))) {
                perror("calloc error ms.matrix_var*");
                exit(1);
            }
            if(!(percentiles[a] = (double **)malloc((unsigned)((totalnloci+2)*sizeof(double *))))) {
                perror("calloc error ms.matrix_var*");
                exit(1);
            }
            for(b=0;b<totalnloci+2;b++) {
                if(!(percentiles[a][b] = (double *)malloc((unsigned)(13*sizeof(double))))) {
                    perror("calloc error ms.matrix_var*");
                    exit(1);
                }
            }
        }
    }
    
    if((*data)->likelihood_line == 0) {
        if((*data)->print_neuttest == 2 || (*data)->print_neuttest == 4) {
            for(a=0;a<c;a++) {
                strcat(dfiles[a],file_out);
                el = strrchr(dfiles[a],'.');
                *el = '\0';
                strcat(dfiles[a],ndfiles[a]);
                if(!(outputind[a] = fopen(dfiles[a],"w"))) return 1;
                for(b=0;b<totalnloci;b++) for(d=0;d<(*data)->max_npop_sampled;d++) fprintf(outputind[a],"%s[locus=%d,pop=%d]\t",names[(int)a],(int)b,(int)d);
                fputs("\n",outputind[a]);
            }
        }
        
        for(x=0;x<(*data)->n_iter;x++) {
            for(a=0;a<c;a++) {
                for(b=0;b<(long int)totalnloci;b++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if((*data)->print_neuttest == 2 || (*data)->print_neuttest == 4) {
                            if(matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] == -10000) fputs("na\t",outputind[a]);
                            else fprintf(outputind[a],"%.6g\t",matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x]);
                        }
                    }
                }
            }
            if((*data)->print_neuttest == 2 || (*data)->print_neuttest == 4)
                for(a=0;a<c;a++) fputs("\n",outputind[a]);
            
            for(a=0;a<c;a++) {
                for(d=0;d<(*data)->max_npop_sampled;d++) {
                    sum = sum2 = (double)0;
                    numloc = 0;
                    for(b=0;b<(long int)totalnloci;b++) {
                        if(matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] != -10000) {
                            sum  += matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x];
                            sum2 += matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] * matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x];
                            numloc++;
                        }
                    }
                    if(numloc > 1) {
                        sum  = sum/(double)numloc;
                        sum2 = (sum2/((double)numloc) - sum*sum)*((double)numloc/((double)numloc-(double)1));
                        if(sum2 <= 1E-37) sum2 = (double)0;
                        fprintf(output,"%.6g\t%.6g\t",sum,sum2);
                        
                        /*observed values*/
                        if(observed == 1) {
                            matrix_avg[a*(*data)->max_npop_sampled+d][x] = sum;
                            matrix_var[a*(*data)->max_npop_sampled+d][x] = sum2;
                        }
                    }
                    else {
                        if(numloc > 0) {
                            sum  = sum/(double)numloc;
                            fprintf(output,"%.6g\t",(double)sum);
                            if(totalnloci > 2) fprintf(output,"na\t");
                            
                            /*observed values*/
                            if(observed == 1) {
                                matrix_avg[a*(*data)->max_npop_sampled+d][x] = sum;
                                matrix_var[a*(*data)->max_npop_sampled+d][x] = (double)-10000;
                            }
                        }
                        else {
                            fputs("na\t",output);
                            if(totalnloci > 2) fputs("na\t",output);
                            
                            /*observed values*/
                            if(observed == 1) {
                                matrix_avg[a*(*data)->max_npop_sampled+d][x] = (double)-10000;
                                matrix_var[a*(*data)->max_npop_sampled+d][x] = (double)-10000;
                            }
                        }
                    }
                }
            }
            fputs("\n",output);
        }
        if((*data)->print_neuttest == 2 || (*data)->print_neuttest == 4) {
            for(a=0;a<c;a++) fclose(outputind[a]);
        }
    }
    
    
    
    /*observed values and percentiles*/
    if(observed == 1) {
        if((*data)->likelihood_line == 0) {
#ifndef __DNASP__
            printf("\n Calculating percentiles and (if required) probabilities for observed values... ");
#endif
            /*sort and calculate probabilities*/
            for(a=0;a<c;a++) {
                for(d=0;d<(*data)->max_npop_sampled;d++) {
                    for(b=0;b<(long int)totalnloci;b++) {
                        qsort(matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a],(long int)(*data)->n_iter,sizeof(double),compare_);
                    }
                    qsort(matrix_avg[a*(*data)->max_npop_sampled+d],(long int)(*data)->n_iter,sizeof(double),compare_);
                    qsort(matrix_var[a*(*data)->max_npop_sampled+d],(long int)(*data)->n_iter,sizeof(double),compare_);
                    /*calculate percentiles and probabilities*/
                    for(b=0;b<(long int)totalnloci;b++) {
                        count = countequal = total = totalobs = (long int)0;
                        for(x=0;x<(*data)->n_iter;x++) {
                            if(matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] != -10000) {
                                total++;
                                if((*data)->obs_statistics_used[a] == 1 && (*data)->obs_statistics[a][b][d+1] != (double)-10000) {
                                    totalobs++;
                                    if((*data)->obs_statistics_used[a] == 1) {
                                        if((*data)->obs_statistics[a][b][d+1] >= matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] + 1e-05 ||
                                           (*data)->obs_statistics[a][b][d+1] >= matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] - 1e-05)
                                            count++;
                                    }
                                    if((*data)->obs_statistics_used[a] == 1) {
                                        if((*data)->obs_statistics[a][b][d+1] <= matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] + 1e-05 &&
                                           (*data)->obs_statistics[a][b][d+1] >= matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] - 1e-05)
                                            countequal++;
                                    }
                                }
                            }
                        }
                        /*probabilitites*/
                        if(totalobs) {
                            matrix_Pvalues[a*(*data)->max_npop_sampled+d][b] = (double)count/(double)totalobs;
                            matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b] = (double)countequal/(double)totalobs;
                        }
                        else {
                            matrix_Pvalues[a*(*data)->max_npop_sampled+d][b] = (double)-10000;
                            matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b] = (double)-10000;
                        }
                        if((*data)->obs_statistics_used[a] == 1 && (*data)->obs_statistics[a][b][d+1] != (double)-10000)
                            validiter[a*(*data)->max_npop_sampled+d][b] = (long int)totalobs;
                        else
                            validiter[a*(*data)->max_npop_sampled+d][b] = (long int)total;
                        /*percentiles*/
                        for(n=0;n<13;n++) {
                            switch(n) {
                                case 0:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * ((double)0.001-(double)0.001/2.0);
                                    break;
                                case 1:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * ((double)0.010-(double)0.010/2.0);
                                    break;
                                case 2:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.025;
                                    break;
                                case 3:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.050;
                                    break;
                                case 4:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.100;
                                    break;
                                case 5:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.250;
                                    break;
                                case 6:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.500;
                                    break;
                                case 7:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.750;
                                    break;
                                case 8:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.900;
                                    break;
                                case 9:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.950;
                                    break;
                                case 10:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * (double)0.975;
                                    break;
                                case 11:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * ((double)0.990+(double)0.010/2.0);
                                    break;
                                case 12:
                                    fcount = (double)validiter[a*(*data)->max_npop_sampled+d][b] * ((double)0.999+(double)0.001/2.0);
                                    break;
                            }
                            if(fcount < (double)1 || fcount > (double)(total-1)) percentiles[a*(*data)->max_npop_sampled+d][b][n] = -10000;
                            else {
                                count=(long int)floor((double)fcount);
                                percentiles[a*(*data)->max_npop_sampled+d][b][n] = calc_quantile(matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a],validiter[a*(*data)->max_npop_sampled+d][b],fcount/(double)validiter[a*(*data)->max_npop_sampled+d][b],(*data)->n_iter-total);
                                
                            }
                        }
                    }
                    /*average*/
                    count = countequal = total = totalobs = (long int)0;
                    numloc = 0;
                    average[a*(*data)->max_npop_sampled+d] = (double)0;
                    for(b=0;b<(long int)totalnloci;b++) {
                        if((*data)->obs_statistics[a][b][d+1] != (double)-10000) {
                            average[a*(*data)->max_npop_sampled+d] += (*data)->obs_statistics[a][b][d+1];
                            numloc++;
                        }
                    }
                    if(numloc) average[a*(*data)->max_npop_sampled+d] /= (double)numloc;
                    else average[a*(*data)->max_npop_sampled+d] = (double)-10000;
                    for(x=0;x<(*data)->n_iter;x++) {
                        if(matrix_avg[a*(*data)->max_npop_sampled+d][x] != -10000) {
                            total++;
                            if((*data)->obs_statistics_used[a] == 1 && average[a*(*data)->max_npop_sampled+d] != (double)-10000) {
                                totalobs++;
                                if((*data)->obs_statistics_used[a] == 1) {
                                    if(average[a*(*data)->max_npop_sampled+d] >= matrix_avg[a*(*data)->max_npop_sampled+d][x] + 1e-05 ||
                                       average[a*(*data)->max_npop_sampled+d] >= matrix_avg[a*(*data)->max_npop_sampled+d][x] - 1e-05)
                                        count++;
                                }
                                if((*data)->obs_statistics_used[a] == 1) {
                                    if(average[a*(*data)->max_npop_sampled+d] <= matrix_avg[a*(*data)->max_npop_sampled+d][x] + 1e-05 &&
                                       average[a*(*data)->max_npop_sampled+d] >= matrix_avg[a*(*data)->max_npop_sampled+d][x] - 1e-05)
                                        countequal++;
                                }
                            }
                        }
                    }
                    /*probabilitites*/
                    if(total) {
                        matrix_Pvalues[a*(*data)->max_npop_sampled+d][totalnloci] = (double)count/(double)total;
                        matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci] = (double)countequal/(double)total;
                    }
                    else {
                        matrix_Pvalues[a*(*data)->max_npop_sampled+d][totalnloci] = (double)-10000;
                        matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci] = (double)-10000;
                    }
                    if((*data)->obs_statistics_used[a] == 1) validiter[a*(*data)->max_npop_sampled+d][totalnloci] = (long int)totalobs;
                    else validiter[a*(*data)->max_npop_sampled+d][totalnloci] = (long int)total;
                    /*percentiles*/
                    for(n=0;n<13;n++) {
                        switch(n) {
                            case 0:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * ((double)0.001-(double)0.001/2.0);
                                break;
                            case 1:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * ((double)0.010-(double)0.010/2.0);
                                break;
                            case 2:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.025;
                                break;
                            case 3:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.050;
                                break;
                            case 4:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.100;
                                break;
                            case 5:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.250;
                                break;
                            case 6:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.500;
                                break;
                            case 7:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.750;
                                break;
                            case 8:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.900;
                                break;
                            case 9:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.950;
                                break;
                            case 10:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * (double)0.975;
                                break;
                            case 11:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * ((double)0.990+(double)0.010/2.0);
                                break;
                            case 12:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci] * ((double)0.999+(double)0.001/2.0);
                                break;
                        }
                        if(fcount < (double)1 || fcount > (double)(total-1)) percentiles[a*(*data)->max_npop_sampled+d][totalnloci][n] = -10000;
                        else {
                            count=(long int)floor((double)fcount);
                            percentiles[a*(*data)->max_npop_sampled+d][totalnloci][n] = calc_quantile(matrix_avg[a*(*data)->max_npop_sampled+d],validiter[a*(*data)->max_npop_sampled+d][totalnloci],fcount/(double)validiter[a*(*data)->max_npop_sampled+d][totalnloci],(*data)->n_iter-total);
                        }
                    }
                    /*variance*/
                    count = countequal = total = totalobs = (long int)0;
                    numloc = 0;
                    average[a*(*data)->max_npop_sampled+d] = variance[a*(*data)->max_npop_sampled+d] = (double)0;
                    for(b=0;b<(long int)totalnloci;b++) {
                        if((*data)->obs_statistics[a][b][d+1] != (double)-10000) {
                            average[a*(*data)->max_npop_sampled+d] += (*data)->obs_statistics[a][b][d+1];
                            variance[a*(*data)->max_npop_sampled+d] += (*data)->obs_statistics[a][b][d+1] * (*data)->obs_statistics[a][b][d+1];
                            numloc++;
                        }
                    }
                    if(numloc>1) {
                        average[a*(*data)->max_npop_sampled+d] /= (double)numloc;
                        variance[a*(*data)->max_npop_sampled+d] = (variance[a*(*data)->max_npop_sampled+d]/((double)numloc) - average[a*(*data)->max_npop_sampled+d]*average[a*(*data)->max_npop_sampled+d])*((double)numloc/((double)numloc-(double)1));
                        if(variance[a*(*data)->max_npop_sampled+d] <= 1E-37) variance[a*(*data)->max_npop_sampled+d] = (double)0;
                    }
                    else {
                        variance[a*(*data)->max_npop_sampled+d] = (double)-10000;
                        if(numloc==0) average[a*(*data)->max_npop_sampled+d] = (double)-10000;
                        else average[a*(*data)->max_npop_sampled+d] /= (double)numloc;
                    }
                    for(x=0;x<(*data)->n_iter;x++) {
                        if(matrix_var[a*(*data)->max_npop_sampled+d][x] != -10000) {
                            total++;
                            if((*data)->obs_statistics_used[a] == 1 && variance[a*(*data)->max_npop_sampled+d] != (double)-10000) {
                                totalobs++;
                                if((*data)->obs_statistics_used[a] == 1) {
                                    if(variance[a*(*data)->max_npop_sampled+d] >= matrix_var[a*(*data)->max_npop_sampled+d][x] + 1e-05 ||
                                       variance[a*(*data)->max_npop_sampled+d] >= matrix_var[a*(*data)->max_npop_sampled+d][x] - 1e-05)
                                        count++;
                                }
                                if((*data)->obs_statistics_used[a] == 1) {
                                    if(variance[a*(*data)->max_npop_sampled+d] <= matrix_var[a*(*data)->max_npop_sampled+d][x] + 1e-05 &&
                                       variance[a*(*data)->max_npop_sampled+d] >= matrix_var[a*(*data)->max_npop_sampled+d][x] - 1e-05)
                                        countequal++;
                                }
                            }
                        }
                    }
                    /*probabilitites*/
                    if(total) {
                        matrix_Pvalues[a*(*data)->max_npop_sampled+d][totalnloci+1] = (double)count/(double)total;
                        matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci+1] = (double)countequal/(double)total;
                    }
                    else {
                        matrix_Pvalues[a*(*data)->max_npop_sampled+d][totalnloci+1] = (double)-10000;
                        matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci+1] = (double)-10000;
                    }
                    if((*data)->obs_statistics_used[a] == 1) validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] = (long int)totalobs;
                    else validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] = (long int)total;
                    /*percentiles*/
                    for(n=0;n<13;n++) {
                        switch(n) {
                            case 0:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * ((double)0.001-(double)0.001/2.0);
                                break;
                            case 1:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * ((double)0.010-(double)0.010/2.0);
                                break;
                            case 2:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.025;
                                break;
                            case 3:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.050;
                                break;
                            case 4:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.100;
                                break;
                            case 5:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.250;
                                break;
                            case 6:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.500;
                                break;
                            case 7:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.750;
                                break;
                            case 8:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.900;
                                break;
                            case 9:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.950;
                                break;
                            case 10:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * (double)0.975;
                                break;
                            case 11:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * ((double)0.990+(double)0.010/2.0);
                                break;
                            case 12:
                                fcount = (double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1] * ((double)0.999+(double)0.001/2.0);
                                break;
                        }
                        if(fcount < (double)1 || fcount > (double)(total-1)) percentiles[a*(*data)->max_npop_sampled+d][totalnloci+1][n] = -10000;
                        else {
                            count=(long int)floor((double)fcount);
                            percentiles[a*(*data)->max_npop_sampled+d][totalnloci+1][n] = calc_quantile(matrix_var[a*(*data)->max_npop_sampled+d],validiter[a*(*data)->max_npop_sampled+d][totalnloci+1],fcount/(double)validiter[a*(*data)->max_npop_sampled+d][totalnloci+1],(*data)->n_iter-total);
                        }
                    }
                }
            }
        }
        else {
            if((*data)->n_iter > 1) {
                printf("\n Calculating likelihoods (2*log(P)) for observed values... ");
                fflush(stdout);
                /*create a new file*/
                file_out1[0] = '\0';
                strncat(file_out1,file_out,420);
                el = strrchr(file_out1,'.');
                *el = '\0';
                strncat(file_out1,"_1.out\0",420);
                if(!(duplik = fopen(file_out1,"w"))) return 1;
                /*init*/
                matrix_Pvaluesequal[0][totalnloci+1] = 0.;
                /*loop*/
                for(a=0;a<c;a++) {
                    if((*data)->obs_statistics_used[a] == 1) {
                        for(d=0;d<(*data)->max_npop_sampled;d++) {
                            matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci] = (double)0;
                            for(b=0;b<(long int)totalnloci;b++) {
                                countequal = total = (long int)0;
                                for(x=0;x<(*data)->n_iter;x++) {
                                    if(matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] != -10000  && (*data)->obs_statistics[a][b][d+1] != (double)-10000) {
                                        total++;
                                        if((*data)->obs_statistics[a][b][d+1] - (*data)->likelihood_error[a] <= matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x]  &&
                                           (*data)->obs_statistics[a][b][d+1] + (*data)->likelihood_error[a] >= matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x])
                                            countequal++;
                                    }
                                }
                                /*probabilitites for each locus, each stistic*/
                                if(total) {
                                    if(countequal) matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b] = (double)2.*(double)log((double)countequal/(double)total);
                                    else matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b] = (double)2.*(double)log((double)1/(double)((*data)->n_iter+1));
                                }
                                else {
                                    matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b] = (double)2.*(double)log((double)1/(double)((*data)->n_iter+1));
                                }
                                /*probabilitites for all loci, each pop, each statistic*/
                                matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci] += matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b];
                                validiter[a*(*data)->max_npop_sampled+d][b] = (long int)total;
                            }
                            /*probabilitites for all loci, all statistics, each pop*/
                            matrix_Pvaluesequal[d][totalnloci+1] += matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci];
                            validiter[a*(*data)->max_npop_sampled+d][totalnloci] = (long int)total;
                        }
                    }
                }
                matrix_Pvaluesequal[0][totalnloci+2] =(double)0;
                for(d=0;d<(*data)->max_npop_sampled;d++) {
                    /*probabilitites for all loci, all statistics, all pops*/
                    matrix_Pvaluesequal[0][totalnloci+2] += matrix_Pvaluesequal[d][totalnloci+1];
                }
                /*PRINT LIKELIHOODS IN A SINGLE LINE in the OUTPUT FILE and in another FILE but FOR A SIGLE LINE*/
                fprintf(output,"%.6g",matrix_Pvaluesequal[0][totalnloci+2]);
                fprintf(duplik,"%.6g",matrix_Pvaluesequal[0][totalnloci+2]);
                for(d=0;d<(*data)->max_npop_sampled;d++) {
                    fprintf(output,"\t%.6g",matrix_Pvaluesequal[d][totalnloci+1]);
                    fprintf(duplik,"\t%.6g",matrix_Pvaluesequal[d][totalnloci+1]);
                    for(a=0;a<c;a++) {
                        if((*data)->obs_statistics_used[a] == 1) {
                            fprintf(output,"\t%.6g",matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci]);
                            fprintf(duplik,"\t%.6g",matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci]);
                        }
                    }
                    for(a=0;a<c;a++) {
                        if((*data)->obs_statistics_used[a] == 1) {
                            for(b=0;b<(long int)totalnloci;b++) {
                                fprintf(output,"\t%.6g",matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b]);
                                fprintf(duplik,"\t%.6g",matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b]);
                            }
                        }
                    }
                }
                fprintf(output,"\n");
                fclose(duplik);
            }
            else {
                /*FOR NITER==1 PRINT THE OBSERVED VALUES FOR EACH LOCUS IN A SINGLE LINE*/
                for(a=0;a<c;a++) {
                    if((*data)->obs_statistics_used[a] == 1) {
                        for(d=0;d<(*data)->max_npop_sampled;d++) {
                            for(b=0;b<(long int)totalnloci;b++) {
                                if(matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] != -10000)
                                    fprintf(output,"%.6g\t",matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x]);
                                else
                                    fprintf(output,"na\t");
                            }
                        }
                    }
                }
                if((*data)->ifgamma == 1 || (*data)->range_thetant) {
                    if((*data)->sfix_allthetas > 0) {
                        thetaptot = (double)1;
                        for(x=0;x<(*data)->n_loci;x++) {
                            thetaptot *= postp[x][(*data)->n_iter].thetap;
                        }
                        fprintf(output,"%G\t",thetaptot);
                    }
                }
                if(((*data)->ifgammar == 1 || (*data)->range_rnt) && !((*data)->ifgamma == 1 || (*data)->range_thetant)){
                    if((*data)->rmfix > 0) {
                        recptot = (double)1;
                        for(x=0;x<(*data)->n_loci;x++) {
                            recptot *= postp[x][(*data)->n_iter].recp;
                        }
                        fprintf(output,"%G\t",recptot);
                    }
                }
                fprintf(output,"\n");
            }
        }
        
        if((*data)->likelihood_line == 0) {
            /*printing values in filePvalue*/
            strcat(dfiles[NEUTVALUES2],file_out);
            el = strrchr(dfiles[NEUTVALUES2],'.');
            *el = '\0';
            strcat(dfiles[NEUTVALUES2],ndfiles[NEUTVALUES2]);
            if(!(filePvalue = fopen(dfiles[NEUTVALUES2],"w"))) return 1;
            
            
            /*DnaSP file*/
            if((*data)->print_neuttest == 5) {
                strcat(filenameDnaSP,file_out);
                el = strrchr(filenameDnaSP,'.');
                *el = '\0';
                strcat(filenameDnaSP,DnaSPtxt);
                if(!(filePvalueDnaSP = fopen(filenameDnaSP,"w"))) return 1;
                /*include input parameters*/
                fprintf(filePvalueDnaSP,MLCOALSIM);
                fprintf(filePvalueDnaSP,"Input file: %s",file_in);
                fprintf(filePvalueDnaSP,"\n");
                
                fprintf(filePvalueDnaSP,"\nNumber of iterations per locus: %ld",(*data)->n_iter);
                fprintf(filePvalueDnaSP,"\nSeed used per locus:");
                for(b=0;b<(long int)totalnloci;b++) fprintf(filePvalueDnaSP, " %ld",(*data)->seed1[b+1]);
                fprintf(filePvalueDnaSP,"\n");
                
                fprintf(filePvalueDnaSP,"\nNumber of loci: %d",totalnloci);
                fprintf(filePvalueDnaSP,"\nNumber of samples per locus:");
                for(b=0;b<(long int)totalnloci;b++) fprintf(filePvalueDnaSP," %d",(*data)->nsam[b+1]);
                fprintf(filePvalueDnaSP,"\nNumber of sites per locus:");
                for(b=0;b<(long int)totalnloci;b++) fprintf(filePvalueDnaSP," %ld",(*data)->nsites[b+1]);
                
                fprintf(filePvalueDnaSP,"\n");
                if((*data)->theta_1[0]) {
                    fprintf(filePvalueDnaSP,"\nTheta (4Nu) per locus:");
                    for(b=0;b<(long int)totalnloci;b++) fprintf(filePvalueDnaSP," %.3f",(*data)->theta_1[b+1] * (*data)->nsites[b+1]);
                }
                else {
                    fprintf(filePvalueDnaSP,"\nFixed mutations per locus:");
                    for(b=0;b<(long int)totalnloci;b++) fprintf(filePvalueDnaSP," %ld",(*data)->mutations[b+1]);
                }
                fprintf(filePvalueDnaSP,"\nRecombination (4Nr) per locus:");
                for(b=0;b<(long int)totalnloci;b++)
                    if((*data)->r[b+1] != 1000000.0)
                        fprintf(filePvalueDnaSP," %f",(*data)->r[b+1] * (*data)->nsites[b+1]);
                    else
                        fprintf(filePvalueDnaSP," infinite");
                
                fprintf(filePvalueDnaSP,"\n");
                if((*data)->npop > 1) {
                    fprintf(filePvalueDnaSP,"\nNumber of populations: %ld\n",(*data)->npop);
                    for(b=0;b<(long int)totalnloci;b++) {
                        fprintf(filePvalueDnaSP,"\nSampled populations in locus %ld = %d:",b,(*data)->npop_sampled[b+1]+1);
                        for(d=1;d<(*data)->npop_sampled[b+1]+1;d++) {
                            fprintf(filePvalueDnaSP," %d",(*data)->ssize_pop[b][d]);
                        }
                    }
                    if((*data)->mig_rate) {
                        fprintf(filePvalueDnaSP,"\nMigration parameter (4Nm): %G",(*data)->mig_rate);
                    }
                    else {
                        if((*data)->mig_rate_matrix[0][0]) {
                            fprintf(filePvalueDnaSP,"\nMigration parameter matrix (4Nm):\n");
                            for(x=0;x<(*data)->npop;x++) {
                                if(x!=0) fprintf(filePvalueDnaSP,"\n");
                                for(y=0;y<(*data)->npop;y++) {
                                    fprintf(filePvalueDnaSP," %G",(*data)->mig_rate_matrix[x][y+1]);
                                }
                            }
                        }
                    }
                }
                if((*data)->npop_events) {
                    fprintf(filePvalueDnaSP,"Number of Demographic Events: %d\n",(int)(*data)->npop_events);
                    fprintf(filePvalueDnaSP,"\tEvent\tPop_I\tPop_II\tNe_I\tTime\tMigration_vectors\n");
                    for(x=0;x<(*data)->npop_events;x++) {
                        if(x!=0) fprintf(filePvalueDnaSP,"\n");
                        for(y=0;y<(*data)->pop_event[x][0];y++) {
                            if(y<3) fprintf(filePvalueDnaSP,"\t%d",(int)(*data)->pop_event[x][y+1]);
                            else {
                                if(y<6) fprintf(filePvalueDnaSP,"\t%G",(*data)->pop_event[x][y+1]);
                                else fprintf(filePvalueDnaSP," %G",(*data)->pop_event[x][y+1]);
                            }
                        }
                    }
                }
                fprintf(filePvalueDnaSP,"\n\n");
                fprintf(filePvalueDnaSP,"STATISTICS:\n");
            }
            /**/
            
            for(a=0;a<c;a++) {
                if((*data)->obs_statistics_used[a] == 1) {
                    pv = 1; /*observed values*/
                    break;
                }
            }
            
            /*DnaSP format*/
            if((*data)->print_neuttest == 5) {
                for(b=0;b<(long int)totalnloci;b++) {
                    fprintf(filePvalueDnaSP,"locus_%ld",b);
                    fprintf(filePvalueDnaSP,"\tNet_replicates");
                    fprintf(filePvalueDnaSP,"\tmean");
                    fprintf(filePvalueDnaSP,"\tmedian");/*6*/
                    fprintf(filePvalueDnaSP,"\t5.0%%");/*3*/
                    fprintf(filePvalueDnaSP,"\t95.0%%");/*9*/
                    fprintf(filePvalueDnaSP,"\t2.5%%");/*2*/
                    fprintf(filePvalueDnaSP,"\t97.5%%");/*10*/
                    fprintf(filePvalueDnaSP,"\t0.5%%");/*1*/
                    fprintf(filePvalueDnaSP,"\t99.5%%");/*11*/
                    fprintf(filePvalueDnaSP,"\t0.05%%");/*0*/
                    fprintf(filePvalueDnaSP,"\t99.95%%");/*12*/
                    fprintf(filePvalueDnaSP,"\tObs");
                    fprintf(filePvalueDnaSP,"\tP(Sim<=Obs)");
                    fprintf(filePvalueDnaSP,"\n");
                    for(a=0;a<c;a++) {
                        if(a!=22 && a!=23 && a!=29 && a!=30 && a!=31 && a!=32 && a!=33 && a!=34) {
                            for(d=0;d<(*data)->max_npop_sampled;d++) {
                                if(totalnloci==1)fprintf(filePvalueDnaSP,"%s",names[a]);
                                else fprintf(filePvalueDnaSP,"%s[pop=%d]",names[a],(int)d);
                                fprintf(filePvalueDnaSP,"\t%ld",validiter[a*(*data)->max_npop_sampled+d][b]);
                                /*mean*/
                                sum = (double)0;
                                total = (long int)0;
                                for(x=0;x<(*data)->n_iter;x++) {
                                    if(matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x] != -10000) {
                                        sum += matrix_test[(b*c*(*data)->max_npop_sampled)+(d*c)+a][x];
                                        total++;
                                    }
                                }
                                if(total > 0)
                                    fprintf(filePvalueDnaSP,"\t%.6g",sum/(double)total);
                                else
                                    fprintf(filePvalueDnaSP,"\tNA");
                                /*probabilitites for all loci, each pop, each statistic*/
                                matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci] += matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b];
                                validiter[a*(*data)->max_npop_sampled+d][b] = (long int)total;
                                
                                /*probabilitites for all loci, all statistics, each pop*/
                                matrix_Pvaluesequal[d][totalnloci+1] += matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci];
                                validiter[a*(*data)->max_npop_sampled+d][totalnloci] = (long int)total;
                                
                                
                                for(y=0;y<9;y++) {
                                    switch (y) {
                                        case 0: n=6;break;
                                        case 1: n=3;break;
                                        case 2: n=9;break;
                                        case 3: n=2;break;
                                        case 4: n=10;break;
                                        case 5: n=1;break;
                                        case 6: n=11;break;
                                        case 7: n=0;break;
                                        case 8: n=12;break;
                                    }
                                    if(percentiles[a*(*data)->max_npop_sampled+d][b][n] != (double)-10000) fprintf(filePvalueDnaSP,"\t%.6g",percentiles[a*(*data)->max_npop_sampled+d][b][n]);
                                    else fprintf(filePvalueDnaSP,"\tNA");
                                }
                                if(pv==1) {
                                    if((*data)->obs_statistics_used[a] == 1) {
                                        if((double)(*data)->obs_statistics[a][b][d+1] == (double)-10000) 
                                            fprintf(filePvalueDnaSP,"\tNA");
                                        else 
                                            fprintf(filePvalueDnaSP,"\t%.6g",(*data)->obs_statistics[a][b][d+1]);
                                        if(matrix_Pvalues[a*(*data)->max_npop_sampled+d][b] != (double)-10000) 
                                            fprintf(filePvalueDnaSP,"\t%.6g",matrix_Pvalues[a*(*data)->max_npop_sampled+d][b]);
                                        else 
                                            fprintf(filePvalueDnaSP,"\tNA");
                                    }
                                    else {
                                        fprintf(filePvalueDnaSP,"\t-");
                                        fprintf(filePvalueDnaSP,"\tNA");
                                    }
                                }
                                fprintf(filePvalueDnaSP,"\n");
                            }
                        }
                    }
                    fprintf(filePvalueDnaSP,"\n");
                }
            }
            /*END DnaSP format*/
            
            if(pv==1) {
                fprintf(filePvalue,"OBSERVED VALUES\n");
                for(a=0;a<c;a++) for(d=0;d<(*data)->max_npop_sampled;d++) fprintf(filePvalue,"\t%s[pop=%d]",names[a],(int)d);
                for(b=0;b<(long int)totalnloci;b++) {
                    fprintf(filePvalue,"\nlocus_%ld",b);
                    for(a=0;a<c;a++) {
                        for(d=0;d<(*data)->max_npop_sampled;d++) {
                            if((*data)->obs_statistics_used[a] == 1) {
                                if((double)(*data)->obs_statistics[a][b][d+1] == (double)-10000) fprintf(filePvalue,"\tna");
                                else fprintf(filePvalue,"\t%.6g",(*data)->obs_statistics[a][b][d+1]);
                            }
                            else fprintf(filePvalue,"\tx");
                        }
                    }
                }
                fprintf(filePvalue,"\naverage");
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if((*data)->obs_statistics_used[a] == 1) {
                            if(average[a*(*data)->max_npop_sampled+d] != (double)-10000) fprintf(filePvalue,"\t%.6g",average[a*(*data)->max_npop_sampled+d]);
                            else fprintf(filePvalue,"\tna");
                        }
                        else fprintf(filePvalue,"\tx");
                    }
                }
                fprintf(filePvalue,"\nvariance");
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if((*data)->obs_statistics_used[a] == 1) {
                            if(variance[a*(*data)->max_npop_sampled+d] != (double)-10000) fprintf(filePvalue,"\t%.6g",variance[a*(*data)->max_npop_sampled+d]);
                            else fprintf(filePvalue,"\tna");
                        }
                        else fprintf(filePvalue,"\tx");
                    }
                }
                fprintf(filePvalue,"\n\nPROBABILITY THAT SIMULATED VALUES BE SMALLER OR EQUAL THAN THE OBSERVED VALUE. P(Sim <= Obs): ACCURACY OF +- 1e-6\n");
                for(a=0;a<c;a++) for(d=0;d<(*data)->max_npop_sampled;d++) fprintf(filePvalue,"\t%s[pop=%d]",names[a],(int)d);
                
                for(b=0;b<(long int)totalnloci;b++) {
                    fprintf(filePvalue,"\nlocus_%ld",b);
                    for(a=0;a<c;a++) {
                        for(d=0;d<(*data)->max_npop_sampled;d++) {
                            if((*data)->obs_statistics_used[a] == 1) {
                                if(matrix_Pvalues[a*(*data)->max_npop_sampled+d][b] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvalues[a*(*data)->max_npop_sampled+d][b]);
                                else fprintf(filePvalue,"\tna");
                            }
                            else fprintf(filePvalue,"\tna");
                        }
                    }
                }
                fprintf(filePvalue,"\naverage");
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if((*data)->obs_statistics_used[a] == 1) {
                            if(matrix_Pvalues[a*(*data)->max_npop_sampled+d][totalnloci] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvalues[a*(*data)->max_npop_sampled+d][totalnloci]);
                            else fprintf(filePvalue,"\tna");
                        }
                        else fprintf(filePvalue,"\tna");
                    }
                }
                fprintf(filePvalue,"\nvariance");
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if((*data)->obs_statistics_used[a] == 1) {
                            if(matrix_Pvalues[a*(*data)->max_npop_sampled+d][totalnloci+1] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvalues[a*(*data)->max_npop_sampled+d][totalnloci+1]);
                            else fprintf(filePvalue,"\tna");
                        }
                        else fprintf(filePvalue,"\tna");
                    }
                }
                
                fprintf(filePvalue,"\n\nPROBABILITY THAT SIMULATED VALUES BE EQUAL THAN THE OBSERVED VALUE. P(Sim = Obs): ACCURACY OF +- 1e-6\n");
                for(a=0;a<c;a++) for(d=0;d<(*data)->max_npop_sampled;d++) fprintf(filePvalue,"\t%s[pop=%d]",names[a],(int)d);
                for(b=0;b<(long int)totalnloci;b++) {
                    fprintf(filePvalue,"\nlocus_%ld",b);
                    for(a=0;a<c;a++) {
                        for(d=0;d<(*data)->max_npop_sampled;d++) {
                            if((*data)->obs_statistics_used[a] == 1) {
                                if(matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][b]);
                                else fprintf(filePvalue,"\tna");
                            }
                            else fprintf(filePvalue,"\tna");
                        }
                    }
                }
                fprintf(filePvalue,"\naverage");
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if((*data)->obs_statistics_used[a] == 1) {
                            if(matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci]);
                            else fprintf(filePvalue,"\tna");
                        }
                        else fprintf(filePvalue,"\tna");
                    }
                }
                fprintf(filePvalue,"\nvariance");
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if((*data)->obs_statistics_used[a] == 1) {
                            if(matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci+1] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvaluesequal[a*(*data)->max_npop_sampled+d][totalnloci+1]);
                            else fprintf(filePvalue,"\tna");
                        }
                        else fprintf(filePvalue,"\tna");
                    }
                }
            }
            fprintf(filePvalue,"\n\nNumber of valid iterations:\n");
            for(a=0;a<c;a++) for(d=0;d<(*data)->max_npop_sampled;d++) fprintf(filePvalue,"\t%s[pop=%d]",names[a],(int)d);
            for(b=0;b<(long int)totalnloci;b++) {
                fprintf(filePvalue,"\nlocus_%ld",b);
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++)
                        fprintf(filePvalue,"\t%ld",validiter[a*(*data)->max_npop_sampled+d][b]);
                }
            }
            fprintf(filePvalue,"\naverage");
            for(a=0;a<c;a++) {
                for(d=0;d<(*data)->max_npop_sampled;d++)
                    fprintf(filePvalue,"\t%ld",validiter[a*(*data)->max_npop_sampled+d][totalnloci]);
            }
            fprintf(filePvalue,"\nvariance");
            for(a=0;a<c;a++) {
                for(d=0;d<(*data)->max_npop_sampled;d++)
                    fprintf(filePvalue,"\t%ld",validiter[a*(*data)->max_npop_sampled+d][totalnloci+1]);
            }
            
            fprintf(filePvalue,"\n\nPERCENTILES:");
            for(b=0;b<(long int)totalnloci;b++) {
                fprintf(filePvalue,"\n\nlocus_%ld\n\n",b);
                for(a=0;a<c;a++) for(d=0;d<(*data)->max_npop_sampled;d++) fprintf(filePvalue,"\t%s[pop=%d]",names[a],(int)d);
                for(n=0;n<13;n++) {
                    switch(n) {
                        case 0:
                            fprintf(filePvalue,"\n 0.05%%");/*fprintf(filePvalue,"\n 0.1%%");*/
                            break;
                        case 1:
                            fprintf(filePvalue,"\n 0.5%%");/*fprintf(filePvalue,"\n 1.0%%");*/
                            break;
                        case 2:
                            fprintf(filePvalue,"\n 2.5%%");
                            break;
                        case 3:
                            fprintf(filePvalue,"\n 5.0%%");
                            break;
                        case 4:
                            fprintf(filePvalue,"\n10.0%%");
                            break;
                        case 5:
                            fprintf(filePvalue,"\n25.0%%");
                            break;
                        case 6:
                            fprintf(filePvalue,"\n50.0%%");
                            break;
                        case 7:
                            fprintf(filePvalue,"\n75.0%%");
                            break;
                        case 8:
                            fprintf(filePvalue,"\n90.0%%");
                            break;
                        case 9:
                            fprintf(filePvalue,"\n95.0%%");
                            break;
                        case 10:
                            fprintf(filePvalue,"\n97.5%%");
                            break;
                        case 11:
                            fprintf(filePvalue,"\n99.5%%");
                            break;
                        case 12:
                            fprintf(filePvalue,"\n99.95%%");
                            break;
                    }
                    for(a=0;a<c;a++) {
                        for(d=0;d<(*data)->max_npop_sampled;d++) {
                            if(percentiles[a*(*data)->max_npop_sampled+d][b][n] != (double)-10000) fprintf(filePvalue,"\t%.6g",percentiles[a*(*data)->max_npop_sampled+d][b][n]);
                            else fprintf(filePvalue,"\tna");
                        }
                    }
                }
            }
            /*average and variance*/
            fprintf(filePvalue,"\n\naverage\n\n");
            for(a=0;a<c;a++) for(d=0;d<(*data)->max_npop_sampled;d++) fprintf(filePvalue,"\t%s[pop=%d]",names[a],(int)d);
            for(n=0;n<13;n++) {
                switch(n) {
                    case 0:
                        fprintf(filePvalue,"\n 0.05%%");
                        break;
                    case 1:
                        fprintf(filePvalue,"\n 0.5%%");
                        break;
                    case 2:
                        fprintf(filePvalue,"\n 2.5%%");
                        break;
                    case 3:
                        fprintf(filePvalue,"\n 5.0%%");
                        break;
                    case 4:
                        fprintf(filePvalue,"\n10.0%%");
                        break;
                    case 5:
                        fprintf(filePvalue,"\n25.0%%");
                        break;
                    case 6:
                        fprintf(filePvalue,"\n50.0%%");
                        break;
                    case 7:
                        fprintf(filePvalue,"\n75.0%%");
                        break;
                    case 8:
                        fprintf(filePvalue,"\n90.0%%");
                        break;
                    case 9:
                        fprintf(filePvalue,"\n95.0%%");
                        break;
                    case 10:
                        fprintf(filePvalue,"\n97.5%%");
                        break;
                    case 11:
                        fprintf(filePvalue,"\n99.5%%");
                        break;
                    case 12:
                        fprintf(filePvalue,"\n99.95%%");
                        break;
                }
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if(percentiles[a*(*data)->max_npop_sampled+d][totalnloci][n] != (double)-10000) fprintf(filePvalue,"\t%.6g",percentiles[a*(*data)->max_npop_sampled+d][totalnloci][n]);
                        else fprintf(filePvalue,"\tna");
                    }
                }
            }
            fprintf(filePvalue,"\n\nvariance\n\n");
            for(a=0;a<c;a++) for(d=0;d<(*data)->max_npop_sampled;d++) fprintf(filePvalue,"\t%s[pop=%d]",names[a],(int)d);
            for(n=0;n<13;n++) {
                switch(n) {
                    case 0:
                        fprintf(filePvalue,"\n 0.05%%");
                        break;
                    case 1:
                        fprintf(filePvalue,"\n 0.5%%");
                        break;
                    case 2:
                        fprintf(filePvalue,"\n 2.5%%");
                        break;
                    case 3:
                        fprintf(filePvalue,"\n 5.0%%");
                        break;
                    case 4:
                        fprintf(filePvalue,"\n10.0%%");
                        break;
                    case 5:
                        fprintf(filePvalue,"\n25.0%%");
                        break;
                    case 6:
                        fprintf(filePvalue,"\n50.0%%");
                        break;
                    case 7:
                        fprintf(filePvalue,"\n75.0%%");
                        break;
                    case 8:
                        fprintf(filePvalue,"\n90.0%%");
                        break;
                    case 9:
                        fprintf(filePvalue,"\n95.0%%");
                        break;
                    case 10:
                        fprintf(filePvalue,"\n97.5%%");
                        break;
                    case 11:
                        fprintf(filePvalue,"\n99.5%%");
                        break;
                    case 12:
                        fprintf(filePvalue,"\n99.95%%");
                        break;
                }
                for(a=0;a<c;a++) {
                    for(d=0;d<(*data)->max_npop_sampled;d++) {
                        if(percentiles[a*(*data)->max_npop_sampled+d][totalnloci+1][n] != (double)-10000) fprintf(filePvalue,"\t%.6g",percentiles[a*(*data)->max_npop_sampled+d][totalnloci+1][n]);
                        else fprintf(filePvalue,"\tna");
                    }
                }
            }
        }
        
        for(a=0;a<c*(*data)->max_npop_sampled;a++) {
            free(validiter[a]);
            free(matrix_avg[a]);
            free(matrix_var[a]);
            free(matrix_Pvalues[a]);
            free(matrix_Pvaluesequal[a]);
            for(n=0;n<(long int)totalnloci+2;n++) 
                free(percentiles[a][n]);
            free(percentiles[a]);
        }
        free(validiter);
        free(average);
        free(variance);
        free(matrix_avg);
        free(matrix_var);
        free(matrix_Pvalues);
        free(matrix_Pvaluesequal);
        free(percentiles);
        if((*data)->likelihood_line == 0) fclose(filePvalue);
        if((*data)->likelihood_line == 0 && (*data)->print_neuttest == 5) fclose(filePvalueDnaSP);
        
    }	
    return 0;
}

double calc_quantile(double *vec, long int lenvec,double quant, long int numit) {
    double res;
    long int pos1,pos2;
    pos1 = (double)0.0  + (double)floor((double)quant * lenvec);/*right position*/
    pos2 = (lenvec-1.0) - (double)floor(((1.0-quant) * lenvec));/*left position*/
    res = (vec[numit + pos1] + vec[numit + pos2])/2.0;
    return(res);
}

void calc_neutparSRH(int valuep,long int segsit,struct var2 **inputp, struct dnapar *ntpar,double valuer,int npopa,int flaghap)
{
    long int S;
    int nhapl;
    int *haplotype = 0;
    int inits;
	int nsam=0;
    long int j;
    int a,b,h,x;
	long int comb;
    char *hapl=0;
	int Min_rec(int,int,int,int,int);
	
	if(npopa == (*inputp)->npop) {
		inits = 0;
		nsam = (*inputp)->nsam;
	}
	else {
		inits = 0;
		for(x=0;x<(*inputp)->npop_sampled;x++) {
			if(x < npopa) inits += (*inputp)->config[x];
			else {
				nsam = (*inputp)->config[x];
				break;
			}
		}
	}

    comb = (int)((double)nsam*((double)nsam-(double)1)/(double)2);
	for(x=0;x<10;x++) (ntpar)->Sanc[x] = 0;

	if(segsit == 0 || nsam < 2) {
		if(nsam == 0) {
			(ntpar)->S = -10000;
			(ntpar)->nhapl = -10000;
			(ntpar)->Rm = (int)-10000;
		}
		else {
			(ntpar)->S = 0;
			(ntpar)->nhapl = 1;
			(ntpar)->Rm = (int)0;
		}
	}
	else {		
		if(valuer) (ntpar)->Rm = Min_rec(0,(int)segsit,nsam,inits,(*inputp)->nsam);
		else (ntpar)->Rm = (int)0;
	}

	if(flaghap > 0) {
		/* calcul de S,nhapl (excluding mhits)*/
		if((hapl = (char *)calloc(nsam*segsit,sizeof(char))) == NULL) {
			perror("calloc error calc_neutpar.0");
			exit(1);
		}
	}

	S = 0;
	nhapl = 0;                
	for(j=0;j<segsit;j++) {
		while(j < segsit) {
			if((h=ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit)) > 0) break; /*h is the frequency of the new mutation*/
			else j++;
		}            
		if(j<segsit) {
			if(flaghap > 0) {
				for(a=inits;a<inits+nsam;a++) {
					hapl[(a-inits)*segsit+S] = list[a][j];
				}
			}
			S++;
		}
	}
	(ntpar)->S = S;
	
	if(flaghap > 0) {
		if((haplotype = (int *)malloc(nsam*sizeof(int))) == NULL) {
			perror("calloc error calc_neutpar.1");
			exit(1);
		}
		(ntpar)->nhapl = 0;
		for(a=0;a<nsam;a++) haplotype[a] = 1;            
		for(a=0;a<nsam-1;a++) {
			if(haplotype[a]) {
				(ntpar)->nhapl += 1;
				for(b=a+1;b<nsam;b++) {
					if(haplotype[b]) {
						if(memcmp(hapl+a*segsit,hapl+b*segsit,S) == 0) { 
							haplotype[a] += 1;
							haplotype[b] = 0;
						}
					}
				}
			}
		}
		if(haplotype[a]) (ntpar)->nhapl += 1;
		free(hapl);
		free(haplotype);
	}
	
	return;
}

void calc_neutpar(int valuep,long int segsit,struct var2 **inputp, struct dnapar *ntpar,double valuer,int npopa)
{
    long int pi;
    double k_,k_e1,k_n1;
    long int S,Se1,Sn1;
    int nhapl,maxhapl;
    int B;
    int A;
    int *haplotype = 0;
    long int *piw = 0;
    long int *pib = 0;
    long int *hapw = 0;
    long int *hapb = 0;
    long int segsitesm1;
	int /*pidcount,*/npw;
	
	#if iHStest == 1
	int *fhap1,*fhap2,*nsam1,*nsam2;
	int n1,n2,flagA,flagD;
	double *EHSA,*EHSD;
	double iHHA,iHHD,uiHS;
	int *nousedA,*nousedD;
	#endif
	
	#if iEStest == 1
	int *fhap_ij;
	double **EHHS,*iES;
	int **nousedM;
	int flagM;
	#endif
	
	#if iEtest == 1 || iHStest == 1
	long int l;
	#endif

    int inits,inits1,inits2,initso/*,initcum*/;
	int nsam=0;
	int *initsq1,*initsq2,*sq;
    int val10,val20,val21;
    long int j,k;
    int a,b,c,d,i,comb2,x; int h=0;
	long int comb;
    char *hapl;
	int **veca;
	
	#if MAXHAP1
	int mh1,mut,maxhapl1,*sshg;
	#endif

	int Min_rec(int,int,int,int,int);
	int pola,polb,polc,pop1,pop2,popo,nt0;
	int *polqa,*polqb;
	int nmh;
	int nsamallpop;
	
	double *freq1,*freq2,*freq3,*freq4,s,s2,g1,g2,moment;
	long int z;
	
	long int *Pwd;
	long int maxpwd;
		
	long int ehh0,ehhs,ehhr;

    if(valuep == 0) {
		if(npopa == (*inputp)->npop) {
			inits = 0;
			nsam = (*inputp)->nsam;
		}
		else {
			inits = 0;
			for(x=0;x<(*inputp)->npop_sampled;x++) {
				if(x < npopa) inits += (*inputp)->config[x];
				else {
					nsam = (*inputp)->config[x];
					break;
				}
			}
		}
	}
    else {
        if(valuep == 1) {
            inits = 0;
            nsam  = (*inputp)->config[0];
        }
        else {
            inits = (*inputp)->config[0];
            nsam  = (*inputp)->config[1];
        } 
    }
    comb = (long int)((double)nsam*((double)nsam-(double)1)/(double)2);
	for(x=0;x<10;x++) (ntpar)->Sanc[x] = 0;
    
	if(segsit == 0 || nsam < 2) {
		if(nsam == 0) {
			(ntpar)->B1 = -10000;
			(ntpar)->Q1 = -10000;
			(ntpar)->k = (double)-10000;
			(ntpar)->piw = (double)-10000;
			(ntpar)->pib = (double)-10000;
			(ntpar)->S = -10000;
			(ntpar)->nhapl = -10000;
			for(a=0;a<2;a++) {
				(ntpar)->freq[a] = -10000;
				(ntpar)->unic[a] = -10000;
				(ntpar)->fhapl[a] = -10000;
			}
			(ntpar)->fhapl[0] = -10000;
			(ntpar)->maxhapl = -10000;
			(ntpar)->maxhapl1 = -10000;
			(ntpar)->Rm = (int)-10000;
			(ntpar)->thetaL = (double)-10000;
			(ntpar)->withinw = (double)-10000;
			(ntpar)->Se1 = -10000;
			(ntpar)->Sn1 = -10000;
			(ntpar)->pie1 = -10000;
			(ntpar)->pin1 = -10000;
			(ntpar)->m_sdev = -10000;
			(ntpar)->m_skew = -10000;
			(ntpar)->m_kurt = -10000;
			(ntpar)->ragg = -10000;
		}
		else {
			(ntpar)->B1 = 0;
			(ntpar)->Q1 = 0;
			(ntpar)->k = (double)0;
			(ntpar)->S = 0;
			(ntpar)->piw = (double)0;
			(ntpar)->pib = (double)0;
			(ntpar)->nhapl = 1;
			for(a=0;a<nsam;a++) {
			   (ntpar)->freq[a] = 0;
			   (ntpar)->unic[a] = 0;
			   (ntpar)->fhapl[a] = 0;
			}
			for(j=0;j<segsit;j++) { /*all valid positions are invariant positions in nsam=1*/
				if((h=ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit)) == 0) (ntpar)->freq[0] += 1;
			}		
			(ntpar)->fhapl[0] = nsam;
			(ntpar)->maxhapl = nsam;
			(ntpar)->maxhapl1 = nsam;
			(ntpar)->Rm = (int)0;
			(ntpar)->thetaL = (double)0;
			(ntpar)->withinw = (double)0;
		}
		(ntpar)->max_iES = (double)-10000;
		(ntpar)->min_uiHS = (double)-10000;
		(ntpar)->max_uiHS = (double)-10000;
		(ntpar)->mhsites = 0;
		(ntpar)->Se1 = 0;
		(ntpar)->Sn1 = 0;
		(ntpar)->pie1 = 0;
		(ntpar)->pin1 = 0;
		(ntpar)->m_sdev = -10000;
		(ntpar)->m_skew = -10000;
		(ntpar)->m_kurt = -10000;
		(ntpar)->ragg = -10000;
	}
	else {
		if((veca = (int **)calloc(segsit,sizeof(int *))) == 0) {
			puts("calloc error veca.1");
			exit(1);
		}
		for(j=0;j<(int)segsit;j++) {
			if((veca[j] = (int *)calloc((*inputp)->nsam,sizeof(int))) == 0) {
				puts("calloc error veca.2");
				exit(1);
			}
		}
		
		/*mismatch distribution*/
		if((freq1 = (double *)calloc(comb,sizeof(double))) == 0) {
			puts("calloc error m.1");
			exit(1);
		}
		if((freq2 = (double *)calloc(comb,sizeof(double))) == 0) {
			puts("calloc error m.2");
			exit(1);
		}
		if((freq3 = (double *)calloc(comb,sizeof(double))) == 0) {
			puts("calloc error m.3");
			exit(1);
		}
		if((freq4 = (double *)calloc(comb,sizeof(double))) == 0) {
			puts("calloc error m.4");
			exit(1);
		}
		   
		if(valuer || (*inputp)->mhits) (ntpar)->Rm = Min_rec(0,(int)segsit,nsam,inits,(*inputp)->nsam);
		else (ntpar)->Rm = (int)0;

        /* calcul B' i Q, pero sense dividir per S (Wall) */
        B = 0;
        A = 0;

        a = 0;
        if((segsitesm1 = segsit-1) < 0) segsitesm1 = 0;
        for(j=0;j<(long int)segsitesm1;) {
            k = j;
            while(k+1 < segsit) { /*calcular k*/
                if((ispolnomhit(k,inits,nsam,(*inputp)->nsam, list, posit)) > 0) break;
                else {
					k++;
				}
            }
            j = k+1;
            while(j < segsit) { /*calcular j*/
                if((ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit)) > 0) break;
                else j++;
            }
            if(j < segsit) {                
                val20 = val21 = -1;
                b = 0;
                for(i=inits;i<inits+nsam;i++) {
                    val10 = (list[i][k] - 48)*4 + (list[i][j]- 48);
                    if(val20 == -1) val20 = val10;
                    else{
                        if(val20 != val10) {
                            if(val21 == -1) val21 = val10;
                            else  if(val21 != val10) {
                                b = 1;
                                break;
                            }
                        }
                    }
                }
				if(!b) {/*congruent*/
					B += 1;
					if(!A) {
						for(i=inits;i<inits+nsam;i++) {
							if(list[i][j] > '0') x = '1';
							else x = '0';
							veca[A][i] = x;
						}
						A += 1;
					}
					else {
						a = 0;
						for(c=0;c<A;c++) {
							d = 0;
							for(i=inits;i<inits+nsam;i++) {
								if(list[i][j] > '0') x = '1';
								else x = '0';
								if(veca[c][i] == x) {
									d += 1;
								}
							}
							if(d == nsam || d == 0) {
								a = 1;
								break;
							}
						}
						if(!a) {
							for(i=inits;i<inits+nsam;i++) {
								if(list[i][j] > '0') x = '1';
								else x = '0';
								veca[A][i] = x;
							}
							A += 1;
						}
					}
				}
			}
		}
		
        (ntpar)->Q1 = B+A;
        (ntpar)->B1 = B;
        
		for(d=0;d<(int)segsit;d++) free(veca[d]);
		free(veca);

        /* calcul de S, pi, fr, nhapl i unics/seq (excluding mhits)*/
        if((hapl = (char *)calloc(nsam*segsit,sizeof(char))) == NULL) {
            perror("calloc error calc_neutpar.0");
            exit(1);
        }
        k_ = (double)0;
        for(a=0;a<nsam;a++) {
            (ntpar)->freq[a] = 0;
            (ntpar)->unic[a] = 0;
        }
		/*
		pidcount = 0;
		for(a=inits;a<inits+nsam-1;a++) {
			for(b=a+1;b<inits+nsam;b++) {
				(ntpar)->pid[pidcount] = (double)0;
				pidcount++;
			}
		}
		*/
        S = 0;
        nhapl = 0;                
		nmh = 0;
		
		Se1 = Sn1 = 0;
		k_e1 = k_n1 = 0.;

        for(j=0;j<segsit;j++) {
            pi = 0;
            while(j < segsit) {
                if((h=ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit)) > 0) break; /*h is the frequency of the new mutation*/
                else {
					j++;
					if(h == -2) nmh += 1;
					if(h ==  0) (ntpar)->freq[0] += 1; /*invariant intrapop*/
				}
            }            
            if(j<segsit) {
                (ntpar)->freq[h] += 1;
                for(a=inits;a<inits+nsam;a++) {
                    hapl[(a-inits)*segsit+S] = list[a][j];
                    /*new two lines: for singleton mutations (no outgroup)*/
                    if(h == 1) if(list[a][j] != '0') (ntpar)->unic[a-inits] += 1;
                    if(h == nsam-1) if(list[a][j] == '0') (ntpar)->unic[a-inits] += 1;
                }
                S++;
				if(h > 1) Se1++;
				if(h > 1 && h < nsam-1) Sn1++;
				/*pidcount = 0;*/
				z = 0;/*mismatch dist*/
                for(a=inits;a<inits+nsam-1;a++) {
                    for(b=a+1;b<inits+nsam;b++) {
						if(list[a][j] != list[b][j]) {
							pi++;
							if(h > 1) k_e1++;
							if(h > 1 && h < nsam-1) k_n1++;
							freq1[z] += 1; /*mismatch dist*/
						}
						z++;/*mismatch dist*/
					}
				}
				k_ += pi;
            }
        }

        /*pi, etc..*/
		(ntpar)->k = k_/(double)comb;
        (ntpar)->S = S;
		(ntpar)->mhsites = nmh;

		(ntpar)->Se1 = Se1;
		(ntpar)->Sn1 = Sn1;
		(ntpar)->pie1 = k_e1/(double)comb;
		(ntpar)->pin1 = k_n1/(double)comb;
		
        if((haplotype = (int *)malloc(nsam*sizeof(int))) == NULL) {
            perror("calloc error calc_neutpar.1");
            exit(1);
        }
        (ntpar)->nhapl = 0;
        for(a=0;a<nsam;a++) haplotype[a] = 1;            
		for(a=0;a<nsam-1;a++) {
            if(haplotype[a]) {
                (ntpar)->nhapl += 1;
                for(b=a+1;b<nsam;b++) {
                    if(haplotype[b]) {
                        if(memcmp(hapl+a*segsit,hapl+b*segsit,S) == 0) { 
                            haplotype[a] += 1;
                            haplotype[b] = 0;
                        }
                    }
                }
            }
        }
        if(haplotype[a]) (ntpar)->nhapl += 1;
        for(a=0;a<nsam;a++) (ntpar)->fhapl[a] = haplotype[a]; /*calcular freq. haplotips*/
        
        maxhapl=0;
        for(a=0;a<nsam;a++) 
			if(haplotype[a]>maxhapl)
				(ntpar)->maxhapl = maxhapl = haplotype[a]; /*calcular maxfreq. haplotips*/
		
		/*mismatch distribution moments*/
		z = 0;
		for(a=inits;a<inits+nsam-1;a++) {
			for(b=a+1;b<inits+nsam;b++) {
				moment = (freq1[z] - (ntpar)->k);/*mismatch dist*/
				freq2[z] = moment * moment;/*mismatch dist*/
				freq3[z] = moment * moment * moment;/*mismatch dist*/
				freq4[z] = moment * moment * moment * moment;/*mismatch dist*/
				z++;/*mismatch dist*/
			}
		}
		/*mismatch distribution moments*/
		s2 = g1 = g2 = 0.;				
		for(z=0;z<comb;z++)  {
			s2 += freq2[z];
			g1 += freq3[z];
			g2 += freq4[z];
		}
		s = sqrt(s2/((double)comb-1));
		(ntpar)->m_sdev = s;
		if(s) {
			(ntpar)->m_skew = g1 * (double)comb/(((double)comb-2)*((double)comb-1)*s*s*s);
			(ntpar)->m_kurt = g2 * ((double)comb-1) * comb /(((double)comb-3)*((double)comb-2)*((double)comb-1)*s*s*s*s)
			- (3. * ((double)comb-1)*((double)comb-1))/(((double)comb-3)*((double)comb-2));
		}
		else {
			(ntpar)->m_skew = -10000.;
			(ntpar)->m_kurt = -10000.;
		}

		/*calculate mismatch distribution to do ragg*/
		maxpwd = 0;
		for(z=0;z<comb;z++) if(freq1[z] > maxpwd) maxpwd = (long int)freq1[z];
		Pwd  = (long int *)calloc((unsigned long)(maxpwd+1),sizeof(long int));
		for(z=0;z<comb;z++) {
			Pwd[(long int)freq1[z]] = Pwd[(long int)freq1[z]] + 1;
		}
		(ntpar)->ragg = raggadeness(Pwd,maxpwd,comb);
		free(Pwd);
		/**/

		free(freq1);
		free(freq2);
		free(freq3);
		free(freq4);
		
		/*EHH statistics*/
		(ntpar)->max_iES  = -10000.;
		(ntpar)->max_uiHS = -10000.;
		(ntpar)->min_uiHS = +10000.;

		#if iHStest == 1
		if((*inputp)->includeEHHrel) {
			if((fhap1 = (int *)malloc(nsam*sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			if((fhap2 = (int *)malloc(nsam*sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			if((nsam1 = (int *)malloc(nsam*sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			if((nsam2 = (int *)malloc(nsam*sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			if((EHSA = (double *)malloc(segsit*sizeof(double))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			if((EHSD = (double *)malloc(segsit*sizeof(double))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			if((nousedA = (int *)calloc(segsit,sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			if((nousedD = (int *)calloc(segsit,sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			/*calculate uiHS for each position; 2 vectors*/
			/*to look at single position, k must be fixed to the closes polymorphism the user is looking at*/
			if((*inputp)->ifehh_fix == 0) {
				ehh0 = 0;
				ehhs = segsit;
				(ntpar)->max_uiHS = (ntpar)->min_uiHS = -10000;
			}
			else {
				/*locate the single segregating site centered in (*inputp)->ehh_fixnt +/- (*inputp)->ehh_margin. if no one skip loop (ehh0=0;ehhs=0;(ntpar)->max_uiHS = (ntpar)->min_uiHS = -10000;) */
				ehhr = (long int)1e7;
				for(j=0;j<segsit;j++) {
					while(j < segsit) {
						if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, list, posit)) > 0) break;
						else j++;
					}                    
					if(j<segsit) {
						if(fabs((float)posit[j] - (float)posit[(*inputp)->ehh_fixnt]) < (float)(*inputp)->ehh_margin) {
							if(fabs((float)posit[j] - (float)posit[(*inputp)->ehh_fixnt]) < ehhr)
								ehhr = j;
						}
						else {
							if(ehhr <  (long int)1e7) 
								break;
						}
					}
				}
				if(ehhr ==  (long int)1e7) {
					ehh0 = ehhs = 0;
				}
				else {
					ehh0 = ehhr;
					ehhs = ehhr + 1;
				}
				(ntpar)->max_uiHS = (ntpar)->min_uiHS = -10000;
			}
			for(k=ehh0;k<ehhs;k++) {
				/*first define the (2) haplotype samples on the SNP core*/
				n1 = n2 = 0;
				for(a=0;a<nsam;a++) {
					if(hapl[a*segsit+k] == '0'/*ancestral*/) {
						nsam1[n1] = a;
						n1++;
					}
					else {
						nsam2[n2] = a;
						n2++;
					}
				}
				if(n1 == nsam || n1 == 0) continue;
				EHSA[k] = EHSD[k] = (double)1;
				nousedA[k] = nousedD[k] = 1;
				
				flagA = flagD = 1;
				for(j=k+1;j<segsit;j++) {
					if(flagA) {
						for(a=0;a<n1;a++) fhap1[a] = 1;
						for(a=0;a<n1-1;a++) {
							if(fhap1[a]) {
								for(b=a+1;b<n1;b++) {
									if(fhap1[b]) {
										if(memcmp(hapl+nsam1[a]*segsit+k,hapl+nsam1[b]*segsit+k,j-k+1) == 0) { 
											fhap1[a] += 1;
											fhap1[b] = 0;
										}
									}
								}
							}
						}
						EHSA[j] = (double)1-(double)testHap(n1,fhap1);
						if(EHSA[j] < 0.05) {
							flagA = 0;
							nousedA[j] = 1;
						}
					}
					else nousedA[j] = 1;
					if(flagD) {
						for(a=0;a<n2;a++) fhap2[a] = 1;
						for(a=0;a<n2-1;a++) {
							if(fhap2[a]) {
								for(b=a+1;b<n2;b++) {
									if(fhap2[b]) {
										if(memcmp(hapl+nsam2[a]*segsit+k,hapl+nsam2[b]*segsit+k,j-k+1) == 0) { 
											fhap2[a] += 1;
											fhap2[b] = 0;
										}
									}
								}
							}
						}
						EHSD[j] = (double)1-(double)testHap(n2,fhap2);
						if(EHSD[j] < 0.05) {
							flagD = 0;
							nousedD[j] = 1;
						}
					}
					else nousedD[j] = 1;
				}
				flagA = flagD = 1;
				for(j=k-1;j>0;j--) {
					if(flagA) {
						for(a=0;a<n1;a++) fhap1[a] = 1;
						for(a=0;a<n1-1;a++) {
							if(fhap1[a]) {
								for(b=a+1;b<n1;b++) {
									if(fhap1[b]) {
										if(memcmp(hapl+nsam1[a]*segsit+j,hapl+nsam1[b]*segsit+j,k-j+1) == 0) { 
											fhap1[a] += 1;
											fhap1[b] = 0;
										}
									}
								}
							}
						}
						EHSA[j] = (double)1-(double)testHap(n1,fhap1);
						if(EHSA[j] < 0.05) {
							flagA = 0;
							nousedA[j] = 1;
						}
					}
					else nousedA[j] = 1;
					if(flagD) {
						for(a=0;a<n2;a++) fhap2[a] = 1;
						for(a=0;a<n2-1;a++) {
							if(fhap2[a]) {
								for(b=a+1;b<n2;b++) {
									if(fhap2[b]) {
										if(memcmp(hapl+nsam2[a]*segsit+j,hapl+nsam2[b]*segsit+j,k-j+1) == 0) { 
											fhap2[a] += 1;
											fhap2[b] = 0;
										}
									}
								}
							}
						}
						EHSD[j] = (double)1-(double)testHap(n2,fhap2);
						if(EHSD[j] < 0.05) {
							flagD = 0;
							nousedD[j] = 1;
						}
					}
					else nousedD[j] = 1;
					/*if(flagA == 0 && flagD == 0) {
						memset(nousedA+0,1,(j)*sizeof(int));
						memset(nousedD+0,1,(j)*sizeof(int));
						break;
					}*/
				}
				/*integrate and obtain uiHS, calculate max value*/
				uiHS = iHHA = iHHD = (double)0;
				for(j=1;j<segsit;j++) {
					if(nousedA[j]) continue;
					l=j-1;
					while(l>= 0 && nousedA[l]) l--;
					if(l==-1) continue;
					iHHA += ((double)(EHSA[l]+EHSA[j])*(double)(posit[j]-posit[l]))/(double)2.0;
				}
				for(j=1;j<segsit;j++) {
					if(nousedD[j]) continue;
					l=j-1;
					while(l>= 0 && nousedD[l]) l--;
					if(l==-1) continue;
					iHHD += ((double)(EHSD[l]+EHSD[j])*(double)(posit[j]-posit[l]))/(double)2.0;
				}
				if(iHHD != (double)0 && (double)((double)iHHA/(double)iHHD) > (double)0) {
					uiHS = (double)log((double)iHHA/(double)iHHD);
					if(k == 0 || (ntpar)->max_uiHS < uiHS) (ntpar)->max_uiHS = uiHS;
					if(k == 0 || (ntpar)->min_uiHS > uiHS) (ntpar)->min_uiHS = uiHS;
				}
				for(j=0;j<segsit;j++) {
					nousedA[j] = 0;
					nousedD[j] = 0;
				}
			}
			free(fhap1);
			free(fhap2);
			free(nsam1);
			free(nsam2);
			free(EHSA);
			free(EHSD);
			free(nousedA);
			free(nousedD);
		}
		#endif
		
		if((ntpar)->min_uiHS == 10000.) (ntpar)->min_uiHS = -10000.;
				
		#if iEStest == 1
		if((*inputp)->includeEHHrel) {
			/*to look at single position, k must be fixed to the closes polymorphism the user is looking at*/
			if((*inputp)->ifehh_fix == 0) {
				ehh0 = 0;
				ehhs = segsit;
				(ntpar)->max_iES = -10000;
			}
			else {
				/*locate the single segregating site centered in (*inputp)->ehh_fixnt +/- (*inputp)->ehh_margin. if no one skip loop (ehh0=0;ehhs=-1;(ntpar)->max_uiHS = (ntpar)->min_uiHS = -10000;) */
				ehhr =  (long int)1e7;
				for(j=0;j<segsit;j++) {
					while(j < segsit) {
						if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, list, posit)) > 0) break;
						else j++;
					}                    
					if(j<segsit) {
						if(fabs((float)posit[j] - (float)posit[(*inputp)->ehh_fixnt]) < (float)(*inputp)->ehh_margin) {
							if(fabs((float)posit[j] - (float)posit[(*inputp)->ehh_fixnt]) < ehhr)
								ehhr = j;
						}
						else {
							if(ehhr <  (long int)1e7) 
								break;
						}
					}
				}
				if(ehhr ==  (long int)1e7) {
					ehh0 = ehhs = 0;
				}
				else {
					ehh0 = ehhr;
					ehhs = ehhr + 1;
				}
				(ntpar)->max_iES = -10000;
			}
			if((EHHS = (double **)calloc(segsit,sizeof(double *))) == 0) {
				puts("calloc error veca.3");
				exit(1);
			}
			for(j=ehh0;j<ehhs;j++) {
				/*the half matrix and the diagonal vector is here included*/
				if((EHHS[j] = (double *)calloc(segsit,sizeof(double))) == 0) {
					puts("calloc error veca.4");
					exit(1);
				}
			}
			if((iES = (double *)calloc(segsit,sizeof(double))) == 0) {
				puts("calloc error veca.4");
				exit(1);
			}
			/*fhap_ij are the haplotype frequencies between 2 positions*/
			if((fhap_ij = (int *)malloc(nsam*sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			if((nousedM = (int **)calloc(segsit,sizeof(int *))) == NULL) {
				perror("calloc error calc_neutpar.1");
				exit(1);
			}
			for(j=ehh0;j<ehhs;j++) {
				if((nousedM[j] = (int *)calloc(segsit,sizeof(int))) == NULL) {
					perror("calloc error calc_neutpar.1");
					exit(1);
				}
			}
						
			/*calculate EHHS for each position; calculate half matrix and duplicate*/
			for(k=ehh0;k<ehhs;k++) {
				flagM = 1;
				for(j=k;j<segsit;j++) {				
					if(flagM) {
						for(a=0;a<nsam;a++) fhap_ij[a] = 1;
						for(a=0;a<nsam-1;a++) {
							if(fhap_ij[a]) {
								for(b=a+1;b<nsam;b++) {
									if(fhap_ij[b]) {
										if(memcmp(hapl+a*segsit+k,hapl+b*segsit+k,j-k+1) == 0) { 
											fhap_ij[a] += 1;
											fhap_ij[b] = 0;
										}
									}
								}
							}
						}
						EHHS[k][j] = EHHS[j][k] = (double)1-(double)testHap(nsam,fhap_ij);
						/*weight for the position 'observed'*/
						if(k!=j) {
							EHHS[k][j] = EHHS[k][j] / EHHS[k][k];
							if(EHHS[k][j] < 0.1) {
								flagM = 0;
								nousedM[k][j] = 1;
							}
						}
					}
					else {
						nousedM[k][j] = 1;
					}
				}
			}
			free(fhap_ij);
			
			/*The other half of matrix must weight for the position 'observed'*/
			/*to look at single position, k must be fixed to the closes polymorphism the user is looking at*/
			for(k=ehh0;k<ehhs;k++) {
				flagM = 1;
				for(j=k+1;j<segsit;j++) {
					if(flagM) {
						EHHS[j][k] = (double)EHHS[j][k] / (double)EHHS[j][j];
						if(EHHS[j][k] < 0.1) {
							flagM = 0;
							nousedM[j][k] = 1;
						}
					}
					else nousedM[j][k] = 1;
				}
			}
			/*integrate by its physical position*/
			(ntpar)->max_iES = 0.;
			/*to look at single position, k must be fixed to the closes polymorphism the user is looking at*/
			for(k=ehh0;k<ehhs;k++) {
				for(j=1;j<segsit;j++) {
					if(nousedM[k][j]) continue;
					l=j-1;
					while(l>=0 && nousedM[k][l]) l--;
					if(l==-1) continue;
					iES[k] += ((double)(EHHS[k][l]+EHHS[k][j])*(double)(posit[j]-posit[l]))/(double)2;
				}
				if((ntpar)->max_iES < log(iES[k])) (ntpar)->max_iES = (double)log((double)iES[k]);
			}
			for(j=ehh0;j<(int)ehhs;j++) free(EHHS[j]);
			free(EHHS);
			free(iES);
			for(j=ehh0;j<(int)ehhs;j++) free(nousedM[j]);
			free(nousedM);
		}
		#endif
				
		#if MAXHAP1 == 1
		if(S>0) {
			if((sshg = (int *)calloc(S,sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.0");
				exit(1);
			}
			maxhapl1 = 0;
			for(a=0;a<nsam;a++) {
				for(j=0;j<S;j++) {
					sshg[j] = 0;
				}
				for(b=0;b<nsam;b++) {
					mut = 0;
					for(j=0;j<S;j++) {
						if(hapl[a*segsit+j] != hapl[b*segsit+j])
							mut += 1;
					}
					if(mut <= 1) {
						for(j=0;j<S;j++) {
							if(hapl[a*segsit+j] != hapl[b*segsit+j]) {
								sshg[j] += 1;
							}
							if(mut == 0) sshg[j] += 1;
						}
					}
				}
				mh1 = 0;
				for(j=0;j<S;j++) {
					if(mh1 < sshg[j]) mh1 = sshg[j];
				}
				if(maxhapl1 < mh1)
					maxhapl1 = mh1;
			}
			(ntpar)->maxhapl1 = maxhapl1;
			free(sshg);
		}
		else {
			(ntpar)->maxhapl1 = nsam;
		}
		#else
		(ntpar)->maxhapl1 = -10000;
		#endif
        free(hapl);
        free(haplotype);

        /*calcular thetaL*/
		(ntpar)->thetaL = 0;
		for(a=1;a<nsam;a++) {
            (ntpar)->thetaL += (double)a * (double)ntpar->freq[a];
		}
		(ntpar)->thetaL *= (double)1/(double)(nsam-1);
	}
	
	if(npopa < (*inputp)->max_npop_sampled) {
		/* Calcular pi_within i pi_between: per poblacions amb config > 0*/
		/* Changed: calculate piwithin(avg for all) and the pib for the current pop vs the rest of populations*/
		
		nsamallpop = 0;
		for(h=0;h<(*inputp)->npop_sampled;h++) nsamallpop += (*inputp)->config[h];

		if((*inputp)->max_npop_sampled > 2 || ((*inputp)->max_npop_sampled == 2 && (*inputp)->pop_outgroup == -1)) {
			if((hapl = (char *)calloc(nsamallpop*segsit,sizeof(char))) == NULL) {
				perror("calloc error calc_neutpar.0");
				exit(1);
			}
			if((piw = (long int *)calloc((*inputp)->max_npop_sampled,sizeof(long int))) == NULL) {
				perror("calloc error calc_neutpar.0");
				exit(1);
			}
			if((pib =(long int *)calloc((*inputp)->max_npop_sampled,sizeof(long int))) ==NULL) {
				perror("calloc error calc_neutpar.0");
				exit(1);
			}
			if((hapw = (long int *)calloc((*inputp)->max_npop_sampled,sizeof(long int))) == NULL) {
				perror("calloc error calc_neutpar.0");
				exit(1);
			}
			if((hapb =(long int *)calloc((*inputp)->max_npop_sampled,sizeof(long int))) ==NULL) {
				perror("calloc error calc_neutpar.0");
				exit(1);
			}
			
			(ntpar)->piw = 0.;
			(ntpar)->withinw = 0.;
			(ntpar)->pib = 0.;
			
			/*calculating piw and pib*/
			S = 0;
			for(j=0;j<segsit;j++) {
				while(j < segsit) {
					if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, list, posit)) > 0) break;
					else j++;
				}                    
				if(j<segsit) {
					inits1 = 0;
					for(h=0;h<(*inputp)->npop_sampled;h++) {
						for(a=inits1;a<inits1+(*inputp)->config[h]-1;a++) {
							hapl[(a)*segsit+S] = list[a][j];
						}
						c = ispolnomhit(j,inits1,(*inputp)->config[h],(*inputp)->nsam, list, posit);
						piw[h] += (c * ((*inputp)->config[h] - c));
						hapl[(inits1+(*inputp)->config[h]-1)*segsit+S] = list[inits1+(*inputp)->config[h]-1][j];
						inits1 += (*inputp)->config[h];
					}
					S++;
					
					inits1 = 0;
					inits2 = inits;
					for(h=0;h<(*inputp)->npop_sampled;h++) {
						if((*inputp)->config[h] > 0 && (*inputp)->config[npopa] > 0 && npopa != h /**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/
						   && ((ispolnomhit(j,inits1,(*inputp)->config[h],(*inputp)->nsam, list, posit)>=0)
						   &&  (ispolnomhit(j,inits2,(*inputp)->config[npopa],(*inputp)->nsam, list, posit)>=0))) {
							for(a=inits1;a<inits1+(*inputp)->config[h];a++)
								for(b=inits2;b<inits2+(*inputp)->config[npopa];b++)
									if(list[a][j] != list[b][j]) pib[h]++;
						}
						inits1 += (*inputp)->config[h];
					}
				}
			}
			/*calculate fixed values in pop npopa*/
			inits2 = inits;
			for(j=0;j<segsit;j++) {
				inits1 = 0;
				for(h=0;h<(*inputp)->pop_outgroup;h++) inits1 += (*inputp)->config[h];
				if(ispolnomhit(j,inits2,(*inputp)->config[npopa],(*inputp)->nsam, list, posit) == 0 &&
				   ispolnomhit(j,inits1,(*inputp)->config[(*inputp)->pop_outgroup],(*inputp)->nsam, list, posit) == 0) {
					if(list[inits2][j] == list[inits1][j]) 
						(ntpar)->freq[0] -= 1;
				}
				if(ispolnomhit(j,inits2,(*inputp)->config[npopa],(*inputp)->nsam, list, posit) == 0 &&
				   ispolnomhit(j,inits1,(*inputp)->config[(*inputp)->pop_outgroup],(*inputp)->nsam, list, posit) < 0) {
					(ntpar)->freq[0] -= 1;
				}
			}			
			/*calculating hapw and hapb*/
			inits1 = 0;
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				for(a=inits1;a<inits1+(*inputp)->config[h]-1;a++) {
					for(b=a+1;b<inits1+(*inputp)->config[h];b++) {
						if(memcmp(hapl+a*segsit,hapl+b*segsit,S) == 0) { 
							hapw[h] += 1;
						}
					}
				}
				inits1 += (*inputp)->config[h];
			}
			inits1 = 0;
			inits2 = inits;
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h] > 0 && (*inputp)->config[npopa] > 0 && npopa != h && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) {
					for(a=inits1;a<inits1+(*inputp)->config[h];a++) {
						for(b=inits2;b<inits2+(*inputp)->config[npopa];b++) {
							if(memcmp(hapl+a*segsit,hapl+b*segsit,S) == 0) { 
								hapb[h] += 1;
							}
						}
					}
				}
				inits1 += (*inputp)->config[h];
			}
			
			/*weighting piw and pib*/
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h]>1) (ntpar)->piwallcomp[h] = piw[h] / (double)(((*inputp)->config[h] * ((*inputp)->config[h] - 1))/2);
			}
            for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h] && (*inputp)->config[npopa] && npopa != h/**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) {
					(ntpar)->piaallcomp[h] = pib[h] / (double)((*inputp)->config[h] * (*inputp)->config[npopa]);/*keeping all comparisons*/
					if((ntpar)->piaallcomp[h]) {
						if((*inputp)->config[h] > 1 && (*inputp)->config[npopa] > 1) 
							(ntpar)->fstallcomp[h] = 1.0 - ((double)((ntpar)->piwallcomp[h]+(ntpar)->piwallcomp[npopa])/2.0)/(double)(ntpar)->piaallcomp[h];
						else {
							if((*inputp)->config[h] == 1 && (*inputp)->config[npopa] > 1)
								(ntpar)->fstallcomp[h] = 1.0 - ((double)((ntpar)->piwallcomp[npopa])/2.0)/(double)(ntpar)->piaallcomp[h];
							if((*inputp)->config[h] > 1 && (*inputp)->config[npopa] == 1)
								(ntpar)->fstallcomp[h] = 1.0 - ((double)((ntpar)->piwallcomp[h])/2.0)/(double)(ntpar)->piaallcomp[h];
							if((*inputp)->config[h] == 1 && (*inputp)->config[npopa] == 1)
								(ntpar)->fstallcomp[h] = 1.0;
						}
					}
					else (ntpar)->fstallcomp[h] = -10000.;
				}
				else {
					(ntpar)->piaallcomp[h] = -10000.;
					(ntpar)->fstallcomp[h] = -10000.;
				}
			}
			
			/*weighting hapw and hapb*/
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h]>1 /**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) (ntpar)->hapwallcomp[h] = 1.0 - hapw[h] / (double)(((*inputp)->config[h] * ((*inputp)->config[h] - 1))/2);
				/*else (ntpar)->hapwallcomp[h] = -10000.;*/
			}
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h] && (*inputp)->config[npopa] && npopa != h /**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) {
					(ntpar)->hapaallcomp[h] = 1.0 - hapb[h] / (double)((*inputp)->config[h] * (*inputp)->config[npopa]);/*keeping all comparisons*/
					if((ntpar)->hapaallcomp[h]) {
						if((*inputp)->config[h] > 1 && (*inputp)->config[npopa] > 1) 
							(ntpar)->fsthapallcomp[h] = 1.0 - ((double)((ntpar)->hapwallcomp[h]+(ntpar)->hapwallcomp[npopa])/2.0)/(double)(ntpar)->hapaallcomp[h];
						else {
							if((*inputp)->config[h] == 1 && (*inputp)->config[npopa] > 1)
								(ntpar)->fsthapallcomp[h] = 1.0 - ((double)((ntpar)->hapwallcomp[npopa])/2.0)/(double)(ntpar)->hapaallcomp[h];
							if((*inputp)->config[h] > 1 && (*inputp)->config[npopa] == 1)
								(ntpar)->fsthapallcomp[h] = 1.0 - ((double)((ntpar)->hapwallcomp[h])/2.0)/(double)(ntpar)->hapaallcomp[h];
							if((*inputp)->config[h] == 1 && (*inputp)->config[npopa] == 1)
								(ntpar)->fsthapallcomp[h] = 1.0;
						}
					}
					else (ntpar)->fsthapallcomp[h] = -10000.;
				}
				else {
					(ntpar)->hapaallcomp[h] = -10000.;
					(ntpar)->fsthapallcomp[h] = -10000.;
				}
			}
			
			comb2 = 0;
			/*es recull el valor a cada subpop individualment*//**/
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h]>1 /**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) {
					comb = (int)((double)(*inputp)->config[h] * ((double)(*inputp)->config[h] - (double)1.0)/(double)2.0);
					(ntpar)->piw += (double)piw[h]/(double)comb; 
					comb2++;
				}
			}			
			/*despres es divideix pel nombre de subpoblacions*/
			if(comb2) {
				(ntpar)->piw = (ntpar)->piw/(double)comb2;
				npw=0;
				for(h=0;h<(*inputp)->npop_sampled;h++) npw += (*inputp)->config[h];
			}
			
			(ntpar)->withinw = -10000;
			comb2 = 0;
			
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h] && (*inputp)->config[npopa] && npopa != h /**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) {
					comb = (*inputp)->config[h] * (*inputp)->config[npopa];
					(ntpar)->pib += (double)pib[h]/(double)comb;
					comb2++;
				}
			}
			if(comb2) (ntpar)->pib = (ntpar)->pib/(double)(comb2);			

			free(piw);
			free(pib);
			free(hapw);
			free(hapb);
			free(hapl);
		}
	}
	/*Calculate ancestral and shared polymorphisms*/
	if((*inputp)->type_ancestral == 2 || (*inputp)->type_ancestral == 3) {
		for(h=0;h<10;h++) {
			(ntpar)->Sanc[h] = 0;
		}
		for(j=0;j<segsit;j++) {
			while(j < segsit) {
				if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, list, posit)) > 0) break;
				else j++;
			}                    
			if(j<segsit) {					
				pola = polb = polc = 0; /*if polymorphic in pop1,pop2,popo*/
				
				/*find the value inits1,inits2 and initso: indicate the position in list[inits12o][j] for comparison of polymorphsms*/
				inits1 = -1;
				for(pop1=1;pop1<=(*inputp)->ancestral_pol[0][0];pop1++) {
					if((*inputp)->config[(*inputp)->ancestral_pol[0][pop1]]) {
						inits1 = 0;
						for(h=0;h<(*inputp)->ancestral_pol[0][pop1];h++) inits1 += (*inputp)->config[h];
						break;
					}
				}
				if(inits1 == -1) pola = -1;
				inits2 = -1;
				for(pop2=1;pop2<=(*inputp)->ancestral_pol[1][0];pop2++) {
					if((*inputp)->config[(*inputp)->ancestral_pol[1][pop2]]) {
						inits2 = 0;
						for(h=0;h<(*inputp)->ancestral_pol[1][pop2];h++) inits2 += (*inputp)->config[h];
						break;
					}
				}
				if(inits2 == -1) polb = -1;

				initso = -1;
				for(popo=1;popo<=(*inputp)->ancestral_pol[2][0];popo++) {
					if((*inputp)->config[(*inputp)->ancestral_pol[2][popo]]) {
						initso = 0;
						for(h=0;h<(*inputp)->ancestral_pol[2][popo];h++) initso += (*inputp)->config[h];
						break;
					}
				}
				
				/*check if each group is polymorphic or not*/
				if(pola != -1 && polb != -1) {
					for(pop1=1;pop1<=(*inputp)->ancestral_pol[0][0];pop1++) {
						nt0 = 0;
						for(h=0;h<(*inputp)->ancestral_pol[0][pop1];h++) nt0 += (*inputp)->config[h];	
						for(h=nt0;h<nt0+(*inputp)->config[(*inputp)->ancestral_pol[0][pop1]];h++) {
							if(list[h][j] != list[inits1][j]) {
								pola = 1;
								break;
							}
						}
						if(pola) break;
					}
					for(pop2=1;pop2<=(*inputp)->ancestral_pol[1][0];pop2++) {
						nt0 = 0;
						for(h=0;h<(*inputp)->ancestral_pol[1][pop2];h++) nt0 += (*inputp)->config[h];	
						for(h=nt0;h<nt0+(*inputp)->config[(*inputp)->ancestral_pol[1][pop2]];h++) {
							if(list[h][j] != list[inits2][j]) {
								polb = 1;
								break;
							}
						}
						if(polb) break;
					}
					for(popo=1;popo<=(*inputp)->ancestral_pol[2][0];popo++) {
						nt0 = 0;
						for(h=0;h<(*inputp)->ancestral_pol[2][popo];h++) nt0 += (*inputp)->config[h];	
						for(h=nt0;h<nt0+(*inputp)->config[(*inputp)->ancestral_pol[2][popo]];h++) {
							if(list[h][j] != list[initso][j]) {
								polc = 1;
								break;
							}
						}
						if(polc) break;
					}
				}
				
				/*define classes*/
				if(pola==1 && polb==0 && polc==0) {
					if((*inputp)->ancestral_pol[2][0]) {
						if(list[inits2][j] == list[initso][j]) (ntpar)->Sanc[0] += 1;/*Sx1*/
						else (ntpar)->Sanc[6] += 1;/*Sx1f2*/
					}
					else {
						(ntpar)->Sanc[0] += 1;/*Sx1*/
					}
				}
				if(pola==0 && polb==1 && polc==0) {
					if((*inputp)->ancestral_pol[2][0]) {
						if(list[inits1][j] == list[initso][j]) (ntpar)->Sanc[1] += 1;/*Sx2*/
						else (ntpar)->Sanc[7] += 1;/*Sx2f1*/
					}
					else {
						(ntpar)->Sanc[1] += 1;/*Sx2*/
					}
				}
				if(polc == 1 && (pola != -1 && polb != -1)) {
					if((pola==0 && polb==0) && (list[inits1][j] == list[inits2][j])) {
						(ntpar)->Sanc[2] += 1;/*Sxo*/
					}
					else {
						(ntpar)->Sanc[9] += 1;/*Sso*/
					}
				}
				if(pola==1 && polb==1 && polc==0) {
					(ntpar)->Sanc[8] += 1;/*Ssh*/
				}
				if(pola==0 && polb==0 && polc==0) {
					if((*inputp)->ancestral_pol[2][0]) {
						if(list[inits2][j] == list[initso][j] && list[inits1][j] != list[initso][j]) {
							(ntpar)->Sanc[3] += 1;/*Sf1*/
						}
						if(list[inits1][j] == list[initso][j] && list[inits2][j] != list[initso][j]) {
							(ntpar)->Sanc[4] += 1;/*Sf2*/
						}
						if(list[inits1][j] == list[inits2][j] && list[inits1][j] != list[initso][j]) {
							(ntpar)->Sanc[5] += 1;/*Sfo*/
						}
					}
					else {
						(ntpar)->Sanc[3] += 1;/*Sf*/
					}
				}
			}
		}
	}
	if((*inputp)->type_ancestral == 1) {
		for(h=0;h<10;h++) {
			(ntpar)->Sanc[h] = 0;
		}
		(ntpar)->Sanc[0] = (int)(ntpar)->S;
	}
	if((*inputp)->type_ancestral > 3) {
		initsq1 = (int *)calloc((*inputp)->type_ancestral,sizeof(int));
		initsq2 = (int *)calloc((*inputp)->type_ancestral,sizeof(int));
		polqa   = (int *)calloc((*inputp)->type_ancestral,sizeof(int));
		polqb   = (int *)calloc((*inputp)->type_ancestral,sizeof(int));
		sq      = (int *)calloc((*inputp)->type_ancestral,sizeof(int));
		
		for(h=0;h<4*(*inputp)->type_ancestral;h++) {
			(ntpar)->Sanc[h] = 0;
		}
		/*find the value initsq1,initsq2,initso: indicate the POSITION in list[inits][j] for comparison of polymorphsms*/
		/*do for the outgroup (the last population is defined as outgroup)*/
		initso = -1;
		for(popo=1;popo<=(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][0];popo++) {
			if((*inputp)->config[(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][popo]]) {
				initso = 0;
				for(h=0;h<(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][popo];h++) initso += (*inputp)->config[h];
				break;
			}
		}
		/*do for all the rest of populations*/
		for(x=0;x<(*inputp)->type_ancestral-1;x++) {
			inits1 = -1;
			for(pop1=1;pop1<=(*inputp)->ancestral_pol[x][0];pop1++) {
				if((*inputp)->config[(*inputp)->ancestral_pol[x][pop1]]) {
					initsq1[x] = 0;
					for(h=0;h<(*inputp)->ancestral_pol[x][pop1];h++) initsq1[x] += (*inputp)->config[h];
					break;
				}
			}
			
			initsq2[x] = sq[x] = 0;
			for(pop1=0;pop1<(*inputp)->type_ancestral-1;pop1++) {
				if(pop1 != x) { 
					for(pop2=1;pop2<=(*inputp)->ancestral_pol[pop1][0];pop2++) {
						if((*inputp)->config[(*inputp)->ancestral_pol[pop1][pop2]]) {
							sq[x] = 1;
							for(h=0;h<(*inputp)->ancestral_pol[pop1][pop2];h++) initsq2[x] += (*inputp)->config[h];
							break;
						}
					}
				}
				if(sq[x] == 1) break;
			}
		}
		for(j=0;j<segsit;j++) {
			while(j < segsit) {
				if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, list, posit)) > 0) break;
				else j++;
			}                    
			if(j<segsit) {
				/*do for each population, first outgroup*/
				polc = 0;
				for(popo=1;popo<=(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][0];popo++) {
					nt0 = 0;
					for(h=0;h<(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][popo];h++) nt0 += (*inputp)->config[h];	
					for(h=nt0;h<nt0+(*inputp)->config[(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][popo]];h++) {
						if(list[h][j] != list[initso][j]) {
							polc = 1;
							break;
						}
					}
					if(polc) break;
				}
				for(x=0;x<(*inputp)->type_ancestral-1;x++) {					
					/*check if each group is polymorphic or not*/
					polqa[x] = polqb[x] = 0;
					if(initsq1[x] != -1 && sq[x] != 0) {
						for(pop1=1;pop1<=(*inputp)->ancestral_pol[x][0];pop1++) {
							nt0 = 0;
							for(h=0;h<(*inputp)->ancestral_pol[x][pop1];h++) nt0 += (*inputp)->config[h];	
							for(h=nt0;h<nt0+(*inputp)->config[(*inputp)->ancestral_pol[x][pop1]];h++) {
								if(list[h][j] != list[initsq1[x]][j]) {
									polqa[x] = 1;
									break;
								}
							}
							if(polqa[x]) break;
						}
						for(pop1=0;pop1<(*inputp)->type_ancestral-1;pop1++) {
							if(pop1 != x) {
								for(pop2=1;pop2<=(*inputp)->ancestral_pol[pop1][0];pop2++) {
									nt0 = 0;
									for(h=0;h<(*inputp)->ancestral_pol[pop1][pop2];h++) nt0 += (*inputp)->config[h];	
									for(h=nt0;h<nt0+(*inputp)->config[(*inputp)->ancestral_pol[pop1][pop2]];h++) {
										if(list[h][j] != list[initsq2[x]][j]) {
											polqb[x] = 1;
											break;
										}
									}
									if(polqb[x]) break;
								}
							}
							if(polqb[x]) break;
						}
					}
					
					/*define classes*/
					if(polqa[x]==1 && polqb[x]==0 && polc==0) {
						if(list[initsq2[x]][j] == list[initso][j]) (ntpar)->Sanc[x*4+0] += 1;/*Sx1*/
						else (ntpar)->Sanc[x*4+2] += 1;/*Sx1f2*/
					}
					
					if(polqa[x]==1 && polqb[x]==1 && polc==0) {
						(ntpar)->Sanc[x*4+3] += 1;/*Ssh*/
					}
					if(polqa[x]==0 && polqb[x]==0 && polc==0) {
						if((list[initsq2[x]][j] == list[initso][j]) && (list[initsq1[x]][j] != list[initso][j])) {
							(ntpar)->Sanc[x*4+1] += 1;/*Sf1*/
						}
					}
				}
				/*for the outgroup is doesn't matter what is the group of populations(x), we chose 0*/
				if(polc == 1 && (polqa[0] != -1 && polqb[0] != -1)) {
					(ntpar)->Sanc[((*inputp)->type_ancestral-1)*4+0] += 1;/*Sxo*/
				}
				if(polqa[0]==0 && polqb[0]==0 && polc==0) {
					if((list[initsq1[0]][j] == list[initsq2[0]][j]) && (list[initsq1[0]][j] != list[initso][j])) {
						(ntpar)->Sanc[((*inputp)->type_ancestral-1)*4+1] += 1;/*Sfo*/
					}
				}	
			}
		}
		(ntpar)->Sanc[((*inputp)->type_ancestral-1)*4+2] = -10000;
		(ntpar)->Sanc[((*inputp)->type_ancestral-1)*4+3] = -10000;
		
		free(initsq1);
		free(initsq2);
		free(polqa);
		free(polqb);
	}
}


double Zns(int valuep,long int segsites,struct var2 **inputp,int npopa)
{
    double ZnS;
    long k,j,comb;
    int i,inits,x;
	int nsam=0;
    int vala,valb,val00,a,b;
	int na=0;
	int nb=0;
    double A,B,C;

    if(valuep == 0) {
		inits = 0;
		for(x=0;x<(*inputp)->npop_sampled;x++) {
			if(x < npopa) inits += (*inputp)->config[x];
			else {
				nsam = (*inputp)->config[x];
				break;
			}
		}
    }
    else {
        if(valuep == 1) {
            inits = 0;
            nsam  = (*inputp)->config[0];
        }
        else {
            inits = (*inputp)->config[0];
            nsam  = (*inputp)->config[1];
        } 
    }
    if(nsam == 0) return(-10000.);
    
    ZnS = (double)0;
    j = 0;
    comb = 0;
    while(j+1 < (long)segsites) {
        while(j < (long)segsites) {
            if(ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)segsites) {
            while(k < (long)segsites) {
                if(ispolnomhit(k,inits,nsam,(*inputp)->nsam, list, posit) > 0) break;
                else k++;
            }
            if(k < (long)segsites) {
                /*calcular freqs p1,q1*/
                vala = valb = 1;
                a = list[inits][j];
                b = list[inits][k];
                for(i=1+inits;i<inits+nsam;i++) {
                    if(list[i][j] == a) vala++;
                    else na = list[i][j];
                    if(list[i][k] == b) valb++;                                
                    else nb = list[i][k];
                }
                if(nsam - vala > vala) {
                    a = na;
                    vala = nsam - vala;
                }
                if(nsam - valb > valb) {
                    b = nb;
                    valb = nsam - valb;
                }
                /*calcular p1q1*/
                val00 = 0;
                for(i=0+inits;i<inits+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
                
                /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
                A = (double)vala/(double)nsam;
                B = (double)valb/(double)nsam;
                C = (double)val00/(double)nsam;
                ZnS += ((C - (A*B)) * (C - (A*B))) / (A*(1.0 - A)*B*(1.0 - B));
                comb++;
            }
            k++;
        }
        j++;
    }
    if(comb) ZnS = ZnS/(double)comb;
    else return(-10000.);
    return(ZnS);
}

double ZnA_(int valuep,long int segsites,struct var2 **inputp,int npopa)
{
    double ZnA;
    long k,j,comb;
    int i,inits,x;
	int nsam=0;
    int vala,valb,val00,a,b;
	int na=0;
	int nb=0;
    double A,B,C;

    if(valuep == 0) {
		inits = 0;
		for(x=0;x<(*inputp)->npop_sampled;x++) {
			if(x < npopa) inits += (*inputp)->config[x];
			else {
				nsam = (*inputp)->config[x];
				break;
			}
		}
    }
    else {
        if(valuep == 1) {
            inits = 0;
            nsam  = (*inputp)->config[0];
        }
        else {
            inits = (*inputp)->config[0];
            nsam  = (*inputp)->config[1];
        } 
    }
    if(nsam == 0) return(-10000.);
    
    ZnA = (double)0;
    j = 0;
    comb = 0;
    while(j+1 < (long)segsites) {
        while(j < (long)segsites) {
            if(ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)segsites) {
            if(ispolnomhit(k,inits,nsam,(*inputp)->nsam, list, posit) > 0) break;
            else k++;
        }
        if(k < (long)segsites) {
            /*calcular freqs p1,q1*/
            vala = valb = 1;
            a = list[inits][j];
            b = list[inits][k];
            for(i=1+inits;i<inits+nsam;i++) {
                if(list[i][j] == a) vala++;
                else na = list[i][j];
                if(list[i][k] == b) valb++;                                
                else nb = list[i][k];
            }
            if(nsam - vala > vala) {
                a = na;
                vala = nsam - vala;
            }
            if(nsam - valb > valb) {
                b = nb;
                valb = nsam - valb;
            }
            /*calcular p1q1*/
            val00 = 0;
            for(i=0+inits;i<inits+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
            
            /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
            A = (double)vala/(double)nsam;
            B = (double)valb/(double)nsam;
            C = (double)val00/(double)nsam;
            /*if(vala != ((double)nsam-1.0) && valb != ((double)nsam-1.0)) if only informative */
            ZnA += ((C - (A*B)) * (C - (A*B))) / (A*(1.0 - A)*B*(1.0 - B));
            comb++;
        }
        j++;
    }
    if(comb) ZnA = ZnA/(double)comb;
    else return(-10000.);
    return(ZnA);
}

double fixoutg(int nsam, int Sout, int freq0)
{
    double div;
    
    if(nsam) {
        div = (double)Sout;
		div += (double)freq0;
		
        return div;
    }
    else return -10000;
}
double koutg(int nsam, int Sout, int *freq)
{
    int x;
    double div;
    
    if(nsam) {
        div = (double)Sout;
		div += (double)freq[0];
		
        for(x=1;x<nsam;x++) 
			div += (double)freq[x]*(double)x/(double)nsam;
        return div;
    }
    else return -10000;
}
double koutgJC(int nsam, int Sout, int *freq, unsigned long nsites)
{
    int x;
    double div;
    
    if(nsam) {
        div = (double)Sout;
		div += (double)freq[0];
		
		for(x=1;x<nsam;x++) 
            div += (double)freq[x]*(double)x/(double)nsam;
		
		if((double)Sout == 0) { /*correction only in case calculating from sequence data*/
			div/= (double)nsites;
			if(div >= (double).75) return (double) -10000;
			else {
				div = (double)-.75*(double)log((double)1-(double)4/(double)3*div);/*Jukes and Cantor correction...*/
				return div*(double)nsites;
			}
		}
		else return div; /*in case calculating directly Sout, we already included the multiple hits in the branch length*/
    }
    else return (double) -10000;
}

void calc_neutpar_windowSRH(struct var2 **inputp,struct dnapar *ntpar,long int  s0,long int  s1,double valuer,int npopa,int flaghap)
{
	long int segsit;    
	long int S;
    int nhapl;
    int *haplotype = 0;
    int inits;
	int nsam=0;
    long int j;
    int a,b,h,x;
	long int comb;
    char *hapl=0;
	int Min_rec(int,int,int,int,int);
	
	if(npopa == (*inputp)->npop) {
		inits = 0;
		nsam = (*inputp)->nsam;
	}
	else {
		inits = 0;
		for(x=0;x<(*inputp)->npop_sampled;x++) {
			if(x < npopa) inits += (*inputp)->config[x];
			else {
				nsam = (*inputp)->config[x];
				break;
			}
		}
	}
	
    comb = (int)((double)nsam*((double)nsam-(double)1)/(double)2);
	for(x=0;x<10;x++) (ntpar)->Sanc[x] = 0;
	
    /*define segsit first*/
    segsit = s1 - s0;

	if(segsit == 0 || nsam < 2) {
		if(nsam == 0) {
			(ntpar)->S = -10000;
			(ntpar)->nhapl = -10000;
			(ntpar)->Rm = (int)-10000;
		}
		else {
			(ntpar)->S = 0;
			(ntpar)->nhapl = 1;
			(ntpar)->Rm = (int)0;
		}
	}
	else {		
		if(valuer || (*inputp)->mhits) (ntpar)->Rm = Min_rec((int)s0,(int)s1,nsam,inits,(*inputp)->nsam);
		else (ntpar)->Rm = (int)0;
	}
	
	/* calcul de S,nhapl (excluding mhits)*/
	if(flaghap > 0) {
		if((hapl = (char *)calloc(nsam*segsit,sizeof(char))) == NULL) {
			perror("calloc error calc_neutpar.0");
			exit(1);
		}
	}
	
	S = 0;
	nhapl = 0;                
	for(j=s0;j<s1;j++) {
		while(j < s1) {
			if((h=ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit)) > 0) break; /*h is the frequency of the new mutation*/
			else j++;
		}            
		if(j<s1) {
			if(flaghap > 0) {
				for(a=inits;a<inits+nsam;a++) {
					hapl[(a-inits)*segsit+S] = list[a][j];
				}
			}
			S++;
		}
	}
	(ntpar)->S = S;
	
	if(flaghap > 0) {
		if((haplotype = (int *)malloc(nsam*sizeof(int))) == NULL) {
			perror("calloc error calc_neutpar.1");
			exit(1);
		}
		(ntpar)->nhapl = 0;
		for(a=0;a<nsam;a++) haplotype[a] = 1;            
		for(a=0;a<nsam-1;a++) {
			if(haplotype[a]) {
				(ntpar)->nhapl += 1;
				for(b=a+1;b<nsam;b++) {
					if(haplotype[b]) {
						if(memcmp(hapl+a*segsit,hapl+b*segsit,S) == 0) { 
							haplotype[a] += 1;
							haplotype[b] = 0;
						}
					}
				}
			}
		}
		if(haplotype[a]) (ntpar)->nhapl += 1;
		free(hapl);
		free(haplotype);
	}
	
	return;
}

double Zns_window(struct var2 **inputp, long int s0, long int s1,int npopa)
{
    double ZnS;
    long k,j,comb;
    int i,nsam,inits,x;
    int vala,valb,val00,a,b;
	int na=0;
	int nb=0;
    double A,B,C;

    nsam  = (*inputp)->nsam;
	inits = 0;
	for(x=0;x<(*inputp)->npop_sampled;x++) {
		if(x < npopa) inits += (*inputp)->config[x];
		else {
			nsam = (*inputp)->config[x];
			break;
		}
	}
    if(nsam == 0) return(-10000.);
	
    ZnS = (double)0;
    j = s0;
    comb = 0;
    while(j+1 < (long)s1) {
        while(j < (long)s1) {
            if(ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)s1) {
            while(k < (long)s1) {
                if(ispolnomhit(k,inits,nsam,(*inputp)->nsam, list, posit) > 0) break;
                else k++;
            }
            if(k < (long)s1) {
                /*calcular freqs p1,q1*/
                vala = valb = 1;
                a = list[inits][j];
                b = list[inits][k];
                for(i=1+inits;i<inits+nsam;i++) {
                    if(list[i][j] == a) vala++;
                    else na = list[i][j];
                    if(list[i][k] == b) valb++;                                
                    else nb = list[i][k];
                }
                if(nsam - vala > vala) {
                    a = na;
                    vala = nsam - vala;
                }
                if(nsam - valb > valb) {
                    b = nb;
                    valb = nsam - valb;
                }
                /*calcular p1q1*/
                val00 = 0;
                for(i=inits+0;i<inits+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
                
                /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
                A = (double)vala/(double)nsam;
                B = (double)valb/(double)nsam;
                C = (double)val00/(double)nsam;
                ZnS += ((C - (A*B)) * (C - (A*B))) / (A*(1.0 - A)*B*(1.0 - B));
                comb++;
            }
            k++;
        }
        j++;
    }
    if(comb) ZnS = ZnS/(double)comb;
    else return(-10000.);
    return(ZnS);
}

double ZnA_window(struct var2 **inputp, long int s0, long int s1,int npopa)
{
    double ZnA;
    long k,j,comb;
    int i,nsam,inits,x;
    int vala,valb,val00,a,b;
	int na=0;
	int nb=0;
    double A,B,C;

    nsam  = (*inputp)->nsam;
	inits = 0;
	for(x=0;x<(*inputp)->npop_sampled;x++) {
		if(x < npopa) inits += (*inputp)->config[x];
		else {
			nsam = (*inputp)->config[x];
			break;
		}
	}
    if(nsam == 0) return(-10000.);
    
    ZnA = (double)0;
    j = s0;
    comb = 0;
    while(j+1 < (long)s1) {
        while(j < (long)s1) {
            if(ispolnomhit(j,inits,nsam,(*inputp)->nsam, list, posit) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)s1) {
            if(ispolnomhit(k,inits,nsam,(*inputp)->nsam, list, posit) > 0) break;
            else k++;
        }
        if(k < (long)s1) {
            /*calcular freqs p1,q1*/
            vala = valb = 1;
            a = list[inits][j];
            b = list[inits][k];
            for(i=1+inits;i<inits+nsam;i++) {
                if(list[i][j] == a) vala++;
                else na = list[i][j];
                if(list[i][k] == b) valb++;                                
                else nb = list[i][k];
            }
            if(nsam - vala > vala) {
                a = na;
                vala = nsam - vala;
            }
            if(nsam - valb > valb) {
                b = nb;
                valb = nsam - valb;
            }
            /*calcular p1q1*/
            val00 = 0;
            for(i=inits+0;i<inits+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
            
            /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
            A = (double)vala/(double)nsam;
            B = (double)valb/(double)nsam;
            C = (double)val00/(double)nsam;
            ZnA += ((C - (A*B)) * (C - (A*B))) / (A*(1.0 - A)*B*(1.0 - B));
            comb++;
        }
        j++;
    }
    if(comb) ZnA = ZnA/(double)comb;
    else return(-10000.);
    return(ZnA);
}

int print_matrix_sfixalltheta(struct var **data,FILE *file_thetap,FILE *file_ttotp,struct prob_par **postp) 
{
	int x;
	long int j;
	
	x=0;
	
    if((*data)->ifgamma == 1 || (*data)->range_thetant) {
		for(x=0;x<(*data)->n_loci;x++) {                                
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_thetap,"thetap[locus=%d]\t",x);
			fprintf(file_thetap,"Prob[locus=%d]\t",x);
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_ttotp,"Ttot[locus=%d]\t",x);
			fprintf(file_ttotp,"Prob[locus=%d]\t",x);
			#endif
		}
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fprintf(file_thetap,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fprintf(file_ttotp,"\n");
		#endif
		for(j=0;j<(*data)->n_iter;j++) {                                
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				if(postp[x][j].thetap > -9999)
					fprintf(file_thetap,"%G\t",postp[x][j].thetap);
				else
					fputs("na\t",file_thetap);
				fprintf(file_thetap,"%f\t",postp[x][j].prob);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				if(postp[x][j].Ttotp > -9999)
					fprintf(file_ttotp,"%G\t",postp[x][j].Ttotp);
				else
					fputs("na\t",file_ttotp);
				fprintf(file_ttotp,"%f\t",postp[x][j].prob);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_thetap,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_ttotp,"\n");
			#endif
		}
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fprintf(file_thetap,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fprintf(file_ttotp,"\n");
		#endif
		if((*data)->sfix_allthetas > 0) {
			if((*data)->mhits == 0) 
				fprintf(file_thetap,"In case mhits = 0, the rejection ratio counts for appropriate trees using a fix S.\n");
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				fprintf(file_thetap,"ratio[locus=%d]\t",x);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				fprintf(file_ttotp,"ratio[locus=%d]\t",x);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_thetap,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_ttotp,"\n");
			#endif
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				if(postp[x][(*data)->n_iter].thetap > -9999)
					fprintf(file_thetap,"%G\t",postp[x][(*data)->n_iter].thetap);
				else
					fputs("na\t",file_thetap);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				if(postp[x][(*data)->n_iter].Ttotp > -9999)
					fprintf(file_ttotp,"%G\t",postp[x][(*data)->n_iter].Ttotp);
				else
					fputs("na\t",file_ttotp);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_thetap,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_ttotp,"\n");
			#endif
		}
	}
	return 0;
}
int print_matrix_rmfix(struct var **data,FILE *file_recp,FILE *file_ttotp,int Sfix,struct prob_par **postp)
{
	int x;
	long int j;
	
    if((*data)->ifgammar == 1 || (*data)->range_rnt) {
		for(x=0;x<(*data)->n_loci;x++) {                                
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_recp,"recp[locus=%d]\t",x);
			fprintf(file_recp,"Prob[locus=%d]\t",x);
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			if(Sfix == 0) {
				fprintf(file_ttotp,"Ttot[locus=%d]\t",x);
				fprintf(file_ttotp,"Prob[locus=%d]\t",x);
			}
			#endif
		}
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fprintf(file_recp,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		if(Sfix == 0) fprintf(file_ttotp,"\n");
		#endif
		for(j=0;j<(*data)->n_iter;j++) {                                
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				if(postp[x][j].recp > -9999)
					fprintf(file_recp,"%G\t",postp[x][j].recp);
				else
					fputs("na\t",file_recp);
				fprintf(file_recp,"%f\t",postp[x][j].prob);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				if(Sfix == 0) {
					if(postp[x][j].Ttotp > -9999) 
						fprintf(file_ttotp,"%G\t",postp[x][j].Ttotp);
					else
						fputs("na\t",file_ttotp);	
					fprintf(file_ttotp,"%f\t",postp[x][j].prob);
				}
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_recp,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			if(Sfix == 0) fprintf(file_ttotp,"\n");
			#endif
		}
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fprintf(file_recp,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		if(Sfix == 0) fprintf(file_ttotp,"\n");
		#endif
		if((*data)->rmfix > 0) {
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				fprintf(file_recp,"ratio[locus=%d]\t",x);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				if(Sfix == 0) fprintf(file_ttotp,"ratio[locus=%d]\t",x);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_recp,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			if(Sfix == 0) fprintf(file_ttotp,"\n");
			#endif
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				if(postp[x][(*data)->n_iter].recp > -9999)
					fprintf(file_recp,"%G\t",postp[x][(*data)->n_iter].recp);
				else
					fputs("na\t",file_recp);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				if(Sfix == 0) {
					if(postp[x][(*data)->n_iter].Ttotp > -9999)
						fprintf(file_ttotp,"%G\t",postp[x][(*data)->n_iter].Ttotp);
					else
						fputs("na\t",file_ttotp);
				}
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_recp,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			if(Sfix == 0) fprintf(file_ttotp,"\n");
			#endif
		}
	}
	return 0;
}

/*compare two double numbers in a long int list*/
int compare_(const void *i,const void *j) 
{
    if(*(double *)i < *(double *)j) return -1;
    if(*(double *)i > *(double *)j) return  1;
    return 0;
}

/*compare two numbers*/
int compare(const void *i,const void *j) 
{
    return (*((int *)i) - (*(int *)j));
}

double logPPoisson2(long int Si, double lambda)
{
    double value;
	double factln(long int);
    
	value = ((double)Si*(double)log((double)lambda) - lambda - factln(Si));
    return value;
}


/*Wall's program for calculating minimum recombination events: modification*/
int Min_rec(int x, int segsit, int nsam, int inits,int totalsam)
{  /* Calculate min # rec. events */
  int a, b, c, e, gtest, flag = 0;
  int h;
  int t11,t12,t21,t22;

	if (segsit<2 || x >= (segsit-1)) return (0);
	
	for (a=x+1; a<segsit; ++a) {
		while(a < segsit) {
			if((h=ispolnomhit((long int)a,inits,nsam,totalsam, list, posit)) > 0) break; /*h is the frequency of the new mutation*/
			else a++;
		}            
		if(a < segsit) {
			for (b=x; b<a; ++b) {
				while(b < a) {
					if((h=ispolnomhit((long int)b,inits,nsam,totalsam, list, posit)) > 0) break; /*h is the frequency of the new mutation*/
					else b++;
				}
				if(b < a) {
					t21 = t22 = list[inits][b];
					t11 = t12 = list[inits][a];
					for (e=inits+1; e<inits+nsam; ++e) {
						if(list[e][b] != t21) t22 = list[e][b];
						if(list[e][a] != t11) t12 = list[e][a];
					}
					
					gtest = 0;
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t21 && list[e][a] == t11) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t21 && list[e][a] == t12) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t22 && list[e][a] == t11) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t22 && list[e][a] == t12) {
							++gtest;
							break;
						}
					}
					if (gtest == 4) {
						flag = 1;
						break;
					}
				}
			}
		}
		if (flag == 1) break;
	}
	if (a >= segsit) return (0);
	else {
		c = Min_rec(a,segsit,nsam,inits,totalsam);
		return (1+c);
	}
	return 0;
}

int do_heter_gamma_sitesrec(double *categories,double gammashape,double poppar,long int nsites) 
{
	/*do a cummulative vector*/
	long int i;
	double gammadist(double);
	
	categories[0] = (double)0;
	
	if(gammashape <= (double)0) {/*we do simply all equal*/
		/*categories[0] = poppar/(double)nsites;*/
		for (i=1;i<nsites;i++) {
			categories[i] = (double)poppar/(double)(nsites-1);
			categories[i] += categories[i-1];
		}
	}
	else {/*gamma distribution*/
		/*categories[0] = gammadist(gammashape)/gammashape;
		categories[0] *= poppar/(double)nsites;*/
		for (i=1;i<nsites;i++) {
			categories[i] = (double)gammadist(gammashape)/(double)gammashape;
			categories[i] *= (double)poppar/(double)(nsites-1);
			categories[i] += categories[i-1];
		}
	}
	return 0;
}

int do_heter_gamma_sites(double *categories,double gammashape,double poppar,long int nsites,long int invariable) 
{
	/*do a cummulative vector*/
	long int i;
	double gammadist(double),newpoppar;
	double ran1(void);
		
	newpoppar = poppar * (double)nsites/(double)(nsites-invariable);
	
	if(gammashape <= (double)0) {/*we do simply all equal*/
		if(invariable > 0) {
			if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[0] = (double)0;
			else categories[0] = (double)newpoppar/(double)nsites;
		}
		else categories[0] = (double)newpoppar/(double)nsites;
		for (i=1;i<nsites;i++) {
			if(invariable > 0) {
				if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[i] = (double)0;
				else categories[i] = (double)newpoppar/(double)nsites;
			}
			else categories[i] = (double)newpoppar/(double)nsites;
			categories[i] += categories[i-1];
		}
	}
	else {/*gamma distribution*/
		if(invariable > 0) {
			if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[0] = (double)0;
			else {
				categories[0] = (double)gammadist(gammashape)/(double)gammashape;
				categories[0] *= (double)newpoppar/(double)nsites;
			}
		}
		else {
			categories[0] = (double)gammadist(gammashape)/(double)gammashape;
			categories[0] *= (double)newpoppar/(double)nsites;
		}
		for (i=1;i<nsites;i++) {
			if(invariable > 0) {
				if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[i] = (double)0;
				else {
					categories[i] = (double)gammadist(gammashape)/(double)gammashape;
					categories[i] *= (double)newpoppar/(double)nsites;
				}
			}
			else {
				categories[i] = (double)gammadist(gammashape)/(double)gammashape;
				categories[i] *= (double)newpoppar/(double)nsites;
			}
			categories[i] += categories[i-1];
		}
	}
	
	return 0;
}

long int localize_positiontop(double *categories,double valuer,long int start,long int end) /*es molt lent probablement*/
{
	long int half;
	
	half = (long int)floor((double)(start+end)/(double)2);
	/**/
	while(half != start) {
		if((double)valuer < categories[half]) end = half;
		else if((double)valuer > categories[half]) start = half;
		half = (long int)floor((double)(start+end)/(double)2);
	}
	/**/
	if(half == start) {
		if((double)valuer < categories[half]) return half;
		else return half+1;
	}	

	return half;
}

double correction_recabs(double f,double sexratio,int m) /*ABSOLUTE value of recombination: This calculation is NOT important*/
{
	#if SEXRATIOA1 == 1 && SEXRATIOX1 == -1
		/*is defined by the value of 4N when the sexratio is 1.0*/
		if(f == (double)1.0) {
			if(m==1) {
				return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));
			}
			else { 
				return((sexratio/(1.+sexratio))/2.);
			}
		}
		if(f == (double)0.75) {
			return((sexratio/(1.+sexratio))/2.);
		}
		if(f == (double)0.25) 
			return(0.);
		if(f == (double)-0.25) /*mythochondrial*/
			return(0.);
	#elif SEXRATIOA1 == 0 && SEXRATIOX1 == -1 
		/*is defined by the value of 4N when the sexratio is whatever defined (factor=1 always for A)*/
		/*I am not sure*/
		if(f == (double)1.0) {
			if(m==1) {
				return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));
			}
			else { 
				return((sexratio/(1.+sexratio))/2.);
			}
		}
		if(f == (double)0.75) {
			return((sexratio/(1.+sexratio))/2.);
		}
		if(f == (double)0.25) 
			return(0.);
		if(f == (double)-0.25) /*mythochondrial*/
			return(0.);
	#elif SEXRATIOX1 == 1 && SEXRATIOA1 == -1
		/*is defined by the value of 3N when the sexratio is 1.0*/
		if(f == (double)1.33) {
			if(m==1) {
				return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));
			}
			else { 
				return((sexratio/(1.+sexratio))/2.);
			}
		}
		if(f == (double)1.00) {
			return((sexratio/(1.+sexratio))/2.);
		}
		if(f == (double)0.33) 
			return(0.);
		if(f == (double)-0.33) /*mythochondrial*/
			return(0.);
	#elif SEXRATIOX1 == 0 && SEXRATIOA1 == -1
		/*is defined by the value of 3N when the sexratio is whatever defined (factor=1 always for X)*/
		/*I am not sure*/
		if(f == (double)1.33) {
			if(m==1) {
				return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));
			}
			else { 
				return((sexratio/(1.+sexratio))/2.);
			}
		}
		if(f == (double)1.00) {
			return((sexratio/(1.+sexratio))/2.);
		}
		if(f == (double)0.33) 
			return(0.);
		if(f == (double)-0.33) /*mythochondrial*/
			return(0.);
	#endif		
	return 1.;
}
