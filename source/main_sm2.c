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
#include <math.h>
#include <time.h>

#if inMPI == 1 
	#include <mpi.h>
#endif

int main (int argc,char *argv[])
{	    
    void input_data(FILE *,struct var **,struct var_priors **);
    void output_data(FILE *, char *,char *, struct var **,struct var_priors *);
    void getpars_fix(struct var **, struct var2 **);
    void getpars_mod(struct var **, struct var2 **,int);
    void free_getpars_fix(struct var **, struct var2 **);
    void free_inputdata(struct var **,struct var_priors *,int);
	int print_matrix_sfixalltheta(struct var **,FILE *,FILE *,struct prob_par **);
	int print_matrix_rmfix(struct var **,FILE *,FILE * /*,FILE * */,int,struct prob_par **);
	int make_priord(struct var_priors *,long int);
	double ran1();
    
    int print_neuttest(struct var **,FILE *,char *,char *,double **,struct prob_par **);
    int ms(struct var2 **,char *,double **,struct prob_par **,long int,long int,long int *,long int *, long int *,struct var_priors *,int,long int);    
	void init_seed1(long int);
	
    struct var *data;
    struct var2 *inputp;
    FILE *file_input;
    FILE *file_output=0;
	
	char file_in[400];
    char file_out[400];
    int x,y;
	
	struct var_priors *priors;
	
	FILE *file_thetap=0;
	FILE *file_ttotp=0;
	FILE *file_recp=0;
	
	char namefile_[420];
	char *el;

#if TIMEDURATION
    time_t start;
    time_t end;
    struct tm *date;
    char s[80];
    int hour,min,sec;
#endif

	struct prob_par **postp = NULL; /*posterior probabilities for theta and rec*/
	double **matrix_test = NULL;/*all results*/
	/*duplicate variables for an easy recover using parallel*/
	struct prob_par **postp2=NULL;
	double **matrix_test2=NULL;
	
	long int *listnumbers=NULL;
	long int jcount,jcount2,mcount2;
    double jc2=0.;
    double mc2=0.;
	int windows;
	/**/
	long int x0,x1;
	int totalnloci,n,m;
	/**/
	/*long int xx;*/

	int a,b,c,d;

	char names[NEUTVALUES2][20]  = {{"TD"},{"Fs"},{"FDn"},{"FFn"},{"FD"},{"FF"},{"H"},
		{"B"},{"Q"},{"ZnA"},{"Fst"},{"Kw"},{"Hw"},{"R2"},
		{"S"},{"piw"},{"pib"},{"ThetaW"},{"ThetaT"},
		{"ThetaFW"},{"D_Dmin"},{"H_norm"},{"maxhap"},{"maxhap1"},{"Rm"},
		{"ThetaFL"},{"ThetaL"},{"ZengE"},{"EW"},{"Fstw"},{"minuiHS"},{"maxuiHS"},{"maxliES"},
		{"fixoutg"},{"KoutgJC"},{"ThetaWA"},{"ThetaWAn"},{"ThetaTA"},{"ThetaTAn"},{"AchY"},{"AchYn"},
        {"m_sdev"},{"m_skew"},{"m_kurt"},{"ragg"},{"ZnS"},{"ZZ"}
	};
	
	/*MPI definitions*/
#if inMPI==1
	int source;
	int nmpiloc=0;
	MPI_Status	status;
	void BuildBcast_derived_data(struct var **,struct var_priors **,int,int *,int **,long int **);
	MPI_Request *recv;
	int *rflag,reqcount,reqmess;
#endif
	void rankprint_inputp(struct var2 **,int,int);
	int *nppr;
	long int *niterpr;
	int maxnppr;
	int npmpi;
	int my_rank;
	
	/***************************** STARTING PROGRAM *******************************************/
#if inMPI == 1	
	/*if MPI is active*/
	MPI_Init(&argc,&argv);/*start up MPI*/
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);/*find out process rank*/
	MPI_Comm_size(MPI_COMM_WORLD,&npmpi);/*find out number of processes*/
#endif

	if(!(nppr = (int *)calloc(2,sizeof(int)))) {
		perror("calloc error ms.097");
		exit(1);
	}
	maxnppr = 1;
	if(!(niterpr = (long int *)calloc(2,sizeof(long int)))) {
		perror("calloc error ms.098");
		exit(1);
	}

#if inMPI == 1
	if(my_rank == 0) { 
		/*WARNING: this is open for inMPI == 1!!*/
#endif		
		/*if inMPI == 1 && my_rank == 0 or if inMPI == 0*/
		#ifndef __DNASP__
		puts(MLCOALSIM);
		#endif
	/******************************* INPUT/OUTPUT FILES ***************************************/
		if(argc != 1 && argc != 3) {
			puts("usage: input_file output_file");
			exit(1);
		}
		if(argc == 3) {
			strcpy(file_in,argv[1]);
			strcpy(file_out,argv[2]);
		}    
		if(argc == 1) {
			puts("\nInput file?");
			scanf("%s",file_in);
		}

		if (!(file_input = fopen (file_in,"r"))) perror("Error in input/output");
		input_data(file_input,&data,&priors);
		fclose(file_input);

		if((*data).n_iter == 0) {
			printf("\nInput file contained apparently NO ERRORS. Exit.\n");
			exit(0);
		}
		
		if(argc == 1) {
			puts("\nOutput file?");
			scanf("%s",file_out);
			if(strrchr(file_out,'.') == 0) {
				printf("\n Sorry. The output file must contain an extension. \n");
				exit(1);
			}
		}

#if TIMEDURATION
		/*Starting time*/
		time(&start);
		date = localtime(&start);
		strftime(s,80,"%c",date);
#endif

		/*TOTALNLOCI*/
		if((*data).linked == 1 && (*data).window && (*data).despl) {
			totalnloci = (int)ceil(((double)(*data).nsites[1] - (double)(*data).window) / ((double)(*data).despl)) + (int)1;
		}	
		else {
			if((*data).linked > 1) totalnloci = (*data).linked;
			else totalnloci = (*data).n_loci;
		}				
		/********************* CALCULATION OF PRIOR DISTRIBUTIONS ***********************/
		for(x=0;x<(*data).npriors;x++) {
			if(make_priord(priors+x,(*data).n_iter) == 0) {
				printf("Error: mlcoalsim was unable to calculate/read prior distributions.\n");
				exit(1);
			}
		}			
		/*********************** CREATING NAMES FOR ADDITIONAL FILES *****************************/
		if((*data).neutral_tests == 0 && (*data).likelihood_line == 0) {
			if(!(file_output = fopen (file_out,"w"))) perror("Error in input/output");
			fputs(MLCOALSIM,file_output);
			output_data(file_output,file_in,file_out,&data,priors);
			fclose(file_output);
		}
		else {
			if((*data).likelihood_line == 0) {
				if((*data).linked == 0 && (*data).print_neuttest > 2 && (*data).print_neuttest <5 && (*data).n_loci >= 1) {
					 namefile_[0] = '\0';
					 strncat(namefile_,file_out,420);
					 el = strrchr(namefile_,'.');
					 *el = '\0';
					 strncat(namefile_,"_summary.out",420);
					 if(!(file_output = fopen(namefile_,"w"))) perror("Error in input/output");
					 fputs(MLCOALSIM,file_output);
					 output_data(file_output,file_in,file_out,&data,priors);
				 }
				 else {
					if(!(file_output = fopen (file_out,"w"))) perror("Error in input/output");
					fputs(MLCOALSIM,file_output);
					output_data(file_output,file_in,file_out,&data,priors);
				 }
			}
			else {
				if(!(file_output = fopen (file_out,"r"))) {
					if(!(file_output = fopen (file_out,"a+"))) perror("Error in input/output");
					c = NOBS_STATISTICS;
					if((*data).n_iter > 1) {
						fprintf(file_output,"Total\t");
						for(b=0;b<(*data).max_npop_sampled;b++) {
							fprintf(file_output,"Total[pop=%d]\t",(int)b);
						}
						for(a=0;a<c;a++) {
							for(b=0;b<(*data).max_npop_sampled;b++) {
								if((*data).obs_statistics_used[a] == 1) {
									fprintf(file_output,"%s[Total(pop=%d)]\t",names[a],(int)b);
								}
							}
						}
					}
					for(a=0;a<c;a++) {
						for(d=0;d<(*data).max_npop_sampled;d++) {
							if((*data).obs_statistics_used[a] == 1) {
								for(b=0;b<(int)totalnloci;b++) {
									fprintf(file_output,"%s[pop=%d,locus=%d]\t",names[a],(int)d,(int)b);
								}
							}
						}
					}
					fprintf(file_output,"\n");
				}
				else {
					fclose(file_output);
					if(!(file_output = fopen (file_out,"a"))) perror("Error in input/output");
				}
				fflush(file_output);
			}
		}
		/******************************** OUTPUT HEADER ********************************************/
		if((*data).likelihood_line == 0) {
			if((*data).neutral_tests && (*data).print_neuttest != 3 && (*data).print_neuttest != 0) {
				fputs("Neutral tests (excluding multiple hits positions)\n",file_output);
				if(totalnloci > 2 || (*data).linked > 2)
					for(a=0;a<NEUTVALUES2;a++) for(b=0;b<(*data).max_npop_sampled;b++) fprintf(file_output,"avg(%s[pop=%d])\tvar(%s[pop=%d])\t",names[a],b,names[a],b);
				else {
					if(totalnloci < 2 && (*data).linked < 2)
						for(a=0;a<NEUTVALUES2;a++) for(b=0;b<(*data).max_npop_sampled;b++) fprintf(file_output,"value(%s[pop=%d])\t",names[a],b);
					else 
						for(a=0;a<NEUTVALUES2;a++) for(b=0;b<(*data).max_npop_sampled;b++) fprintf(file_output,"avg(%s[pop=%d])\t",names[a],b);
				}
				fputs("\n",file_output); 
			}
		}
		/******************************** OUTPUT Sfixall_thetas ***************************************/
		if(((*data).ifgamma == 1 && (*data).likelihood_line == 0) || ((*data).range_thetant && (*data).likelihood_line == 0)) {
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			namefile_[0] = '\0';
			strncat(namefile_,file_out,420);
			el = strrchr(namefile_,'.');
			*el = '\0';
			strncat(namefile_,"_thetapost.out",420);
			if (!(file_thetap = fopen(namefile_,"w"))) perror("Error in input/output");
			#endif
		}
		/******************************** OUTPUT rmfix ***************************************/
		if(((*data).ifgammar == 1 && (*data).likelihood_line == 0) || ((*data).range_rnt && (*data).likelihood_line == 0)) {
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			namefile_[0] = '\0';
			strncat(namefile_,file_out,420);
			el = strrchr(namefile_,'.');
			*el = '\0';
			strncat(namefile_,"_recpost.out",420);
			if (!(file_recp = fopen(namefile_,"w"))) perror("Error in input/output");
			#endif
		}
		if((((*data).ifgamma == 1  && (*data).likelihood_line == 0) || (((*data).range_thetant)  && (*data).likelihood_line == 0)) ||
		   (((*data).ifgammar == 1  && (*data).likelihood_line == 0) || ((*data).range_rnt && (*data).likelihood_line == 0))) {
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			namefile_[0] = '\0';
			strncat(namefile_,file_out,420);
			el = strrchr(namefile_,'.');
			*el = '\0';
			strncat(namefile_,"_treelength.out",420);
			if (!(file_ttotp = fopen(namefile_,"w"))) perror("Error in input/output");
			#endif
		}
		/**************** COALESCENT SIMULATIONS FOR EACH INDEPENDENT LOCI **************************/
		#if SHOWPROGRESS == 1
		#ifndef __DNASP__
		printf("\n Starting coalescent simulations...");
		/*printf("\n Warning: in case conditioning for mutation parameter, the first dot may take long time.");*/
		printf("\n Each dot indicate that aproximately 2%% of the simulation is done.");
		printf("\n         1    2    3    4    5    6    7    8    9  100%%");
		printf("\n RUN ");
		#endif
		fflush(stdout);
		#endif
		
#if inMPI == 1
	}	
	/*send/receive &data and &priors to/from the rest/root of process/es*/
	BuildBcast_derived_data(&data,&priors,my_rank,&npmpi,&nppr,&niterpr);
	for(n=0;n<npmpi;n++) {
		maxnppr = ((maxnppr > nppr[n])? maxnppr:nppr[n]); 
	}
	MPI_Bcast(file_out,400,MPI_CHAR,0,MPI_COMM_WORLD);
#endif
	/***************** init seed **********************/
	/*init_seed1((*data).seed1[1]);*//*it is already defined in input_data()!*/
	/*********** DEFINING VARIABLES FOR MS **************************************/	
	if(!(inputp = (struct var2 *)calloc(1,sizeof(struct var2)))) perror("calloc error.main.0");		
	getpars_fix(&data,&inputp);
	/*********** define the number of regions *******************/
	if((*inputp).linked == 0)
		windows = (*inputp).tloci;
	else if((*inputp).linked == 1 && (*data).window && (*data).despl) /*sliding windows*/
		windows = (int)ceil(((double)(*data).nsites[1] - (double)(*inputp).window) / ((double)(*inputp).despl)) + 1;
	else if((*inputp).linked > 1)/*separated but linked regions*/
		windows = (*inputp).linked;
	else windows = 1; /*single locus*/
	/************** define vectors for posterior prob in theta/R ***************/
	if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
#if inMPI == 1
		if(my_rank == 0) {
			if(!(postp = (struct prob_par **)calloc((*inputp).tloci,sizeof(struct prob_par *)))) {
				perror("calloc error ms.09");
				exit(1);
			}
			for(n=0;n<(*inputp).tloci;n++) {
				if(!(postp[n] = (struct prob_par *)calloc(((*inputp).howmany + (*inputp).mc_jump + 1),sizeof(struct prob_par)))) {
					perror("calloc error ms.09");
					exit(1);
				}
			}
		}
		else {
			if(!(postp = (struct prob_par **)calloc((unsigned)maxnppr,sizeof(struct prob_par *)))) {
			    perror("calloc error ms.09");
			    exit(1);
			}
			for(n=0;n<maxnppr;n++) {
				if(!(postp[n] = (struct prob_par *)calloc(((*inputp).howmany + (*inputp).mc_jump + 1),sizeof(struct prob_par)))) {
					perror("calloc error ms.09");
					exit(1);
				}
			}
		}
		if(!(postp2 = (struct prob_par **)calloc(1,sizeof(struct prob_par *)))) {
			perror("calloc error ms.09");
			exit(1);
		}
		if(!(postp2[0] = (struct prob_par *)calloc(((*inputp).howmany + (*inputp).mc_jump + 1),sizeof(struct prob_par)))) {
			perror("calloc error ms.09");
			exit(1);
		}
#else
		if(!(postp = (struct prob_par **)calloc((*inputp).tloci,sizeof(struct prob_par *)))) {
			perror("calloc error ms.09");
			exit(1);
		}
		if(!(postp2 = (struct prob_par **)calloc(1,sizeof(struct prob_par *)))) {
			perror("calloc error ms.09");
			exit(1);
		}
		for(n=0;n<(*inputp).tloci;n++) {
			if(!(postp[n] = (struct prob_par *)calloc(((*inputp).howmany + (*inputp).mc_jump + 1),sizeof(struct prob_par)))) {
				perror("calloc error ms.09");
				exit(1);
			}
		}
		if(!(postp2[0] = (struct prob_par *)calloc(((*inputp).howmany + (*inputp).mc_jump + 1),sizeof(struct prob_par)))) {
			perror("calloc error ms.09");
			exit(1);
		}
#endif
	}
	/************** define more structs and matrix ***************/
	if((*inputp).neutral_tests) {		
		/*define matrix*/
#if inMPI == 1
		if(my_rank == 0) {
			if((*data).print_neuttest == 3) {
				matrix_test = NULL;
			}
			else {
				if((matrix_test = (double **)calloc(windows*NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
					perror("calloc error ms.1b");
					exit(1);
				}
			}
			if((*inputp).linked) {
				if((matrix_test2 = (double **)calloc(windows*NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
					perror("calloc error ms.1b");
					exit(1);
				}
			}
			else {
				if((matrix_test2 = (double **)calloc(NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
					perror("calloc error ms.1b");
					exit(1);
				}
			}
			for(n=0;n<windows*NEUTVALUES2*(*data).max_npop_sampled;n++) {
				if((*data).print_neuttest == 3) {
					if((*inputp).linked) {
						if((matrix_test2[n] = (double *)calloc(1,sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
					}
				}
				else {
					if((matrix_test[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
						if((*inputp).linked == 0) {
							printf("Not enough memory to calculate multilocus data: information for each region is shown separately.\n");
							for(x=0;x<n;x++) free(matrix_test[x]);
							free(matrix_test);
							matrix_test = NULL;
							break;
						}
						perror("calloc error ms.1c");
						exit(1);
					}
					if((*inputp).linked) {
						if((matrix_test2[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
					}
				}
			}
			if((*inputp).linked == 0) {
				for(n=0;n<NEUTVALUES2*(*data).max_npop_sampled;n++) {
					if((*data).print_neuttest == 3) {
						if((matrix_test2[n] = (double *)calloc(1,sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
					}
					else {
						if((matrix_test2[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
					}
				}
			}
		}
		else {
			if((*inputp).linked) {
				if((*data).print_neuttest == 3) {
					matrix_test = NULL;
				}
				else {
					if((matrix_test = (double **)calloc(windows*NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
						perror("calloc error ms.1b");
						exit(1);
					}
				}
				if((matrix_test2 = (double **)calloc(windows*NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
					perror("calloc error ms.1b");
					exit(1);
				}
				for(n=0;n<windows*NEUTVALUES2*(*data).max_npop_sampled;n++) {
					if((*data).print_neuttest == 3) {
						if((matrix_test2[n] = (double *)calloc(1,sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
					}
					else {
						if((matrix_test[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
						if((matrix_test2[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
					}
				}
			}
			else {
				if((*data).print_neuttest == 3) {
					matrix_test = NULL;
				}
				else {
					if((matrix_test = (double **)calloc(maxnppr*NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
						perror("calloc error ms.1b");
						exit(1);
					}
				}
				if((matrix_test2 = (double **)calloc(NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
					perror("calloc error ms.1b");
					exit(1);
				}
				for(n=0;n<maxnppr*NEUTVALUES2*(*data).max_npop_sampled;n++) {
					if(!((*data).print_neuttest == 3)) {
						if((matrix_test[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
							if((*inputp).linked == 0) {
								printf("Not enough memory to calculate multilocus data: information for each region is shown separately.\n");
								for(x=0;x<n;x++) free(matrix_test[x]);
								free(matrix_test);
								matrix_test = NULL;
								break;
							}
							perror("calloc error ms.1c");
							exit(1);
						}
					}
				}
				for(n=0;n<NEUTVALUES2*(*data).max_npop_sampled;n++) {
					if((*data).print_neuttest == 3) {
						if((matrix_test2[n] = (double *)calloc(1,sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
					}
					else {
						if((matrix_test2[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
							perror("calloc error ms.1c");
							exit(1);
						}
					}
				}
			}		
		}
#else
		if((*data).print_neuttest == 3) {
			matrix_test = NULL;
		}
		else {
			if((matrix_test = (double **)calloc(windows*NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
				perror("calloc error ms.1b");
				exit(1);
			}
		}
		if((*inputp).linked) {
			if((matrix_test2 = (double **)calloc(windows*NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
				perror("calloc error ms.1b");
				exit(1);
			}
		}
		else {
			if((matrix_test2 = (double **)calloc(NEUTVALUES2*(*data).max_npop_sampled,sizeof(double *))) == NULL) {
				perror("calloc error ms.1b");
				exit(1);
			}
		}		
		for(n=0;n<windows*NEUTVALUES2*(*data).max_npop_sampled;n++) {
			if((*data).print_neuttest == 3) {
				if((*inputp).linked) {
					if((matrix_test2[n] = (double *)calloc(1,sizeof(double))) == NULL) {
						perror("calloc error ms.1c");
						exit(1);
					}
				}
			}
			else {
				if((matrix_test[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
					if((*inputp).linked == 0) {
						printf("Not enough memory to calculate multilocus data: information for each region is shown separately.\n");
						for(x=0;x<n;x++) free(matrix_test[x]);
						free(matrix_test);
						matrix_test = NULL;
						break;
					}
					perror("calloc error ms.1c");
					exit(1);
				}
				if((*inputp).linked) {
					if((matrix_test2[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
						perror("calloc error ms.1c");
						exit(1);
					}
				}
			}
		}
		if((*inputp).linked == 0) {
			for(n=0;n<NEUTVALUES2*(*data).max_npop_sampled;n++) {
				if((*data).print_neuttest == 3) {
					if((matrix_test2[n] = (double *)calloc(1,sizeof(double))) == NULL) {
						perror("calloc error ms.1c");
						exit(1);
					}
				}
				else {
					if((matrix_test2[n] = (double *)calloc(((*inputp).howmany + (*inputp).mc_jump),sizeof(double))) == NULL) {
						perror("calloc error ms.1c");
						exit(1);
					}
				}
			}
		}
#endif		
	}
	/******************************** ******************** **************************************/	
	/******************************** MAIN LOOP (GO TO MS) **************************************/	
	/******************************** ******************** **************************************/	
#if inMPI == 0
	npmpi = 1;
	my_rank = 0;
	nppr[my_rank+1] = (*data).n_loci;
	niterpr[my_rank+1] = (*inputp).howmany + (*inputp).mc_jump;
#endif
	y = 0;
	if((*data).n_loci > 1) {		
		/*In case MHMCMC/MCMCRA, we need to mix values to avoid high correlations. ELIMINATED*/
		if(!(listnumbers = (long int *)malloc((unsigned)((*inputp).howmany + (*inputp).mc_jump)*sizeof(long int)))) {
			perror("calloc error ms.theta_acc");
			exit(1);
		}

		for(x=nppr[my_rank];x<nppr[my_rank+1];x++) {
			/*change the seed for each loci:*/
			/*if(x)*/ init_seed1((*data).seed1[x+1]);/*init_seed1((long int)-(21474*x+(*data).seed1));*//*modify seeds. the user has to give the seed for each locus*/
			getpars_mod(&data,&inputp,x);			

			/*rankprint_inputp(&inputp,my_rank,x);
			fflush(stdout);*/
			/*exit(1);*/
			
			/*****In case MHMCMC/MCMCRA (inactivated).*******/
			for(jcount=0;jcount<(*inputp).howmany + (*inputp).mc_jump;jcount++) {
				if((*data).print_neuttest == 3) listnumbers[jcount] = 0;
				else listnumbers[jcount] = jcount;
			}
			jcount2 = mcount2 = 0;

			/* DO NEUTRALITY TEST, MORE THAN ONE LOCUS*/
			if(ms(&inputp,file_out,matrix_test2,postp2,0,(*inputp).howmany + (*inputp).mc_jump,listnumbers,&jcount2,&mcount2,priors,my_rank,(*data).seed1[x+1])) {
				y = 1;
				break;
			}			 
			if((*data).neutral_tests && matrix_test)
				for(n=0;n<NEUTVALUES2*(*data).max_npop_sampled;n++)
					memcpy(matrix_test[(x-nppr[my_rank])*(NEUTVALUES2*(*data).max_npop_sampled)+n],matrix_test2[0*(NEUTVALUES2*(*data).max_npop_sampled)+n],sizeof(double)*((*inputp).howmany + (*inputp).mc_jump));
			if((*inputp).ifgamma == 1 || (*inputp).range_thetant  || (*inputp).ifgammar == 1 || (*inputp).range_rnt)
				memcpy(postp[x-nppr[my_rank]],postp2[0],sizeof(struct prob_par)*((*inputp).howmany + (*inputp).mc_jump + 1));				

			
			if((*inputp).Sfix_alltheta == 1) 
				postp[(x-nppr[my_rank])][(*inputp).howmany + (*inputp).mc_jump].thetap = postp[(x-nppr[my_rank])][(*inputp).howmany + (*inputp).mc_jump].Ttotp = (double)/*(*inputp).howmany*/mcount2/(double)jcount2;
			if((*inputp).rmfix == 1) 
				postp[(x-nppr[my_rank])][(*inputp).howmany + (*inputp).mc_jump].recp = postp[(x-nppr[my_rank])][(*inputp).howmany + (*inputp).mc_jump].Ttotp = (double)/*(*inputp).howmany*/mcount2/(double)jcount2;
		}
		free(listnumbers);
#if inMPI == 1
		/*COMMUNICATION BETWEEN PROCESSES*/
		if(my_rank == 0) {
			if((*data).print_neuttest == 3) {
				MPI_Finalize(); /*shut down MPI*/
			}
			else {
				/*DO REQUEST FOR NONBLOCKING MPI*/
				if(!(recv = (MPI_Request *)malloc(npmpi*sizeof(MPI_Request)))) {
					perror("calloc error ms.MPI(01)");
					exit(1);
				}
				if(!(rflag = (int *)calloc(npmpi,sizeof(int)))) {
					perror("calloc error ms.MPI(02)");
					exit(1);
				}
				/**/
				for(source=1;source<npmpi;source++) {
					MPI_Irecv(&reqmess,1,MPI_INT,source,12345+source,MPI_COMM_WORLD,&recv[source]);
				}
				reqcount = 0;
				do {
					for(source=1;source<npmpi;source++) {
						if(rflag[source] == 0) {
							MPI_Test(&recv[source],&rflag[source],&status);
							if(rflag[source]) {
								reqcount++;
								nmpiloc = nppr[source+1] - nppr[source];
								if((*data).neutral_tests && matrix_test) {
									for(x=nppr[source];x<nppr[source+1];x++) {
										for(n=0;n<NEUTVALUES2*(*data).max_npop_sampled;n++) {
											MPI_Recv(matrix_test[x*(NEUTVALUES2*(*data).max_npop_sampled)+n],((*inputp).howmany + (*inputp).mc_jump),MPI_DOUBLE,source,x*(NEUTVALUES2*(*data).max_npop_sampled)+n,MPI_COMM_WORLD,&status);
										}
									}
								}
								if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
									for(x=nppr[source];x<nppr[source+1];x++) {
										MPI_Recv(postp[x],4*((*inputp).howmany + (*inputp).mc_jump + 1),MPI_DOUBLE,source,nppr[npmpi]*NEUTVALUES2*(*data).max_npop_sampled+source,MPI_COMM_WORLD,&status);
									}
								}
							}
						}
					}
				}while(reqcount < npmpi-1);
				
				free(recv);
				MPI_Finalize(); /*shut down MPI*/
			}
		}
		else {
			if((*data).print_neuttest == 3) {
				if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
					for(x=0;x<nmpiloc;x++) {
						MPI_Send(postp[x],4*((*inputp).howmany + (*inputp).mc_jump + 1),MPI_DOUBLE,0,nppr[npmpi]*NEUTVALUES2*(*data).max_npop_sampled+my_rank,MPI_COMM_WORLD);
					}
				}
				/*EXIT of the process 'my_rank'*/
				if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
					for(x=0;x<maxnppr;x++) free(postp[x]);
					free(postp);
					free(postp2[0]);
					free(postp2);
				}
				for(x=0;x<NEUTVALUES2*(*data).max_npop_sampled;x++) 
					free(matrix_test2[x]);
				free(matrix_test2);
				free_getpars_fix(&data,&inputp);
				free(inputp);
				free_inputdata(&data,priors,my_rank);
				free(data);
				free(priors);
				free(nppr);
				free(niterpr);
				MPI_Finalize(); /*shut down MPI*/
				return 0; /*exit(0);*/
			}
			else {
				/*flag to rank=0 to indicate the process has finished*/
				MPI_Send(&my_rank,1,MPI_INT,0,12345+my_rank,MPI_COMM_WORLD);
				
				nmpiloc = nppr[my_rank+1]-nppr[my_rank];
				if((*data).neutral_tests && matrix_test) {
					for(x=0;x<nmpiloc;x++) {
						for(n=0;n<NEUTVALUES2*(*data).max_npop_sampled;n++) {
							MPI_Send(matrix_test[x*(NEUTVALUES2*(*data).max_npop_sampled)+n],((*inputp).howmany + (*inputp).mc_jump),MPI_DOUBLE,0,(x+nppr[my_rank])*(NEUTVALUES2*(*data).max_npop_sampled)+n,MPI_COMM_WORLD);			
						}
					}
				}
				if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
					for(x=0;x<nmpiloc;x++) {
						MPI_Send(postp[x],4*((*inputp).howmany + (*inputp).mc_jump + 1),MPI_DOUBLE,0,nppr[npmpi]*NEUTVALUES2*(*data).max_npop_sampled+my_rank,MPI_COMM_WORLD);
					}
				}
				/*EXIT of the process 'my_rank'*/
				if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
					for(x=0;x<maxnppr;x++) free(postp[x]);
					free(postp);
					free(postp2[0]);
					free(postp2);
				}
				
				if((*inputp).neutral_tests && matrix_test) {
					for(x=0;x<maxnppr*NEUTVALUES2*(*data).max_npop_sampled;x++)
						free(matrix_test[x]);
					for(x=0;x<NEUTVALUES2*(*data).max_npop_sampled;x++) 
						free(matrix_test2[x]);
					free(matrix_test2);
					free(matrix_test);
				}
				free_getpars_fix(&data,&inputp);
				free(inputp);
				free_inputdata(&data,priors,my_rank);
				free(data);
				free(priors);
				free(nppr);
				free(niterpr);
				MPI_Finalize(); /*shut down MPI*/
				return 0; /*exit(0);*/
			}
		}
#endif
	}
	else {
		init_seed1((*data).seed1[1]);
        
        /* SINGLE LOCUS*/
        getpars_mod(&data,&inputp,0);
		
		/*****In case MHMCMC/MCMCRA, (inactivated) *******/
		if(!(listnumbers = (long int *)malloc((unsigned)((*inputp).howmany + (*inputp).mc_jump)*sizeof(long int)))) {
			perror("calloc error ms.theta_acc");
			exit(1);
		}
		for(jcount=0;jcount<(*inputp).howmany + (*inputp).mc_jump;jcount++) {
			if((*data).print_neuttest == 3) listnumbers[jcount] = 0;
			else listnumbers[jcount] = jcount;
		}
		
		x1 = niterpr[my_rank+1]-niterpr[my_rank];
		for(x0=niterpr[my_rank];x0<niterpr[my_rank+1];x0+=x1) {
			/*debug*//*
			printf("\nx0=%ld\tx1=%ld\tmy_rank=%d",niterpr[my_rank],niterpr[my_rank+1],my_rank);
			fflush(stdout);
			*/
			jcount2 = mcount2 = 0;
			
			/*rankprint_inputp(&inputp,my_rank,0);*/
			/*printf("\nmy_rank: %d, x0: %ld",my_rank,x0);*/
			/*fflush(stdout);*/
			/*exit(1);*/

			if(ms(&inputp,file_out,matrix_test2,postp2,x0,x0+x1,listnumbers,&jcount2,&mcount2,priors,my_rank,(*data).seed1[1])) {
				y = 1;
				break;
			}
			if(y==1) break;
			/* Assuming that 'listnumbers' are consecutive: !!*/
			if((*data).neutral_tests && matrix_test)
				for(n=0;n<windows;n++)
					for(m=0;m<NEUTVALUES2*(*data).max_npop_sampled;m++)
						memcpy(&(matrix_test[n*(NEUTVALUES2*(*data).max_npop_sampled)+m][listnumbers[x0]]),&(matrix_test2[n*(NEUTVALUES2*(*data).max_npop_sampled)+m][listnumbers[x0]]),sizeof(double)*x1);
			if((*inputp).ifgamma == 1 || (*inputp).range_thetant  || (*inputp).ifgammar == 1 || (*inputp).range_rnt)
				memcpy(&(postp[0][listnumbers[x0]]),&(postp2[0][listnumbers[x0]]),sizeof(struct prob_par)*x1);				
			/*
			for(xx=x0;xx<x0+x1;xx++) {
				printf("\nORIGINAL: %f\t%d",postp[0][xx].thetap,my_rank);
			}
			*/
			/*
			for(xx=x0;xx<x0+x1;xx++) {
				if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
					postp[0][listnumbers[xx]].thetap = postp2[0][listnumbers[xx]].thetap;
					postp[0][listnumbers[xx]].recp = postp2[0][listnumbers[xx]].recp;
					postp[0][listnumbers[xx]].Ttotp = postp2[0][listnumbers[xx]].Ttotp;
				}
				if((*inputp).neutral_tests && matrix_test)
					for(n=0;n<windows;n++)
						for(m=0;m<NEUTVALUES2*(*data).max_npop_sampled;m++)
							matrix_test[n*(NEUTVALUES2*(*data).max_npop_sampled)+m][listnumbers[xx]] = matrix_test2[n*(NEUTVALUES2*(*data).max_npop_sampled)+m][listnumbers[xx]];
			}
			*/
			if((*inputp).Sfix_alltheta == 1 || (*inputp).rmfix == 1) {
				postp[0][(*inputp).howmany + (*inputp).mc_jump].thetap += (double)mcount2; 
				postp[0][(*inputp).howmany + (*inputp).mc_jump].recp  += (double)mcount2;
				postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp  += (double)jcount2;
			}
		}
		free(listnumbers);

#if inMPI == 0
		if((*inputp).ifgamma == 1 || (*inputp).range_thetant  || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
			mc2 = postp[0][(*inputp).howmany + (*inputp).mc_jump].thetap;
			jc2 = postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp;
		}
		if((*inputp).Sfix_alltheta == 1 || (*inputp).rmfix == 1) {
			postp[0][(*inputp).howmany + (*inputp).mc_jump].thetap = (double)/*((*inputp).howmany + (*inputp).mc_jump)*/mc2/(double)jc2;
			postp[0][(*inputp).howmany + (*inputp).mc_jump].recp  = (double)/*((*inputp).howmany + (*inputp).mc_jump)*/mc2/(double)jc2;
			postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp  = (double)/*((*inputp).howmany + (*inputp).mc_jump)*/mc2/(double)jc2;
		}
#else
		/*COMMUNICATION BETWEEN PROCESSES*/		
		/* Assuming that 'listnumbers' are consecutive: !!*/
		if(my_rank == 0) {
			if((*data).print_neuttest == 3) {
				MPI_Finalize(); /*shut down MPI*/
			}
			else {
				/*DO REQUEST FOR NONBLOCKING MPI*/
				if(!(recv = (MPI_Request *)malloc(npmpi*sizeof(MPI_Request)))) {
					perror("calloc error ms.MPI");
					exit(1);
				}
				if(!(rflag = (int *)calloc(npmpi,sizeof(int)))) {
					perror("calloc error ms.MPI(02)");
					exit(1);
				}
				/**/
				for(source=1;source<npmpi;source++) {
					MPI_Irecv(&reqmess,1,MPI_INT,source,12345+source,MPI_COMM_WORLD,&recv[source]);
				}
				reqcount = 0;
				do {
					for(source=1;source<npmpi;source++) {
						if(rflag[source] == 0) {
							MPI_Test(&recv[source],&rflag[source],&status);
							if(rflag[source]) {
								reqcount++;
								x0 = niterpr[source];
								x1 = niterpr[source+1]-niterpr[source];
								if((*data).neutral_tests && matrix_test) {
									for(x=0;x<windows;x++) {
										for(n=0;n<NEUTVALUES2*(*data).max_npop_sampled;n++) {
											MPI_Recv(&(matrix_test[x*(NEUTVALUES2*(*data).max_npop_sampled)+n][listnumbers[x0]]),x1,MPI_DOUBLE,source,x*(NEUTVALUES2*(*data).max_npop_sampled)+n,MPI_COMM_WORLD,&status);
										}
									}
								}
								if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
									MPI_Recv(&(postp[0][listnumbers[x0]]),4*x1,MPI_DOUBLE,source,windows*NEUTVALUES2*(*data).max_npop_sampled,MPI_COMM_WORLD,&status);
									/* 
									 for(xx=x0;xx<x0+x1;xx++) {
									 printf("\nRECEIVED: %f\t%d\tx0:%ld\tx1:%ld\txx:%ld",postp[0][xx].thetap,source,x0,x1,xx);
									 }
									 */
									MPI_Recv(&(jc2),1,MPI_DOUBLE,source,windows*NEUTVALUES2*(*data).max_npop_sampled+1,MPI_COMM_WORLD,&status);
									postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp += jc2;
									if((*inputp).Sfix_alltheta == 1 || (*inputp).rmfix == 1) {
										MPI_Recv(&(mc2),1,MPI_DOUBLE,source,windows*NEUTVALUES2*(*data).max_npop_sampled+2,MPI_COMM_WORLD,&status);
										postp[0][(*inputp).howmany + (*inputp).mc_jump].thetap += mc2;
									}
								}
							}
							if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
								jc2 = postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp;
								if((*inputp).Sfix_alltheta == 1 || (*inputp).rmfix == 1) 
									mc2 = postp[0][(*inputp).howmany + (*inputp).mc_jump].thetap;
								postp[0][(*inputp).howmany + (*inputp).mc_jump].thetap = postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp = (double)/*((*inputp).howmany + (*inputp).mc_jump)*/mc2/(double)jc2;
								postp[0][(*inputp).howmany + (*inputp).mc_jump].recp = postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp = (double)/*((*inputp).howmany + (*inputp).mc_jump)*/mc2/(double)jc2;
							}
						}
					}
				}while(reqcount < npmpi-1);
				
				free(recv);
				MPI_Finalize(); /*shut down MPI*/
			}
		}
		else {
			if((*data).print_neuttest == 3) {
				if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
					MPI_Send(&(postp[0][listnumbers[x0]]),4*x1,MPI_DOUBLE,0,windows*NEUTVALUES2*(*data).max_npop_sampled,MPI_COMM_WORLD);
					MPI_Send(&(postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp),1,MPI_DOUBLE,0,windows*NEUTVALUES2*(*data).max_npop_sampled+1,MPI_COMM_WORLD);
					if((*inputp).Sfix_alltheta == 1 || (*inputp).rmfix == 1)
						MPI_Send(&(postp[0][(*inputp).howmany + (*inputp).mc_jump].thetap),1,MPI_DOUBLE,0,windows*NEUTVALUES2*(*data).max_npop_sampled+2,MPI_COMM_WORLD);
					/*EXIT of the process 'my_rank'*/
					for(x=0;x<windows;x++) free(postp[x]);
					free(postp);
					free(postp2[0]);
					free(postp2);
				}
				if((*inputp).linked) {
					for(x=0;x<windows*NEUTVALUES2*(*data).max_npop_sampled;x++) {
						free(matrix_test2[x]);
					}
				}
				else {
					free(matrix_test2[0]);
				}
				free(matrix_test2);
				free_getpars_fix(&data,&inputp);
				free(inputp);
				free_inputdata(&data,priors,my_rank);
				free(data);
				free(priors);
				free(nppr);
				free(niterpr);
				MPI_Finalize(); /*shut down MPI*/
				return 0; /*exit(0);*/
			}
			else {
				/*flag to rank=0 to indicate the process has finished*/
				MPI_Send(&my_rank,1,MPI_INT,0,12345+my_rank,MPI_COMM_WORLD);

				x0 = niterpr[my_rank];
				x1 = niterpr[my_rank+1]-niterpr[my_rank];
				if((*data).neutral_tests && matrix_test)	{
					for(x=0;x<windows;x++) {
						for(n=0;n<NEUTVALUES2*(*data).max_npop_sampled;n++) {
							MPI_Send(&(matrix_test[x*(NEUTVALUES2*(*data).max_npop_sampled)+n][listnumbers[x0]]),x1,MPI_DOUBLE,0,x*(NEUTVALUES2*(*data).max_npop_sampled)+n,MPI_COMM_WORLD);
						}
					}
				}
				if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
					MPI_Send(&(postp[0][listnumbers[x0]]),4*x1,MPI_DOUBLE,0,windows*NEUTVALUES2*(*data).max_npop_sampled,MPI_COMM_WORLD);
					MPI_Send(&(postp[0][(*inputp).howmany + (*inputp).mc_jump].Ttotp),1,MPI_DOUBLE,0,windows*NEUTVALUES2*(*data).max_npop_sampled+1,MPI_COMM_WORLD);
					if((*inputp).Sfix_alltheta == 1 || (*inputp).rmfix == 1)
						MPI_Send(&(postp[0][(*inputp).howmany + (*inputp).mc_jump].thetap),1,MPI_DOUBLE,0,windows*NEUTVALUES2*(*data).max_npop_sampled+2,MPI_COMM_WORLD);
					/*
					for(xx=x0;xx<x0+x1;xx++) {
						printf("\nSENT: %f\t%d\tx0:%ld\tx1:%ld\txx:%ld",postp[0][xx].thetap,my_rank,x0,x1,xx);
					}
					*/
					
					/*EXIT of the process 'my_rank'*/
					for(x=0;x<windows;x++) free(postp[x]);
					free(postp);
					free(postp2[0]);
					free(postp2);
				}
				if((*inputp).neutral_tests && matrix_test) {
					if((*inputp).linked) {
						for(x=0;x<windows*NEUTVALUES2*(*data).max_npop_sampled;x++) {
							free(matrix_test[x]);
							free(matrix_test2[x]);
						}
					}
					else {
						for(x=0;x<windows*NEUTVALUES2*(*data).max_npop_sampled;x++) 
							free(matrix_test[x]);
						free(matrix_test2[0]);
					}
					free(matrix_test);
					free(matrix_test2);
				}
				free_getpars_fix(&data,&inputp);
				free(inputp);
				free_inputdata(&data,priors,my_rank);
				free(data);
				free(priors);
				free(nppr);
				free(niterpr);
				MPI_Finalize(); /*shut down MPI*/
				return 0; /*exit(0);*/
			}
		}
#endif		
	}

	/*************************** in case inMPI, only the process my_rank = 0 should arrive here! ************************************/
	/******************************** END MAIN LOOP **************************************/
	#if SHOWPROGRESS == 1
	#ifndef __DNASP__
	printf(" Simulation finished.");
	#endif
	fflush(stdout);
	#endif

	#if SHOWPROGRESS == 1
	#ifndef __DNASP__
	printf("\n Saving results in output file/s...");
	#endif
	fflush(stdout);
	#endif

	/******************************** OUTPUT Sfixall_thetas  ************************************/
	if(((*data).ifgamma == 1  && (*data).likelihood_line == 0) || ((*data).range_thetant && (*data).likelihood_line == 0)) {
		if(print_matrix_sfixalltheta(&data,file_thetap,file_ttotp,postp)) return 1;  /*exit*/
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fclose(file_thetap);
		#endif
	}
	/******************************** OUTPUT 4Nr  ************************************/
	if(((*data).ifgammar == 1  && (*data).likelihood_line == 0) || ((*data).range_rnt && (*data).likelihood_line == 0)) {
		if(print_matrix_rmfix(&data,file_recp,file_ttotp,((*data).ifgamma == 1 || (*data).range_thetant),postp)) return 1;  /*exit*/
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fclose(file_recp);
		#endif
	}
	/******************************** OUTPUT trees  ************************************/
	if((((*data).ifgamma == 1 || (*data).range_thetant)  && (*data).likelihood_line == 0) || (((*data).ifgammar == 1 || (*data).range_rnt) && (*data).likelihood_line == 0)) {
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fclose(file_ttotp);
		#endif
	}

	/****************** EXIT IF ERROR, AFTER PRINTING ADDITIONAL VALUES  ************************/
    if(y == 1) return 1; /*exit*/
	/****************************** PRINT STATS IN OUTPUT FILE **********************************/
	if((*data).print_neuttest && (*data).neutral_tests && matrix_test) 
		if(print_neuttest(&data,file_output,file_out,file_in,matrix_test,postp)) return 1; /*exit*/

    /*************************************** END FREE VECTORS ***************************************/
#if TIMEDURATION
    /*time and duration*/
    if((*data).likelihood_line == 0) {
        time(&end);
        date = localtime(&end);
        strftime(s,80,"%c",date);
        fprintf(file_output,"\n\nDate of completion: %s\n",s);
        fprintf(file_output,"\nDuration of the process: %.0f seconds.",difftime(end,start));
        hour =  difftime(end,start)/3600.;
        min  = (difftime(end,start) - (double)hour*3600.)/60.;
        sec  =  difftime(end,start) - (double)hour*3600. - (double)min*60.;
        fprintf(file_output," (%dh:%dm:%ds)\n",hour,min,sec);
    }
#endif
    
    /*************************************** CLOSE FILE ***************************************/
    if(!((*data).neutral_tests == 0 && (*data).likelihood_line == 0))
        fclose(file_output);
	
	/*************************************** FREE VECTORS ***************************************/
	/*in case inMPI, only the process my_rank = 0 arrive here*/	
	if((*inputp).ifgamma == 1 || (*inputp).range_thetant || (*inputp).ifgammar == 1 || (*inputp).range_rnt) {
		for(x=0;x<(*data).n_loci;x++)
			free(postp[x]);
		free(postp);
		free(postp2[0]);
		free(postp2);
	}
	if((*inputp).neutral_tests) {
		for(x=0;x<windows*NEUTVALUES2*(*data).max_npop_sampled;x++) {
			if(matrix_test) free(matrix_test[x]);
			if((*inputp).linked) free(matrix_test2[x]);
		}
		if((*inputp).linked==0) {
			for(x=0;x<NEUTVALUES2*(*data).max_npop_sampled;x++)
				free(matrix_test2[x]);
		}
		if(matrix_test) free(matrix_test);
		free(matrix_test2);
	}
    free_getpars_fix(&data,&inputp);
    free(inputp);
    free_inputdata(&data,priors,my_rank);
    free(data);
	free(priors);
	free(nppr);
	free(niterpr);
	
    /*************************************** EXIT ***************************************/
    #ifndef __DNASP__
	printf("\n %s exited succesfully.\n\n",MLCOALSIM);
	#endif

	return 0;
}
