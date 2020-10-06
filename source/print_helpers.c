#include <stdio.h>
#include "mlsp_sm.h"

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
