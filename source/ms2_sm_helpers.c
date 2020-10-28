#include <stdio.h>
#include <math.h>
#include "ran1.h"
#include "distrib.h"

// Expand matrix with mutations
void biggerlist(int nsam, char **mutations_matrix, long int maxsites)
{
    int i;
    for(i=0;i<nsam;i++) {
		if(!(mutations_matrix[i] = realloc(mutations_matrix[i],sizeof *mutations_matrix[i] *maxsites)))
			perror("realloc error. biggerlist");
    }
}

/*Function used to compare numbers in `qsort`*/
int compare_(const void *i,const void *j)
{
    if(*(double *)i < *(double *)j) return -1;
    if(*(double *)i > *(double *)j) return  1;
    return 0;
}

int do_heter_gamma_sites(double *categories,double gammashape,double poppar,long int nsites,long int invariable) 
{
	/*do a cummulative vector*/
	long int i;
	double newpoppar;
		
	newpoppar = poppar * (double)nsites/(double)(nsites-invariable);
	
	if(gammashape <= 0.0) {/*we do simply all equal*/
		if(invariable > 0) {
			if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[0] = 0.0;
			else categories[0] = newpoppar/(double)nsites;
		}
		else categories[0] = newpoppar/(double)nsites;
		for (i=1;i<nsites;i++) {
			if(invariable > 0) {
				if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[i] = 0.0;
				else categories[i] = newpoppar/(double)nsites;
			}
			else categories[i] = newpoppar/(double)nsites;
			categories[i] += categories[i-1];
		}
	}
	else {/*gamma distribution*/
		if(invariable > 0) {
			if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[0] = 0.0;
			else {
				categories[0] = gammadist(gammashape)/gammashape;
				categories[0] *= newpoppar/(double)nsites;
			}
		}
		else {
			categories[0] = gammadist(gammashape)/gammashape;
			categories[0] *= newpoppar/(double)nsites;
		}
		for (i=1;i<nsites;i++) {
			if(invariable > 0) {
				if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[i] = 0.0;
				else {
					categories[i] = gammadist(gammashape)/gammashape;
					categories[i] *= newpoppar/(double)nsites;
				}
			}
			else {
				categories[i] = gammadist(gammashape)/gammashape;
				categories[i] *= newpoppar/(double)nsites;
			}
			categories[i] += categories[i-1];
		}
	}
	
	return 0;
}

int do_heter_gamma_sitesrec(double *categories,double gammashape,double poppar,long int nsites)
{
	/*do a cummulative vector*/
	long int i;
	
	categories[0] = 0.0;
	
	if(gammashape <= 0.0) {/*we do simply all equal*/
		for (i=1;i<nsites;i++) {
			categories[i] = (double)poppar/(double)(nsites-1);
			categories[i] += categories[i-1];
		}
	}
	else {/*gamma distribution*/
		for (i=1;i<nsites;i++) {
			categories[i] = gammadist(gammashape)/gammashape;
			categories[i] *= (double)poppar/(double)(nsites-1);
			categories[i] += categories[i-1];
		}
	}
	return 0;
}

double fixoutg(int nsam, int Sout, int freq0)
{   
    if(nsam) {
        return Sout + freq0;
    } else {
        return -10000;
    }
}

int ispolnomhit(long int j,int init, int nsam, int totalsam, char **mutations_matrix, long int *positions)
{
    int i,h,g;
    int a0,a1;
    
    if(j && positions[j-1] == positions[j]) 
		return(-1);/*estem a la segona mutació o més*/
    
	/*mhit, including outgroup*//**/
    h = g = a1 = 0;
	a0 = '0';
    for(i=0;i<totalsam;i++) {
        if((int)mutations_matrix[i][j] != a0) {
            if(!a1)
				a1 = (int)mutations_matrix[i][j];
            if((int)mutations_matrix[i][j] != a1) 
				return(-2);
        }
    }
	/**/
    h = g = a1 = 0;
	a0 = '0';
    for(i=init;i<init+nsam;i++) {
        if((int)mutations_matrix[i][j] == a0) h++;
        else {
            if(!a1)
				a1 = (int)mutations_matrix[i][j];
            if((int)mutations_matrix[i][j] == a1)
				g++;
            else 
				return(-2);/*en el cas de la primera mutacio hagin almenys tres variants*/
        }
    }
    if(h==nsam || h == 0)
		return(0); /*invariant intrapop*/
    return(nsam-h); /*gives the frequency of the variant that is different of '0'*/
}

double koutgJC(int nsam, int Sout, int *freq, unsigned long nsites)
{
    int x;
    double div;
    
    if(nsam) {
        div = Sout + freq[0];
		
		for(x=1;x<nsam;x++) 
            div += freq[x]*(double)x/(double)nsam;
		
		if(Sout == 0) { /*correction only in case calculating from sequence data*/
			div/= nsites;
			if(div >= .75) return -10000.0;
			else {
				div = -.75*log(1.0-4.0/3.0*div);/*Jukes and Cantor correction...*/
				return div*nsites;
			}
		}
		else return div; /*in case calculating directly Sout, we already included the multiple hits in the branch length*/
    }
    else return -10000.0;
}

long int localize_positiontop(const double *categories,double valuer,long int start,long int end) /*es molt lent probablement*/
{
	long int half;
	
	half = (long int)floor((double)(start+end)/2.0);
	while(half != start) {
		if((double)valuer < categories[half]) end = half;
		else if((double)valuer > categories[half]) start = half;
		half = (long int)floor((double)(start+end)/2.0);
	}

	if(half == start) {
		if((double)valuer < categories[half]) return half;
		else return half+1;
	}	

	return half;
}

double logPPoisson2(long int Si, double lambda)
{
    double value;
    
	value = Si*log(lambda) - lambda - factln(Si);
    return value;
}
