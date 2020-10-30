#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mlsp_sm.h"

void calc_neutparSRH(long int segsit, struct var2 **inputp, struct dnapar *ntpar, double valuer, int npopa, int flaghap, char **mutations_matrix, long int *positions)
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

    comb = (int)((double)nsam*((double)nsam-1.0)/(double)2);
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
		if(valuer) (ntpar)->Rm = Min_rec(0,(int)segsit,nsam,inits,(*inputp)->nsam, mutations_matrix, positions);
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
			if((h=ispolnomhit(j,inits,nsam,(*inputp)->nsam, mutations_matrix, positions)) > 0) break; /*h is the frequency of the new mutation*/
			else j++;
		}            
		if(j<segsit) {
			if(flaghap > 0) {
				for(a=inits;a<inits+nsam;a++) {
					hapl[(a-inits)*segsit+S] = mutations_matrix[a][j];
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

void calc_neutpar_windowSRH(struct var2 **inputp, struct dnapar *ntpar, long int s0, long int s1, double valuer, int npopa, int flaghap, char **mutations_matrix, long int *positions)
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
	
    comb = (int)((double)nsam*((double)nsam-1.0)/2.0);
	for(x=0;x<10;x++) (ntpar)->Sanc[x] = 0;
	
    /*define segsit first*/
    segsit = s1 - s0;

	if(segsit == 0 || nsam < 2) {
		if(nsam == 0) {
			(ntpar)->S = -10000;
			(ntpar)->nhapl = -10000;
			(ntpar)->Rm = -10000;
		}
		else {
			(ntpar)->S = 0;
			(ntpar)->nhapl = 1;
			(ntpar)->Rm = 0;
		}
	}
	else {		
		if(valuer || (*inputp)->mhits) (ntpar)->Rm = Min_rec((int)s0,(int)s1,nsam,inits,(*inputp)->nsam, mutations_matrix, positions);
		else (ntpar)->Rm = 0;
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
			if((h=ispolnomhit(j,inits,nsam,(*inputp)->nsam, mutations_matrix, positions)) > 0) break; /*h is the frequency of the new mutation*/
			else j++;
		}            
		if(j<s1) {
			if(flaghap > 0) {
				for(a=inits;a<inits+nsam;a++) {
					hapl[(a-inits)*segsit+S] = mutations_matrix[a][j];
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

/*Wall's program for calculating minimum recombination events: modification*/
int Min_rec(int x, int segsit, int nsam, int inits, int totalsam, char **mutations_matrix, long int *positions)
{  /* Calculate min # rec. events */
  int a, b, c, e, gtest, flag = 0;
  int h;
  int t11,t12,t21,t22;

	if (segsit<2 || x >= (segsit-1)) return (0);
	
	for (a=x+1; a<segsit; ++a) {
		while(a < segsit) {
			if((h=ispolnomhit((long int)a,inits,nsam,totalsam, mutations_matrix, positions)) > 0) break; /*h is the frequency of the new mutation*/
			else a++;
		}            
		if(a < segsit) {
			for (b=x; b<a; ++b) {
				while(b < a) {
					if((h=ispolnomhit((long int)b,inits,nsam,totalsam, mutations_matrix, positions)) > 0) break; /*h is the frequency of the new mutation*/
					else b++;
				}
				if(b < a) {
					t21 = t22 = mutations_matrix[inits][b];
					t11 = t12 = mutations_matrix[inits][a];
					for (e=inits+1; e<inits+nsam; ++e) {
						if(mutations_matrix[e][b] != t21) t22 = mutations_matrix[e][b];
						if(mutations_matrix[e][a] != t11) t12 = mutations_matrix[e][a];
					}
					
					gtest = 0;
					for (e=inits; e<inits+nsam; ++e) {
						if (mutations_matrix[e][b] == t21 && mutations_matrix[e][a] == t11) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (mutations_matrix[e][b] == t21 && mutations_matrix[e][a] == t12) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (mutations_matrix[e][b] == t22 && mutations_matrix[e][a] == t11) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (mutations_matrix[e][b] == t22 && mutations_matrix[e][a] == t12) {
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
		c = Min_rec(a,segsit,nsam,inits,totalsam, mutations_matrix, positions);
		return (1+c);
	}
	return 0;
}
