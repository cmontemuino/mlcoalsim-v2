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
		if(valuer < categories[half]) end = half;
		else if(valuer > categories[half]) start = half;
		half = (long int)floor((double)(start+end)/2.0);
	}

	if(half == start) {
		if(valuer < categories[half]) return half;
		else return half+1;
	}	

	return half;
}

/* localitza les mutacions en un fragment */
void locate(long int n, long int beg, long int *ptr, int mhits, double *weightmut, long int end)
{
     ordran(n, ptr, mhits, beg, weightmut, end);	/* mutacions en [0,len) ordenades de major a menor */
}

/* localitza les mutacions en un fragment */
void locate2(long int n, long int beg, long int *ptr, int mhits, int nsam, double *weightmut, long int end, char **mutations_matrix)
{    
     ordran2(n,ptr,mhits,nsam,beg,weightmut,end, mutations_matrix);	/* mutacions en [0,len) ordenades de major a menor */
}

double logPPoisson2(long int Si, double lambda)
{
    double value;
    
	value = Si*log(lambda) - lambda - factln(Si);
    return value;
}

/* localitza les mutacions en un fragment */
void locate_psel(long int n, long int beg, long int *ptr, int mhits, long int sel_nt, double *weightmut, long int end)
{    
     ordran_psel(n,ptr,mhits,sel_nt,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
}

/* localitza les mutacions en un fragment */
void locate2_psel(long int n, long int beg, long int *ptr, int mhits, int nsam, long int sel_nt, double *weightmut, long int end, char **mutations_matrix)
{     
     ordran2_psel(n,ptr,mhits,nsam,sel_nt,beg,weightmut,end, mutations_matrix);	/* mutacions en [0,len) ordenades de major a menor */
}

void order(long int n, long int *pbuf)/* ordena els valors: es molt lent per un gran nombre de polimofismes! */
{
    long int gap,i,j;
    long int temp;
    
    for(gap= n/2; gap>0;gap /= 2)
        for(i=gap;i<n;i++)
            for(j=i-gap;j>=0 && pbuf[j]>pbuf[j+gap];j -= gap) {
                temp = pbuf[j];
                pbuf[j] = pbuf[j+gap];
                pbuf[j+gap] = temp;
            }
}

/* posa un nombre entre [0,len) */
void ordran(long int n, long int *pbuf, int mhits, long int beg, double *weightmut, long int end)
{
    ranvec(n,pbuf,mhits,beg,weightmut,end);
    order(n,pbuf);
}

void ordran2(long int n, long int *pbuf, int mhits, int nsam, long int beg, double *weightmut, long int end, char **mutations_matrix)
{
    /* posa un nombre entre [0,len) */

    ranvec2(n,pbuf,mhits,nsam,beg,weightmut,end, mutations_matrix);
    order(n,pbuf);
}


void ordran_psel(long int n, long int *pbuf, int mhits, long int sel_nt, long int beg, double *weightmut, long int end)
{
    ranvec_psel(n,pbuf,mhits,sel_nt,beg,weightmut,end);
    order(n,pbuf);
}

void ordran2_psel(long int n,long int *pbuf,int mhits,int nsam,long int sel_nt,long int beg,double *weightmut,long int end, char **mutations_matrix)
{
    ranvec2_psel(n,pbuf,mhits,nsam,sel_nt,beg,weightmut,end, mutations_matrix);
    order(n,pbuf);
}


void ranvec(long int n, long int *pbuf, int mhits, long int beg, double *weightmut, long int end) /* posa un nombre entre [0,len) */
{
    long int i,x;
	double valuer,wstartm1;
    
    
	if(beg==0) wstartm1 = 0.0;
	else wstartm1 = weightmut[beg-1];
    for(i=0;i<n;i++) {
        valuer = ran1()*(weightmut[end] - wstartm1) + wstartm1;
		pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
        if(!mhits) {	/* bucle per no mhits (mhits=0) */
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

/* posa un nombre entre [0,len) */
void ranvec2(long int n, long int *pbuf, int mhits, int nsam, long int beg, double *weightmut, long int end, char **mutations_matrix)
{
    long int i,x;
    int y;
    int z,f;
	
    double dlen;
    double a;
	double valuer,wstartm1;
    
	if(beg==0) wstartm1 = 0.0;
	else wstartm1 = weightmut[beg-1];
    if(mhits) {
        dlen = (weightmut[end] - wstartm1);
        for(i=0;i<n;i++) {
            a = ran1()*dlen + wstartm1;
            pbuf[i] = localize_positiontop(weightmut,a,beg,end+1);   
            x = i-1;
            y = 1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    y++;
                    if(y == 3 || y == nsam) {/*no mes de 4 nt per posicio*/
                        a = ran1()*dlen + wstartm1;
                        pbuf[i] = localize_positiontop(weightmut,a,beg,end+1);   
                        x = i;
                        y = 1;
                    }
					f=0;
					for(z=0;z<nsam;z++) /*all mutations should be observed. equal pattern in different position*/
						if((mutations_matrix[z][i] == '0' && mutations_matrix[z][x] == '0') ||
						   (mutations_matrix[z][i] != '0' && mutations_matrix[z][x] != '0')) 
							f++;
					if(f==nsam) {
						a = ran1()*dlen + wstartm1;
						pbuf[i] = localize_positiontop(weightmut,a,beg,end+1);   
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
        for(i=0;i<n;i++) {
			valuer = ran1()*(weightmut[end] - wstartm1) + wstartm1;
			pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

/* posa un nombre entre [0,len) */
void ranvec_psel(long int n, long int *pbuf, int mhits, long int sel_nt, long int beg, double *weightmut, long int end)
{
    long int i,x;
	double valuer,wstartm1;
    
	/*include nsites*/
    pbuf[0] = sel_nt;
	if(beg==0) wstartm1 = 0.0;
	else wstartm1 = weightmut[beg-1];
    for(i=1;i<n;i++) {
        valuer = ran1()*(weightmut[end] - wstartm1) + wstartm1;
		pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
		if(pbuf[i] == pbuf[0]) {
			i--;
			continue;
		}   
        if(!mhits) {	/* bucle per no mhits (mhits=0) */
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

/* posa un nombre entre [0,len) */
void ranvec2_psel(long int n,long int *pbuf,int mhits,int nsam,long int sel_nt,long int beg,double *weightmut,long int end, char **mutations_matrix)
{
    long int i,x;
    int y;
    int z,f;
    double dlen;
    double a;
	double valuer,wstartm1;
    
	if(beg==0) wstartm1 = 0.0;
	else wstartm1 = weightmut[beg-1];
    if(mhits) {
        dlen = weightmut[end] - wstartm1;
        for(i=0;i<n;i++) {
            if(i) {
				a = ran1()*dlen + wstartm1;
				pbuf[i] = localize_positiontop(weightmut,a,beg,end+1);
			}
			else pbuf[0] = sel_nt;
			
            x = i-1;
            y = 1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    y++;
                    if(y == 3 || y == nsam || x == 0) {/*no mes de 4 nt per posicio*/
                        a = ran1()*dlen + wstartm1;
                        pbuf[i] = localize_positiontop(weightmut,a,beg,end+1);   
                        x = i;
                        y = 1;
                    }
					f=0;
					for(z=0;z<nsam;z++) /*all mutations should be observed. equal pattern in different position*/
						if(mutations_matrix[z][pbuf[i]] == mutations_matrix[z][pbuf[x]]) f++;
					if(f==nsam) {
						a = ran1()*dlen + wstartm1;
						pbuf[i] = localize_positiontop(weightmut,a,beg,end+1);  
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
		pbuf[0] = sel_nt;
        for(i=1;i<n;i++) {
			valuer = ran1()*(weightmut[end] - wstartm1) + wstartm1;
			pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}
