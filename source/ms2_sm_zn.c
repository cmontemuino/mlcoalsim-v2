#include "neutpar.h"
#include "mlsp_sm.h"
#include "ms2_sm_helpers.h"

double ZnA_(int valuep,long int segsites,struct var2 **inputp, int npopa, char **mutations_matrix, long int *positions)
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
    
    ZnA = 0.0;
    j = 0;
    comb = 0;
    while(j+1 < (long)segsites) {
        while(j < (long)segsites) {
            if(ispolnomhit(j,inits,nsam,(*inputp)->nsam, mutations_matrix, positions) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)segsites) {
            if(ispolnomhit(k,inits,nsam,(*inputp)->nsam, mutations_matrix, positions) > 0) break;
            else k++;
        }
        if(k < (long)segsites) {
            /*calcular freqs p1,q1*/
            vala = valb = 1;
            a = mutations_matrix[inits][j];
            b = mutations_matrix[inits][k];
            for(i=1+inits;i<inits+nsam;i++) {
                if(mutations_matrix[i][j] == a) vala++;
                else na = mutations_matrix[i][j];
                if(mutations_matrix[i][k] == b) valb++;                                
                else nb = mutations_matrix[i][k];
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
            for(i=0+inits;i<inits+nsam;i++) if(mutations_matrix[i][j] == a && mutations_matrix[i][k] == b) val00++;
            
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

double ZnA_window(struct var2 **inputp, long int s0, long int s1, int npopa, char **mutations_matrix, long int *positions)
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
    
    ZnA = 0.0;
    j = s0;
    comb = 0;
    while(j+1 < (long)s1) {
        while(j < (long)s1) {
            if(ispolnomhit(j,inits,nsam,(*inputp)->nsam, mutations_matrix, positions) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)s1) {
            if(ispolnomhit(k,inits,nsam,(*inputp)->nsam, mutations_matrix, positions) > 0) break;
            else k++;
        }
        if(k < (long)s1) {
            /*calcular freqs p1,q1*/
            vala = valb = 1;
            a = mutations_matrix[inits][j];
            b = mutations_matrix[inits][k];
            for(i=1+inits;i<inits+nsam;i++) {
                if(mutations_matrix[i][j] == a) vala++;
                else na = mutations_matrix[i][j];
                if(mutations_matrix[i][k] == b) valb++;                                
                else nb = mutations_matrix[i][k];
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
            for(i=inits+0;i<inits+nsam;i++) if(mutations_matrix[i][j] == a && mutations_matrix[i][k] == b) val00++;
            
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

double Zns(int valuep,long int segsites,struct var2 **inputp,int npopa, char **mutations_matrix, long int *positions)
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
    
    ZnS = 0.0;
    j = 0;
    comb = 0;
    while(j+1 < segsites) {
        while(j < segsites) {
            if(ispolnomhit(j,inits,nsam,(*inputp)->nsam, mutations_matrix, positions) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < segsites) {
            while(k < segsites) {
                if(ispolnomhit(k,inits,nsam,(*inputp)->nsam, mutations_matrix, positions) > 0) break;
                else k++;
            }
            if(k < segsites) {
                /*calcular freqs p1,q1*/
                vala = valb = 1;
                a = mutations_matrix[inits][j];
                b = mutations_matrix[inits][k];
                for(i=1+inits;i<inits+nsam;i++) {
                    if(mutations_matrix[i][j] == a) vala++;
                    else na = mutations_matrix[i][j];
                    if(mutations_matrix[i][k] == b) valb++;                                
                    else nb = mutations_matrix[i][k];
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
                for(i=0+inits;i<inits+nsam;i++) if(mutations_matrix[i][j] == a && mutations_matrix[i][k] == b) val00++;
                
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

double Zns_window(struct var2 **inputp, long int s0, long int s1, int npopa, char **mutations_matrix, long int *positions)
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
	
    ZnS = 0.0;
    j = s0;
    comb = 0;
    while(j+1 < (long)s1) {
        while(j < (long)s1) {
            if(ispolnomhit(j,inits,nsam,(*inputp)->nsam, mutations_matrix, positions) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)s1) {
            while(k < (long)s1) {
                if(ispolnomhit(k,inits,nsam,(*inputp)->nsam, mutations_matrix, positions) > 0) break;
                else k++;
            }
            if(k < (long)s1) {
                /*calcular freqs p1,q1*/
                vala = valb = 1;
                a = mutations_matrix[inits][j];
                b = mutations_matrix[inits][k];
                for(i=1+inits;i<inits+nsam;i++) {
                    if(mutations_matrix[i][j] == a) vala++;
                    else na = mutations_matrix[i][j];
                    if(mutations_matrix[i][k] == b) valb++;                                
                    else nb = mutations_matrix[i][k];
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
                for(i=0+inits;i<inits+nsam;i++) if(mutations_matrix[i][j] == a && mutations_matrix[i][k] == b) val00++;
                
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
