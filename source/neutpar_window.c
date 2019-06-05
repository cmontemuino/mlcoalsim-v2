#include "mlsp_sm.h"
#include "neutpar_common.h"
#include "neut_tests.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "zf_log.h"

void calc_neutpar_window(struct var2 **inputp,struct dnapar *ntpar,long int  s0,long int  s1,double valuer,int npopa, char **dna_matrix, long int *posit)
{
    long int segsit;    
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
    
    int inits,inits1,inits2,initso;
	int nsam=0;
	int *initsq1,*initsq2,*sq;
    int val10,val20,val21;
    long int j,k;
    int a,b,c,d,i,comb,comb2,x; int h=0;
    char *hapl;

	int **veca;
	
	#if MAXHAP1
	int mut,mh1,maxhapl1,*sshg;
	#endif
	
	int Min_rec(int,int,int,int,int);
	int pola,polb,polc,pop1,pop2,popo,nt0;
	int *polqa,*polqb;
	int nmh;
	int nsamallpop;
	
	double *freq1,*freq2,*freq3,*freq4,s,s2,g1,g2,moment;
	long int z;
	
	long int ehh0,ehhs,ehhr;

	int npw;
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
	
	#if iEStest == 1 || iHStest == 1
	double testHap(int, int *);
	long int l;
	#endif
	
	long int *Pwd;
	long int maxpwd;
	
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
    /*nsam  = (*inputp)->nsam;*/
    comb = (int)((double)nsam*((double)nsam-1.0)/2.0);
    
    /*define segsit first*/
    segsit = s1 - s0;
    
	if(segsit == 0 || nsam < 2) {
		if(nsam == 0) {
			(ntpar)->B1 = -10000;
			(ntpar)->Q1 = -10000;
			(ntpar)->k = (double)-10000;
			(ntpar)->S = -10000;
			(ntpar)->piw = (double)-10000;
			(ntpar)->pib = (double)-10000;
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
			(ntpar)->thetaL = -10000.0;
			(ntpar)->withinw = -10000.0;
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
			(ntpar)->nhapl = 1;
			for(a=0;a<nsam;a++) {
				(ntpar)->freq[a] = 0;
				(ntpar)->unic[a] = 0;
				(ntpar)->fhapl[a] = 0;
			}
			for(j=0;j<segsit;j++) { /*all valid positions are invariant positions in nsam=1*/
				if((h=ispolnomhit(j,inits,nsam,(*inputp)->nsam, dna_matrix, posit)) == 0) (ntpar)->freq[0] += 1;
			}		
			(ntpar)->fhapl[0] = nsam;
			(ntpar)->piw = (double)0;
			(ntpar)->pib = (double)0;
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
		ZF_LOGD("--> Min_rec(s0=%ld, s1=%ld, nsam=%d, inits=%d, totalsam=%d)", s0, s1, nsam, inits, (*inputp)->nsam);
		if(valuer) (ntpar)->Rm = Min_rec((int)s0,(int)s1,nsam,inits,(*inputp)->nsam);
		else (ntpar)->Rm = (int)0;

        /* calcul B' i Q, pero sense dividir per S (Wall) */
        B = 0;
        A = 0;

        a = 0;
        for(j=s0;j<s1;) {
            k = j;
            while(k+1 < s1) { /*calcular k*/
                if((ispolnomhit(k, inits, nsam, (*inputp)->nsam, dna_matrix, posit)) > 0) break;
                else {
					k++;
				}
            }
            j = k+1;
            while(j < s1) { /*calcular j*/
                if((ispolnomhit(j,inits,nsam,(*inputp)->nsam, dna_matrix, posit)) > 0) break;
                else j++;
            }
            if(j < s1) {                
                val20 = val21 = -1;
                b = 0;
                for(i=inits;i<inits+nsam;i++) {
                    val10 = (dna_matrix[i][k] - 48)*4 + (dna_matrix[i][j]- 48);
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
							if(dna_matrix[i][j] > '0') x = '1';
							else x = '0';
							veca[A][i] = x;
						}
						A += 1;
					}
					else {
						a = 0;
						for(c=0;c<A;c++) {
							d = 0;
							for(i=0;i<inits+nsam;i++) {
								if(dna_matrix[i][j] > '0') x = '1';
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
								if(dna_matrix[i][j] > '0') x = '1';
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
		S = 0;
        nhapl = 0;                
		nmh = 0;
		Se1 = Sn1 = 0;
		k_e1 = k_n1 = 0.;
		
        for(j=s0;j<s1;j++) {
            pi = 0;
            while(j < s1) {
                if((h=ispolnomhit(j,inits,nsam,(*inputp)->nsam, dna_matrix, posit)) > 0) break; /*h is the frequency of the new mutation*/
                else {
					j++;
					if(h == -2) nmh += 1;
					if(h ==  0) (ntpar)->freq[0] += 1; /*invariant intrapop*/
				}
            }            
            if(j<s1) {
                (ntpar)->freq[h] += 1;
                for(a=inits;a<inits+nsam;a++) {
                    hapl[(a-inits)*segsit+S] = dna_matrix[a][j];
                    /*new two lines: for singleton mutations (no outgroup)*/
                    if(h == 1) if(dna_matrix[a][j] != '0') (ntpar)->unic[a-inits] += 1;
                    if(h == nsam-1) if(dna_matrix[a][j] == '0') (ntpar)->unic[a-inits] += 1;
                }
                S++;
				if(h > 1) Se1++;
				if(h > 1 && h < nsam-1) Sn1++;
 				/*pidcount = 0;*/
				z = 0;/*mismatch dist*/
                for(a=inits;a<inits+nsam-1;a++) {
                    for(b=a+1;b<inits+nsam;b++) {
						if(dna_matrix[a][j] != dna_matrix[b][j]) {
							pi++;
							if(h > 1) k_e1++;
							if(h > 1 && h < nsam-1) k_n1++;
							/*(ntpar)->pid[pidcount] += (double)1;*/
							freq1[z] += 1; /*mismatch dist*/
						}
						z++;/*mismatch dist*/
						/*pidcount++;*/
					}
				}
                k_ += pi;
            }
        }
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
			/*Based on Voight et al PLoS 2006 4(3) e72. Calculate hapl div for the SNP (ancestral and derived separately) 
			 across all the region until the value reach 0.05 and sum all together. Divide ancestral by derived 
			 and we take the maximum and the minimum value. Large negative, unusually long derived haplotypes.
			 Positive values indicate unusual large ancestral haplotypes*/
			
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
				ehh0 = s0;
				ehhs = s1;
				(ntpar)->max_uiHS = (ntpar)->min_uiHS = -10000;
			}
			else {
				/*locate the single segregating site centered in (*inputp)->ehh_fixnt +/- (*inputp)->ehh_margin. if no one skip loop (ehh0=0;ehhs=0;(ntpar)->max_uiHS = (ntpar)->min_uiHS = -10000;) */
				ehhr =  (long int)1e7;
				for(j=s0;j<s1;j++) {
					while(j < s1) {
						if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, dna_matrix, posit)) > 0) break;
						else j++;
					}                    
					if(j<s1) {
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
					/*if(flagA == 0 && flagD == 0) {
					 memset(nousedA+j,1,(segsit-j)*sizeof(int));
					 memset(nousedD+j,1,(segsit-j)*sizeof(int));
					 break;
					 }*/
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
				}
				/*integrate and obtain uiHS, calculate max value*/
				uiHS = iHHA = iHHD = (double)0;
				for(j=1;j<segsit;j++) {
					if(nousedA[j]) continue;
					l=j-1;
					while(l>= 0 && nousedA[l]) l--;
					if(l==-1) continue;
					iHHA += ((double)(EHSA[l]+EHSA[j])*(double)(posit[s0+j]-posit[s0+l]))/(double)2.0;
				}
				for(j=1;j<segsit;j++) {
					if(nousedD[j]) continue;
					l=j-1;
					while(l>= 0 && nousedD[l]) l--;
					if(l==-1) continue;
					iHHD += ((double)(EHSD[l]+EHSD[j])*(double)(posit[s0+j]-posit[s0+l]))/(double)2.0;
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
			/*Based on Tang,Thornton,Stoneking PLoS B 2007 5(7)e171 hapl homozigosity between i and j weighted by the homozigosity at i
			  Summarize all values for each single position and calulate the log value. take the max and the min.
			  It is made to see differences between populations*/
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
						if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, dna_matrix, posit)) > 0) break;
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
			/*to look at single position, k must be fixed to the closes polymorphism the user is looking at*/
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
						/*memset(nousedM+k*segsit+j,1,(segsit-j)*sizeof(int));
						 break;*/
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
					iES[k] += ((double)(EHHS[k][l]+EHHS[k][j])*(double)(posit[s0+j]-posit[s0+l]))/(double)2;
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
		
		#if MAXHAP1
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
			for(j=s0;j<s1;j++) {
				while(j < s1) {
					if((c=ispolnomhit(j,0,nsamallpop,nsamallpop, dna_matrix, posit)) > 0) break;
					else j++;
				}                    
				if(j<s1) {
					inits1 = 0;
					for(h=0;h<(*inputp)->npop_sampled;h++) {
						for(a=inits1;a<inits1+(*inputp)->config[h]-1;a++) {
							hapl[(a)*segsit+S] = dna_matrix[a][j];
						}
						/**/
						c = ispolnomhit(j,inits1,(*inputp)->config[h],(*inputp)->nsam, dna_matrix, posit);
						piw[h] += (c * ((*inputp)->config[h] - c));
						/**/
						hapl[(inits1+(*inputp)->config[h]-1)*segsit+S] = dna_matrix[inits1+(*inputp)->config[h]-1][j];
						inits1 += (*inputp)->config[h];
					}
					S++;
					
					inits1 = 0;
					inits2 = inits;
					for(h=0;h<(*inputp)->npop_sampled;h++) {
						if((*inputp)->config[h] > 0 && (*inputp)->config[npopa] > 0 && npopa != h /**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/
						   && ((ispolnomhit(j,inits1,(*inputp)->config[h],(*inputp)->nsam, dna_matrix, posit)>=0) 
						   &&  (ispolnomhit(j,inits2,(*inputp)->config[npopa],(*inputp)->nsam, dna_matrix, posit)>=0))) {
						   for(a=inits1;a<inits1+(*inputp)->config[h];a++)
							   for(b=inits2;b<inits2+(*inputp)->config[npopa];b++)
								   if(dna_matrix[a][j] != dna_matrix[b][j]) pib[h]++;
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
				if(ispolnomhit(j,inits2,(*inputp)->config[npopa],(*inputp)->nsam, dna_matrix, posit) == 0 &&
				   ispolnomhit(j,inits1,(*inputp)->config[(*inputp)->pop_outgroup],(*inputp)->nsam, dna_matrix, posit) == 0) {
					if(dna_matrix[inits2][j] == dna_matrix[inits1][j]) 
						(ntpar)->freq[0] -= 1;
				}
				if(ispolnomhit(j,inits2,(*inputp)->config[npopa],(*inputp)->nsam, dna_matrix, posit) == 0 && 
				   ispolnomhit(j,inits1,(*inputp)->config[(*inputp)->pop_outgroup],(*inputp)->nsam, dna_matrix, posit) < 0) {
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
				if((*inputp)->config[h] > 0 && (*inputp)->config[npopa] > 0 && npopa != h /**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) {
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
				/*else (ntpar)->piwallcomp[h] = -10000.;*/
			}
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h] && (*inputp)->config[npopa] && npopa != h /**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) {
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
				if((*inputp)->config[h]>1) (ntpar)->hapwallcomp[h] = 1.0 - hapw[h] / (double)(((*inputp)->config[h] * ((*inputp)->config[h] - 1))/2);
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
				if((*inputp)->config[h]>1) {
					comb = (int)((double)(*inputp)->config[h] * ((double)(*inputp)->config[h] - (double)1.0)/(double)2.0);
					(ntpar)->piw += (double)piw[h]/(double)comb; 
					comb2++;
				}
			}
			/*despres es divideix pel nombre de subpoblacions*/
			if(comb2) {
				(ntpar)->piw = (ntpar)->piw/(double)comb2;
				npw=0;
				for(h=0;h<(*inputp)->npop_sampled;h++) /*if((*inputp)->config[h]>1) npw +=1;*/ npw +=(*inputp)->config[h]>1;
			}
			
			comb2 = 0;
			
			for(h=0;h<(*inputp)->npop_sampled;h++) {
				if((*inputp)->config[h] && (*inputp)->config[npopa] && npopa != h/**/ && h != (*inputp)->pop_outgroup/*Included to eliminate the outgroup in the analysis*/) {
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
		/*find the value inits1,inits2 and initso: indicate the position in list[inits12o][j] for comparison of polymorphsms*/
		pola = polb = polc = 0; /*if polymorphic in pop1,pop2,popo*/
		
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
		
		for(j=s0;j<s1;j++) {
			while(j < s1) {
				if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, dna_matrix, posit)) > 0) break;
				else j++;
			}                    
			if(j<s1) {					
				if(pola != -1) pola = 0;
				if(polb != -1) polb = 0;
				if(polc != -1) polc = 0; 
				/*if polymorphic in pop1,pop2,popo*/
				
				/*check if each group is polymorphic or not*/
				if(pola != -1 && polb != -1) {
					for(pop1=1;pop1<=(*inputp)->ancestral_pol[0][0];pop1++) {
						nt0 = 0;
						for(h=0;h<(*inputp)->ancestral_pol[0][pop1];h++) nt0 += (*inputp)->config[h];	
						for(h=nt0;h<nt0+(*inputp)->config[(*inputp)->ancestral_pol[0][pop1]];h++) {
							if(dna_matrix[h][j] != dna_matrix[inits1][j]) {
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
							if(dna_matrix[h][j] != dna_matrix[inits2][j]) {
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
							if(dna_matrix[h][j] != dna_matrix[initso][j]) {
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
						if(dna_matrix[inits2][j] == dna_matrix[initso][j]) (ntpar)->Sanc[0] += 1;/*Sx1*/
						else (ntpar)->Sanc[6] += 1;/*Sx1f2*/
					}
					else {
						(ntpar)->Sanc[0] += 1;/*Sx1*/
					}
				}
				if(pola==0 && polb==1 && polc==0) {
					if((*inputp)->ancestral_pol[2][0]) {
						if(dna_matrix[inits1][j] == dna_matrix[initso][j]) (ntpar)->Sanc[1] += 1;/*Sx2*/
						else (ntpar)->Sanc[7] += 1;/*Sx2f1*/
					}
					else {
						(ntpar)->Sanc[1] += 1;/*Sx2*/
					}
				}
				if(polc == 1 && (pola != -1 && polb != -1))/*(*inputp)->ancestral_pol[2][0] == 1*/ {
					if((pola==0 && polb==0) && (dna_matrix[inits1][j] == dna_matrix[inits2][j])) {
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
						if(dna_matrix[inits2][j] == dna_matrix[initso][j] && dna_matrix[inits1][j] != dna_matrix[initso][j]) {
							(ntpar)->Sanc[3] += 1;/*Sf1*/
						}
						if(dna_matrix[inits1][j] == dna_matrix[initso][j] && dna_matrix[inits2][j] != dna_matrix[initso][j]) {
							(ntpar)->Sanc[4] += 1;/*Sf2*/
						}
						if(dna_matrix[inits1][j] == dna_matrix[inits2][j] && dna_matrix[inits1][j] != dna_matrix[initso][j]) {
							(ntpar)->Sanc[5] += 1;/*Sfo*/
						}
					}
					else {
						if(dna_matrix[inits1][j] != dna_matrix[inits2][j])
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
		for(j=s0;j<s1;j++) {
			while(j < s1) {
				if((c=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam, dna_matrix, posit)) > 0) break;
				else j++;
			}                    
			if(j<s1) {
				/*do for each population, first outgroup*/
				polc = 0;
				for(popo=1;popo<=(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][0];popo++) {
					nt0 = 0;
					for(h=0;h<(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][popo];h++) nt0 += (*inputp)->config[h];	
					for(h=nt0;h<nt0+(*inputp)->config[(*inputp)->ancestral_pol[(*inputp)->type_ancestral-1][popo]];h++) {
						if(dna_matrix[h][j] != dna_matrix[initso][j]) {
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
								if(dna_matrix[h][j] != dna_matrix[initsq1[x]][j]) {
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
										if(dna_matrix[h][j] != dna_matrix[initsq2[x]][j]) {
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
						if(dna_matrix[initsq2[x]][j] == dna_matrix[initso][j]) (ntpar)->Sanc[x*4+0] += 1;/*Sx1*/
						else (ntpar)->Sanc[x*4+2] += 1;/*Sx1f2*/
					}
					
					if(polqa[x]==1 && polqb[x]==1 && polc==0) {
						(ntpar)->Sanc[x*4+3] += 1;/*Ssh*/
					}
					if(polqa[x]==0 && polqb[x]==0 && polc==0) {
						if((dna_matrix[initsq2[x]][j] == dna_matrix[initso][j]) && (dna_matrix[initsq1[x]][j] != dna_matrix[initso][j])) {
							(ntpar)->Sanc[x*4+1] += 1;/*Sf1*/
						}
					}
				}
				/*for the outgroup is doesn't matter what is the group of populations(x), we chose 0*/
				if(polc == 1 && (polqa[0] != -1 && polqb[0] != -1)) {
					(ntpar)->Sanc[((*inputp)->type_ancestral-1)*4+0] += 1;/*Sxo*/
				}
				if(polqa[0]==0 && polqb[0]==0 && polc==0) {
					if((dna_matrix[initsq1[0]][j] == dna_matrix[initsq2[0]][j]) && (dna_matrix[initsq1[0]][j] != dna_matrix[initso][j])) {
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
