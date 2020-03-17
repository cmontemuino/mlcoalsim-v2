#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mlsp_sm.h"

void print_matrix(long int segsites,struct var2 **inputp,FILE *output,int x,long int count0, struct var_priors *priors,int mhits, char **dna_matrix, long int *posit)
{
    long j,k,Ss,i;
    int h,y;
    int ispolnomhit(long int,int,int,int);
    void function_atcg(int,long int,char **,const double *);

    int xx;
    long int kk,ll,pnsites;
    long int s0,s1,psegsites=0;
    char **list2;
    char ss[1];

    /*x=1 pr_matrix, x=2 print Hudson format (modified for mhits), x=3 print matrix excluding positions with mhits*/

    if(x==3) fputs("\nMatrix of dna excluding all positions with multiple hits\n",output);

    if((*inputp)->linked == 1) {
        s0=s1=0;
        psegsites = 0;
        pnsites = (*inputp)->nsites;
        for(kk=0,ll=(*inputp)->window;kk<pnsites;kk += (long int)(*inputp)->despl) {
            while(s0 < segsites && posit[s0]<kk) s0++;
            while(s1 < segsites && posit[s1]<ll) s1++;

            for(i=s0;i<s1;i++) {
                if(i==s0) psegsites += 1;
                else if(posit[i-1] != posit[i])
                    psegsites ++;
            }
            if(ll == pnsites) break;
            else ll += (*inputp)->despl;
            if(ll > pnsites) ll = pnsites;
            s0=s1;
        }
    }
    if((*inputp)->linked == 2) {
        psegsites = 0;
        pnsites = 0;
        s0 = s1 = 0;
        for(xx=0;xx<(*inputp)->linked;xx++) {
            kk=(*inputp)->loci_linked[xx][0];
            ll=(*inputp)->loci_linked[xx][1]+1;
            while(s0 < segsites && posit[s0]<kk) s0++;
            while(s1 < segsites && posit[s1]<ll) s1++;

            psegsites += 1;
            for(i=s0;i<s1;i++) {
                if(i==s0) psegsites += 1;
                else if(posit[i-1] != posit[i])
                    psegsites ++;
            }
            pnsites += ll - kk;

            s0=s1;
        }
    }

    if(x==2) {
        fprintf(output,"\n//");
        /*include here the prior values!*/
        for(y=0;y<(*inputp)->npriors;y++) {
            if((*inputp)->pointtoprior[y])
                fprintf(output,"\t%f",(float)priors[y].priordist[count0]);
        }
        if((*inputp)->linked == 0) {
            Ss=0;
            for(i=0;i<segsites;i++)
                if(posit[i-1] != posit[i]) Ss++;
            fprintf(output,"\nsegsites: %ld",Ss);
            if(Ss > 0) fprintf(output,"\npositions: ");
            else fprintf(output,"\n ");
            for(i=0;i<(int)segsites;i++)
                if(posit[i-1] != posit[i])
                    fprintf(output,"%g ",(double)posit[i]/(double)(*inputp)->nsites);
            fprintf(output,"\n");
            if(Ss > 0) {
                for(i=0;i<(*inputp)->nsam;i++) dna_matrix[i][segsites] = '\0';
                for(i=0;i<(*inputp)->nsam;i++) {
                    if(mhits){
                        for(j=0;j<(int)segsites;j++)
                            if(posit[j-1] != posit[j])
                                fprintf(output,"%c",dna_matrix[i][j]);
                    } else {
                        for(j=0;j<(int)segsites;j++) {
                            if(dna_matrix[i][j] != '0')
                                ss[0] = '1';
                            else
                                ss[0] = '0';
                            fprintf(output,"%c",ss[0]);
                        }
                    }
                    fprintf(output,"\n");
                }
            }
        } else {
            if((*inputp)->linked == 1) {
                fprintf(output,"\nsegsites: %ld",psegsites);
                if(psegsites > 0) {
                    fprintf(output,"\npositions: ");
                    s0=s1=0;
                    for(kk=0,ll=(*inputp)->window;kk<(*inputp)->nsites;kk += (long int)(*inputp)->despl) {
                        while(s0 < segsites && posit[s0]<kk) s0++;
                        while(s1 < segsites && posit[s1]<ll) s1++;

                        for(i=s0;i<(int)s1;i++)
                            if(i>0 && posit[i-1] != posit[i])
                                fprintf(output,"%g ",(double)posit[i]/(double)(*inputp)->nsites);
                            else /*we cannot avoid those positions because segsites is already defined and printed!!. RE-CODE*/
                                fprintf(output,"%g ",(double)posit[i]/(double)(*inputp)->nsites);/*error: segsites already printed!!*/
                        if(ll == (*inputp)->nsites) break;
                        else ll += (*inputp)->despl;
                        if(ll > (*inputp)->nsites) ll = (*inputp)->nsites;
                        s1=s0;
                    }
                    fprintf(output,"\n");
                    for(i=0;i<(*inputp)->nsam;i++) {
                        s0=s1=0;
                        for(kk=0,ll=(*inputp)->window;kk<(*inputp)->nsites;kk += (long int)(*inputp)->despl) {
                            while(s0 < segsites && posit[s0]<kk) s0++;
                            while(s1 < segsites && posit[s1]<ll) s1++;

                            for(j=s0;j<s1;j++) {
                                if(mhits) {
                                    if(posit[j-1] != posit[j])
                                        fprintf(output,"%c",dna_matrix[i][j]);
                                    else/*we cannot avoid those positions because segsites is already defined and printed!!. RE-CODE*/
                                        fprintf(output,"%c",dna_matrix[i][j]);/*error: segsites already printed!!*/
                                }
                                else {
                                    if(dna_matrix[i][j] != '0')
                                        ss[0] = '1';
                                    else
                                        ss[0] = '0';
                                    fprintf(output,"%c",ss[0]);
                                }
                            }

                            if(ll == (*inputp)->nsites) break;
                            else ll += (*inputp)->despl;
                            if(ll > (*inputp)->nsites) ll = (*inputp)->nsites;
                            s1=s0;
                        }
                        fprintf(output,"\n");
                    }
                }
                else fprintf(output,"\n ");
            }
            if((*inputp)->linked == 2) {
                fprintf(output,"\nsegsites: %ld",psegsites);
                if(psegsites > 0) {
                    fprintf(output,"\npositions: ");
                }
                else fprintf(output,"\n ");
                s0=s1=0;
                for(xx=0;xx<(*inputp)->linked;xx++) {
                    kk=(*inputp)->loci_linked[xx][0];
                    ll=(*inputp)->loci_linked[xx][1]+1;
                    while(s0 < segsites && posit[s0] < kk) s0++;
                    while(s1 < segsites && posit[s1] < ll) s1++;

                    for(i=s0;i<(int)s1;i++)
                        if(posit[i-1] != posit[i])
                            fprintf(output,"%g ",(double)posit[i]/(double)(*inputp)->nsites);
                        else /*we cannot avoid those positions because segsites is already defined and printed!!. RE-CODE*/
                            fprintf(output,"%g ",(double)posit[i]/(double)(*inputp)->nsites);/*error: segsites already printed!!*/
                    s0 = s1;
                }
                fprintf(output,"\n");
                for(i=0;i<(*inputp)->nsam;i++) {
                    s0=s1=0;
                    for(xx=0;xx<(*inputp)->linked;xx++) {
                        kk=(*inputp)->loci_linked[xx][0];
                        ll=(*inputp)->loci_linked[xx][1]+1;
                        while(s0 < segsites && posit[s0] < kk) s0++;
                        while(s1 < segsites && posit[s1] < ll) s1++;

                        for(j=s0;j<s1;j++) {
                            if(mhits) {
                                if(posit[j-1] != posit[j])
                                    fprintf(output,"%c",dna_matrix[i][j]);
                                else/*we cannot avoid those positions because segsites is already defined and printed!!. RE-CODE*/
                                    fprintf(output,"%c",dna_matrix[i][j]);/*error: segsites already printed!!*/
                            }
                            else {
                                if(dna_matrix[i][j] != '0') ss[0] = '1';
                                else ss[0] = '0';
                                fprintf(output,"%c",ss[0]);
                            }
                        }

                        s0 = s1;
                    }
                    fprintf(output,"\n");
                }
            }
        }
    }
    if(x==1) {
        if(!(list2 = (char **)malloc((unsigned)((*inputp)->nsam)*sizeof(char *)))) perror("malloc error in ms.3a");
        for(i=0;i<(*inputp)->nsam;i++) {
            if(!(list2[i] = (char *)malloc((long int)((*inputp)->nsites+1)*sizeof(char))))
                perror("malloc error in ms.3");
            memset(list2[i],'A',(*inputp)->nsites);
            list2[i][(*inputp)->nsites] = '\0';
        }
        for(i=0;i<(*inputp)->nsam;i++) {
            if(dna_matrix[i][0] == '1') list2[i][posit[0]] = 'G';
            if(dna_matrix[i][0] == '2') list2[i][posit[0]] = 'C';
            if(dna_matrix[i][0] == '3') list2[i][posit[0]] = 'T';
        }
        for(j=1;j<(int)segsites;j++) {
            if(posit[j-1] != posit[j]) {
                for(i=0;i<(*inputp)->nsam;i++) {
                    if(dna_matrix[i][j] == '1') list2[i][posit[j]] = 'G';
                    if(dna_matrix[i][j] == '2') list2[i][posit[j]] = 'C';
                    if(dna_matrix[i][j] == '3') list2[i][posit[j]] = 'T';
                }
            }
        }

        if((*inputp)->patcg[0] != -1.)
            function_atcg((*inputp)->nsam,(*inputp)->nsites,list2,(*inputp)->patcg);

        fprintf(output,"\nFASTA file: locus %d, nsam %d, total sites %ld, mutations %ld, iteration %ld\n",(*inputp)->nloci,(*inputp)->nsam,(*inputp)->nsites,segsites,count0);
        j=h=0;
        for(i=0;i<(*inputp)->nsam;i++) {
            while(i-h >= (*inputp)->config[j]) {
                h += (*inputp)->config[j];
                j++;
            }
            fprintf(output,">L%ldPop%ld\n%s\n",i,j,list2[i]);
        }
        fprintf(output,"\n");
        for(i=0;i<(*inputp)->nsam;i++) free(list2[i]);
        free(list2);
    }
    if(x==3) { /* exclude mhits from the matrix, but also the outgroup sequence in speciation with mhits */
        if(!(list2 = (char **)malloc((unsigned)((*inputp)->nsam)*sizeof(char *)))) perror("malloc error in ms.3a");
        for(i=0;i<(*inputp)->nsam;i++) {
            if(!(list2[i] = (char *)malloc((long int)((*inputp)->nsites+1)*sizeof(char))))
                perror("malloc error in ms.3");
            memset(list2[i],'A',(*inputp)->nsites);
            list2[i][(*inputp)->nsites] = '\0';
        }
        k = 0;
        if(segsites > 0) {
            if((h=ispolnomhit(0,0,(*inputp)->nsam,(*inputp)->nsam)) > 0) {
                for(i=0;i<(*inputp)->nsam;i++) {
                    if(dna_matrix[i][0] == '1') list2[i][posit[0]-k] = 'G';
                    if(dna_matrix[i][0] == '2') list2[i][posit[0]-k] = 'C';
                    if(dna_matrix[i][0] == '3') list2[i][posit[0]-k] = 'T';
                }
            }
            else if(h == -2) k++;/*k=0; check positions*/
            for(j=1;j<(int)segsites;j++) {
                if((h=ispolnomhit(j,0,(*inputp)->nsam,(*inputp)->nsam)) > 0) {
                    for(i=0;i<(*inputp)->nsam;i++) {
                        if(dna_matrix[i][j] == '1') list2[i][posit[j]-k] = 'G';
                        if(dna_matrix[i][j] == '2') list2[i][posit[j]-k] = 'C';
                        if(dna_matrix[i][j] == '3') list2[i][posit[j]-k] = 'T';
                    }
                }
                else if(h == -2) k++;/*k=0; check positions*/
            }
        }
        for(i=0;i<(*inputp)->nsam;i++) {
            list2[i][((*inputp)->nsites)-k] = '\0';
        }

        fprintf(output,"\nFASTA file: locus %d, nsam %d, nsites %ld, mutations %ld, iteration %ld\n",(*inputp)->nloci,(*inputp)->nsam,((*inputp)->nsites)-k,segsites,count0);
        h = j = 0;
        for(i=0;i<(*inputp)->nsam;i++) {
            while(i-h >= (*inputp)->config[j]) {
                h += (*inputp)->config[j];
                j++;
            }
            fprintf(output,">L%ldPop%ld\n%s\n",i,j,list2[i]);
        }
        fprintf(output,"\n");
        for(i=0;i<(*inputp)->nsam;i++) free(list2[i]);
        free(list2);
    }
}

void function_atcg(int nsam,long int nsites,char **list2,const double *patcg)
{
    double _random(void);

    int i,flagp,s;
    long int j;
    char p[4],n[4];
    double r;

    for(j=0;j<nsites;j++) {
        flagp = 0;
        p[0] = p[1] = p[2] = p[3] = (char)0;

        if(list2[0][j] > 64 && list2[0][j] < 85)
            p[0] = list2[0][j];
        else continue;
        /*defining positions, mono and polymorphic*/
        for(i=1;i<nsam;i++) {
            if(list2[i][j] > 64 && list2[i][j] < 85) {
                if(p[0] != list2[i][j]) {
                    if(p[1] == (char)0) p[1] = list2[i][j];
                    else {
                        if(p[1] != list2[i][j]) {
                            if(p[2] == (char)0) p[2] = list2[i][j];
                            else {
                                if(p[2] != list2[i][j]) {
                                    if(p[3] == (char)0) p[3] = list2[i][j];
                                    else {
                                        printf("\nError printing matrix.\n");
                                        exit(1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(p[1])
            flagp = 1;
        /*monomorphic*/
        if(flagp == 0) {
            r = _random();
            if(r<patcg[0]) n[0] = 'A';
            else if(r<patcg[0]+patcg[1]) n[0] = 'T';
            else if(r<patcg[0]+patcg[1]+patcg[2]) n[0] = 'C';
            else n[0] = 'G';
        }
        else {
            /*biallelic*/
            if(p[2] == 0) {
                /*transitions*/
                if((p[0] == 'A' && p[1] == 'G') ||
                   (p[0] == 'G' && p[1] == 'A') ||
                   (p[0] == 'C' && p[1] == 'T') ||
                   (p[0] == 'T' && p[1] == 'C')) {
                    r = _random();
                    if(r<patcg[0]) {n[0] = 'A'; n[1] = 'G';}
                    else if(r<patcg[0]+patcg[1]) {n[0] = 'T'; n[1] = 'C';}
                    else if(r<patcg[0]+patcg[1]+patcg[2]) {n[0] = 'C'; n[1] = 'T';}
                    else {n[0] = 'G'; n[1] = 'A';}
                }
                else {
                    /*transversions*/
                    r = _random();
                    if(r<patcg[0]) {
                        n[0] = 'A';
                        r = _random();
                        if(r < (patcg[1]/(patcg[1]+patcg[2]))) n[1] = 'T';
                        else n[1] = 'C';
                    }
                    else
                    if(r<patcg[0]+patcg[1]) {
                        n[0] = 'T';
                        r = _random();
                        if(r < (patcg[0]/(patcg[0]+patcg[3]))) n[1] = 'A';
                        else n[1] = 'G';
                    }
                    else
                    if(r<patcg[0]+patcg[1]+patcg[2]) {
                        n[0] = 'C';
                        r = _random();
                        if(r < (patcg[0]/(patcg[0]+patcg[3]))) n[1] = 'A';
                        else n[1] = 'G';
                    }
                    else {
                        n[0] = 'G';
                        r = (double)rand()/(double)RAND_MAX;
                        if(r < (patcg[1]/(patcg[1]+patcg[2]))) n[1] = 'T';
                        else n[1] = 'C';
                    }
                }
            }
            else {
                /*triallelic or tetra*/
                if((p[0] == 'A' && p[1] == 'G') ||
                   (p[0] == 'G' && p[1] == 'A') ||
                   (p[0] == 'C' && p[1] == 'T') ||
                   (p[0] == 'T' && p[1] == 'C')) s = 1;
                else if((p[0] == 'A' && p[2] == 'G') ||
                        (p[0] == 'G' && p[2] == 'A') ||
                        (p[0] == 'C' && p[2] == 'T') ||
                        (p[0] == 'T' && p[2] == 'C')) s = 2;
                else if((p[0] == 'A' && p[3] == 'G') ||
                        (p[0] == 'G' && p[3] == 'A') ||
                        (p[0] == 'C' && p[3] == 'T') ||
                        (p[0] == 'T' && p[3] == 'C')) s = 3;
                else s = 0;

                if(p[3] == 0) {/*tri*/
                    r = _random();
                    if(r<patcg[0]) {
                        n[0] = 'A';
                        if(s) n[s] = 'G';
                        r = _random();
                        if(s==1) {
                            if(r<0.5) {n[2] = 'T';}
                            else {n[2] = 'C';}
                        }
                        else
                        if(s==2){
                            if(r<0.5) {n[1] = 'T';}
                            else {n[1] = 'C';}
                        }
                        else {
                            if(r<0.5) {n[1] = 'T'; n[2] = 'C';}
                            else {n[1] = 'C'; n[2] = 'T';}
                        }
                    }
                    else
                    if(r<patcg[0]+patcg[1]) {
                        n[0] = 'T';
                        if(s) n[s] = 'C';
                        r = _random();
                        if(s==1) {
                            if(r<0.5) {n[2] = 'A';}
                            else {n[2] = 'G';}
                        }
                        else
                        if(s==2){
                            if(r<0.5) {n[1] = 'A';}
                            else {n[1] = 'G';}
                        }
                        else {
                            if(r<0.5) {n[1] = 'A'; n[2] = 'G';}
                            else {n[1] = 'G'; n[2] = 'A';}
                        }
                    }
                    else
                    if(r<patcg[0]+patcg[1]+patcg[2]) {
                        n[0] = 'C';
                        if(s) n[s] = 'T';
                        r = _random();
                        if(s==1) {
                            if(r<0.5) {n[2] = 'A';}
                            else {n[2] = 'G';}
                        }
                        else
                        if(s==2){
                            if(r<0.5) {n[1] = 'A';}
                            else {n[1] = 'G';}
                        }
                        else {
                            if(r<0.5) {n[1] = 'A'; n[2] = 'G';}
                            else {n[1] = 'G'; n[2] = 'a';}
                        }
                    }
                    else {
                        n[0] = 'G';
                        if(s) n[s] = 'A';
                        r = _random();
                        if(s==1) {
                            if(r<0.5) {n[2] = 'T';}
                            else {n[2] = 'C';}
                        }
                        else
                        if(s==2){
                            if(r<0.5) {n[1] = 'T';}
                            else {n[1] = 'C';}
                        }
                        else {
                            if(r<0.5) {n[1] = 'T'; n[2] = 'C';}
                            else {n[1] = 'C'; n[2] = 'T';}
                        }
                    }
                }
                else {/*tetra*/
                    r = _random();
                    if(r<patcg[0]) {
                        n[0] = 'A';
                        n[s] = 'G';
                        r = _random();
                        if(s==1) {
                            if(r<0.5) {n[2] = 'T'; n[3] = 'C';}
                            else {n[2] = 'C'; n[3] = 'T';}
                        }
                        else if(s==2){
                            if(r<0.5) {n[1] = 'T'; n[3] = 'C';}
                            else {n[1] = 'C'; n[3] = 'T';}
                        }
                        else {
                            if(r<0.5) {n[2] = 'T'; n[1] = 'C';}
                            else {n[2] = 'C'; n[1] = 'T';}
                        }
                    }
                    else
                    if(r<patcg[0]+patcg[1]) {
                        n[0] = 'T';
                        n[s] = 'C';
                        r = _random();
                        if(s==1) {
                            if(r<0.5) {n[2] = 'A'; n[3] = 'G';}
                            else {n[2] = 'G'; n[3] = 'A';}
                        }
                        else if(s==2){
                            if(r<0.5) {n[1] = 'A'; n[3] = 'G';}
                            else {n[1] = 'G'; n[3] = 'A';}
                        }
                        else {
                            if(r<0.5) {n[2] = 'A'; n[1] = 'G';}
                            else {n[2] = 'G'; n[1] = 'A';}
                        }
                    }
                    else
                    if(r<patcg[0]+patcg[1]+patcg[2]) {
                        n[0] = 'C';
                        n[s] = 'T';
                        r = _random();
                        if(s==1) {
                            if(r<0.5) {n[2] = 'A'; n[3] = 'G';}
                            else {n[2] = 'G'; n[3] = 'A';}
                        }
                        else if(s==2){
                            if(r<0.5) {n[1] = 'A'; n[3] = 'G';}
                            else {n[1] = 'G'; n[3] = 'A';}
                        }
                        else {
                            if(r<0.5) {n[2] = 'A'; n[1] = 'G';}
                            else {n[2] = 'G'; n[1] = 'A';}
                        }
                    }
                    else {
                        n[0] = 'G';
                        n[s] = 'A';
                        r = _random();
                        if(s==1) {
                            if(r<0.5) {n[2] = 'T'; n[3] = 'C';}
                            else {n[2] = 'C'; n[3] = 'T';}
                        }
                        else if(s==2){
                            if(r<0.5) {n[1] = 'T'; n[3] = 'C';}
                            else {n[1] = 'C'; n[3] = 'T';}
                        }
                        else {
                            if(r<0.5) {n[2] = 'T'; n[1] = 'C';}
                            else {n[2] = 'C'; n[1] = 'T';}
                        }
                    }
                }
            }
        }
        for(i=0;i<nsam;i++) {
            if(list2[i][j] == p[0]) list2[i][j] = n[0];
            else if(list2[i][j] == p[1]) list2[i][j] = n[1];
            else if(list2[i][j] == p[2]) list2[i][j] = n[2];
            else if(list2[i][j] == p[3]) list2[i][j] = n[3];
        }
    }
}

void show_progress(double counterp, double *restp, long int *p)
{
    if((double)*p + *restp >= counterp) {
        *restp += (double)*p - counterp;
        if(*restp/counterp > 1.0) {
            for(int i=0;i<(int)floor(*restp/counterp);i++) printf(".");
            *restp -= floor(*restp/counterp) * counterp;
        }
        else printf(".");
        fflush(stdout);
        *p = 1;
    } else {
        *p += 1;
    }
}

double _random(void)
{
  return rand()/(double)RAND_MAX;
}
