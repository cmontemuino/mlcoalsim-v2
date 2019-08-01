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

/*Hudson streec file modified*/

#include "mlsp_sm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SEGINC 1000
#define SELECTION_RECW_ALLOWED	1
#define SELECTION_CLASSIC 2 /*0:calculates the time directly using deterministic equation (NOT WORKING). 1: given dt and the deterministic equation, calculates events (sfreqinit option not included!) 2: given dt, DRIFT and selection, calculates events*/
#define FSTFROMNODIVISION 0
#define PRINT_TRACE_SELECTION 2
/*PRINT_TRACE_SELECTION=0 means standard for soft sweeps: many alleles are possible, like Petrov defines (i.e., the mutation may come from different lineages)*/
/*PRINT_TRACE_SELECTION=1: a single allele for soft sweep and print trajectory; 2. This is strict softsweep from a single lineag. Conditioning on fixing (losing in coalescence) the derived allele: change the effect of selection for neutral scenario: pop_sel = 2.0/(1.0-xdt); Solution from Zhao,Lascoux,Overall and Waxman, Genetics 2013, from paragraph after Equation 8. pop_sel is 4Ns and s is 1/(2Nx) following the conditional fixing, so pop_sel is 2/x */
/*PRINT_TRACE_SELECTION=2: a single allele for soft sweep but no trajectory printed. This is strict softsweep from a single lineag. Conditioning on fixing (losing in coalescence) the derived allele: change the effect of selection for neutral scenario: pop_sel = 2.0/(1.0-xdt); Solution from Zhao,Lascoux,Overall and Waxman, Genetics 2013, from paragraph after Equation 8. pop_sel is 4Ns and s is 1/(2Nx) following the conditional fixing, so pop_sel is 2/x */
#define PARAMETER_FREQINIT 2
/*PARAMETER_FREQINIT=2 means the parameter used for the frequency of current selected alleles at past is the SAMPLE     frequency*/
/*PARAMETER_FREQINIT=1 means the parameter used for the frequency of current selected alleles at past is the POPULATION frequency and not the sampled*/

long int nchrom;
long int begs;
long int nsegs;
long int nlinks;
double t, cleft,pc,lnpc;

long int total_nts,sel_nts_glob,new_chrom;	/*in case selection*/
int ifsel_glob;					/*in case selection*/

double nlinksr,total_ntsr;

static long int seglimit = SEGINC;
static long int maxchr;

/*Cada chrom té un nombre de segments, els quals es troben a una població. La direccio dels segments es indicada*/
/*Cada chrom té els seus propis segments. S'indica l'inici i final de cada segment, així com el nombre del segment amb el que continua aquest chrom*/
struct seg {
    long int beg;
    long int end;
    int desc;
};
struct chromo {
    long int nseg;
    long int pop;
    struct seg *pseg;
};

struct node *ptree1,*ptree2;

/* El resultat és la matriu de segments-length, que té associat tots els arbres per a cada segment. */
/* Les matrius chrom i pseg només són neccesàries per fer la coalescència, però després ja no són necessàries. */

/* Hudson routine*/
/**********  segtre_mig.c **********************************
*
*	This subroutine uses a Monte Carlo algorithm described in
*	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
*	a history of a random sample of gametes under a neutral
*	Wright-Fisher model with recombination and geographic structure. 
*	Input parameters
*	are the sample size (nsam), the number of sites between
*	which recombination can occur (nsites), and the recombination
*	rate between the ends of the gametes (r). The function returns
*	nsegs, the number of segments the gametes were broken into
*	in tracing back the history of the gametes.  The histories of
*	these segments are passed back to the calling function in the
*	array of structures seglst[]. An element of this array,  seglst[i],
* 	consists of three parts: (1) beg, the starting point of
*	of segment i, (2) ptree, which points to the first node of the
*	tree representing the history of the segment, (3) next, which
*	is the index number of the next segment.
*	     A tree is a contiguous set of 2*nsam nodes. The first nsam
*	nodes are the tips of the tree, the sampled gametes.  The other
*	nodes are the nodes ancestral to the sampled gametes. Each node
*	consists of an "abv" which is the number of the node ancestral to
*	that node, an "ndes", which is not used or assigned to in this routine,
*	and a "time", which is the time (in units of 4N generations) of the
*	node. For the tips, time equals zero.
*	Returns a pointer to an array of segments, seglst.

**************************************************************************/

struct segl *segtre_mig(long int npop,int nsam,int *inconfig,long int nsites,double r,double f,double track_len,
    double mig_rate,long int *pnsegs,long int iteration, double *factor,int ifselection, double pop_sel, double sinit,
    double pop_size,long int sel_nt,int *all_sel,int *selnsam, int *nintn,double **nrec,double **npast, double **tpast,
    int split_pop, double time_split, double time_scoal, double factor_anc, double *freq, double tlimit,int iflogistic,
    double *Tts,double factor_chrn,double *weightrec, double **migrate_matrix,int my_rank, int npop_events,
    struct events *pop_event,int event_forsexratio, double event_sexratio, double sex_ratio, int no_rec_males, 
    double sendt, double sfreqend,double sfreqinit)
{
	long int *nnodes = NULL;
	struct segl *seglst = NULL;
	struct chromo *chrom = NULL;

    long int j,dec,c1,c2,ind,rchrom;
    long int migrant,*config;
    long int pop,source_pop;/*modificat a long int*/
    
    double trm,ttemp,rft,clefta;
	double tcoal = 0.;
    double prec,cin,prect,mig,ran,coal_prob,rdum; 
    int re(int,double *,long int,double, long int *, struct segl *, struct chromo *);
	int ca(int,long int,long int,long int,double *,double,long int, long int *, struct segl *, struct chromo *);
    void pick2_chrom(long int,long int *,long int *,long int *, struct chromo *);
    int cinr(int,long int,double *,double,long int, long int *, struct segl *, struct chromo *);
	int cleftr(int,double *,long int,double, long int *, struct segl *, struct chromo *chrom); 

    double ran1(void);
    double binomialdist(double,int);
	double largebinomialdist(double,double);
    
    double xdt,eps,tf,ts,tot_rec,freqend,freqinit;	/*for selection*/
	#if SELECTION_CLASSIC == 0
	double coal_prob0,coal_prob1,tcoal0,tcoal1;
	#endif
	#if SELECTION_CLASSIC
	double dt;
	double rec_prob;
	double prob;
	#endif
	#if SELECTION_CLASSIC == 2
	double fitness;
	double avgftns;
	double coal_prob0;
	double coal_prob1;
    int flag_freq=0;
	#endif
	#if SELECTION_CLASSIC == 1 || SELECTION_CLASSIC == 0
    double eps2n,freqpa;
    long int timeg;
    #endif
	
    double poissondist(double xm); /*finding the time xdt reach to eps (before the deterministic selective event starts)*/
	
	double *nref;
	double *alphag;
	double *Ts;
	int *intn;
	double tpop;
	int coalpop,flagintn;
	
	double xacc,x1,x2;
	int zbracn(double(*)(double,double,double,double,double,double,double,double,double,double),double *,double *,double,double,double,double,double,double,double,double,double);
	double zriddrn(double(*)(double,double,double,double,double,double,double,double,double,double),double,double,double,double,double,double,double,double,double,double,double,double);
	double functiont_freqp_sel(double,double,double,double,double,double,double,double,double,double);
	double functiont_freqq_nsel(double,double,double,double,double,double,double,double,double,double);
	double functiont_logistic(double,double,double,double,double,double,double,double,double,double);
	
    double correction_rec(double,double,int);
    double correction_theta(double,double);
	double factor_chrnall;
	/*long int lirec;*/
	/*double *weightrec;*/
	double sumfreq,pip,spip;			/*refugia*/
    
    #if FSTFROMNODIVISION == 1
		int cf=0; /*per Fst*/
	#endif
	
	int a,b;
	double migran,migsum,migppop;
	
	double **migm;
	double *factpop;
	double time_first_event;
	int nevent;
	double sex_ratio_original;
	
	int po1,po2;
	double pof1,pof2;
	
	long int jj;/*for infinite recombination*/
	long int indtrees = 1;
	long int tr = 0;
    
    #if PRINT_TRACE_SELECTION == 1
    FILE *filexdt;/*debugging t and xdt*/
    #endif
    
	#if PARAMETER_FREQINIT == 2
	double trapzd_inv(double (*distp_binom)(int,int,double),double,double,int,int,int,double *,double *);
	double distp_binom(int n, int k, double p);
	double distp_binom(int,int,double);
	double cummtot;
	static double *xm = 0;
	static double *psum = 0;
	static int itm;
	int ni;
	#endif
	
	sex_ratio_original = sex_ratio;
	if(iteration) j = 0;
	else j = 0;
	
	/*correcting the tree depending on the sexratio and the chromosome*/
	factor_chrnall = 1.0/(double)correction_theta(factor_chrn,sex_ratio);
	/*correcting recombination according the sexratio and the chromosome*/
	r *= correction_rec(factor_chrn,sex_ratio,no_rec_males);
    
	sel_nts_glob = sel_nt;			/*for selection*/
    ifsel_glob = ifselection;			/*for selection*/

    
	maxchr = nsam + 50;	/* + 50 ... */
	if((chrom = (struct chromo *)malloc((unsigned)(maxchr*sizeof(struct chromo)))) == NULL)
		perror("malloc error. segtre_mig.1");

	if((nnodes = (long int *)malloc((unsigned)(seglimit*sizeof(long int)))) == NULL)
		perror("malloc error. segtre_mig.2");

	if((seglst = (struct segl *)malloc((unsigned)(seglimit*sizeof(struct segl)))) == NULL)
		perror("malloc error. segtre_mig.3");

	if((intn = (int *)malloc((unsigned)(npop*sizeof(int)))) == NULL)
		perror("malloc error. segtre_mig.32");
	if((nref = (double *)malloc((unsigned)(npop*sizeof(double)))) == NULL)
		perror("malloc error. segtre_mig.34");
	if((alphag = (double *)malloc((unsigned)(npop*sizeof(double)))) == NULL)
		perror("malloc error. segtre_mig.35");
	if((Ts = (double *)malloc((unsigned)(npop*sizeof(double)))) == NULL)
		perror("malloc error. segtre_mig.35");
		
	config = (long int *)malloc((unsigned)((npop+1)*sizeof(long int)));	/* nou vector de mida de poblacions */
	if(config==NULL) perror("malloc error segtre_mig.4");

	if(r/nsites < 1e6) {
		indtrees = 1; 

		seglst[0].beg = 0;
		if(!(seglst[0].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node))))
			perror("calloc error segtre_mig.6");
		nnodes[0] = nsam -1;	/* nombre de fragments per a cada segment .. */
		
		/*new values to take into account heterogeneity in recombination*/
		nlinksr = (double)nsam * weightrec[nsites-1] * r;
		total_ntsr = nlinksr;
		
		nsegs = 1;
	}
	else {
		/* for infinite recombination.*/
		indtrees = nsites; 

		if((seglst = (struct segl *)realloc(seglst,(indtrees*sizeof(struct segl)))) == NULL)
			perror("malloc error. segtre_mig.1nsites");
		if((nnodes = (long int *)realloc(nnodes,(indtrees*sizeof(long int)))) == NULL)
			perror("malloc error. segtre_mig.2");
		for(jj=0;jj<indtrees;jj++) {
			seglst[jj].beg = jj;
			if(!(seglst[jj].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node))))
				perror("calloc error segtre_mig.6nsites");
			nnodes[jj] = nsam - 1;	/* nombre de fragments per a cada segment .. */
			if(jj<indtrees-1) seglst[jj].next = jj+1;
			else seglst[jj].next = jj;
		}
		
		/*new values to take into account heterogeneity in recombination*/
		r = 0.;
		nlinksr = (double)nsam * weightrec[nsites-1] * r; /*<- recombination now eliminated. Every position is already separated. Only coalescence*/
		total_ntsr = nlinksr;
	}

	tr = 0;
	while(tr < indtrees) {
		nlinks = ((long int)(nsam))*(nsites-1);
		total_nts = nlinks; 		/*included for selection*/
		if(ifselection == 0) sel_nt = 0;	/*included for selection*/
		

		if((int)maxchr < nsam) {
			maxchr = (long int)nsam + 50;
			if((chrom = (struct chromo *)realloc(chrom,(maxchr*sizeof(struct chromo)))) == NULL)
				perror("malloc error. segtre_mig.1");
		} else {
			if((chrom = (struct chromo *)realloc(chrom,(maxchr*sizeof(struct chromo)))) == NULL)
				perror("malloc error. segtre_mig.1");
		}
		
		for(pop=0;pop<npop;pop++) config[pop] = (long int)inconfig[pop];	/* inicialitzar amb valors del input */
		for(pop=ind=0;pop<npop;pop++)	/* matriu chrom, indica els individus de cada població i els segments/ind */
			for(j=0;j<inconfig[pop];j++,ind++) {
				chrom[ind].nseg = 1;
				if(!(chrom[ind].pseg = (struct seg *)malloc((unsigned)sizeof(struct seg))))
					perror("malloc error segtre_mig.5");
				(chrom[ind].pseg)->beg = 0;
				(chrom[ind].pseg)->end = nsites - 1;
				(chrom[ind].pseg)->desc = (int)ind;
				chrom[ind].pop = pop;
			}
		
		nchrom = nsam; 		/* nombre de fragments que hi han a l'arbre .. */
		nsegs = tr+1; /*in case independent positions, 1 segment each time starting from indtrees, otherwise 1*/

		t = (double)0;
		trm = 0.;

		if(ifselection == 0 || (ifselection==1 && pop_sel <= (double)-10000)) {
			
			if(f > 0.0) pc = (track_len - (double)1)/track_len;	/* tot això conversió gènica */
			else pc = (double)1;
			lnpc = (double)log((double)pc);
			cleft = nsam * ((double)1 - pow(pc,(double)(weightrec[nsites-2]*r/*nsites-1*/)));
			rft = /*r**/f*track_len;
			
			if(split_pop == 0) {
				/********************* ROUTINE modified FROM HUDSON *************************/
				for(a=0;a<npop;a++) {
					intn[a] = 1;
					alphag[a] = (double)2*(double)10/(tpast[a][intn[a]]-tpast[a][intn[a]-1]);/*logistic*/
					if(iflogistic) {
						Ts[a] = (double)Tts[a];
						nref[a] = ((double)nrec[a][1] + ((double)npast[a][1]-(double)nrec[a][1])/((double)1 + (double)exp((double)(-(double)alphag[a]*(Ts[a]-((double)tpast[a][0]+(double)tpast[a][1])/(double)2)))));
					}
				}
				/*new*/
				migm = (double **)malloc((unsigned)((npop)*sizeof(double *)));	/* new migration matrix */
				if(mig_rate > 0.) {
					for(a=0;a<npop;a++) {
						migm[a] = (double *)calloc((unsigned)(npop),sizeof(double));
						for(b=0;b<npop;b++) 
							if(a != b) 
								migm[a][b] = mig_rate/(((double)npop)*((double)npop-1.));
					}
				}
				else {
					for(a=0;a<npop;a++) {
						migm[a] = (double *)malloc((unsigned)((npop)*sizeof(double)));
						for(b=0;b<npop;b++) migm[a][b] = migrate_matrix[a][b];
					}
				}
				factpop = (double *)malloc((unsigned)((npop)*sizeof(double)));	/* new factor_pop matrix */
				for(a=0;a<npop;a++) {
					factpop[a] = factor[a];
				}
				if(npop_events) time_first_event = pop_event[0].time_event;
				else time_first_event = 1E07;
				nevent = 0;
				/******************************** Main loop *********************************/
				while(nchrom > 1) {
									 
					/*modification to increase the speed with high Recombination ..., if t > tlimit r = 0.;*/
					if(t > tlimit) {
						r = (double)0; 
						nlinksr = (double)0;
					}
					
					prec = (double) nlinksr;	/* prob. recombinació. Cada posicio entre 2nt te una prob. de recombinar */
					cin  = (double) nlinksr*f;  /* prob. cin event*/
					clefta = cleft * rft;       /* prob. cleft event*/
					prect = prec + cin + clefta;/* prob. recombinació + cin + cleft */
					
					/*mig = nchrom * mig_rate;*/	/* prob. migració. Cada ind. té una probabilitat de migrar */
					mig = 0.;
					/*if(t < time_scoal) {*/
						for(a=0;a<npop;a++) {
							for(b=0;b<npop;b++) {
								mig += config[a] * migm[a][b] /* * factpop[b] weighting for the popsize of the source?NO*/;
							}
						}
					/*}*/
					
					if(prect+mig > (double)0) {
						while((rdum = (double)ran1()) == 1.0);/* de fet, mai pot ser 1, es [0.1) en ran1 ... */
						trm = -(double)log((double)1-(double)rdum)/(prect+mig);/* temps a la recombinació, (conversió), migració o extincio*/
					}/* no es troben relacionats amb No, per tant no afecta canvi demogràfic */
					else trm = (double)0;
					
					coalpop = 0;
					tcoal = tpop = trm + 999999.;
					for(pop=0;pop<npop;pop++) {
						tpop = 999999.;
						coal_prob = ((double)config[pop])*(config[pop]-(double)1)*factor_chrnall/factpop[pop]; /* arbres en funció de 4No in a diploid sp*/
						if(coal_prob > (double)0) {
							while((rdum = (double)ran1()) == 1.0);
							/* ponderar canvi demogràfic: ACTIVAT. totes les poblacions canvien a la vegada*//*5.5.03*/
							if(t >= time_first_event || iflogistic==0) {
								/*instantaneous*/
								tpop = -(double)log((double)1.0 - (double)rdum)/coal_prob;
								tpop *= npast[pop][intn[pop]];
							}
							else {
								/*logistic*/
								tpop = -(double)log((double)1.0-(double)rdum)/coal_prob;
								/*calculate the range*/
								if(intn[pop] > 1) Ts[pop] = (double)0;
								x1 = (double)0.;
								x2 = (double)0.5;
								if(zbracn(functiont_logistic,&x1,&x2,(double)tpop,(double)t,(double)tpast[pop][intn[pop]-1],(double)tpast[pop][intn[pop]],(double)alphag[pop],(double)nrec[pop][intn[pop]],(double)npast[pop][intn[pop]],nref[pop],Ts[pop]) == 1) {
									/*estimate the value of T*/
									xacc = (double)1e-6*x2; /* accuracy*/
									tpop = (double)zriddrn(functiont_logistic,x1,x2,xacc,(double)tpop,(double)t,(double)tpast[pop][intn[pop]-1],(double)tpast[pop][intn[pop]],(double)alphag[pop],(double)nrec[pop][intn[pop]],(double)npast[pop][intn[pop]],nref[pop],Ts[pop]);
								}
							}
						}
						if(tpop < tcoal) {
							tcoal = tpop;
							coalpop = (int)pop;
						}
					}
					
					if((prect + mig) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
					else ttemp = tcoal;                
					
					flagintn = 0;
					if(nevent < npop_events) {
						if((t <= pop_event[nevent].time_event) && (t+ttemp > pop_event[nevent].time_event)) { /*CHANGES IN COALPOP,MIG,FACTOR_POP!!*/
							t = pop_event[nevent].time_event /*+1E-17*//*precission error*/;/*modify time*/
							/*modify factor_pop,mig matrix*/
							for(pop=0;pop<npop;pop++) {
								migm[pop_event[nevent].npop1][pop] = pop_event[nevent].mig_fw[pop];
								migm[pop][pop_event[nevent].npop1] = pop_event[nevent].mig_rv[pop];
							}
							factpop[pop_event[nevent].npop1] = pop_event[nevent].factor_pop;
							if(pop_event[nevent].npop1 != pop_event[nevent].npop2) {
								if(pop_event[nevent].npop2 >= 0) { /*in case npop2 is a new pop, pop is negative: split option*/
									/* two population coalesce */
									for(j=0;j<nchrom;j++) 
										if(chrom[j].pop == pop_event[nevent].npop2) 
											chrom[j].pop = pop_event[nevent].npop1;
									config[pop_event[nevent].npop1] += config[pop_event[nevent].npop2];
									config[pop_event[nevent].npop2] = 0;
									for(pop=0;pop<npop;pop++) {
										migm[pop_event[nevent].npop2][pop] = 0.;
										migm[pop][pop_event[nevent].npop2] = 0.;
									}
								}
								else {
									/* one population split in two: SIMPLIFIED: the split asigns the proportion factor_pop to each population */
									pof1 = pop_event[nevent].factor_pop;
									po1  = pop_event[nevent].npop1;
									nevent += 1;
									pof2 = pop_event[nevent].factor_pop;
									po2  = pop_event[nevent].npop1;
									for(pop=0;pop<npop;pop++) {
										migm[pop_event[nevent].npop1][pop] = pop_event[nevent].mig_fw[pop];
										migm[pop][pop_event[nevent].npop1] = pop_event[nevent].mig_rv[pop];
									}
									factpop[po2] = pof2;
									/* From one to 2 pops */
									sumfreq = pof1 + pof2;								
									for(j=0;j<nchrom;j++) {
										if((chrom[j].pop) == po1) {
											if((double)ran1() <= ((double)1.0 - pof1/sumfreq)) {
												chrom[j].pop = po2;
												config[po2] += 1;
												config[po1] -= 1;
											}
										}
									}
								}
							}
							/*CHANGES IN SEX_RATIO*/
							if(pop_event[nevent].sex_ratio > 0) {
								/*calculate all corrections, factor_chrnall,matrix theta and rec */
								factor_chrnall = 1.0/correction_theta(factor_chrn,pop_event[nevent].sex_ratio);
								r *= correction_rec(factor_chrn,pop_event[nevent].sex_ratio,no_rec_males)/correction_rec(factor_chrn,sex_ratio,no_rec_males);
								nlinksr *= correction_rec(factor_chrn,pop_event[nevent].sex_ratio,no_rec_males)/correction_rec(factor_chrn,sex_ratio,no_rec_males);
								sex_ratio = pop_event[nevent].sex_ratio;
							}
							else {
								if((nevent > 0) && (pop_event[nevent].sex_ratio !=  pop_event[nevent-1].sex_ratio)) { /*the sex ratio becomes again the original. SR changes only in one of the steps*/
									/*calculate all corrections, factor_chrnall,matrix theta and rec again*/
									factor_chrnall = 1.0/correction_theta(factor_chrn,sex_ratio_original);
									r *= correction_rec(factor_chrn,sex_ratio_original,no_rec_males)/correction_rec(factor_chrn,pop_event[nevent].sex_ratio,no_rec_males);
									nlinksr *= correction_rec(factor_chrn,sex_ratio_original,no_rec_males)/correction_rec(factor_chrn,pop_event[nevent].sex_ratio,no_rec_males);
									sex_ratio = sex_ratio_original;
								}
							}
							flagintn = 1;
							nevent += 1;
						}
					}
					if(nevent == 0) {				
						/*en cas un canvi No. Si ttemp > tpast aleshores tornem a calcular probs. i t=ttemp : Ara ACTIVAT*//*5.5.03*/
						for(pop=0;pop<npop;pop++) {
							if((intn[pop] < nintn[pop]) && (t + ttemp > tpast[pop][intn[pop]])) {
								t = tpast[pop][intn[pop]++];
								alphag[pop] = (double)2*(double)10/(tpast[pop][intn[pop]] - tpast[pop][intn[pop]-1]);/*logistic*/
								flagintn = 1;
							}
						}
					}
					if(flagintn == 0) {/*same process of intn*/
						t += ttemp;
						if(((prect+mig ) > 0.0) && (trm < tcoal)) {
							if((ran = ran1()) < (prec/(prect + mig))) {
								/* recombination */
								rchrom = re(nsam,weightrec,nsites,r, nnodes, seglst, chrom);	/* rchrom és l'individu a on té lloc la recombinació */		
								config[chrom[rchrom].pop] += 1;  /* aumenta 1 la població on pertany l'individu */
							}
							else {
								/**/
								if(ran  < ((prec + clefta)/(prect + mig))) {
									/* cleft event*/
									rchrom = cleftr(nsam,weightrec,nsites,r, nnodes, seglst, chrom);
									config[chrom[rchrom].pop] += 1;
								}
								else {
									if(ran < (prect/(prect + mig))) {
										/* cin event*/
										rchrom = cinr(nsam,nsites,weightrec,r,tr, nnodes, seglst, chrom);
										if(rchrom >= 0) config[chrom[rchrom].pop] += 1;
									}
									else {
										migran = mig*ran1();
										migsum = 0.;
										a = b = 0;
										migppop = config[a] * migm[a][b] /* * factpop[b] weighting for the popsize of the source?NO*/;
										while(migsum + migppop < migran) {
											migsum += migppop;
											b++;
											if(b==npop) {
												b=0;
												a++;
											}
											migppop = config[a] * migm[a][b] /* * factpop[b] weighting for the popsize of the source?NO*/;
										}
										migppop = migm[a][b];
										for(j=0;j<nchrom;j++) {
											if(chrom[j].pop == a) {
												if(migsum + migppop < migran) migsum += migppop;
												else break;
											}
										}
										migrant = j;
										source_pop = b;
										/*end new code*/
										config[chrom[migrant].pop] -= 1;
										config[source_pop] += 1;
										chrom[migrant].pop = source_pop;
									}
								}
							}
						}
						else {
							/* coalescent event */
							/* pick the two, c1, c2 */						
							pick2_chrom(coalpop,config,&c1,&c2, chrom);	/* escull c1 i c2 */
							dec = ca(nsam,nsites,c1,c2,weightrec,r,tr, nnodes, seglst, chrom);	/* dec és el nombre de fragments a restar */
							config[coalpop] = config[coalpop] - dec;			/* si hi ha MRCA aleshores és més d'un */
						}/*coal event*/
					}/*event*/
				}/*en cas chrom > 1*/
				for(a=0;a<npop;a++) free(migm[a]);
				free(migm);
				free(factpop);
			} else {
				if(split_pop == 1) {
					/*********************** REFUGIA **************************/
					sumfreq = (double)0;
					for(pop=0;pop<npop;pop++) sumfreq += freq[pop];
					for(pop=0;pop<npop;pop++) freq[pop] /= sumfreq;
					/* Main loop */                    
					while(nchrom > 1) {

						/*modification to increase the speed with high R, if t > tlimit r = 0.;*/
						if(t > tlimit) {
							r = (double)0; 
							nlinksr = (double)0;
						}

						prec = (double) nlinksr;	/* prob. recombinació. Cada posicio entre 2nt te una prob. de recombinar */
						cin  = (double) nlinksr*f;  /* prob. cin event*/
						clefta = cleft * rft;       /* prob. cleft event*/
						prect = prec + cin + clefta;/* prob. recombinació + cin + cleft */
						
						coalpop = 0;
						if((t >= time_split) && (t < time_scoal)) {
							coal_prob = (double)0;
							tcoal = 999999.;
							for(pop=0;pop<npop;pop++) {
								tpop = 999999.;
								coal_prob = ((double)config[pop])*(config[pop]-(double)1.0)*factor_chrnall/factor[pop];
								if(coal_prob > (double)0) {
									while((rdum = (double)ran1()) == 1.0);
									tpop = -(double)log((double)1.0 - (double)rdum)/coal_prob;								
								}
								if(tpop < tcoal) {
									tcoal = tpop;
									coalpop = (int)pop;
								}
							}
							coal_prob = coalpop;
						}
						else {
							if(t < time_split)
								coal_prob = ((double)config[0])*((double)config[0]-(double)1.0)*factor_chrnall;
							else /*t >= time_scoal*/
								coal_prob = ((double)config[0])*((double)config[0]-(double)1.0)*factor_chrnall/factor_anc;

							if(coal_prob > 0.0) {
								while((rdum = (double)ran1()) == (double)1.0);
								tcoal = -(double)log((double)1.0 - (double)rdum)/coal_prob;
							}
						}
						if(coal_prob == 0.0)
							tcoal = 999999.;
					   
						if(t<time_split || t>=time_scoal) mig = 0.;
						else mig = nchrom * mig_rate;	/* prob. migració. Cada ind. té una probabilitat de migrar */
						
						if(prect+mig > (double)0) {
							while((rdum = (double)ran1()) == (double)1.0);	
							trm = -(double)log((double)1.0-(double)rdum)/(prect+mig);	
						}
						
						if((prect) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
						else ttemp = tcoal;
						if((t < time_split) && (t+ttemp >= time_split)) {
							t = time_split /*+1E-17*//*precission error*/;
							/* From one to npop pops */
							for(j=0;j<nchrom;j++) {
								if((double)ran1() <= ((double)1.0/*sumfreq*/ - freq[0])) {
									pip = (double)ran1()*((double)1.0/*sumfreq*/ - freq[0]);
									spip = (double)0;
									for(pop=1;pop<npop;pop++)
										if((spip += freq[pop]) > pip) break;
									chrom[j].pop = pop;
									config[pop] += 1;
									config[0] -= 1;
								}
							}
						}
						else {
							if((t >= time_split) && (t < time_scoal) && (t+ttemp >= time_scoal)) {
								t = time_scoal /*+1E-17*//*precission error*/;
								/* From npop to one pop */
								for(j=0;j<nchrom;j++) chrom[j].pop = 0;
								config[0] = nchrom;
								for(pop=1;pop<npop;pop++) config[pop] = 0;
								/*printf("config[0]=%d\tconfig[1]=%d\n",config[0],config[1]);*/
							}
							else {
								t += ttemp;
															
								if(((prect+mig) > 0.0) && (trm < tcoal)) {
									if((ran = ran1()) < (prect/(prect + mig))) {
										/* recombination */
										rchrom = re(nsam,weightrec,nsites,r, nnodes, seglst, chrom);		
										config[chrom[rchrom].pop] += 1;
									}
									else {
										if(ran  < ((prec + clefta)/(prect + mig))) {
											/* cleft event*/
											rchrom = cleftr(nsam,weightrec,nsites,r, nnodes, seglst, chrom);
											config[chrom[rchrom].pop] += 1;
										}
										else {
											if(ran < (prect/(prect + mig))) {
												/* cin event*/
												rchrom = cinr(nsam,nsites,weightrec,r,tr, nnodes, seglst, chrom);
												if(rchrom >= 0) config[chrom[rchrom].pop] += 1;
											}
											else {
												/* migration event */
												migrant = (int)(nchrom*ran1());
												while((source_pop = (int)(npop*ran1())) == chrom[migrant].pop);
												config[chrom[migrant].pop] -= 1;
												config[source_pop] += 1;
												chrom[migrant].pop = source_pop;
											}
										}
									}
								}
								else {
									/* coalescent event */
									/* pick the two, c1, c2 */
									pick2_chrom(coalpop,config,&c1,&c2, chrom);
									dec = ca(nsam,nsites,c1,c2,weightrec,r,tr,nnodes, seglst, chrom);
									config[coalpop] -= dec;
								}
							}
						}
					}
				}
			}
		} else {/******************************SELECTION MODULE*********************************************/
			/*TRY TO INTRODUCE INSTANTANEOUS CHANGES IN Ne: IN CONSTRUCTION...*/
			factpop = (double *)malloc((unsigned)((1)*sizeof(double)));	/* new factor_pop matrix */
			for(a=0;a<1;a++) {
				factpop[a] = factor[a];
			}
			if(npop_events) time_first_event = pop_event[0].time_event;
			else time_first_event = 1E07;
			nevent = 0;
			flagintn = 0;

			/*EQUATION 3A and 3B IN STEPHAN 1992.*/
			/*selected pop is number 0. Non-selected pop is number 1*/
			eps = 1./(pop_sel/4.);/*100./(2.*(double)pop_size);*//*1./(2.*(double)pop_size);*//*5./(pop_sel/2.);*/
            /*frequency of selective allele when selective phase finish (1-eps when starts) */
			if(eps > 0.1) eps = 0.1;
			
            #if SELECTION_CLASSIC == 1 || SELECTION_CLASSIC == 0
            /*look for the generation when the deterministic event starts (ts). */
            /*Binomial, perhaps using a poisson given p << 0.1 and n >> 30*/
            eps2n = eps * (double)2* (double)pop_size;
			do{
				timeg = 0;
				freqpa = (double)1;/*one individual has the mutation*/
				while(freqpa > (double)0 && freqpa < eps2n) {
                    freqpa = binomialdist(freqpa/pop_size,pop_size); /*poissondist((double)freqpa);*/
					timeg++;
				}
			}while(freqpa < eps2n);
			ts = (double)timeg/((double)4*(double)pop_size);
            /*that is, we consider selection starts when ts starts*/

            if((sfreqend > 0.0) && (sendt != 1E06)) {
                /*incompatible parameters for calssic tye 0 and 1! we should calculate what is first!*/
                sfreqend  = 0.0;
            }
            #else
			ts = 0.0;
            /*in case selection type 2 it is considered always under selection but drift is more powerful*/
            #endif
            
            if(sfreqend > /*0.*/eps) freqend = sfreqend;/*modification for stopping the selective process at freqend given!!!!*/
            else freqend = eps;/*in case selective mutation is a new mutation*/
            if(sfreqinit < 1.0) freqinit = sfreqinit;/*modification for starting the selective process at freqinit given!!!!*/
            else freqinit = 1.0;/*in case selective mutation finishing normally*/
            if(sfreqinit<sfreqend) {
                printf("\nError: sfreqinit (%f) must not be smaller than sfreqend (%f). Only positive (or zero) selection is allowed.\n",sfreqinit,sfreqend);
                exit(1);
            }
            /*for selection type 0 or 1:*/
            if(sendt != 1E06) tf = sendt-sinit;/*modification for stopping the selective process at time sendt given!!!!*/
            else tf = (pop_sel*(ts)+(double)log((double)((1.-freqend-eps+eps*freqend)/(eps*freqend))))/(pop_sel);/*??*/
            /*time starting/ending selective phase*//*the strength of selection is 2Nes because affects the individual*/
            
			*all_sel = (int)config[0];
			if(sinit < 0.0) {/*selection has not fixed yet the selective allele */
                xdt = 1.0 - eps/(eps+(1-eps)*exp((-pop_sel * (0.0 -sinit-ts))));
                /*freq sel. allele in N individuals*/
				if(-sinit < ts) {
					*all_sel = (int)config[0];/*all alleles have the favorable mutation at time 0*/
					for(ind=0;ind<*all_sel;ind++) selnsam[ind] = (int)ind;/*fav lines in a vector (will be useful in mig+sel)*/
					for(ind=*all_sel;ind<(int)config[0];ind++) chrom[ind].pop = 1; /*chrom with unfav allele*/
					config[0] = nchrom;
					config[1] = 0;
				} 
				else {
					if(-sinit >= tf) {
						*all_sel = 0;/*no allele has the favorable mutation at time 0*/
						for(ind=0;ind<*all_sel;ind++) selnsam[ind] = (int)ind;/*fav lines in a vector (will be useful in mig+sel)*/
						for(ind=*all_sel;ind<(int)config[0];ind++) chrom[ind].pop = 1; /*chrom with unfav allele*/
						config[0] = 0;
						config[1] = nchrom;
					}
					else {
						*all_sel = (int)binomialdist((float)xdt,(int)config[0]);
                        /*how many alleles have the favorable mutation at time 0*/
						for(ind=0;ind<*all_sel;ind++) selnsam[ind] = (int)ind;/*fav lines in a vector (will be useful in mig+sel)*/
						for(ind=*all_sel;ind<(int)config[0];ind++) chrom[ind].pop = 1; /*chrom with unfav allele*/
						config[1] = config[0] - *all_sel;
						config[0] = *all_sel;
					}
				}
			}
			else xdt = 1.;
            /*sinit plus sfreqinit give the initial time where selection starts with a sfreqinit frequency!*/
			/* selective position is at 'sel_nt' starting from the left of the studied sequence (it can be neg)*/
			/* recombination includes all the region until the selected position */
			/* recombination is only performed when nt from the studied region are involved. */
			/* For the rest we only assign the group given xdt and r=r/(nsites-1). sites are the sites in the region.*/
			total_nts = ((((long int)nsites-1 > sel_nt-1) ? (long int)nsites-1 : sel_nt-1) - ((sel_nt < 0) ? sel_nt : 0));
			/*total_ntsr*/
			if((long int)nsites-1 > sel_nt-1) total_ntsr = weightrec[nsites-1]*r;
			else total_ntsr = ((double)(sel_nt-1)-(nsites-1))*r/(double)(nsites) + weightrec[nsites-1]*r;
			if(sel_nt < 0) total_ntsr -= sel_nt*r/(double)(nsites);
			/*it is necessary to know the total positions that are able to recombine (until the selected pos.)*/
			/*r *= (double)(nsites-1)/(double)total_nts;*/
			/*in selection module, R=4Nr for the fragment studied + all until the selective position*/
			total_nts *= nchrom; /*like nlinks but for all the region*/
			total_ntsr *= (double)nchrom; /*like nlinksr but for all the region*/
			/*************************************** Main loop for selection *****************************************/
			#if SELECTION_CLASSIC == 1
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            /*NOT COMPLETELY FUNCTIONAL: NOT INCLUDED THE OPTIONS WITH sfreqinit!!!*/
            #if PRINT_TRACE_SELECTION == 1
            if(iteration == 1){
                filexdt = fopen("./mlcoalsim_output_sel2_xdt_all.txt","w");
                fprintf(filexdt,"\nTime\txdt\tnsam[0]\tnsam[1]");
            }
            else filexdt = fopen("./mlcoalsim_output_sel2_xdt_all.txt","a");
                /**/
            #endif
			while(nchrom > 1) {
				/*modification to increase the speed with high R, if t > tlimit r = 0.;*/
				if(t > tlimit) {
					r = (double)0; 
					nlinksr = (double)0;
				}

				prect = (double)/*nlinks*/nlinksr /* * r*/;        /* 'tot_rec' is r in all region.*/
				tot_rec = (double)/*total_nts*/total_ntsr /* * r*/;   /*like prect but for all the region*/      

				/************ non-selective phase ****************************/
				if((t < sinit + ts) || (t >= sinit + tf)) { 
					coal_prob = ((double)nchrom)*(nchrom-(double)1)*factor_chrnall/factpop[0];
					if(prect > (double)0) {
						while((rdum = (double)ran1()) == (double)1);	
						trm = -(double)log((double)1.0-(double)rdum)/(prect);	
					}
					while((rdum = (double)ran1()) == 1.0);
					tcoal = -(double)log((double)1.0 - (double)rdum)/coal_prob;
					if((prect) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
					else ttemp = tcoal;

					if(nevent < npop_events) {
						if((t <= pop_event[nevent].time_event) && (t+ttemp > pop_event[nevent].time_event)) { /*CHANGES IN POPULATION SIZE: TO CHECK!!*/
							t = pop_event[nevent].time_event /*+1E-17*//*precission error*/;/*modify time*/
							/*modify factor_pop*/
							factpop[pop_event[nevent].npop1/*0*/] = pop_event[nevent].factor_pop;
							/*modification to consider the selective event: under construction: ONLY FOR pop_sel*factor[0] > 10 ALWAYS*/
							if(t < sinit + ts) { 
								eps = 1./(pop_sel*factpop[0]);
								if(eps > 0.1) eps = 0.1;
								freqend = eps; ts = t;
								tf = (pop_sel*factpop[0]*(ts)+(double)log((double)((1.-freqend-eps+eps*freqend)/(eps*freqend))))/(pop_sel*factpop[0]);
							}
							if((t <= sinit + ts) && (t < sinit + tf)) {
								eps = 1./(pop_sel*factpop[0]);
								if(eps > 0.1) eps = 0.1;
								freqend = xdt; ts = t;/*WORKS??*/
								/*INCORRECT?*/tf = (pop_sel*factpop[0]*(ts)+(double)log((double)((1.-freqend-eps+eps*freqend)/(eps*freqend))))/(pop_sel*factpop[0]);
							}
							/*if(t >= sinit + tf) do nothing, the selective event already happened*/
							/*end*/
							flagintn = 1;
							nevent += 1;
						}
					}
					if(flagintn == 0) {/*no change in Ne*/
						t += ttemp;
						if((t < sinit + ts && t - ttemp < sinit + ts) || (t >= sinit + tf && t - ttemp >= sinit + tf)) {                                            
							if(((prect) > (double)0) && (trm < tcoal)) {
								#if SELECTION_RECW_ALLOWED
								/* recombination */
								rchrom = re(nsam,weightrec,nsites,r);		
								config[chrom[rchrom].pop] += 1;
								#endif
							}
							else {
								/* coalescent event */
								if(config[0] == 0) pop = 1;
								else pop = 0;
								/* pick the two, c1, c2 */
								pick2_chrom(pop,config,&c1,&c2);	
								dec = ca(nsam,nsites,c1,c2,weightrec,r,tr);		
								config[pop] -= dec;
							}
						}
						else {
							if((t >= sinit + ts) && (t - ttemp < sinit + ts)) {
								/*starting selective process. We assume no recombination between selected-deselected occured before*/
								t = sinit + ts;
								xdt = 1. - eps;
							}
						}
					}
					flagintn = 0;
				}
				/************************* selective phase ***************************/
				else {
					/**/dt = (tf-ts)/((tf-ts)*(4.*(double)pop_size*factpop[0]));/**/
                    /*do each generation individually in discrete steps*/
					/*dt = 1./(100.*pop_sel/2.);*//*braverman*/
					ttemp = dt;
					
					xdt = (double)1 - eps/(eps+((double)1-eps)*(double)exp((double)(-pop_sel*factpop[0] * (t+dt-sinit-ts)))); 
					/*xdt = here add the contribution of drift. Assumed pop_size*factpop[0] is always large*/
					coal_prob = dt * factor_chrnall * (((double)config[0])*(config[0]-(double)1)/((double)   xdt) + 
													   ((double)config[1])*(config[1]-(double)1)/((double)1- xdt));
					rec_prob  = dt * tot_rec;
					while(0.1 < (coal_prob+rec_prob)) {/*in case two or more events can happen at the same time...*/
						dt =  dt * 0.1;
						xdt = (double)1 - eps/(eps+((double)1-eps)*(double)exp((double)(-pop_sel*factpop[0] * (t+dt-sinit-ts))));
						/*xdt = here add the contribution of drift. Assumed pop_size*factpop[0] is always large*/
						coal_prob = dt*factor_chrnall * (((double)config[0])*(config[0]-(double)1)/((double)   xdt) + 
														 ((double)config[1])*(config[1]-(double)1)/((double)1- xdt));
						rec_prob  = dt*tot_rec;
					}
					
					if(nevent < npop_events) {
						if((t <= pop_event[nevent].time_event) && (t+ttemp > pop_event[nevent].time_event)) { /*CHANGES IN POPULATION SIZE: TO CHECK!!*/
							t = pop_event[nevent].time_event /*+1E-17*//*precission error*/;/*modify time*/
							/*modify factor_pop*/
							factpop[pop_event[nevent].npop1/*0*/] = pop_event[nevent].factor_pop;
							/*modification to consider the selective event: under construction: ONLY FOR pop_sel*factor[0] > 10 ALWAYS*/
							if(t < sinit + ts) { 
								eps = 1./(pop_sel*factpop[0]);
								if(eps > 0.1) eps = 0.1;
								freqend = eps; ts = t;/*WORKS??*/
								/*INCORRECT?*/tf = (pop_sel*factpop[0]*(ts)+(double)log((double)((1.-freqend-eps+eps*freqend)/(eps*freqend))))/(pop_sel*factpop[0]);
							}
							if((t <= sinit + ts) && (t < sinit + tf)) {
								eps = 1./(pop_sel*factpop[0]);
								if(eps > 0.1) eps = 0.1;
								freqend = xdt; ts = t;/*WORKS??*/
								/*INCORRECT?*/tf = (pop_sel*factpop[0]*(ts)+(double)log((double)((1.-freqend-eps+eps*freqend)/(eps*freqend))))/(pop_sel*factpop[0]);
							}
							/*if(t >= sinit + tf) do nothing, the selective event already happened*/
							/*end*/
							flagintn = 1;
							nevent += 1;
						}
					}
					if(flagintn == 0) {/*no change in Ne*/
						if(t+ttemp < sinit+tf) {/*in case xdt <= eps and config[0] > 0, we join ind pop 0 to 1*/
							t += ttemp;
							/*the program would be faster if we keep the product of coal_prob+rec_prob until having an event*/
							if((rdum = (double)ran1()) < (coal_prob + rec_prob)) { /*then we have one event, otherwise nothing*/
								if((rdum = (double)ran1()) >= coal_prob/(coal_prob + rec_prob)) {/*rec or coal?*/
									/* recombination */
									while((rdum = (double)ran1()) == 1.0);
									if(rdum < (double)/*nlinks*/nlinksr/(double)/*total_nts*/total_ntsr) { /*recombination is only effective within nlinks*/
										#if SELECTION_RECW_ALLOWED
										rchrom = re(nsam,weightrec,nsites,r);
										while((rdum = (double)ran1()) == 1.0);
										rdum = freqend+rdum*((double)1-(double)eps-freqend);/*rdum goes from freqend to 1-eps*/
										if(rdum < xdt) chrom[new_chrom].pop = 0;
										else chrom[new_chrom].pop = 1;
										config[chrom[new_chrom].pop] += 1;
										#endif
									}
									else {/*recombination in one chrom outside studied region.We do changes only some times.Slow*/
										while((rdum = (double)ran1()) == 1.0);
										new_chrom = rdum * nchrom;	/*find the chrom*/
										config[chrom[new_chrom].pop] -= 1;/*we check the pop. We erase that from the former pop*/
										while((rdum = (double)ran1()) == 1.0);
										rdum = freqend+rdum*((double)1-(double)eps-freqend);/*rdum goes from freqend to 1-eps*/
										if(rdum < xdt) chrom[new_chrom].pop = 0;
										else chrom[new_chrom].pop = 1;
										config[chrom[new_chrom].pop] += 1;
										/*we assign the new pop. many times is the same than before*/
									}
								}
								else { /*coalescent*/
									while((ran = (double)ran1()) == 1.0);
									prob = (ttemp*((double)config[0])*(config[0]-(double)1)*factor_chrnall/xdt)/coal_prob;
									if(ran < prob) pop = 0;
									else pop = 1; 
									/* pick the two, c1, c2 */
									pick2_chrom(pop,config,&c1,&c2);	
									dec = ca(nsam,nsites,c1,c2,weightrec,r,tr);	
									config[pop] -= dec;
								}
							}
						}
						else {/*End selective phase. we cut at 1-eps.*/
							t = sinit + tf;
                            /*coalescent: in case xdt==eps, we force coalescent to have the mutation event at time sinit+tf*/
                            if(xdt <= eps) {
                                while(config[0] > 1) {
                                    pick2_chrom(pop=0,config,&c1,&c2);
                                    dec = ca(nsam,nsites,c1,c2,weightrec,r,tr);	
                                    config[pop] -= dec;
                                    t += 1E-17;
                                }
                            }
							for(new_chrom=0;new_chrom<nchrom;new_chrom++) {
								if(chrom[new_chrom].pop == 0) {
									config[0] -= 1;
									chrom[new_chrom].pop = 1;
									config[1] += 1;
									/*break; only in case forcing coalescent, we have only one ind with favorable mut*/
								}
							}
						}
					}
					flagintn = 0;
				}
                #if PRINT_TRACE_SELECTION == 1
                    fprintf(filexdt,"\n%f\t%f\t%ld\t%ld",t,xdt,config[0],config[1]);
                #endif
			}
            #if PRINT_TRACE_SELECTION == 1
               fprintf(filexdt,"\nNA\tNA\t%ld\t%ld",config[0],config[1]);
                fclose(filexdt);
            #endif
			/****************************************************************************************************/
			#elif SELECTION_CLASSIC == 2
            #if PRINT_TRACE_SELECTION == 1
                /*filexdt = fopen("/Users/sramos/Documents/mlcoalsim_project/mlcoalsim_output_sel2_xdt_all.txt","a");*/
                if(iteration == 1){
                    filexdt = fopen("./mlcoalsim_output_sel2_xdt_all.txt","w");
                    fprintf(filexdt,"\nTime\txdt\tnsam[0]\tnsam[1]");
                }
                else filexdt = fopen("./mlcoalsim_output_sel2_xdt_all.txt","a");
            #endif
            flag_freq = 0;
            while(nchrom > 1) {
				/*modification to increase the speed with high R, if t > tlimit r = 0.;*/
				if(t > tlimit) {
					r = (double)0; 
					nlinksr = (double)0;
				}				
				prect = (double)/*nlinks*/nlinksr /* * r*/;        /* 'tot_rec' is r in all region.*/
				tot_rec = (double)/*total_nts*/total_ntsr /* * r*/;   /*like prect but for all the region*/      
				
                /*selection process finish for a given custom freqend*/
                /*for freqend == eps continue the selection with drift*/
                if (flag_freq== 0 && ((xdt < freqend && freqend > eps) || t >= sendt)) {
                    flag_freq = 1;
                    #if PRINT_TRACE_SELECTION == 0
                    /*collect all lineages to non-selective_pop: join ind pop 0 to 1*/
                     for(new_chrom=0;new_chrom<nchrom;new_chrom++) {
                        if(chrom[new_chrom].pop == 0) {
                            config[0] -= 1;
                            chrom[new_chrom].pop = 1;
                            config[1] += 1;
                        }
                    }
                    xdt = 0.0;
                    /*#else*/
                    /*pop_sel = 1e-3;*//*like erasing selection force*/
                    #endif
                }
                
				/************ non-selective phase ****************************/
				if(xdt == 0.0 || t < sinit) {
					coal_prob = ((double)nchrom)*(nchrom-(double)1)*factor_chrnall/factpop[0];
					if(prect > (double)0) {
						while((rdum = (double)ran1()) == (double)1);	
						trm = -(double)log((double)1.0-(double)rdum)/(prect);	
					}
					while((rdum = (double)ran1()) == 1.0);
					tcoal = -(double)log((double)1.0 - (double)rdum)/coal_prob;
					if((prect) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
					else ttemp = tcoal;
					
					if(nevent < npop_events) {
						if((t <= pop_event[nevent].time_event) && (t+ttemp > pop_event[nevent].time_event)) { /*CHANGES IN POPULATION SIZE: TO CHECK!!*/
							t = pop_event[nevent].time_event /*+1E-17*//*precission error*/;/*modify time*/
							factpop[pop_event[nevent].npop1/*0*/] = pop_event[nevent].factor_pop;/*modify factor_pop*/
							flagintn = 1;
							nevent += 1;
						}
					}
					if(flagintn == 0) {/*no change in Ne*/
						t += ttemp;
						if(xdt == 0 || (t < sinit && t - ttemp < sinit)) {                                            
							if(((prect) > (double)0) && (trm < tcoal)) {
								#if SELECTION_RECW_ALLOWED
								/* recombination */
								rchrom = re(nsam,weightrec,nsites,r, nnodes, seglst, chrom);	
								config[chrom[rchrom].pop] += 1;
								#endif
							}
							else {
								/* coalescent event */
								if(config[0] == 0) pop = 1;
								else pop = 0;
								/* pick the two, c1, c2 */
								pick2_chrom(pop,config,&c1,&c2, chrom);	
								dec = ca(nsam,nsites,c1,c2,weightrec,r,tr, nnodes, seglst, chrom);
								config[pop] -= dec;
							}
						}
						else {
							if(t >= sinit && t - ttemp < sinit) 
								t = sinit;
						}
					}
					flagintn = 0;
				}
				/************************* selective phase ***************************/
				else {
                    /*added for considering sfreqinit!*/
                    if(t == sinit) {
						#if PARAMETER_FREQINIT == 1
						xdt = freqinit;
                        *all_sel = (int)binomialdist((float)xdt,(int)config[0]);
						#endif						
						#if PARAMETER_FREQINIT == 2
						/*in case using *all_sel as the 'parameter' then we have to estimate xdt from a CDF binomial*/
						*all_sel = (int)(freqinit * config[0]);
						if(xm==0) {
							itm = pow(2.0,20-2.0);
							xm=(double *)calloc(itm,sizeof(double));
							psum=(double *)calloc(itm,sizeof(double));
						}
						cummtot=trapzd_inv(distp_binom,0.0,1.0,20,(int)config[0],*all_sel,xm,psum);
						ran=ran1();
						ni=0;while(psum[ni]<ran)ni++;
						xdt=xm[ni];
						#endif
                        /*how many alleles have the favorable mutation at time sinit with freqinit*/
                        /*for(ind=0;ind<*all_sel;ind++) selnsam[ind] = (int)ind;*//*fav lines in a vector (will be useful in mig+sel)*/
                        for(ind=*all_sel;ind<(int)config[0];ind++) chrom[ind].pop = 1; /*chrom with unfav allele*/
                        config[1] = config[0] - *all_sel;
                        config[0] = *all_sel;
                    }
                    if(flag_freq == 1) {
                        /*Conditioning on fixing (losing in coalescence) the derived allele: change the effect of selection for neutral scenario: */
                        pop_sel = 2.0/(1.0-xdt); /*Solution from Zhao,Lascoux,Overall and Waxman, Genetics 2013, from paragraph after Equation 8*/
                        /*pop_sel is 4Ns and s is 1/(2Nx) following the conditional fixing, so pop_sel is 2/x*/
                    }
                    /**/

                    dt = 1./(4.*pop_size*factpop[0]);/*do each generation individually in discrete steps*/
					ttemp = dt;
					
					fitness = -pop_sel/(2.*pop_size); /*fitness: negative because it goes back in time*/
					avgftns = xdt * (xdt*(1.+fitness)+(1.-xdt)*(1.+fitness/2.)) + (1.-xdt) * (xdt*(1.+fitness/2.)+(1.-xdt)*1.);
                    /*average fitness*/
                    
					if(xdt > 1. - dt)
						xdt = 1. - dt;
 					xdt = xdt * (xdt*(1.0+fitness)+(1.0-xdt)*(1.0+fitness/2.0))/avgftns;
					xdt = largebinomialdist(xdt,pop_size*factpop[0])/(pop_size*factpop[0]);/*DRIFT*/
					
					while(xdt == 0. && config[0]>1) {
                        /*do all coalescences at the same time?*/
                        ttemp = ttemp/(config[0]-1);
                        t += ttemp;
                        while((ran = ran1()) == 1.0);
                        prob = coal_prob0/coal_prob;
                        if(ran < prob) pop = 0;
                        else pop = 1;
                        /* pick the two, c1, c2 */
                        pick2_chrom(pop,config,&c1,&c2, chrom);
                        dec = ca(nsam,nsites,c1,c2,weightrec,r,tr, nnodes, seglst, chrom);
                        config[pop] -= dec;
                        #if PRINT_TRACE_SELECTION == 1
                            fprintf(filexdt,"\n%.7f\t%f\t%ld\t%ld",t,xdt,config[0],config[1]);
                        #endif
                    }
                    if(xdt == 0 && config[0]==1) {
                        new_chrom = 0;
                        while(chrom[new_chrom].pop == 1) new_chrom++;
                        config[chrom[new_chrom].pop] -= 1;
                        chrom[new_chrom].pop = 1;
                        config[chrom[new_chrom].pop] += 1;
                        #if PRINT_TRACE_SELECTION == 1
                            fprintf(filexdt,"\n%.7f\t%f\t%ld\t%ld",t,xdt,config[0],config[1]);
                            fflush(filexdt);
                        #endif
                    }
                    
                    if(xdt)
                        coal_prob0 = dt * factor_chrnall * ((double)config[0])*(config[0]-(double)1)/((double)   xdt);
                    else coal_prob0 = (config[0]-(double)1) * 99999.;
                    if(1.-xdt)
                        coal_prob1 = dt * factor_chrnall * ((double)config[1])*(config[1]-(double)1)/((double)1- xdt);
                    else coal_prob1 = (config[1]-(double)1) * 99999.;
                    coal_prob  = coal_prob0 + coal_prob1;
                    rec_prob  = dt * tot_rec;
                    
					if(nevent < npop_events) {
						if((t <= pop_event[nevent].time_event) && (t+ttemp > pop_event[nevent].time_event)) { /*CHANGES IN POPULATION SIZE: TO CHECK!!*/
							t = pop_event[nevent].time_event /*+1E-17*//*precission error*/;/*modify time*/
							factpop[pop_event[nevent].npop1/*0*/] = pop_event[nevent].factor_pop;/*modify factor_pop*/
							flagintn = 1;
							nevent += 1;
						}
					}
					if(flagintn == 0) {/*no change in Ne*/
						t += ttemp;
						if((rdum = (double)ran1()) < (coal_prob + rec_prob)) { /*then we have one event, otherwise nothing*/
							if((rdum = (double)ran1()) >= coal_prob/(coal_prob + rec_prob)) {/*rec or coal?*/
								/* recombination */
								while((rdum = ran1()) == 1.0);
								if(rdum < nlinksr/total_ntsr) { /*recombination is only effective within nlinks*/
									#if SELECTION_RECW_ALLOWED
									rchrom = re(nsam,weightrec,nsites,r, nnodes, seglst, chrom);
									while((rdum = ran1()) == 1.0);
									rdum = freqend+rdum*(1-eps-freqend);/*rdum goes from freqend to 1-eps*/
									if(rdum < xdt) chrom[new_chrom].pop = 0;
									else chrom[new_chrom].pop = 1;
									config[chrom[new_chrom].pop] += 1;
									#endif
								}
								else {/*recombination in one chrom outside studied region.We do changes only some times.Slow*/
									while((rdum = ran1()) == 1.0);
									new_chrom = (long int)(rdum * nchrom);	/*find the chrom*/
									config[chrom[new_chrom].pop] -= 1;/*we check the pop. We erase that from the former pop*/
									while((rdum = ran1()) == 1.0);
									rdum = freqend+rdum*(1-eps-freqend);/*rdum goes from freqend to 1-eps*/
									if(rdum < xdt) chrom[new_chrom].pop = 0;
									else chrom[new_chrom].pop = 1;
									config[chrom[new_chrom].pop] += 1;
									/*we assign the new pop. many times is the same than before*/
								}
							}
							else { /*coalescent*/
								while((ran = ran1()) == 1.0);
								/*prob = (ttemp*((double)config[0])*(config[0]-(double)1)*factor_chrnall/xdt)/coal_prob;*/
								prob = coal_prob0/coal_prob;
								if(ran < prob) pop = 0;
								else pop = 1; 
								/* pick the two, c1, c2 */
								pick2_chrom(pop,config,&c1,&c2, chrom);	
								dec = ca(nsam,nsites,c1,c2,weightrec,r,tr, nnodes, seglst, chrom);
								config[pop] -= dec;
							}
						}
					}
					flagintn = 0;
				}
                #if PRINT_TRACE_SELECTION == 1
                    fprintf(filexdt,"\n%.7f\t%f\t%ld\t%ld",t,xdt,config[0],config[1]);
                #endif
            }
            #if PRINT_TRACE_SELECTION == 1
                fprintf(filexdt,"\nNA\tNA\t%ld\t%ld",config[0],config[1]);
                fclose(filexdt);
            #endif
			/************ IN CASE USING SELECTIVE APPROACH BY INTEGRATING THE DETERMINISTIC EQUATION *************/
            /*********************         !!!!!!!  NOT WORKING  !!!!!!!!!!!    *******************/
            /*********************         !!!!!!!  NOT WORKING  !!!!!!!!!!!    *******************/
            /*********************         !!!!!!!  NOT WORKING  !!!!!!!!!!!    *******************/
            /*********************         !!!!!!!  NOT WORKING  !!!!!!!!!!!    *******************/
            /*********************         !!!!!!!  NOT WORKING  !!!!!!!!!!!    *******************/
            /*********************         !!!!!!!  NOT WORKING  !!!!!!!!!!!    *******************/
            /*********************         !!!!!!!  NOT WORKING  !!!!!!!!!!!    *******************/
			#else
            #if PRINT_TRACE_SELECTION == 1
                if(iteration == 1){
                    filexdt = fopen("./mlcoalsim_output_sel2_xdt_all.txt","w");
                    fprintf(filexdt,"\nTime\txdt\tnsam[0]\tnsam[1]");
                }
            else filexdt = fopen("./mlcoalsim_output_sel2_xdt_all.txt","a");
            #endif
			while(nchrom > 1) {
				/*modification to increase the speed with high R, if t > tlimit, r = 0.;*/
				if(t > tlimit) {
					r = (double)0; 
					nlinksr = (double)0;
				}

				prect = (double)/*nlinks*/nlinksr /* * (double)r*/;        /* 'tot_rec' is r in all studied region.*/
				tot_rec = (double)/*total_nts*/total_ntsr /* * (double)r*/;   /*like prect but for all the region (in case selection is out from studied region is bigger)*/      

				/*recombination*/
				if(t < sinit + ts || t >= sinit + tf) 
					tot_rec = prect;/*in case no selection, we do not check recombinants in the non-studied region*/
				if(tot_rec > (double)0) {
					while((rdum = (double)ran1()) == 1.0);	
					trm = -(double)log((double)1-(double)rdum)/tot_rec;
				}
				/*selected alleles*/
				if(t < sinit + tf) 
					coal_prob0 = ((double)config[0])*(config[0]-(double)1)*factor_chrnall/factpop[0];
				else coal_prob0 = (double)0;
				if(coal_prob0 > (double)0) {
					while((rdum = (double)ran1()) == 1.0);				
					tcoal0 = -(double)log((double)1 - (double)rdum)/coal_prob0;
					/*integrating xdt in function of t, we calculate the time of coalescent directly*/
					if((t >= sinit + ts) && (t < sinit + tf)) {
						/*calculate the range*/
						x1 = 0;
						x2 = /*(double)sinit + */(double)tf;
						if(zbracn(functiont_freqp_sel,&x1,&x2,(double)tcoal0,(double)t,(double)ts,(double)sinit,(double)eps,(double)pop_sel*factpop[0],(double)0,(double)0,(double)0) == 1) {
							/*estimate the value of T*/
							xacc = (double)1e-6*x2; /* accuracy*/
							tcoal0 = (double)zriddrn(functiont_freqp_sel,x1,x2,xacc,(double)tcoal0,(double)t,(double)ts,(double)sinit,(double)eps,(double)pop_sel*factpop[0],(double)0,(double)0,(double)0);
                            /*possibly an error in functiont_freqp_sel ??*/
						}
					}
				}
				else tcoal0 = trm + 99999.;
				/*non-selected alleles*/
				if(t >= sinit + ts)
					coal_prob1 = ((double)config[1])*(config[1]-(double)1)*factor_chrnall/factpop[0];
				else coal_prob1 = (double)0;
				if(coal_prob1 > (double)0) {
					while((rdum = (double)ran1()) == 1.0);				
					tcoal1 = -(double)log((double)1 - (double)rdum)/coal_prob1;
					/*integrating 1-xdt in function of t, we calculate the time of coalescent directly: INCORRECT!!*/
					if((t >= sinit + ts) && (t < sinit + tf)) {
						/*calculate the range*/
						x1 = 0;
						x2 = /*(double)sinit + */(double)tf;
						if(zbracn(functiont_freqq_nsel,&x1,&x2,(double)tcoal1,(double)t,(double)ts,(double)sinit,(double)eps,pop_sel*factpop[0],0,0,0) == 1) {
							/*estimate the value of T*/
							xacc = (double)1e-6*x2; /* accuracy*/
							tcoal1 = (double)zriddrn(functiont_freqq_nsel,x1,x2,xacc,(double)tcoal1,(double)t,(double)ts,(double)sinit,(double)eps,pop_sel*factpop[0],(double)0,(double)0,(double)0);
                            /*possibly an error in functiont_freqq_nsel ??*/
						}
					}
				}
				else tcoal1 = trm + 99999.;
				
				/*decision*/
				tcoal = ((tcoal0 < tcoal1) ? tcoal0 : tcoal1); 
				if((tot_rec) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
				else ttemp = tcoal;
					
				if(nevent < npop_events) {
					if((t <= pop_event[nevent].time_event) && (t+ttemp > pop_event[nevent].time_event)) { /*CHANGES IN POPULATION SIZE: TO CHECK!!*/
						t = pop_event[nevent].time_event /*+1E-17*//*precission error*/;/*modify time*/
						/*modify factor_pop*/
						factpop[pop_event[nevent].npop1/*0*/] = pop_event[nevent].factor_pop;
						/*modification to consider the selective event: under construction: ONLY FOR pop_sel*factor[0] > 10 ALWAYS*/
						if(t < sinit + ts) { 
							eps = 1./(pop_sel*factpop[0]);
							if(eps > 0.1) eps = 0.1;
							freqend = eps; ts = t;
							tf = (pop_sel*factpop[0]*(ts)+(double)log((double)((1.-freqend-eps+eps*freqend)/(eps*freqend))))/(pop_sel*factpop[0]);
						}
						if((t <= sinit + ts) && (t < sinit + tf)) {
							eps = 1./(pop_sel*factpop[0]);
							if(eps > 0.1) eps = 0.1;
							freqend = xdt; ts = t;/*WORKS??*/
							/*INCORRECT?*/tf = (pop_sel*factpop[0]*(ts)+(double)log((double)((1.-freqend-eps+eps*freqend)/(eps*freqend))))/(pop_sel*factpop[0]);
						}
						/*if(t >= sinit + tf) do nothing, the selective event already happened*/
						/*end*/
						flagintn = 1;
						nevent += 1;
					}
				}
				if(flagintn == 0) {/*no change in Ne*/
					if(t+ttemp >= sinit + ts && t < sinit + ts) {
						/*starting selective process. We assume no recombination between selected-deselected occured before*/
						ttemp = sinit + ts - t + 1E-17/*precission error*/;
					}
					if((t+ttemp >= sinit + tf && t < sinit + tf && t >= sinit + ts) || (coal_prob0 == (double)0 && coal_prob1 == (double)0 && tot_rec == (double)0)) {
						/*finishing selective process at 1-eps...*/
						ttemp = sinit + tf - t + 1E-17/*precission error*/;
                        /*coalescent: in case xdt==eps, we force coalescent to have the mutation event at time sinit+tf*/
                        if(t+ttemp >= sinit + (pop_sel*(ts)+(double)log((double)((1.-freqend-eps+eps*eps)/(eps*eps))))/(pop_sel)) {
                            while(config[0] > 1) {
                                pick2_chrom(pop=0,config,&c1,&c2);
                                dec = ca(nsam,nsites,c1,c2,weightrec,r,tr);
                                config[pop] -= dec;
                                t += 1E-17;
                            }
                        }
						for(new_chrom=0;new_chrom<nchrom;new_chrom++) {
							if(chrom[new_chrom].pop == 0) {
								config[0] -= 1;
								chrom[new_chrom].pop = 1;
								config[1] += 1;
								/*break; only in case forcing coalescent, we have only one ind with favorable mut*/
							}
						}
					}
					
					/*modify the frequency of selective alleles: xdt changes from ts to tf*/
					if(t < sinit + tf) xdt = (double)1 - eps/(eps+((double)1-eps)*(double)exp((double)(-pop_sel*factpop[0] * (t+ttemp-sinit-ts))));
					else xdt = 0.;
					t += ttemp;  
					
					if(ttemp == trm){
						/* recombination */
						while((rdum = (double)ran1()) == 1.0);
						if(rdum < (double)/*nlinks*/nlinksr/(double)/*total_nts*/total_ntsr || (t < sinit + ts || t >= sinit + tf)) { 
							/*recombination is only effective within nlinks*/
							#if SELECTION_RECW_ALLOWED
							rchrom = re(nsam,weightrec,nsites,r);
							while((rdum = (double)ran1()) == 1.0);
							rdum = freqend+rdum*((double)1-(double)eps-freqend);/*rdum goes from freqend to 1-eps*/
							if(rdum < xdt) chrom[new_chrom].pop = 0;
							else chrom[new_chrom].pop = 1;
							config[chrom[new_chrom].pop] += 1;
							#endif
						}
						else {/*recombination in one chrom outside studied region.We do changes only some times.Slow*/
							while((rdum = (double)ran1()) == 1.0);
							new_chrom = (long int)((double)rdum * (double)nchrom);	/*find the chrom*/
							config[chrom[new_chrom].pop] -= 1;/*we check the pop. We erase that from the former pop*/
							while((rdum = (double)ran1()) == 1.0);
							rdum = freqend+rdum*((double)1-(double)eps-freqend);/*rdum goes from freqend to 1-eps*/
							if(rdum < xdt) chrom[new_chrom].pop = 0;
							else chrom[new_chrom].pop = 1;
							config[chrom[new_chrom].pop] += 1;
							/*we assign the new pop. many times is the same than before*/
						}
					}
					if(ttemp == tcoal) { /*coalescent*/
						pop = ((tcoal0 < tcoal1) ? 0 : 1);
						/* pick the two, c1, c2 */
						pick2_chrom(pop,config,&c1,&c2);	
						dec = ca(nsam,nsites,c1,c2,weightrec,r,tr);	
						config[pop] -= dec;
					}
				}
				flagintn = 0;
                #if PRINT_TRACE_SELECTION == 1
                    fprintf(filexdt,"\n%f\t%f\t%ld\t%ld",t,xdt,config[0],config[1]);
                #endif
            }
            #if PRINT_TRACE_SELECTION == 1
                ffprintf(filexdt,"\nNA\tNA\t%ld\t%ld",config[0],config[1]);
                fclose(filexdt);
            #endif
			#endif
			free(factpop);
		}    
		/*printf("\n%f\t%f",-log((double)ran1())/2.,t);*//*time under coalescent 2 indiv and total time under selection*/
		/************************ End of making tree *********************************/
		tr++;
	}
	*pnsegs = nsegs;
    
    free(config);
	free(nnodes);
	free(chrom);
	free(nref);
	free(intn);
	free(alphag);
	free(Ts);

    return seglst;
}
/*************************************************** Hudson routine *******************************************/
/* recombination */
int re(int nsam,double *weightrec,long int nsites,double r, long int *nnodes, struct segl *seglst, struct chromo *chrom)
{
    struct seg *pseg=0;
    long int /*el,spot,*/is;
    int lsg,ic;
	int lsgm1=0;
	double elr,isr,spotr;
    
    double ran1(void);
	void xover(int, int, long int,double *,long int,double, long int **, struct segl **, struct chromo **);
	long int localize_positionrec(double *,double,long int,long int,double);
    
    /* First generate a random x-over spot, then locate it as to chrom and seg */
	spotr = ((double)nlinksr * ran1());
	/* get chromosome number (ic) */
    for(ic=0;ic<nchrom;ic++) { 	/* Busca els valors MÀXIM i MÍNIM de l'individu escollit. (Què nt. segment i individu) */
        lsg = (int)chrom[ic].nseg;	/* el nombre de segments a l'individu ic */
        lsgm1 = lsg - 1;	/* de fet lsg-1 conté la informació de l'ultim segment */
        pseg = chrom[ic].pseg;	/* punter al primer segment de l'individu ic */
		elr = (weightrec[((pseg+lsgm1)->end)] - weightrec[(pseg->beg)])*r;/* mida del segments (max-min) a l'individu ic*/
        if(spotr <= elr) break; 	/* anem restant 'el' fins trobar l'individu */
		if(ic==nchrom-1) {/*precission problem?*/
			nlinksr -= spotr-elr;
			spotr = elr;
			break;
		}
		spotr -= elr;
    }
    isr = weightrec[pseg->beg]*r + spotr; 	/* posició dins l'individu ic */
	is  = localize_positionrec(weightrec,(double)isr,pseg->beg,(pseg+lsgm1)->end,r);
    xover(nsam,ic,is,weightrec,nsites,r, &nnodes, &seglst, &chrom);
    return ic;
}

// ic: the chromosome position in chrom
void xover(int nsam,int ic,long int is,double *weightrec,long int nsites,double r, long int **nnodes, struct segl **seglst, struct chromo **chrom)
{
    struct seg *pseg,*pseg2;
    long int i,lsg,lsgm1,newsg,jseg,k,in;
    double len,lenr;
    double ran(void);
    struct chromo chromosome = (*chrom)[ic];

    pseg = chromosome.pseg;	/* punter al primer segment de l'individu ic */
    lsg  = chromosome.nseg;	/* nombre de segments de ic */
    len  = (double)(pseg + lsg-1)->end - (double)pseg->beg;	/* max-min de l'individu ic */
	lenr = (weightrec[(pseg + lsg-1)->end] - weightrec[pseg->beg])*r;	/* weighted max-min de l'individu ic */

    cleft -= 1.0 - pow(pc, lenr);    /* per conversió */
	if(cleft < 0.) cleft = 0.;
    /* get segment number (jseg)*/
    for(jseg=0;is >= (pseg+jseg)->end;jseg++);	/* Busca el segment a on es troba la recombinació */
    if(is >= (pseg+jseg)->beg) in = (long int)1;		/* mira si la recombinació es troba entre els segments o dins el segment */
    else in = (long int)0;				/* in=0 entre segments, separació d'individus, pero no es fa un nou arbre */
    newsg = lsg - jseg; 			/* Això indica el desplaçament dels fragments de l'individu */
    
    /* copy LAST part of chrom to nchrom */
    nchrom++; 					/* fem un individu més, variable externa, val per tot el fitxer */
    if((long int)nchrom >= (long int)maxchr) {
        maxchr += 50;				/* afegeix 50 independent individus cada vegada que hem d'ampliar */
        if(!(*chrom = (struct chromo *)realloc(chrom,(long int)(maxchr*sizeof(struct chromo)))))
            perror("realloc error. xover.1");
    }
    if(!(pseg2 = (*chrom)[nchrom-1].pseg = (struct seg *)calloc((unsigned)newsg,sizeof(struct seg))))/* pseg2 apunta a chrom */
        perror("calloc error. xover.2");	/* pseg2 apunta a chrom[nchrom-1].pseg. Són newsg segments nous */
    (*chrom)[nchrom-1].nseg = newsg; 		/* el nou individu te newsg segments, els de la dreta */	
    (*chrom)[nchrom-1].pop = chromosome.pop;	/* la mateixa població de la qual prové, es clar *//*pero no en cas de seleccio*/
    pseg2->end = (pseg+jseg)->end;		/* el primer segment pseg2->end apunta a final del segment de jseg, ok */
    
    if(in) {					/* només al cas que haguem de fer nous arbres */
        pseg2->beg = is + (long int)1;			/* pseg2->beg és a on hi ha hagut la rec. + 1 */
        (pseg+jseg)->end = is;			/* aleshores (pseg+jseg) acaba a is, el punt de rec., clar */
    }
    else pseg2->beg = (pseg+jseg)->beg;		/* si no es fan nous arbres, l'inici de pseg2 és a l'inici del segment */
    
    pseg2->desc = (pseg+jseg)->desc;		/* el nombre del desc és el mateix que el de pseg+jseg. */
    for(k=1;k<newsg;k++) {
        (pseg2+k)->beg = (pseg+jseg+k)->beg;	/* creant tots els segments per la dreta del nou individu */
        (pseg2+k)->end = (pseg+jseg+k)->end;
        (pseg2+k)->desc = (pseg+jseg+k)->desc;	/* el descendent de seg indica l'individu de tree. IMPORTANT ! */
    }
    lsg = chromosome.nseg = lsg - newsg + in;	/* el nombre de segments de ic és jseg més in (1 o 0) */
    lsgm1 = lsg - 1;				/* l'últim fragment és lsgm1 */
    nlinksr/*nlinks*/ -= (weightrec[pseg2->beg] - weightrec[(pseg+lsgm1)->end])*r;
	if(nlinksr < (double)1E-07)nlinksr =0.0;
    /* posicions a recombinar: restem el principi d'un segment amb el final de l'altre: normalment resta només un nt. pero 
    els max i min canvien depenent dels segments, aixi que es poden restar molts més nt. si cau entre segments */
    
     /*in case SELECTION*/
	if(ifsel_glob) {
        total_ntsr/*total_nts*/ -= (weightrec[(long int)pseg2->beg] - weightrec[(long int)(pseg+lsgm1)->end])*r;/*restem entre zones recombinants*/
        if(sel_nts_glob > (long int)(pseg+lsgm1)->end) {
			if(sel_nts_glob >= (long int)nsites) {
				total_ntsr += (weightrec[nsites-1] - weightrec[(pseg+lsgm1)->end])*r;
				total_ntsr += (double)((sel_nts_glob-1) - (long int)(nsites-1))*r/(double)nsites;
			}
			else total_ntsr += (weightrec[sel_nts_glob] - weightrec[(pseg+lsgm1)->end])*r;
		}
        /*sumem si sel_nts dreta del chrom esquerra*/
        if(sel_nts_glob < (long int)pseg2->beg) {
			if(sel_nts_glob < 0) {
				total_ntsr += (weightrec[pseg2->beg] - weightrec[0])*r;
				total_ntsr += 0.0 - (double)sel_nts_glob*r/(double)nsites;
			}
			else total_ntsr += (weightrec[pseg2->beg] - weightrec[sel_nts_glob])*r;
		}
		if(total_ntsr < (double)1E-07)total_ntsr = 0.0;
        /*sumem si sel_nts esquerra del chrom dreta*/
        if(sel_nts_glob < (long int)is) 
            new_chrom = nchrom-1;/*quin es el chrom separat de sel_nts? assignar a 'new_chrom'*/
        else 
            new_chrom = ic;
    }
    
    lenr = ((double)weightrec[(pseg+lsgm1)->end] - weightrec[pseg->beg])*r;	/* conversió */
    cleft += 1.0 - pow(pc,lenr);    		/* conversió */
	lenr = (weightrec[(pseg2 + newsg - 1)->end] - weightrec[pseg2->beg])*r;/* llargada del nou individu */
    cleft += 1.0 - pow(pc,lenr);   		/* conversió */
    
	if(!(chromosome.pseg = (struct seg *)realloc(chromosome.pseg,(long int)(lsg*sizeof(struct seg)))))
        perror("realloc error. xover.3");	/* només es deixen lsg(=jseg+in) segments a ic */
    if(in) {
        begs = (long int)pseg2->beg;	/* inici del primer segment */
        for(i=0,k=0;(k < (long int)nsegs-1) && (begs > (long int)(*seglst)[(*seglst)[i].next].beg-1);i= (*seglst)[i].next, k++);
        /* arribem a on es troba begs, i tenim el nombre del segment (i)-> o arribem al final, l'ultim segment */
        if(begs != (long int)(*seglst)[i].beg) {	/* en cas que no hi hagi hagut recombinació al mateix lloc */
            /* new tree */
            if(nsegs >= seglimit) {	
                seglimit += SEGINC;
                if(!(*nnodes = (long int *)realloc(*nnodes,(unsigned)(sizeof(long int)*seglimit))))
                    perror("realloc error. xover.4");
                if(!(*seglst = (struct segl *)realloc(seglst,(unsigned)(sizeof(struct segl)*seglimit))))
                    perror("realloc error. xover.5");
            }
            (*seglst)[nsegs].next = (*seglst)[i].next;	/* Crear un segment entre els altres segments */
            (*seglst)[i].next = nsegs;			/* MOLT BO ! */
            (*seglst)[nsegs].beg = (long int)begs;
            if(!((*seglst)[nsegs].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node))))
                perror("calloc error. xover.6");	/* crear els nodes de l'arbre del nou segment */
            *nnodes[nsegs] = *nnodes[i];			/* el nombre de nodes de nsegs és el mateix que i */
            ptree1 = (*seglst)[i].ptree;			/* punter a l'arbre de i */
            ptree2 = (*seglst)[nsegs].ptree;		/* punter a l'arbre de nsegs */
            nsegs++;					/* nsegs és un segment més gran */
            for(k=0;k<=*nnodes[i];k++) {			/* donem els mateixos valors a i que a nsegs */
                (ptree2+k)->abv = (ptree1+k)->abv;
                (ptree2+k)->time = (ptree1+k)->time;
            }
        }
    }
}
/* coalescent functions */
void pick2_chrom(long int pop,long int *config,long int *pc1,long int *pc2, struct chromo *chrom)
{
    long int c1,c2,cs,cb,i,count;
    void pick2(long int,long int *,long int *);
    double ran1(void);
    
    pick2(config[pop],&c1,&c2);	/* trobar els dos individus de la població pop que tindran coalescència, nombre c1 i c2 */
    cs = (c1 > c2) ? c2 : c1;	/* Ara buscar quins individus pertanyen a pop, i trobar dins pop quin es c1 i c2 */
    cb = (c1 > c2) ? c1 : c2;	/* cs és el mínim i cb el màxim entre c1 i c2 */
    i = count = 0;
    for(;;) {
        while(chrom[i].pop != pop) i++;
        if(count == cs) break;	/* anem buscant individus de pop, quan trobem, count++ fins cs, i ja tenim l'individu */
        count++;
        i++;
    }
    *pc1 = i;
    i++;
    count++;
    for(;;) {
        while(chrom[i].pop != pop) i++;
        if(count == cb) break;	/* igual per cb */
        count++;
        i++;
    }
    *pc2 = i;
}
void pick2(long int n, long int *i,long int *j)
{
    double ran1(void);
   
    *i = (long int)floor((double)n*ran1());
    while((*j = (long int)floor((double)n * ran1())) == *i); /* dos valors aleatoris entre 0 i n-1 */
}
/* coalescent */
int ca(int nsam,long int nsites,long int c1, long int c2,double *weightrec,double r,long int tr, long int *nnodes, struct segl *seglst, struct chromo *chrom)
{
    int yes1,yes2,seg1,seg2;
	long int seg;
    long int start,end;
    long int tseg,desc,k;
    struct seg *pseg;
    struct node *ptree;
    int isseg(long int,long int,int *, struct chromo *);
    double linksr(long int,double *,double, struct chromo *chrom);
    double calc_total_ntsr(long int,double *,long int,double, struct chromo *);
     
	seg1=0;	/* valor de la primera posició del segment actiu de l'individu c1 */
    seg2=0;	/* valor de la primera posició del segment actiu de l'individu c2 */
    
    /*if(c1==12 && c2==14 && nchrom==18 && nsam==28) {printf("\nExited in ca: nsegs: %ld\n"); exit(1);}*/	
    if(!(pseg=(struct seg *)calloc((unsigned)nsegs,sizeof(struct seg)))) {	/* vector de segments pel nou node */
        printf("\nCalloc error in ca(pseg)");
		fflush(stdout);
		perror("calloc error. ca.1");
		exit(1);
	}
    tseg = -1;						/* nombre de segments del nou node */
    
	for(seg=tr,k=tr;k<(long int)nsegs;seg=seglst[seg].next,k++) {	/* mirem tots els segments */
        start = seglst[seg].beg;				/* 1a posició del segment que mirem */
        yes1  = isseg(start,c1,&seg1, chrom);	/* yes1=1 si el segment es troba dins c1 */
        yes2  = isseg(start,c2,&seg2, chrom);	/* yes2=1 si el segment es troba dins c2 */
        if(yes1 || yes2) {				/* si un dels dos té el segment */
            tseg++;					/* sumem el segment */
            (pseg+tseg)->beg = seglst[seg].beg;		/* l'inici del segment pel nou node és l'inici del segment actiu */
            end = (k < (long int)nsegs-1 ? (long int)seglst[seglst[seg].next].beg-(long int)1 : (long int)nsites-(long int)1);
            (pseg+tseg)->end = end;			/* i el final és el final del segment */
            
            if(yes1 && yes2) {				/* si els dos tenen el node hi ha COALESCÈNCIA */
                nnodes[seg]++;				/* per aquell segment sumem el nombre de nodes de l'arbre */
                if(nnodes[seg] >= (2*nsam-2)) tseg--;	/* en cas sigui el MRCA, aleshores restem el segment */
                else (pseg+tseg)->desc = (int)nnodes[seg];	/*si no,el nombre del node del segment de l'individu és nnodes[seg]*/
                ptree = seglst[seg].ptree;		/* punter a l'arbre de seglst[seg] */
                desc = (chrom[c1].pseg+seg1)->desc;	/* el node de c1 és desc */
                (ptree+desc)->abv = (int)nnodes[seg];	/* el node de desc apunta a nnodes[seg], el nou node */
                desc = (chrom[c2].pseg + seg2)->desc;	/* el node de c2 és desc */
                (ptree+desc)->abv = (int)nnodes[seg];	/* el node de desc apunta a nnodes[seg], el nou node */
                (ptree+nnodes[seg])->time = (double)t;		/* per últim, indicar el temps de la coalescència */
            }
            else (pseg+tseg)->desc = (yes1 ? (chrom[c1].pseg + seg1)->desc : (chrom[c2].pseg + seg2)->desc);
        }	/* en cas només un dels individus té el segment, indicar el desc de l'individu que el té */
    }

	nlinksr -= linksr(c1,weightrec,r, chrom);	/* mida de posicions de recombinació. restar els individus que tenen coal i sumar el nou */
	if(nlinksr < (double)1E-07)nlinksr =(double)0;
    if(ifsel_glob) total_ntsr -= calc_total_ntsr(c1,weightrec,nsites,r, chrom); /*in case SELECTION*/
	if(total_ntsr < (double)1E-07)total_ntsr=(double)0;
    cleft  -= (double)1 - pow(pc,(double)linksr(c1,weightrec,r, chrom));	/* conversió */
	if(cleft < 0.) cleft = 0.;
    free(chrom[c1].pseg);	/* els segments de c1 ara ja no interessen */
    if(tseg < 0) {		/* en cas el nou node sigui MRCA per TOTS els segments que contenien els descendents */
        free(pseg);				/* en aquest cas tampoc interessa res de pseg */
        chrom[c1].pseg = chrom[nchrom-1].pseg;	/* assignem les direccions de c1 a l'últim cromosoma, c1 queda lliure */
        chrom[c1].nseg = chrom[nchrom-1].nseg;
        chrom[c1].pop  = chrom[nchrom-1].pop;
        if(c2==nchrom-1) c2 = c1;		/* c2 també queda liure, si c2 és l'últim, c2=c1 */
        nchrom--;				/* un individu menys! */
    }
    else {			/* en cas no hi hagi MRCA per tots segments */
        if(!(pseg = (struct seg *)realloc(pseg,(unsigned)((tseg+1)*sizeof(struct seg)))))
            perror("realloc error. ca.2");
        chrom[c1].pseg = pseg;			/* ara c1 passa a tenir els valors de pseg */
        chrom[c1].nseg = tseg + 1;
        nlinksr += linksr(c1,weightrec,r, chrom);
		if(nlinksr < (double)1E-07)nlinksr =(double)0;
        if(ifsel_glob) total_ntsr += calc_total_ntsr(c1,weightrec,nsites,r, chrom); /*in case SELECTION*/
		if(total_ntsr < (double)1E-07)total_ntsr=(double)0;
        cleft += 1.0 - pow(pc,(double)linksr(c1,weightrec,r, chrom)); /* conversió */
    }

    nlinksr -= linksr(c2,weightrec,r, chrom);	
	if(nlinksr < (double)1E-07)nlinksr = 0.0;
    if(ifsel_glob) total_ntsr -= calc_total_ntsr(c2,weightrec,nsites,r, chrom); 	/*in case SELECTION*/
	if(total_ntsr < (double)1E-07) total_ntsr=0.0;
    cleft  -= 1.0 - pow(pc,(double)linksr(c2,weightrec,r, chrom)); /* conversió */
	if(cleft < 0.) cleft = 0.;
    free(chrom[c2].pseg);			/* eliminem c2 i apuntem a l'últim individu */
    chrom[c2].pseg = chrom[nchrom-1].pseg;
    chrom[c2].nseg = chrom[nchrom-1].nseg;
    chrom[c2].pop  = chrom[nchrom-1].pop;
    nchrom--;					/* un individu menys! */
    if(tseg<0)
		return 2; /* decrease of nchrom is two */
    else
		return 1;
}

/* isseg: does chromosome c contain the segment on seglst which starts at start? *psg is the segment of chrom[c] at wich one is to begin looking. */
int isseg(long int start, long int c, int *psg, struct chromo *chrom) /**psg és inicialment 0, però varia quan anem avançant pels segments */
{
    int ns;
    struct seg *pseg;
    
    ns = (int)chrom[c].nseg;		/* nombre de segments de l'individu c */
    pseg = chrom[c].pseg;	/* punter al primer segment de l'individu c */
    
    /*Sylvain diu que es incorrecte, aleshores: for(;((*psg) < ns) && ((pseg+(*psg))->beg <= start);  ++(*psg)) */
    /*for(;((pseg+(*psg))->beg <= start) && ((*psg) < ns); ++(*psg))*//*des del segment psg fins que sigui més gran de start*/
    for(;((*psg) < ns) && ((pseg+(*psg))->beg <= start);  ++(*psg)) /*Sylvain modification*/
        if((pseg+(*psg))->end >= start) return(1);	/* psg ja és 1+ perque ++(*psg), i no (*psg)++ ?? NO POT SER */
        /* en cas el final del segment sigui més gran o igual start, tenim el segment buscat a l'individu c */
        /* com c pot tenir segments més grans, start ha d'estar entre o = begin i end, aleshores el segment és inclós */
    return(0);
}
double calc_total_ntsr(long int c,double *weightrec,long int nsites,double r, struct chromo *chrom)	/*IN CASE SELECTION*/
{
    double linksr(long int,double *,double, struct chromo *);
    double lenr;
    int ns;
	
    lenr = linksr(c,weightrec,r, chrom);					/*afegim regio estudiada*/
    ns = (int)chrom[c].nseg - 1;

	if(sel_nts_glob > (long int) (chrom[c].pseg + ns)->end) {
		if(sel_nts_glob >= (long int)nsites) {
			lenr += (weightrec[nsites-1] - weightrec[(chrom[c].pseg + ns)->end])*r;
			lenr += (double)((sel_nts_glob-1) - (long int)(nsites-1))*r;
		}
		else lenr += (weightrec[sel_nts_glob] - weightrec[(chrom[c].pseg + ns)->end])*r;
	}
	/*sumem si sel_nts dreta del chrom esquerra*/
	if(sel_nts_glob < (long int)(chrom[c].pseg)->beg) {
		if(sel_nts_glob < 0) {
			lenr += (weightrec[(chrom[c].pseg)->beg] - weightrec[0])*r;
			lenr += (double)0 - (double)sel_nts_glob*r/(double)nsites;
		}
		else lenr += (weightrec[(chrom[c].pseg)->beg] - weightrec[sel_nts_glob])*r;
	}	
    return lenr;
}

double linksr(long int c,double *weightrec,double r, struct chromo *chrom)
{
    int ns;
    ns = (int)chrom[c].nseg - 1;
    return((weightrec[(chrom[c].pseg + ns)->end] - weightrec[(chrom[c].pseg)->beg])*r);	/* max - min de l'individu c in recombinational units*/
}

/* No pot encara treballar amb seleccio ni amb heterogeneuous mutation/recombination rates.... */
int cleftr(int nsam,double *weightrec,long int nsites,double r, long int *nnodes, struct segl *seglst, struct chromo *chrom)
{
	/*conversion coming from the left region, the effect is like recombination: Warning: if no recombination rate, is not working*/
    struct seg *pseg;
    int ic,lsgm1;
    double x, sum;
    long int /*len,*/is;
    void xover(int, int, long int,double *,long int,double, long int **, struct segl **, struct chromo **);
    double linksr(long int,double *,double, struct chromo *);
    double ran1(void);
	double isr,lenr;
	long int localize_positionrec(double *,double,long int,long int,double);
	
    while((x = (double)cleft*ran1()) == 0.0);
    sum = 0.0;
    ic = -1;

    while(sum < x) {
        sum += (double)1.0 - (double)pow((double)pc,(double)linksr((long int)++ic,weightrec,r, chrom)); /*trobem l'individu ic*/
    }
    pseg = chrom[ic].pseg; 	/*pseg apunta a l'individu ic*/
	lsgm1 = (int)chrom[ic].nseg - 1;
    lenr/*len*/ = linksr((long int)ic,weightrec,r, chrom);		/*mirem la llargada del max al min de ic*/
    isr = (weightrec[pseg->beg] + weightrec[(long int)floor((double)(1.0 + log((1.0 - (1.0- pow( pc, lenr))*ran1()))/lnpc))-1])*r; /*localitzem el punt de rec a is*/
    is  = localize_positionrec(weightrec,(double)isr,pseg->beg,(pseg+lsgm1)->end,r);
	xover(nsam,ic,is,weightrec,nsites,r, &nnodes, &seglst, &chrom);		/*recombinacio*/
    return ic;			/*a ic*/
}
/* No pot encara treballar amb seleccio ni amb heterogeneuous mutation/recombination rates........ */
int cinr(int nsam, long int nsites,double *weightrec,double r,long int tr, long int *nnodes, struct segl *seglst, struct chromo *chrom)
{
    /*conversion from a point to right (when is not finish inside the studied region is like recombination). Warning: if no recombination rate, is not working*/
	struct seg *pseg=0;
    int lsg, ic;
	int lsgm1 = 0;
    long int /*spot,el,*/len,is,endic;
	double spotr,lenr,elr,isr;
    
    double ran1(void);
	long int localize_positionrec(double *,double,long int,long int,double);

    spotr = /*(long int)floor*/((double)(/*nlinks*/nlinksr * ran1()));
    /* get chromosome number (ic) */
    for(ic=0;ic<nchrom;ic++) { 	/* Busca els valors MÀXIM i MÍNIM de l'individu escollit. (Què nt. segment i individu) */
        lsg = (int)chrom[ic].nseg; /* el nombre de segments a l'individu ic */
        lsgm1 = lsg - 1;/* de fet lsg-1 conté la informació de l'ultim segment */
        pseg = chrom[ic].pseg;/* punter al primer segment de l'individu ic */
        /*el = ((pseg+lsgm1)->end) - (pseg->beg);*//* mida dels segments (max-min) a l'individu ic*/
		elr = (weightrec[((pseg+lsgm1)->end)] - weightrec[(pseg->beg)])*r;/* mida de tots els segments (max-min) a l'individu ic*/
        if(spotr <= elr) break;/* anem restant 'el' fins trobar l'individu */
		if(ic==nchrom-1) {/*precission problem?*/
			/*printf("%f\n",spotr-elr);*/
			nlinksr -= spotr-elr;
			spotr = elr;
			break;
		}
        spotr -= elr;
    }
    isr = weightrec[pseg->beg]*r + spotr; 	/* posició dins l'individu ic */
	is  = localize_positionrec(weightrec,(double)isr,pseg->beg,(pseg+lsgm1)->end,r);
	endic = (pseg+lsgm1)->end; 			/*posicio final de l'individu ic*/
    xover(nsam,ic,is,weightrec,nsites,r, &nnodes, &seglst, &chrom);

    lenr = /*(long int)floor*/((double)(1.0 + (double)log((double)ran1())/lnpc));	/*llargada del evente de conversio...*/
	len = localize_positionrec(weightrec,(double)lenr,(chrom[ic].pseg)->beg,((chrom[ic].pseg)+(chrom[ic].nseg - 1))->end,r);
    if(is+len >= endic) return(ic);  		/*si es mes llarg, es igual que rec, acabem*/
    if(is+len < (chrom[nchrom-1].pseg)->beg){	/*si es mes curt que l'inici del nou chrom*/
        ca(nsam,nsites,(long int)ic,(long int)nchrom-1,weightrec,r,tr, nnodes, seglst, chrom);			/*aleshores eliminem el nou chrom, coalescencia*/
        return(-1);					/*i no ha passat res (...?)*/
    }
    xover(nsam,(int)nchrom-1,is+len,weightrec,nsites,r, &nnodes, &seglst, &chrom);		/*tornem a recombinar el fragment a llargada is+len*/
    ca(nsam,nsites,(long int)ic,(long int)nchrom-1,weightrec,r,tr, nnodes, seglst, chrom);		/*... i fem coalescencia a l'extrem. conversio finalitzada*/
    return ic;					/*i tenim un mes*/
}

/*MORE FUNCTIONS*/

double functiont_freqp_sel(double x,double tcoal,double tt,double ts,double sinit,double eps,double pop_sel,double no1,double no2,double no3)
/*include no1 and no2 to have all functions with the same number of variables*/
{
	double fdT;
	no1=no2=no3;

	if(pop_sel*(tt+(double)x-sinit-(ts)) - (pop_sel*(tt+(double)0-sinit-(ts))) < 700.) {
		fdT = (eps/pop_sel)/((double)1-eps)*(double)exp(pop_sel*(tt+(double)x-sinit-(ts)) - (pop_sel*(tt+(double)0-sinit-(ts)))) + x - tcoal;
	}
	else fdT = 1E07;

	/*
	fdT  = ((eps*(double)exp((double)(pop_sel*(tt+(double)x-sinit-(ts))))/pop_sel)+(tt+(double)x)*((double)1-eps))/((double)1-eps);
	fdT -= ((eps*(double)exp((double)(pop_sel*(tt+(double)0-sinit-(ts))))/pop_sel)+(tt+(double)0)*((double)1-eps))/((double)1-eps);
	fdT -= tcoal;
	*/
	return fdT;
}

double functiont_freqq_nsel(double x,double tcoal,double tt,double ts,double sinit,double eps,double pop_sel,double no1,double no2,double no3)
/*include no1 and no2 to have all functions with the same number of variables*/
{
	double fdT;
	no1=no2=no3;

	if((-pop_sel*(tt+(double)x-sinit-(ts))) - (-pop_sel*(tt+(double)0-sinit-(ts))) < 700.) {
		fdT = ((eps-(double)1)/pop_sel)/eps *(double)exp((-pop_sel*(tt+(double)x-sinit-(ts))) - (-pop_sel*(tt+(double)0-sinit-(ts)))) + x - tcoal;
	}
	else fdT = 1E07;
	/*
	fdT  = (((eps-(double)1)*(double)exp((double)(-pop_sel*(tt+(double)x-sinit-(ts))))/pop_sel)+(tt+(double)x)*eps)/eps;
	fdT -= (((eps-(double)1)*(double)exp((double)(-pop_sel*(tt+(double)0-sinit-(ts))))/pop_sel)+(tt+(double)0)*eps)/eps;
	fdT -= tcoal;
	*/
	return fdT;
}

double functiont_logistic(double x,double tcoal,double t,double tpast0,double tpast1,double alphag,double nrec1,double npast,double nrecf,double Ts)
{
	double fdT,expw0,expw1,logo0,logo1;
	
	if(fabs(nrec1-npast) < (double)1E-08) {
		fdT = nrecf * (((t+(double)x)/npast)) - nrecf * (((t+(double)0)/npast)) - tcoal;
	}
	else {
		expw0 = -alphag*((t+Ts+(double)x)-(tpast0+tpast1)/(double)2);
		expw1 = -alphag*((t+Ts+(double)0)-(tpast0+tpast1)/(double)2);
		if(expw0 > (double)60) {
			logo0 = expw0;
			logo1 = expw1;
		}
		else {
			logo0 = (double)log((double)(npast+(double)exp((double)expw0)*nrec1));
			logo1 = (double)log((double)(npast+(double)exp((double)expw1)*nrec1));
		}
	
		fdT  =  (nrecf * (((t+(double)x)/npast) + ((nrec1-npast)*logo0)/(alphag*npast*nrec1)));
		fdT -=  (nrecf * (((t+(double)0)/npast) + ((nrec1-npast)*logo1)/(alphag*npast*nrec1)));
		fdT -=  tcoal;
	}
		
	return fdT;
}

long int localize_positionrec(double *categories,double valuer,long int start,long int end,double r)
{
	long int half;
	
	half = (long int)floor((double)(start+end)/(double)2);
	if(half == start) return half;
	
	if((double)valuer < categories[half]*r) half = localize_positionrec(categories,valuer,start,half,r);
	else if((double)valuer > categories[half]*r) half = localize_positionrec(categories,valuer,half,end,r);
	
	return half;
}

double correction_theta(double f,double sexratio) { 
	/*from classical equations of calculation of Ne given different ratios male/female*/	
	#if SEXRATIOA1 == 1 && SEXRATIOX1 == -1
		/*theta is defined by the value of 4N when the sexratio is 1.0*/
		if(f == (double)1.0)
			return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));
		if(f == (double)0.75) 
			return(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio)));
		if(f == (double)0.25) 
			return((1./(1.+sexratio))/2.);
		if(f == (double)-0.25) 
			return((sexratio/(1.+sexratio))/2.);/*mytochondrial*/
	#elif SEXRATIOA1 == 0 && SEXRATIOX1 == -1
		/*theta is defined by the value of 4N when the sexratio is whatever defined (factor=1 always for A)*/
		if(f == (double)1.0)
			return(1.0/(4.*sexratio/((1.+sexratio)*(1.+sexratio))) * 4.*sexratio/((1.+sexratio)*(1.+sexratio)));
		if(f == (double)0.75) 
			return(1.0/(4.*sexratio/((1.+sexratio)*(1.+sexratio))) * 9.*sexratio/((2.*sexratio+4.)*(1.+sexratio)));
		if(f == (double)0.25) 
			return(1.0/(4.*sexratio/((1.+sexratio)*(1.+sexratio))) * (1./(1.+sexratio))/2.);
		if(f == (double)-0.25) 
			return(1.0/(4.*sexratio/((1.+sexratio)*(1.+sexratio))) * (sexratio/(1.+sexratio))/2.);/*mytochondrial*/
	#elif SEXRATIOX1 == 1 && SEXRATIOA1 == -1
			/*theta is defined by the value of 3N when the sexratio is 1.0*/
			if(f == (double)1.33)
				return(4./3./**//*1./(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio)))*/ * 4.*sexratio/((1.+sexratio)*(1.+sexratio)));
			if(f == (double)1.00) 
				return(4./3./**//*1./(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio)))*/ * 9.*sexratio/((2.*sexratio+4.)*(1.+sexratio)));
			if(f == (double)0.33) 
				return(4./3./**//*1./(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio)))*/ * (1./(1.+sexratio))/2.);
			if(f == (double)-0.33) 
				return(4./3./**//*1./(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio)))*/ * (sexratio/(1.+sexratio))/2.);/*mytochondrial*/
	#elif SEXRATIOX1 == 0 && SEXRATIOA1 == -1
			/*theta is defined by the value of 3N when the sexratio is whatever defined (factor=1 always for X)*/
			if(f == (double)1.33)
				return(1./(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio))) * 4.*sexratio/((1.+sexratio)*(1.+sexratio)));
			if(f == (double)1.00) 
				return(1./(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio))) * 9.*sexratio/((2.*sexratio+4.)*(1.+sexratio)));
			if(f == (double)0.33) 
				return(1./(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio))) * (1./(1.+sexratio))/2.);
			if(f == (double)-0.33) 
				return(1./(9.*sexratio/((2.*sexratio+4.)*(1.+sexratio))) * (sexratio/(1.+sexratio))/2.);/*mytochondrial*/
	#endif
	
	return 1.;
}

double correction_rec(double f,double sexratio,int m) /*relative to theta, because theta correction is already changing the tree size*/
{
	#if SEXRATIOA1 == 1 && SEXRATIOX1 == -1
		if(f == (double)1.0) {
			if(m==1) {
				/*return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));*/
				return((1.+sexratio)/4.);
			}
			else { 
				/*return((sexratio/(1.+sexratio))/2.);*/
				return(1.);
			}
		}
		if(f == (double)0.75) {
			/*return((sexratio/(1.+sexratio))/2.);*/
			return((2.*sexratio+4)/9.);
		}
		if(f == (double)0.25) 
			return(0.);
		if(f == (double)-0.25) /*mytochondrial*/
			return(0.);
	#elif SEXRATIOA1 == 0 && SEXRATIOX1 == -1 
	/*I am not sure*/
		if(f == (double)1.0) {
			if(m==1) {
				/*return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));*/
				return((1.+sexratio)/4.);
			}
			else { 
				/*return((sexratio/(1.+sexratio))/2.);*/
				return(1.);
			}
		}
		if(f == (double)0.75) {
			/*return((sexratio/(1.+sexratio))/2.);*/
			return((2.*sexratio+4)/9.);
		}
		if(f == (double)0.25) 
			return(0.);
		if(f == (double)-0.25) /*mytochondrial*/
			return(0.);
	#elif SEXRATIOX1 == 1 && SEXRATIOA1 == -1
		if(f == (double)1.33) {
			if(m==1) {
				/*return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));*/
				return((1.+sexratio)/4.);
			}
			else { 
				/*return((sexratio/(1.+sexratio))/2.);*/
				return(1.);
			}
		}
		if(f == (double)1.00) {
			/*return((sexratio/(1.+sexratio))/2.);*/
			return((2.*sexratio+4)/9.);
		}
		if(f == (double)0.33) 
			return(0.);
		if(f == (double)-0.33) /*mytochondrial*/
			return(0.);
	#elif SEXRATIOX1 == 0 && SEXRATIOA1 == -1
	/*I am not sure*/
		if(f == (double)1.33) {
			if(m==1) {
				/*return(4.*sexratio/((1.+sexratio)*(1.+sexratio)));*/
				return((1.+sexratio)/4.);
			}
			else { 
				/*return((sexratio/(1.+sexratio))/2.);*/
				return(1.);
			}
		}
		if(f == (double)1.00) {
			/*return((sexratio/(1.+sexratio))/2.);*/
			return((2.*sexratio+4)/9.);
		}
		if(f == (double)0.33) 
			return(0.);
		if(f == (double)-0.33) /*mytochondrial*/
			return(0.);
	#endif	
	
	return 1.;
}

