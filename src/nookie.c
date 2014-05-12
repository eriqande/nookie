/*!
\file nookie.c
\brief Source and main for combining genotypes of parents.
\author Eric C. Anderson
\date Created 2005-11-08
*/

/*!
\mainpage Nookie Program Source Code Documentation

Nookie is a simple program.  These writings provide the documentation for the 
source code.   You can use the tabs above to navigate to see the structures and 
functions defined in different files, but you are probably going to be best served by
just going straight to the <a href="./globals_func.html" list of functions </a> in the 
source distributed with Nookie.  Note that not all those functions are used in Nookie.  Many 
of the functions are just ones that appear in my various libraries.  

\section copy Copyright

I hope to ultimately be able to release it under the GPL.  But I 
have to figure out if that is possible as a Fed employee.  Also, I
must check with authors of some borrowed code for rational functions.


*/


/*!
This macro must be defined in exactly one source file.  I usually define it
in the source file that contains main().
*/
#define UN_EXTERN




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "ECA_MemAlloc.h"
#include "MathStatRand.h"
#include "ranlib.h"
#include "ECA_Opt3.h"

#define MAX_POPS 1000
#define MAX_ALLELES 1000
#define MAX_LOCI 1000
#define MAX_SPEC_MALES 300000
#define MAX_MALES_TOTAL 500000
#define MAX_SPEC_FEMALES 300000
#define MAX_FEMALES_TOTAL 3000000
#define MAX_CLUTCHES 50000
#define MAX_ADMIXED_INDS 10000
#define MAX_OUTPUT_STORE_STRING 10000


/*!
\brief  Stores the allele frequecies of a single population in the simulation

Not too much discussion of this thing yet.  
*/
typedef struct {
	
	/*! The index of the population.  Starting with 0, the population gets this index 
		according to the order in which it is read in on the command line. */
	int Idx;
	
	/*! The name of the population.  Maximum name length is 999 characters. */
	char Name[1000];
	
	/*! x The allele counts.  Subscripted by [l][k] for the k-th allele at the l-th locus */
	double **x;
	
	/*! dirichlet pars for simulating via a CDM dsn */
	double **xBP;
	
	/*! The allele frequencies.  Subscripted by [l][k] for the k-th allele at the l-th locus */
	double **p;
	
	/*! The sum of x's at each locus. Subscripted by [l]. */
	double *xsum;
	
	/*! The number of gene copies sampled at each locus.  (The rounded version of xsum.) Subscripted by [l]. */
	int *NumGC;
	
	/*! The number of loci */
	int NumLoc;
	
	
	/*! The number of alleles at each locus.  Subscripted by [l]. */
	int *NumAlle;
	
	/*! Dirichlet parameter factor */
	double Factor;
	
			
	int ScalePrior;
	
		
} pop_struct;

/*!
\brief Enum for whether a parents genotype was specified, or is to be drawn randomly from allele freqs

*/
typedef enum {
	SPEC,  /*!< The genotype was specified on the command line   */
	RANDO /*!< The genotype is to be randomly drawn, given allele freqs.   */
} IndGeno;



/*!
\brief Stores information about and genotypes of the parents to be combined
*/
typedef struct {
	int L; /*!< Number of loci---used to check for correct number in #SPEC males and females */
	int P; /*!< For RANDO individuals, to denote their population of origin */
	IndGeno G; /*!< Is the genotype pre-specified or not.  */
	int **Y; /*!< The actual genotype subscripted by [locus][gene copy] */
} parent_struct;


/*! 
\brief Stores information about how many offspring to make from each male and female
*/
typedef struct {
	int Mom; /*!< The number that identifies the female parent */
	int Dad; /*!< The number that identifies the the male parent */
	int NumKids; /*!< Number of offspring to produce in this clutch */
} clutch_struct;

/*!
\brief Stores information for the admixed-ind option.
*/
typedef struct {
	int Num;  /*!< Number of individuals with this Q-vector to make */
	double *Q; /*!< The "structure Q-vector" for each individual */ 
} admixed_ind_struct;


/* prototypes */
int GetNookie_Options(pop_struct **P, int argc, char *argv[], parent_struct *Dads, parent_struct *Moms, int *NumDads, int *NumMoms,
						clutch_struct *Clutches, int *NumClutches, int *NumReps, double *MissingRates,
						admixed_ind_struct **AdmInds, int *NumAdmInd);
void ReadLocFile(pop_struct **P, int L, char *FileName);
pop_struct *InitPopStruct(int I);
pop_struct *AllocPopStructForLoc(int NumLoc, int *NumAlle, char *PopName, int Idx);



/* Some globals */
int gAlleleNameIncrementer = 0;  /* can be used to add 1 to all allele names so that we get 1 instead of 0 for the 
									lowest index on an allele.  This helps if 0 should denote missing data in downstream
									analyses */









int main(int argc, char *argv[]) 
{
	int i,j,k,r,NumPops=0,NumDads=0,NumMoms=0, NumReps=1, NumLoc=0, NumClutches=0;
	parent_struct *Dads,*Moms;
	clutch_struct *Clutches;
	pop_struct **Pops=(pop_struct **)ECA_CALLOC(MAX_POPS, sizeof(pop_struct *));
	double *MissingRates=(double *)ECA_CALLOC(MAX_LOCI,sizeof(double));
	admixed_ind_struct **AdmInds = (admixed_ind_struct **)ECA_CALLOC(MAX_ADMIXED_INDS,sizeof(admixed_ind_struct *));
	int NumAdmInds=0;
	int Bail=0;
	char *ZoutStr = (char *)ECA_CALLOC(MAX_OUTPUT_STORE_STRING,sizeof(char));
	char *YoutStr = (char *)ECA_CALLOC(MAX_OUTPUT_STORE_STRING,sizeof(char));
	




	/**************************************************
	*
	*	Memory Allocation and Initialization
	*
	**************************************************/
	/* allocate memory and initialize space for the parents */
	Dads = (parent_struct *)ECA_CALLOC(MAX_MALES_TOTAL,sizeof(parent_struct));
	Moms = (parent_struct *)ECA_CALLOC(MAX_FEMALES_TOTAL,sizeof(parent_struct));
	for(i=0;i<MAX_MALES_TOTAL;i++)  {
		Dads[i].Y=NULL;
		Dads[i].L=0;
		Dads[i].G=RANDO;
	}
	for(i=0;i<MAX_FEMALES_TOTAL;i++)  {
		Moms[i].Y=NULL;
		Moms[i].L=0;
		Moms[i].G=RANDO;
	}
	
	/* allocate memory and initialize space for the clutch designations */
	Clutches = (clutch_struct *)ECA_CALLOC(MAX_CLUTCHES, sizeof(clutch_struct));
	for(i=0;i<MAX_CLUTCHES;i++)  {
		Clutches[i].Mom = -1;
		Clutches[i].Dad = -1;
		Clutches[i].NumKids = 0;
	}
	
	
	/**********************************************************
	*
	*	Getting all the user options
	*
	***********************************************************/
	NumPops = GetNookie_Options(Pops, argc, argv, Dads, Moms, &NumDads, &NumMoms, Clutches, &NumClutches, &NumReps, MissingRates, AdmInds, &NumAdmInds);
	
		
	/* set the RNG */
	SeedFromFile("nookie_seeds");
	
	
	
	
	/*************************************************************
	*
	*	Check the inputs for errors in number of loci and indexes
	*	of populations of origin and indexes of parents of clutches
	*
	***************************************************************/
	/* first check if there are any alle-freq-files specified.  If so,
	then they should all have the right number of loci. */
	if(NumPops>0) {
		NumLoc = Pops[0]->NumLoc;
	}
	/* then cycle over all moms, and inspect the SPEC ones */
	for(i=0;i<NumMoms;i++) {
		if(Moms[i].G == SPEC)  {
			if(NumLoc==0) 
				NumLoc = Moms[i].L;
			else if(NumLoc != Moms[i].L) {
				fprintf(stderr, "Error! Mom %d has mismatching locus number -- %d as opposed to %d\n",i+1,Moms[i].L,NumLoc);
				Bail=1;
			}
			printf("CHECKING: Mom %d NumLoc %d\n",i+1,Moms[i].L);
		}
		else if(Moms[i].P > NumPops) {
			fprintf(stderr, "Error! Mom %d has RANDO origin in population %d which is greater than NumPops=%d\n",i,Moms[i].P,NumPops);
			Bail=1;
		}
	}
	/* and do the same with the dads */
	for(i=0;i<NumDads;i++) {
		if(Dads[i].G == SPEC)  {
			if(NumLoc==0) 
				NumLoc = Dads[i].L;
			else if(NumLoc != Dads[i].L) {
				fprintf(stderr, "Error! Dad %d has mismatching locus number -- %d as opposed to %d\n",i+1,Dads[i].L,NumLoc);
				Bail=1;
			}
			printf("CHECKING: Dad %d NumLoc %d\n",i+1,Dads[i].L);
		}
		else if(Dads[i].P > NumPops) {
			fprintf(stderr, "Error! Dad %d has RANDO origin in population %d which is greater than NumPops=%d\n",i+1,Dads[i].P,NumPops);
			Bail=1;
		}
	}
	/* then cycle over all the clutches and make sure that the indexes of the parents of the 
	clutches are within the range given NumDads and NumMoms */
	for(i=0;i<NumClutches;i++)  {
		if(Clutches[i].Mom > NumMoms) {
			fprintf(stderr,"Error! Clutch %d calls for mother number %d, but there are only %d mothers\n",i+1,Clutches[i].Mom,NumMoms);
			Bail=1;
		}
		if(Clutches[i].Dad > NumDads) {
			fprintf(stderr,"Error! Clutch %d calls for father number %d, but there are only %d fathers\n",i+1,Clutches[i].Dad,NumDads);
			Bail=1;
		}
	}
	/* and while we are at it, we should print out the Missing data rates, and check to make sure they are all less than 1 */
	printf("MISSING_DATA_RATES: ");
	for(i=0;i<NumLoc;i++)  {
		printf(" %.6f",MissingRates[i]);
		if(MissingRates[i] >= 1.0 || MissingRates[i] < 0.0 ) {
			fprintf(stderr,"Error! Missing data rate at locus %d is %f which is >1 of <0.\n",i+1,MissingRates[i]);
			Bail=1;
		}
	}
	printf("\n");
	
	if(Bail>0) {
		fprintf(stderr,"Exiting with errors\n");
		exit(1);
	}
	
	/*************************************
	*
	*	Assemble the clutches...
	*
	*************************************/
	/* cycle over all the reps */
	for(r=0;r<NumReps;r++)  { char Str[1000];
		printf("REP_NUMBER: %d\n",r+1);
		
		/* then, go through the individuals and simulate the genotypes for the RANDO ones, and print them all */
		for(i=0;i<NumMoms;i++) {
			if(Moms[i].G==RANDO) {
				if(Moms[i].Y==NULL) { /* if it needs to have some memory allocated to it */
					Moms[i].Y = (int **)ECA_CALLOC(NumLoc,sizeof(int *));
					for(j=0;j<NumLoc;j++)  {
						Moms[i].Y[j] = (int *)ECA_CALLOC(2,sizeof(int));
					}
				}
				/* and then assign allelic types */
				for(j=0;j<NumLoc;j++)  {
					for(k=0;k<2;k++) {
						Moms[i].Y[j][k] = IntFromProbsRV(Pops[Moms[i].P-1]->p[j], 0, Pops[Moms[i].P-1]->NumAlle[j]);
					}
				}
			}
			/* and then print out their genotypes */
			if(Moms[i].G==RANDO) {
				sprintf(Str,"%d",Moms[i].P-1);
			}
			else {
				sprintf(Str,"S");
			}
			
			printf("Indiv_Fem_%d_%s   ",i+1,Str);
			for(j=0;j<NumLoc;j++)  {
				printf("   %d %d",Moms[i].Y[j][0] + gAlleleNameIncrementer, Moms[i].Y[j][1] + gAlleleNameIncrementer);
			}
			printf("\n");

		}
		
		
		/* go through the dad's like that too. */
		for(i=0;i<NumDads;i++) {
			if(Dads[i].G==RANDO) {
				if(Dads[i].Y==NULL) { /* if it needs to have some memory allocated to it */
					Dads[i].Y = (int **)ECA_CALLOC(NumLoc,sizeof(int *));
					for(j=0;j<NumLoc;j++)  {
						Dads[i].Y[j] = (int *)ECA_CALLOC(2,sizeof(int));
					}
				}
				/* and then assign allelic types */
				for(j=0;j<NumLoc;j++)  {
					for(k=0;k<2;k++) {
						Dads[i].Y[j][k] = IntFromProbsRV(Pops[Dads[i].P-1]->p[j], 0, Pops[Dads[i].P-1]->NumAlle[j]);
					}
				}
			}
			
			/* and then print out their genotypes */
			if(Dads[i].G==RANDO) {
				sprintf(Str,"%d",Dads[i].P-1);
			}
			else {
				sprintf(Str,"S");
			}
			
			printf("Indiv_Male_%d_%s   ",i+1,Str);
			for(j=0;j<NumLoc;j++)  {
				printf("   %d %d",Dads[i].Y[j][0] + gAlleleNameIncrementer, Dads[i].Y[j][1] + gAlleleNameIncrementer);
			}
			printf("\n");

		}
		
		/* and now we just cycle over the clutches and produce the offspring as needed */
		for(i=0;i<NumClutches;i++)  { int size; int ma; int pa; int a; int b;
			size = Clutches[i].NumKids;
			ma = Clutches[i].Mom - 1;
			pa = Clutches[i].Dad - 1;
			printf("CLUTCH: %d Size %d\n",i+1, size);
			
			/* print the mother's genotype */
			printf("MOTHER_Fem_%d  ",ma+1);
			for(j=0;j<NumLoc;j++)  {
				printf("   %d %d",Moms[ma].Y[j][0]+gAlleleNameIncrementer, Moms[ma].Y[j][1]+gAlleleNameIncrementer);
			}
			printf("\n");
			
			/* print the father's genotype */
			printf("FATHER_Male_%d  ",pa+1);
			for(j=0;j<NumLoc;j++)  {
				printf("   %d %d",Dads[pa].Y[j][0]+gAlleleNameIncrementer, Dads[pa].Y[j][1]+gAlleleNameIncrementer);
			}
			printf("\n");
			
			/* and then we assemble the kids */
			for(k=0;k<size;k++)  {
				printf("KID_%d  ",k+1); 
				for(j=0;j<NumLoc;j++)  {
					if(ranf()<MissingRates[j]) { /* make the genotypes missing at random according to MissingRates */
						a=-1+gAlleleNameIncrementer;
						b=-1+gAlleleNameIncrementer;
					}
					else {
						a = Moms[ma].Y[j][ UniformRV(0,1) ];  /* draw one of the alleles from ma */
						b = Dads[pa].Y[j][ UniformRV(0,1) ];  /* draw one of the alleles from pa */
					}
					printf("   %d %d",a+gAlleleNameIncrementer,b+gAlleleNameIncrementer);
				}
				printf("\n");
			}
			
			
		}
		
		/*************************************
		*
		*	Do the structure-like output.  Note
		*   that this is also inside of the Rep
		*	loop.
		*
		*************************************/
		for(i=0;i<NumAdmInds;i++)  {
			for(j=0;j<AdmInds[i]->Num;j++)  {  int a; int l; int Z; int Y; char sep;
				ZoutStr[0]='\0';  /* initialize to the empty string */
				YoutStr[0]='\0';
				for(l=0;l<NumLoc;l++)  {  /* cycle over loci */
					for(a=0;a<2;a++)  {  /* cycle over gene copies */
						/* choose a population of origin for the gene copy */
						Z = IntFromProbsRV(AdmInds[i]->Q,0,NumPops);  
						/* then choose an allelic type for the gene copy */
						Y = IntFromProbsRV(Pops[Z]->p[l],0,Pops[Z]->NumAlle[l]);
						/* then record those on the output strings */
						if(a==0) {
							sep=' ';
						}
						else {
							sep='/';
						}
							
						sprintf(ZoutStr,"%s%c%d",ZoutStr,sep,Z+1);
						sprintf(YoutStr,"%s%c%d",YoutStr,sep,Y+1);
					}
					/* add extra space between loci on the output strings */
					sprintf(ZoutStr,"%s  ",ZoutStr);
					sprintf(YoutStr,"%s  ",YoutStr);
				}
				
				/* output the individual in three lines.  The first shows his Q, the second the Z and the third the Y */
				printf("ADMIX_IND_Q: ");
				for(l=0;l<NumPops;l++)  {
					printf(" %.5f",AdmInds[i]->Q[l]);
				}
				printf("\n");
				printf("ADMIX_IND_Z: %s\n",ZoutStr);
				printf("ADMIX_IND_Y: %s\n",YoutStr);
				
			}
		}

	}
	
			
	SeedToFile("nookie_seeds");

	return(0);
}


/*! 
\brief  Allocate memory to and initialize a pop_struct

This function allocates memory to a pointer to pop_struct.  It
also sets all values within it to default values, and pointers within it to
NULL.
\param I The value that will be assigned to the Idx field of the created pop_struct.
*/
pop_struct *InitPopStruct(int I) 
{
	pop_struct *temp;
	
	temp = (pop_struct *)ECA_MALLOC(sizeof(pop_struct));
	
	temp->Idx = I;
	sprintf(temp->Name,"NoName");
	temp->NumLoc = 0;
	temp->x = NULL;
	temp->NumAlle = NULL;
	temp->Factor = 1.0;
	temp->ScalePrior = 0;
	

	
	return(temp);
}


/*! 
\brief Allocate memory to a population struct and the arrays for specific alleles.

This function allocates more memory than InitPopStruct.  It fills out all the arrays
using calloc (thus setting their values to zero.

\param NumLoc  the number of loci
\param NumAlle array of the number of alleles at all the loci
\param PopName name to be given the population
\param Idx population index to be given to the population

*/
pop_struct *AllocPopStructForLoc(int NumLoc, int *NumAlle, char *PopName, int Idx)
{
	int j;
	pop_struct *temp;
	
	temp = (pop_struct *)ECA_MALLOC(sizeof(pop_struct));
	
	temp->Idx = Idx;
	sprintf(temp->Name,"%s",PopName);
	temp->NumLoc = NumLoc;
	temp->NumAlle = (int *)ECA_CALLOC(NumLoc, sizeof(int));
	temp->NumGC = (int *)ECA_CALLOC(NumLoc, sizeof(int));

	
	temp->x = (double **)ECA_CALLOC(temp->NumLoc,sizeof(double*));
	temp->p = (double **)ECA_CALLOC(temp->NumLoc,sizeof(double*));
	temp->xBP = (double **)ECA_CALLOC(temp->NumLoc,sizeof(double*));
	temp->xsum = (double *)ECA_CALLOC(temp->NumLoc,sizeof(double));
	
	for(j=0;j<temp->NumLoc;j++)  {
		temp->NumAlle[j] = NumAlle[j];
		temp->x[j] = (double *)ECA_CALLOC(temp->NumAlle[j],sizeof(double));
		temp->p[j] = (double *)ECA_CALLOC(temp->NumAlle[j],sizeof(double));
		temp->xBP[j] = (double *)ECA_CALLOC(temp->NumAlle[j],sizeof(double));
	}
	
	
	temp->Factor = 1.0;
	temp->ScalePrior = 0;
	
	return(temp);
}


/*!
\brief Read allele count data out of file FileName

This function reads allele count data out of file 
FileName and inserts it into the appropriate fields in the P[L] which
is a pointer to a pop_struct.  In the process, it does many of the 
calculations necessary to fill fields such as "stp" and "xsum" in P[L]
from the raw allele count data in FileName.

In the process, this function prints a lot of information out to stdout informing
the user of the allele counts and functions thereof encountered in the file.

\param P an array of pointers to the pop_structs.  This function assumes that P[L]
has memory allocated to it and that it has been initialized.
\param L the element of P in which the allele count data should be inserted.
\param FileName the name of the file that contains the allele count data.
*/
void ReadLocFile(pop_struct **P, int L, char *FileName)
{
	int j,k;
	FILE *in;
	double SterSum;
	
	/* now get all the locus pars */
	if( (in=fopen(FileName,"r"))==NULL ) {
		fprintf(stderr,"\n\nCouldn't open file \"%s\" for locus information.\nExiting...\n\n",FileName);
		exit(1);
	}
	
	printf("ALLECOUNTS : opened file \"%s\" to get genetic data for population %d\n",FileName,L+1);
	while(eat_comments(in,'&')) ;
	/* get the number of loci */
	fscanf(in,"%d",&(P[L]->NumLoc) );
	
	if(L>0) {
		if(P[L]->NumLoc != P[L-1]->NumLoc) {
			fprintf(stderr,"Number of loci in file %s is %d, but for previous population was %d.  Exiting...\n",FileName,P[L]->NumLoc,P[L-1]->NumLoc);
			exit(1);
		}
	}
	
	printf("ALLEFREQS : number of loci in file %s is %d\n",FileName,P[L]->NumLoc);
	
	while(eat_comments(in,'&')) ;
	/* then get all the rest of the locus parameters */
	P[L]->NumAlle = (int *)ECA_CALLOC(P[L]->NumLoc, sizeof(int));
	
	while(eat_comments(in,'&')) ;
	P[L]->x = (double **)ECA_CALLOC(P[L]->NumLoc, sizeof(double *));
	P[L]->xBP = (double **)ECA_CALLOC(P[L]->NumLoc, sizeof(double *));
	P[L]->p = (double **)ECA_CALLOC(P[L]->NumLoc, sizeof(double *));
	P[L]->xsum = (double *)ECA_CALLOC(P[L]->NumLoc, sizeof(double));
	P[L]->NumGC = (int *)ECA_CALLOC(P[L]->NumLoc, sizeof(int));
	
	for(j=0;j<P[L]->NumLoc;j++)  {
		while(eat_comments(in,'&')) ;
		fscanf(in,"%d",&(P[L]->NumAlle[j]));
		
		if(L>0) {
			if(P[L]->NumAlle[j] != P[L-1]->NumAlle[j]) {
				fprintf(stderr,"Number of alleles in file %s at locus %d is %d, but for previous population was %d.  Exiting...\n",FileName,j+1,P[L]->NumAlle[j],P[L-1]->NumAlle[j]);
				exit(1);
			}
		}
		
		P[L]->x[j] = (double *)ECA_CALLOC(P[L]->NumAlle[j], sizeof(double));
		P[L]->xBP[j] = (double *)ECA_CALLOC(P[L]->NumAlle[j], sizeof(double));
		P[L]->p[j] = (double *)ECA_CALLOC(P[L]->NumAlle[j], sizeof(double));
		
		
		for(k=0;k<P[L]->NumAlle[j];k++)  { /* in this loop, we collect all of the allele counts and compute their sum. */
			while(eat_comments(in,'&')) ;
			fscanf(in,"%lf",&(P[L]->x[j][k]));
			
			/* add the x to the xsum */
			P[L]->xsum[j] += P[L]->x[j][k];
		}
		
		
		/* print the allele counts */
		printf("ALLE_COUNTS : Locus %d : %d Alleles : ", j+1,P[L]->NumAlle[j]);
		for(k=0;k<P[L]->NumAlle[j];k++)  {
			printf("%f ",P[L]->x[j][k]);
		}
		printf("\n");
		
		/* compute and print the allele frequencies and compute a sum for the stereographic projection */
		SterSum = 0.0;  /* initialize to get the sum of (p/k)^(1/2)  */
		printf("ALLE_FREQUENCIES : Locus %d : %d Alleles : ", j+1,P[L]->NumAlle[j]);
		for(k=0;k<P[L]->NumAlle[j];k++)  {
			P[L]->p[j][k] = P[L]->x[j][k] / P[L]->xsum[j];  /* compute the allele frequencies */
			printf("%f ",P[L]->p[j][k]);
			SterSum += sqrt( P[L]->p[j][k] / P[L]->NumAlle[j]);
		}
		printf("\n");
		
		/* make the int version of xsum */
		P[L]->NumGC[j] = (int)(P[L]->xsum[j] + .49999);	
		printf("NUMBER_OF_GENE_COPIES : Locus %d : %d  :  AsDouble= %f\n",j+1,P[L]->NumGC[j],P[L]->xsum[j]);
		
	}
	fclose(in);
}




 
/*!
\brief  Read the command line options and get input data.

This function reads the command line option in and does considerable amounts of
memory allocation, etc.  It uses the ECA_Opt options processing package.  
*/
int GetNookie_Options(pop_struct **P, int argc, char *argv[], parent_struct *Dads, parent_struct *Moms, int *NumDads, int *NumMoms,
						clutch_struct *Clutches, int *NumClutches, int *NumReps, double *MissingRates,
						admixed_ind_struct **AdmInds, int *NumAdmInd) 
{
	char locfilename[10000];
	int SpecFem_F=0,
		SpecMale_F=0,
		RandoFem_F=0,
		RandoMale_F=0,
		LocFileF=0,
		clutch_F=0,
		NumReps_F = 0,
		Missing_F=0,
		admixed_indF=0,
		FactorF=0,
		ScaleF=0,
		alle_incr_F=0;
		int CurrentPop=0; 
		int CurAdmInd=0;
	DECLARE_ECA_OPT_VARS;


	SET_OPT_WIDTH(25);
	SET_ARG_WIDTH(34);
	SET_PROGRAM_NAME("nookie");
	SET_PROGRAM_SHORT_DESCRIPTION("a program for making offpspring genotypes");
	SET_PROGRAM_LONG_DESCRIPTION(
		This is a silly little program that does Mendelian inheritance of genes from a mother
		and a father (or several mothers and fathers) into collections of offspring.  You can either
		specify the exact genotypes of males and females (where alleles are denoted by integers) and then
		mate specific males and females so that the produce a specified number of offspring. Or you can 
		simulate the genotypes of males and females from population allele frequencies as specified in one 
		or several allele frequency files (one for each population).  In this case\054 alleles are specified by
		the integers 0\0541\054...\054K-1 where K is the number of alleles at the locus. 
		
		\n\nOne major limitation of the program is that when you specify allele frequencies in an allele frequence
		file\054 you cannot specify the integer identifiers of those alleles\054 i.e.\054 they are numbered 
		0\0541\054...\054K-1.  This makes it complicated to mix spec-females and spec-males from a real population
		with individuals simulated from the estimated allele frequencies from that real population.  It can be done\054
		but it involves recoding allele lengths\054 etc.  I will fix this eventually.
		
		If the file nookie_seeds is  present in the current working directory then its contents are used to seed the random
		number generator.  If not\054 then seeds are generated from the current time.  The current state of the random number generator
		is written to nookie_seeds at the end of execution so it may be used over and over again in a sane fashion with 
		respect to the random number generator. 
		
		This has been an incomplete description.
	)
	SET_VERSION("VERSION: 1.0 Beta. 18 December 2005")
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson (eric.anderson@noaa.gov)");
	SET_VERSION_HISTORY(" No Version History ")
	
	/* allocate memory to the first population */
	P[0] = InitPopStruct(CurrentPop);

	BEGIN_OPT_LOOP 	 
		
		OPEN_SUBSET(Nookie Input Options, Nookie Input Options, These options control what sorts of data sets nookie will simulate.)
		
		if(MULT_USE_OPTION(
			Specific Female Parent,
			SpecFem_F,
			F,
			spec-fem,
			J1\0541 J1\0542 ... JL\0541 JL\0542,
			Create a female parent with genotype as specified,
			This command creates a female parent with genotype given by
			J1\0541 J1\0542 ... JL\0541 JL\0542.  i.e. the first allele at the first locus is
			has the index J1\0541; the second allele at the first locus has the index J1\0542; and so on up to 
			the L-th locus.  If there are K alleles at a locus in the total population then the alleles 
			will be specified as an integer between 0 and K-1 inclusive.  -1 is used for missing data., 
			MAX_SPEC_FEMALES) ) { int i; int temp;
				Moms[*NumMoms].G = SPEC;
				temp = COUNT_ARGS;
				if(temp<2 || temp%2==1) {
					fprintf(stderr, "Error! Female Number %d has %d args on the --spec-fem or -F option.\n",*NumMoms+1,temp);
					OPT_ERROR;
				}
				else {
					Moms[*NumMoms].L = temp/2;
					Moms[*NumMoms].Y = (int **)ECA_CALLOC(Moms[*NumMoms].L,sizeof(int *));
					for(i=0;i<Moms[*NumMoms].L;i++)  {
						Moms[*NumMoms].Y[i] = (int *)ECA_CALLOC(2,sizeof(int));
						Moms[*NumMoms].Y[i][0] = GET_INT;
						Moms[*NumMoms].Y[i][1] = GET_INT;
					}
					*NumMoms += 1;
				}
		}
		
		if(MULT_USE_OPTION(
			Specific Male Parent,
			SpecMale_F,
			M,
			spec-male,
			J1\0541 J1\0542 ... JL\0541 JL\0542 ,
			Create a male parent with genotype as specified,
			See documentation for --spec-fem option.  This is parallel to that., 
			MAX_SPEC_MALES) ) { int i; int temp;
				Dads[*NumDads].G = SPEC;
				temp = COUNT_ARGS;
				if(temp<2 || temp%2==1) {
					fprintf(stderr, "Error! Male Number %d has %d args on the --spec-fem or -F option.\n",*NumDads+1,temp);
					OPT_ERROR;
				}
				else {
					Dads[*NumDads].L = temp/2;
					Dads[*NumDads].Y = (int **)ECA_CALLOC(Dads[*NumDads].L,sizeof(int *));
					for(i=0;i<Dads[*NumDads].L;i++)  {
						Dads[*NumDads].Y[i] = (int *)ECA_CALLOC(2,sizeof(int));
						Dads[*NumDads].Y[i][0] = GET_INT;
						Dads[*NumDads].Y[i][1] = GET_INT;
					}
					*NumDads += 1;
				}
		}
		
		if(MULT_USE_OPTION(
			Random Female Parent,
			RandoFem_F,
			f,
			rand-fem,
			J1 J2,
			create J1 female parents randomly drawn from allele-freq-file J2,
			This will put J1 female parents on the female parent list and each one of them
			will be slated to get genotypes by randomly drawing from allele freqs as given
			in allele freq file J2.  Note that allele freq file 1 is the allele freq file
			designated in the first invocation of the -a or --allele-freq-file option.  Allele freq
			file 2 is the file designated in the second invocation of the -a or --allele-freq-file option
			and so on. This option may be issued 50 times., 
			50))
		{
			if(ARGS_EQ(2)) { int J1; int J2; int i;
				J1 = GET_INT;
				J2 = GET_INT;
				for(i=0;i<J1;i++)  {
					Moms[*NumMoms].G = RANDO;
					Moms[*NumMoms].P = J2;
					*NumMoms += 1;
				}
			}
		}
		
		
		
		if(MULT_USE_OPTION(
			Random Male Parent,
			RandoMale_F,
			m,
			rand-male,
			J1 J2,
			create J1 male parents randomly drawn from allele-freq-file J2,
			This will put J1 male parents on the male parent list and each one of them
			will be slated to get genotypes by randomly drawing from allele freqs as given
			in allele freq file J2.  Note that allele freq file 1 is the allele freq file
			designated in the first invocation of the -a or --allele-freq-file option.  Allele freq
			file 2 is the file designated in the second invocation of the -a or --allele-freq-file option
			and so on. This option may be issued 50 times., 
			50))
		{
			if(ARGS_EQ(2)) { int J1; int J2; int i;
				J1 = GET_INT;
				J2 = GET_INT;
				for(i=0;i<J1;i++)  {
					Dads[*NumDads].G = RANDO;
					Dads[*NumDads].P = J2;
					*NumDads += 1;
				}
			}
		}

		
		if ( MULT_USE_OPTION(
			Allele Freq File,
			LocFileF,
			a,
			allele-freq-file,
			F,
			S=pathname to file with locus information,
			S=pathname to file with locus information.  The format of this allele-frequency information
					file is simple.  It includes only integers or real numbers (or comments enclosed by a pair of & characters).
					The first integer is the number of loci in the file.  Then for each locus you give the
					number of alleles at that locus followed by the counts (possibly real-valued) of each allele at the locus.
					All these things should be separated by whitespace.  No punctuation!  If the number of
					allele counts listed does not match the number of alleles given for the locus then the program
					may fail without warning or give otherwise unexpected results.  Each time this command is executed
					the population counter is incremented by one.  So any commands previous to the current --freq-file option
					but after the last --freq-file option apply to the population with allele counts given in the current
					--freq-file option., 
					MAX_POPS) ) {
			if(  ARGS_EQ(1) ) {
				if(admixed_indF) {
					fprintf(stderr,"Error!  The -a/--allele-freq-file option was given after the --admixed-inds option was given.  This cannot be done.  Exiting\n");
					OPT_ERROR;
					exit(1);  /* I'm going to completely bail at this point */
				}
				GET_STR(locfilename);
				ReadLocFile(P,CurrentPop,locfilename);
				
				/* increment the population counter, and allocate to the next population */
				CurrentPop++;
				P[CurrentPop] = InitPopStruct(CurrentPop);
			}
		}
		
		
		if(MULT_USE_OPTION(
			Create Clutch,
			clutch_F,
			c,
			clutch,
			J1 J2 J3,
			Create a clutch of J1 individuals with mother J2 and father J3,
			This will cause the program to output a clutch of J1 individuals who are formed by
			Mendelian segregation of the genes in the female who is number J2 on the list of
			female parents and the male who is number J3 on the list of male parents.  May be used 
			1000 times., 
			MAX_CLUTCHES))
		{
			if(ARGS_EQ(3)) {
				Clutches[*NumClutches].NumKids = GET_INT;
				Clutches[*NumClutches].Mom = GET_INT;
				Clutches[*NumClutches].Dad = GET_INT;
				*NumClutches += 1;
			}
		}
		
		
		if(OPTION(
			Number of Replicates,
			NumReps_F,
			n,
			num-reps,
			J1,
			Number of reps of all the clutches to do,
			The default is 1.  But if you want to do the whole thing over---i.e. sample new genotypes
			for the random male and female parents and draw all the clutches over again---then you can 
			use this option.)) {
			if(ARGS_EQ(1)) {
				(*NumReps) = GET_INT;
			}
		}
		
		
		if(OPTION(
			Missing Data Rate,
			Missing_F,
			,
			miss-data-rate,
			R1 ... RL,
			per-locus rates of missing data in the offspring genotypes,
			R1 is the rate at which data will be simulated to be missing---completely at random---from
			individuals in the clutches at locus 1.  R2 is the same for locus 2.  Etc.  The default missing
			rate is 0.0.  Missing data at a genotype is expressed in the output as a -1 at each allele of the genotype. Currently the program
			does not check to make sure that the number of loci expressed on this line is consistent with the number of loci
			in other individuals in the input.  Be aware!)) { int j; int NumArgs;
			NumArgs = COUNT_ARGS;
			for(j=0;j<NumArgs;j++)  {
				MissingRates[j] = GET_DUB;
			}
		}
	
	
	
	if(OPTION(
			  Allele Incrementer,
			  alle_incr_F,
			  i,
			  allele-addition,
			  J,
			  Add this amount to the allele names of the non-admixed output of the program,
			  By default\054 nookie names alleles starting at 0 and going up from there.  This is a major 
			  hassle if you are feeding it into a program that thinks that missing data should be denoted by a 0.  Using this option\054
			  J will be added to every allele name (for the the non-admixed individual outputs).  It will also be added to the -1 that
			  typically denotes missing data from the miss-data-rate option\054 turning them into 0s.  This option will not affect the 
			  admixed-inds output as that is hard-wired to name alleles starting with 1\054 and the missing data does not affect the
			  admixed-inds output anyway.)) {
		if(ARGS_EQ(1)) {
			gAlleleNameIncrementer = GET_INT;
		}
	}
		
		if(MULT_USE_OPTION(
			Make Admixed Individuals,
			admixed_indF,
			,
			admixed-inds,
			J R1 ... RK,
			make J admixed indivs using structure model with q=R1...RK ,
			This command is separate from all the mom and dad and clutch stuff that comes
			before.  I am using it here because I needed this functionality and nookie was the 
			easiest of my programs to add it into.  This command must be given only after all uses of
			the -a/--allele-freq-file command have been given because it needs to know how many populations there
			are.  R1 is the expected proportion of genes in an individual coming from the first population specified
			with the -a option.  R2 is the expected proportion from the second population and so on.  If -a/--allele-freq-file
			was given K times then there should be K+1 arguments to the option---the first is the number of such individuals
			to make and the remaining K are the R1 ... RK values that define the q-vector of the individual.  May be used
			MAX_ADMIXED_INDS times.  Currrently that is set to 10000. ,
			MAX_ADMIXED_INDS))
		{
			if(ALREADY_HAS(LocFileF,"-a/--allele-freq-file")) {
				if(ARGS_EQ(CurrentPop + 1) ) { int i; double normo;
					AdmInds[CurAdmInd] = ECA_MALLOC(sizeof(admixed_ind_struct));
					AdmInds[CurAdmInd]->Num = GET_INT;  /* get the number */
					/* allocate to the Q vector */
					AdmInds[CurAdmInd]->Q = (double *)ECA_CALLOC(CurrentPop,sizeof(double));
					/* get the Q's */
					normo = 0.0;
					for(i=0;i<CurrentPop;i++)  {
						AdmInds[CurAdmInd]->Q[i] = GET_DUB;
						normo += AdmInds[CurAdmInd]->Q[i];
					}
					/* then normalize the Q's */
					for(i=0;i<CurrentPop;i++)  {
						AdmInds[CurAdmInd]->Q[i] /= normo;
					}
					/* increment the CurAdmInd counter */
					CurAdmInd++;
				}
			}
		}
		
		CLOSE_SUBSET;
		OPEN_SUBSET(Unimplemented Options, Unimplemented Options, Things I may or may not get around to doing...);
		
		if(MULT_USE_OPTION(
			Sample Size Factor,
			FactorF,
			,
			sample-size-factor,
			R,
			Set Dirichlet pars to allele counts times [R],
			Allele frequencies for each population are simulated from something
				akin to the predictive posterior for them given the baselines.  Using this
				option causes the Dirichlet parameters defining that predictive-posterior-like
				distribution to be multiplied by [R] for the current population.  Default is 1.0. NOT IMPLEMENTED,
			MAX_POPS)) {
			
			if(ARGS_EQ(1)) {
				P[CurrentPop]->Factor = GET_DUB;
				printf("SAMPLE_SIZE_FACTOR: Population number %d gets sample size factor %f\n",CurrentPop+1,P[CurrentPop]->Factor);
			}
		}
		
		
		if(MULT_USE_OPTION(
			Scale Prior By...,
			ScaleF,
			s,
			scale-prior,
			,
			Make the prior for allele freqs p/K where K is # of alleles,
			Issuing this option makes the prior for allele freqs p/K for the current 
				population where
				 p is the argument to the --prior or -p option and K is the number
				 of alleles at the locus.  The default is for the prior to be p for
				 each allele. NOT IMPLEMENTED,
			MAX_POPS)) {
			
			if(ARGS_EQ(0)) {
				P[CurrentPop]->ScalePrior = 1;
				printf("SCALE_PRIOR: Population number %d given the scale-prior option.\n",CurrentPop+1);
			}
		}
		
		
		CLOSE_SUBSET;
		
		
				
			
	
	END_OPT_LOOP
	
	*NumAdmInd = CurAdmInd;

	return(CurrentPop);
	
}


