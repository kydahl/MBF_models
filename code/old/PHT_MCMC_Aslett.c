#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "PHT_MCMC_Aslett.h"
#include "utility.h"
#include "Simulate_AbsCTMC_gt_Aslett_DCS.h"
#include "Simulate_AbsCTMC_eq_Aslett_DCS.h"
#include "Simulate_AbsCTMC_gt_Bladt_MHRS.h"
#include "Simulate_AbsCTMC_eq_Bladt_MHRS.h"
#include "Simulate_AbsCTMC_eq_Aslett_ECS.h"
#include "Simulate_AbsCTMC_gt_Hobolth_DCS.h"
#include "Simulate_AbsCTMC_eq_AslettHobolth_DCS.h"

//#define DEBUG
/*

===> TODO <===
-> improve matrix exponential stuff
-> fix so that pi is part of sampling scheme
-> censoring (done?)


===> WHAT'S NEW RECENTLY? <===

-> Bladt no longer allocates any dynamic memory
-> Aslett no longer allocates dynamic memory
-> ARMS no longer allocates any dynamic memory
-> fixed callocing in LAPACK calls (computes optimal at start and uses that on global scope)
-> partial ARMS/Exp sampling
-> Gibbs is done in C fully
-> start at the mode rather than a random prior draw

*/

//#define DEBUG

typedef struct p_into_double_ {
	double *data;
	double dataFixed;
	struct p_into_double_ *next;
} p_into_double;
typedef struct p_into_int_ {
	int *data;
	struct p_into_int_ *next;
} p_into_int;

p_into_int *listAddI(p_into_int *p, int *elt) {
	p_into_int *newNode = (p_into_int *) R_alloc(1, sizeof(p_into_int));
	newNode->next = p;
	p = newNode;
	newNode->data = elt;
	return p;
}
p_into_double *listAddDD(p_into_double *p, double *elt, double eltFixed) {
	p_into_double *newNode = (p_into_double *) R_alloc(1, sizeof(p_into_double));
	newNode->next = p;
	p = newNode;
	newNode->data = elt;
	newNode->dataFixed = eltFixed;
	return p;
}
p_into_double *listAddD(p_into_double *p, double *elt) {
	return listAddDD(p, elt, 1.0);
}

// Bit masks for sampling method
static int METHOD_MHRS = 0x1; // MHRS (Bladt et al.)
static int METHOD_ECS = 0x2; // ECS (Aslett & Wilson)
static int METHOD_DCS = 0x4; // DCS (Aslett changes to Hobolth)
//// Performs the Gibbs sampling (calling MH within)
// it (input)
//     number of Gibbs iterates to perform
// mhit (input)
//     number of MH within Gibbs iterates to perform
// method (input)
//     Binary mask: 0x1 = Forward Uncond, 0x2 = Reverse Uncond, 0x4 = Forward Cond, 0x8 = Reverse Cond
// n (input)
//     dimension (of non-absorbing part)
// m (input)
//     number of parameters in the inference
// nu (input)
//     m dimensional vector of nu prior parameters
// zeta (input)
//     m dimensional vector of zeta prior parameters
// T (input)
//     form of the matrix to infer, coded by integer (0=fix at zero, 1+=variable index in res)
// C (input)
//     matrix of constant multiples of the parameters in T (should just be 1.0 in every element if unsure)
// y (input)
//     l dimensional array of absorption times to be conditioned on
// l (input)
//     number of observations
// censored (input)
//     l dimensional array of 1/0 for is/not a censored observation
// start (input)
//     m dimensional array of starting points of chain for each variable ... if first entry is -1, draw our own starting values
// silent (input)
//     0 = verbose progress output, 1 = quiet progress output
// res (output)
//     it x m matrix to hold the samples for each variable
//
void LJMA_Gibbs(int *it, int *mhit, int *method, int *n, int *m, double *nu, double *zeta, int *T, double *C, double *y, int *l, int *censored, double *start, int *silent, double *res) {
	LJMA_GetRNGstate();

	//static const char bss[] = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
	//char b[250];
	//int bs=0;

	// Setup general storage space
	double rsum, *pi, *TT, *S, *s, *z, *zsum, *P, *Pfull;
	pi = (double *) R_alloc(*n, sizeof(double));
	TT = (double *) R_alloc((*n+1)*(*n+1), sizeof(double));
	S = (double *) R_alloc(*n * *n, sizeof(double));
	s = (double *) R_alloc(*n, sizeof(double));
	z = (double *) R_alloc(*n, sizeof(double));
	zsum = (double *) R_alloc(*m, sizeof(double));
	P = (double *) R_alloc(*n * *n, sizeof(double));
	Pfull = (double *) R_alloc(*n * (*n+1), sizeof(double));
	//stationary = (double *) R_alloc(*n, sizeof(double));
	int *N, *B, *Nsum;
	N = (int *) R_alloc(*n * *n, sizeof(int));
	B = (int *) R_alloc(*n, sizeof(int));
	Nsum = (int *) R_alloc(*m, sizeof(int));

	// Setup high-speed linked list pointer space
	p_into_int **NbyVar;
	NbyVar = (p_into_int **) R_alloc(*m, sizeof(p_into_int*)); for(int i=0; i<*m; i++) NbyVar[i] = NULL;
	p_into_double **zbyVar, **TTbyVar, **SbyVar, **sbyVar, **TTDiag;
	zbyVar = (p_into_double **) R_alloc(*m, sizeof(p_into_double*)); for(int i=0; i<*m; i++) zbyVar[i] = NULL;
	TTbyVar = (p_into_double **) R_alloc(*m, sizeof(p_into_double*)); for(int i=0; i<*m; i++) TTbyVar[i] = NULL;
	SbyVar = (p_into_double **) R_alloc(*m, sizeof(p_into_double*)); for(int i=0; i<*m; i++) SbyVar[i] = NULL;
	sbyVar = (p_into_double **) R_alloc(*m, sizeof(p_into_double*)); for(int i=0; i<*m; i++) sbyVar[i] = NULL;
	TTDiag = (p_into_double **) R_alloc(*n+1, sizeof(p_into_double*)); for(int i=0; i<*n+1; i++) TTDiag[i] = NULL;

	// Extra sampling specific storage;
	double *Q=NULL, *Qinv=NULL, *evals=NULL, *e=NULL, *Qinv_1=NULL, *Qinv_s=NULL, *b=NULL, *Qinv_b=NULL;
	double oneD = 1.0, zeroD = 0.0;
	int oneI = 1;
	char transN = 'N';
	if(((*method & METHOD_ECS) | (*method & METHOD_DCS)) > 0) {
		Q = (double *) R_alloc(*n * *n, sizeof(double));
		Qinv = (double *) R_alloc(*n * *n, sizeof(double));
		evals = (double *) R_alloc(*n, sizeof(double));
	}
	if((*method & METHOD_ECS) > 0) {
		Qinv_1 = (double *) R_alloc(*n, sizeof(double));
		e = (double *) R_alloc(*n, sizeof(double));
		for(int i=0; i<*n; i++) e[i] = 1.0;
		Qinv_s = (double *) R_alloc(*n, sizeof(double));
	}
	if((*method & METHOD_DCS) > 0) {
		Qinv_b = (double *) R_alloc(*n, sizeof(double)); // Do we need these really still?
		b = (double *) R_alloc(*n, sizeof(double));
	}

	// Setup workspace storage
	double *workD;
	int *workI, workIsize = 0, workDsize = 0;
	// ---> FIX THIS BY RECOUNTING MEMORY USAGE BY SAMPLING METHOD <--
	if((*method & METHOD_MHRS) > 0) {
		workDsize = (100000 < workDsize) ? workDsize : 100000;
		workIsize = (100000 < workIsize) ? workIsize : 100000;
	}
	if((*method & METHOD_ECS) > 0) {
		workDsize = (100000 < workDsize) ? workDsize : 100000;
		workIsize = (100000 < workIsize) ? workIsize : 100000;
	}
	if((*method & METHOD_DCS) > 0) {
		workDsize = (100000 < workDsize) ? workDsize : 100000;
		workIsize = (100000 < workIsize) ? workIsize : 100000;
	}
	workD = (double *) R_alloc(workDsize, sizeof(double));
	workI = (int *) R_alloc(workIsize, sizeof(int));

	if(((*method & METHOD_ECS) | (*method & METHOD_DCS)) > 0) {
		// LAPACK workspace ... make NULL calls to figure out optimal workspace size for speed
		char balanc = 'B', jobvl = 'V', jobvr = 'V', sense = 'B'; int lwork = -1, info, ilo, ihi; double work, abnrm;
		F77_CALL(dgeevx)(&balanc, &jobvl, &jobvr, &sense, n, Q, n, evals, workD, workD, n, Q, n, &ilo, &ihi, workD, &abnrm, workD, workD, &work, &lwork, NULL, &info FCONE FCONE FCONE FCONE);
		LJMA_LAPACK_lwork = (int) work;
		F77_CALL(dgetri)(n, NULL, n, NULL, &work, &lwork, &info);
		if((int) work > LJMA_LAPACK_lwork) LJMA_LAPACK_lwork = (int) work;
		LJMA_LAPACK_work = (double *) R_alloc(LJMA_LAPACK_lwork, sizeof(double));
	}

	// Step 0 - draw from priors
	Rprintf("Setting up Gibbs run ...\n"); LJMA_GUI();

	// ******************** FIX ME ********************
	*pi = 1.0;
	for(int i=1; i<*n; i++) pi[i] = 0.0;
	// ******************** FIX ME ********************

	if(*start < 0) {
		for(int i=0; i<*m; i++) {
			if(nu[i] > 1) { // if k>1 can choose modal value of prior
				res[0 + i * *it] = (nu[i]-1.0) / zeta[i];
			} else { // else just do a random draw
				res[0 + i * *it] = rgamma(nu[i], 1.0/zeta[i]);
			}
		}
	} else {
		for(int i=0; i<*m; i++) {
			res[0 + i * *it] = start[i];
		}
	}
	//for(int i=0; i<*m; i++) { Rprintf("%lf ", res[0 + i * *it]); } Rprintf("\n\n"); LJMA_GUI(); // Print starting values
	rsum = 0.0;
	for(int i=0; i<*n+1; i++) {
		for(int j=0; j<*n+1; j++) {
			if(T[i + j * (*n+1)] == 0) {
				TT[i + j * (*n+1)] = 0.0;
			} else {
				rsum -= TT[i + j * (*n+1)] = res[0 + (T[i + j * (*n+1)]-1) * *it] * C[i + j * (*n+1)];
				if(j == *n) {
					NbyVar[T[i + j * (*n+1)]-1] = listAddI(NbyVar[T[i + j * (*n+1)]-1], &N[i + i * *n]);
					sbyVar[T[i + j * (*n+1)]-1] = listAddDD(sbyVar[T[i + j * (*n+1)]-1], &s[i], C[i + j * (*n+1)]);
				} else {
					NbyVar[T[i + j * (*n+1)]-1] = listAddI(NbyVar[T[i + j * (*n+1)]-1], &N[i + j * *n]);
					SbyVar[T[i + j * (*n+1)]-1] = listAddDD(SbyVar[T[i + j * (*n+1)]-1], &S[i + j * *n], C[i + j * (*n+1)]);
				}
				zbyVar[T[i + j * (*n+1)]-1] = listAddDD(zbyVar[T[i + j * (*n+1)]-1], &z[i], C[i + j * (*n+1)]); // For computing sum of z's
				TTbyVar[T[i + j * (*n+1)]-1] = listAddDD(TTbyVar[T[i + j * (*n+1)]-1], &TT[i + j * (*n+1)], C[i + j * (*n+1)]); // For filling TT matrix
//				TTbyVar[T[i + j * (*n+1)]-1] = listAddD(TTbyVar[T[i + j * (*n+1)]-1], &TT[i + j * (*n+1)]); // For filling TT matrix
				TTDiag[i] = listAddD(TTDiag[i], &TT[i + j * (*n+1)]);
			}
		}
		TT[i + i * (*n+1)] = rsum;
		rsum = 0.0;
	}
		/*Rprintf("\n");
		for(int i=0; i<*n+1; i++) {
			for(int j=0; j<*n+1; j++) {
				Rprintf(" %.3lf ", TT[i + j * (*n+1)]);
			}
			Rprintf("\n");
		}*/
	for(int i=0; i<*n; i++) { // get S from TT
		for(int j=0; j<*n; j++) {
			S[i + j * *n] = TT[i + j * (*n+1)];
		}
	}
	for(int i=0; i<*n; i++) { // get s from TT
		s[i] = TT[i + (*n) * (*n+1)];
	}
	if((*method & METHOD_DCS) > 0) {
		for(int i=0; i<*n; i++) {
			if(s[i] > 0.0) b[i] = 1.0; else b[i] = 0.0;
		}
	}
		/*Rprintf("\n");
		for(int i=0; i<*n; i++) {
			for(int j=0; j<*n; j++) {
				Rprintf(" %.3lf ", S[i + j * *n]);
			}
			Rprintf("\n");
		}
		Rprintf("\n");
		for(int i=0; i<*n; i++) {
			Rprintf(" %.3lf ", s[i]);
		}*/

	Rprintf("Starting phase-type MCMC sampler ...\n\nBegining processing ..."); LJMA_GUI();
	if(*silent) {
		Rprintf(" silent processing selected, there will be no further feedback until MCMC run complete"); LJMA_GUI();
	}
	for(int iter=1; iter<*it; iter++) {
		int i, j;
		//bs = -bs - 1 + sprintf(b, "%.*sProcessing iteration %d of %d (%.1lf\%%%%)", bs, bss, iter+1, *it, (100.0*(iter+1)) / *it);
		//Rprintf(b);
		if(!*silent) {
			Rprintf("\rProcessing iteration %d of %d (%.1lf%%)\r", iter+1, *it, (100.0*(iter+1)) / *it); LJMA_GUI();
		}
		Rprintf("I'm on line 275"); LJMA_GUI();
		// Step 1 - simulate the sufficient statistics of the unobserved processes
		// First: compute the embedded chain transition probabilities
		// --> FIX: this isn't needed for DCS <--
		double rsum, rsumfull;
		for(i=0; i<*n; i++) {
		  Rprintf("I'm on line 281"); LJMA_GUI();
			rsumfull = 0.0;
			for(j=0; j<*n; j++) {
				rsumfull += Pfull[i + j * *n] = P[i + j * *n] = -S[i + j * *n] / S[i + i * *n];
			}
			rsum = rsumfull - P[i + i * *n];
			rsumfull += Pfull[i + *n * *n] = -s[i] / S[i + i * *n];
			rsumfull -= Pfull[i + i * *n];
			Pfull[i + i * *n] = P[i + i * *n] = 0.0;
			for(j=0; j<*n; j++) {
				P[i + j * *n] = P[i + j * *n] / rsum;
				Pfull[i + j * *n] = Pfull[i + j * *n] / rsumfull;
				//Rprintf("%lf ", Pfull[i + j * *n]);
			}
			Pfull[i + *n * *n] = Pfull[i + *n * *n] / rsumfull;
			//Rprintf("%lf ", Pfull[i + *n * *n]);
			//Rprintf("\n");
		}

		// Second: check if we should allow reverse simulation (if requested) [NB: starting dist varys by time, so can't precompute!]
		// Get stationary dist
		/*
		allowReverse = 0;
		if((*method & METHOD_REV_UNC) | (*method & METHOD_REV_CON) > 0) {
			allowReverse = 1;
			LJMA_stationary(S, n, stationary, workD, workI);
			// Check the thing is even reversible!
			for(int i=0; i<*n; i++) {
				for(int j=0; j<*n; j++) {
					if(i!=j && fabs(stationary[i]*S[i + j * *n] - stationary[j]*S[j + i * *n]) > 1e-9) {
						allowReverse = 0;
						goto out; // Recommended by K & R §3.8
					}
				}
			}
		}
		out:
		*/

		// Third: eigen-decompose if necessary
		if(((*method & METHOD_ECS) | (*method & METHOD_DCS)) > 0) {
			LJMA_eigen(n, S, evals, Q, Qinv, workD, workI);
		}

		// Then do the simulation
		if((*method & METHOD_MHRS) > 0) {
			LJMA_MHsample_Bladt(y, censored, l, pi, S, s, Pfull, n, mhit, z, B, N, workD, workI);
		} else if((*method & METHOD_DCS) > 0) {
			F77_CALL(dgemv)(&transN, n, n, &oneD, Qinv, n, b, &oneI, &zeroD, Qinv_b, &oneI FCONE);
			LJMA_MHsample_Hobolth2(y, censored, l, pi, S, s, Q, evals, Qinv_b, b, Qinv, n, mhit, z, B, N, workD, workI);
		} else if((*method & METHOD_ECS) > 0) {
			F77_CALL(dgemv)(&transN, n, n, &oneD, Qinv, n, s, &oneI, &zeroD, Qinv_s, &oneI FCONE);
			F77_CALL(dgemv)(&transN, n, n, &oneD, Qinv, n, e, &oneI, &zeroD, Qinv_1, &oneI FCONE);
			LJMA_MHsample_Aslett2(y, censored, l, pi, S, s, Q, evals, Qinv_s, Qinv_1, P, Pfull, n, z, B, N, workD, workI);
		} else {
			Rprintf("CRITICAL ERROR: Unknown sampling method (code = %d)\n\n", *method);
			continue;
		}

		// Fourth: compile the stuff together
		for(i=0; i<*m; i++) {
			Nsum[i] = 0;
			zsum[i] = 0.0;
			p_into_int *pI;
			pI = NbyVar[i];
			while(pI!=NULL) {
				Nsum[i] += *(pI->data);
				pI = pI->next;
			}
			p_into_double *pD;
			pD = zbyVar[i];
			while(pD!=NULL) {
				zsum[i] += *(pD->data) / pD->dataFixed;
				pD = pD->next;
			}
		}
//		for(i=0; i<*n; i++) Rprintf("%lf ",z[i]);
//		Rprintf("\n");
//		for(i=0; i<*m; i++) Rprintf("%lf ",zsum[i]);
//		Rprintf("\n");
//		for(i=0; i<*m; i++) Rprintf("%d ",Nsum[i]);
//		Rprintf("\n");

		// Fifth: draw from posteriors
		double tmp;
		for(i=0; i<*m; i++) {
			tmp = res[iter + i * *it] = rgamma(nu[i]+Nsum[i], 1.0/(zeta[i]+zsum[i]));

			p_into_double *pD;
			pD = TTbyVar[i];
			while(pD!=NULL) {
				*(pD->data) = tmp * pD->dataFixed;
				pD = pD->next;
			}
			pD = SbyVar[i];
			while(pD!=NULL) {
				*(pD->data) = tmp * pD->dataFixed;
				pD = pD->next;
			}
			pD = sbyVar[i];
			while(pD!=NULL) {
				*(pD->data) = tmp * pD->dataFixed;
				pD = pD->next;
			}
		}
		for(i=0; i<*n; i++) {
			tmp = 0.0;

			p_into_double *pD;
			pD = TTDiag[i];
			while(pD!=NULL) {
				tmp -= *(pD->data);
				pD = pD->next;
			}

			TT[i + i * (*n+1)] = tmp;
			S[i + i * *n] = tmp;
		}
//		for(i=0; i<(*n+1); i++) {
//			for(j=0; j<(*n+1); j++) {
//				Rprintf("%lf ", TT[i + j * (*n+1)]);
//			}
//			Rprintf("\n");
//		}
//		Rprintf("\n");
	}

	Rprintf("\n\nCompleted MCMC run, returning results ...\n"); LJMA_GUI();

	LJMA_PutRNGstate();
}



