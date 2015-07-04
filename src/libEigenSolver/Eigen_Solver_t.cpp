
#include <iostream>
#include <math.h>

using namespace std;

#include <stdio.h>
#include <stdarg.h>

#include <blitz/array.h>
using namespace blitz;

#include "Utils.h"
#include "lapack_mallocs.h"
#include "Eigen_Solution_t.h"
#include "Eigen_Solver_t.h"

/* =====================================================================   
 -- LAPACK driver routine (version 3.0) --   
 Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
 Courant Institute, Argonne National Lab, and Rice University   
 June 30, 1999   
 
 
 Purpose   
 =======   
 
 ZGGEVX computes for a pair of N-by-N complex nonsymmetric matrices   
 (A,B) the generalized eigenvalues, and optionally, the left and/or   
 right generalized eigenvectors.   
 
 Optionally, it also computes a balancing transformation to improve   
 the conditioning of the eigenvalues and eigenvectors (ILO, IHI,   
 LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for   
 the eigenvalues (RCONDE), and reciprocal condition numbers for the   
 right eigenvectors (RCONDV).   
 
 A generalized eigenvalue for a pair of matrices (A,B) is a scalar   
 lambda or a ratio alpha/beta = lambda, such that A - lambda*B is   
 singular. It is usually represented as the pair (alpha,beta), as   
 there is a reasonable interpretation for beta=0, and even for both   
 being zero.   
 
 The right eigenvector v(j) corresponding to the eigenvalue lambda(j)   
 of (A,B) satisfies   
 A * v(j) = lambda(j) * B * v(j) .   
 The left eigenvector u(j) corresponding to the eigenvalue lambda(j)   
 of (A,B) satisfies   
 u(j)**H * A  = lambda(j) * u(j)**H * B.   
 where u(j)**H is the conjugate-transpose of u(j).   
 
 
 Arguments   
 =========   
 
 BALANC  (input) CHARACTER*1   
 Specifies the balance option to be performed:   
 = 'N':  do not diagonally scale or permute;   
 = 'P':  permute only;   
 = 'S':  scale only;   
 = 'B':  both permute and scale.   
 Computed reciprocal condition numbers will be for the   
 matrices after permuting and/or balancing. Permuting does   
 not change condition numbers (in exact arithmetic), but   
 balancing does.   
 
 JOBVL   (input) CHARACTER*1   
 = 'N':  do not compute the left generalized eigenvectors;   
 = 'V':  compute the left generalized eigenvectors.   
 
 JOBVR   (input) CHARACTER*1   
 = 'N':  do not compute the right generalized eigenvectors;   
 = 'V':  compute the right generalized eigenvectors.   
 
 SENSE   (input) CHARACTER*1   
 Determines which reciprocal condition numbers are computed.   
 = 'N': none are computed;   
 = 'E': computed for eigenvalues only;   
 = 'V': computed for eigenvectors only;   
 = 'B': computed for eigenvalues and eigenvectors.   
 
 N       (input) INTEGER   
 The order of the matrices A, B, VL, and VR.  N >= 0.   
 
 A       (input/output) COMPLEX*16 array, dimension (LDA, N)   
 On entry, the matrix A in the pair (A,B).   
 On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'   
 or both, then A contains the first part of the complex Schur   
 form of the "balanced" versions of the input A and B.   
 
 LDA     (input) INTEGER   
 The leading dimension of A.  LDA >= max(1,N).   
 
 B       (input/output) COMPLEX*16 array, dimension (LDB, N)   
 On entry, the matrix B in the pair (A,B).   
 On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'   
 or both, then B contains the second part of the complex   
 Schur form of the "balanced" versions of the input A and B.   
 
 LDB     (input) INTEGER   
 The leading dimension of B.  LDB >= max(1,N).   
 
 ALPHA   (output) COMPLEX*16 array, dimension (N)   
 BETA    (output) COMPLEX*16 array, dimension (N)   
 On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the generalized   
 eigenvalues.   
 
 Note: the quotient ALPHA(j)/BETA(j) ) may easily over- or   
 underflow, and BETA(j) may even be zero.  Thus, the user   
 should avoid naively computing the ratio ALPHA/BETA.   
 However, ALPHA will be always less than and usually   
 comparable with norm(A) in magnitude, and BETA always less   
 than and usually comparable with norm(B).   
 
 VL      (output) COMPLEX*16 array, dimension (LDVL,N)   
 If JOBVL = 'V', the left generalized eigenvectors u(j) are   
 stored one after another in the columns of VL, in the same   
 order as their eigenvalues.   
 Each eigenvector will be scaled so the largest component   
 will have abs(real part) + abs(imag. part) = 1.   
 Not referenced if JOBVL = 'N'.   
 
 LDVL    (input) INTEGER   
 The leading dimension of the matrix VL. LDVL >= 1, and   
 if JOBVL = 'V', LDVL >= N.   
 
 VR      (output) COMPLEX*16 array, dimension (LDVR,N)   
 If JOBVR = 'V', the right generalized eigenvectors v(j) are   
 stored one after another in the columns of VR, in the same   
 order as their eigenvalues.   
 Each eigenvector will be scaled so the largest component   
 will have abs(real part) + abs(imag. part) = 1.   
 Not referenced if JOBVR = 'N'.   
 
 LDVR    (input) INTEGER   
 The leading dimension of the matrix VR. LDVR >= 1, and   
 if JOBVR = 'V', LDVR >= N.   
 
 ILO,IHI (output) INTEGER   
 ILO and IHI are integer values such that on exit   
 A(i,j) = 0 and B(i,j) = 0 if i > j and   
 j = 1,...,ILO-1 or i = IHI+1,...,N.   
 If BALANC = 'N' or 'S', ILO = 1 and IHI = N.   
 
 LSCALE  (output) DOUBLE PRECISION array, dimension (N)   
 Details of the permutations and scaling factors applied   
 to the left side of A and B.  If PL(j) is the index of the   
 row interchanged with row j, and DL(j) is the scaling   
 factor applied to row j, then   
 LSCALE(j) = PL(j)  for j = 1,...,ILO-1   
 = DL(j)  for j = ILO,...,IHI   
 = PL(j)  for j = IHI+1,...,N.   
 The order in which the interchanges are made is N to IHI+1,   
 then 1 to ILO-1.   
 
 RSCALE  (output) DOUBLE PRECISION array, dimension (N)   
 Details of the permutations and scaling factors applied   
 to the right side of A and B.  If PR(j) is the index of the   
 column interchanged with column j, and DR(j) is the scaling   
 factor applied to column j, then   
 RSCALE(j) = PR(j)  for j = 1,...,ILO-1   
 = DR(j)  for j = ILO,...,IHI   
 = PR(j)  for j = IHI+1,...,N   
 The order in which the interchanges are made is N to IHI+1,   
 then 1 to ILO-1.   
 
 ABNRM   (output) DOUBLE PRECISION   
 The one-norm of the balanced matrix A.   
 
 BBNRM   (output) DOUBLE PRECISION   
 The one-norm of the balanced matrix B.   
 
 RCONDE  (output) DOUBLE PRECISION array, dimension (N)   
 If SENSE = 'E' or 'B', the reciprocal condition numbers of   
 the selected eigenvalues, stored in consecutive elements of   
 the array.   
 If SENSE = 'V', RCONDE is not referenced.   
 
 RCONDV  (output) DOUBLE PRECISION array, dimension (N)   
 If JOB = 'V' or 'B', the estimated reciprocal condition   
 numbers of the selected eigenvectors, stored in consecutive   
 elements of the array. If the eigenvalues cannot be reordered   
 to compute RCONDV(j), RCONDV(j) is set to 0; this can only   
 occur when the true value would be very small anyway.   
 If SENSE = 'E', RCONDV is not referenced.   
 Not referenced if JOB = 'E'.   
 
 WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)   
 On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   
 
 LWORK   (input) INTEGER   
 The dimension of the array WORK. LWORK >= max(1,2*N).   
 If SENSE = 'N' or 'E', LWORK >= 2*N.   
 If SENSE = 'V' or 'B', LWORK >= 2*N*N+2*N.   
 
 If LWORK = -1, then a workspace query is assumed; the routine   
 only calculates the optimal size of the WORK array, returns   
 this value as the first entry of the WORK array, and no error   
 message related to LWORK is issued by XERBLA.   
 
 RWORK   (workspace) DOUBLE PRECISION array, dimension (6*N)   
 Real workspace.   
 
 IWORK   (workspace) INTEGER array, dimension (N+2)   
 If SENSE = 'E', IWORK is not referenced.   
 
 BWORK   (workspace) LOGICAL array, dimension (N)   
 If SENSE = 'N', BWORK is not referenced.   
 
 INFO    (output) INTEGER   
 = 0:  successful exit   
 < 0:  if INFO = -i, the i-th argument had an illegal value.   
 = 1,...,N:   
 The QZ iteration failed.  No eigenvectors have been   
 calculated, but ALPHA(j) and BETA(j) should be correct   
 for j=INFO+1,...,N.   
 > N:  =N+1: other than QZ iteration failed in ZHGEQZ.   
 =N+2: error return from ZTGEVC.   
 
 Further Details   
 ===============   
 
 Balancing a matrix pair (A,B) includes, first, permuting rows and   
 columns to isolate eigenvalues, second, applying diagonal similarity   
 transformation to the rows and columns to make the rows and columns   
 as close in norm as possible. The computed reciprocal condition   
 numbers correspond to the balanced matrix. Permuting rows and columns   
 will not change the condition numbers (in exact arithmetic) but   
 diagonal scaling will.  For further explanation of balancing, see   
 section 4.11.1.2 of LAPACK Users' Guide.   
 
 An approximate error bound on the chordal distance between the i-th   
 computed generalized eigenvalue w and the corresponding exact   
 eigenvalue lambda is   
 
 chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)   
 
 An approximate error bound for the angle between the i-th computed   
 eigenvector VL(i) or VR(i) is given by   
 
 EPS * norm(ABNRM, BBNRM) / DIF(i).   
 
 For further explanation of the reciprocal condition numbers RCONDE   
 and RCONDV, see section 4.11 of LAPACK User's Guide.   
 
 =====================================================================   */

extern "C" int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
					   sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *b,
					   integer *ldb, doublecomplex *alpha, doublecomplex *beta, 
					   doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
					   integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale,
					   doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
					   rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
					   integer *iwork, logical *bwork, integer *info);

//----------------------------------------------------------------------------
Eigen_Solver_t::Eigen_Solver_t()
{
}

//----------------------------------------------------------------------------
Eigen_Solver_t::~Eigen_Solver_t()
{
	free_Eigen_Solution_t(eigen_solution,n);
	free_clapack_variables();
}

//----------------------------------------------------------------------------
void Eigen_Solver_t::allocate_memory(int n, int n_block)
{
	// Set flags for clapack routine
	balanc = 'N';
	jobvl = 'V';
	jobvr = 'V';
	sense = 'E';										// Something wrong with the eigenvector reciprocal condition
	
	// Assign variables
	this->n = n;
	this->n_block = n_block;
	lda = n;	ldvl = n;
	ldb = n;	ldvr = n;
	ilo = 0;	ihi = 0;
	abnrm = 0;	bbnrm = 0;								// The one norm of the matrics 'a' and 'b'

//	if ((&sense=="N") || (&sense=="E")) {
	if ( (strcmp(&sense,"N")==0) || (strcmp(&sense,"E")==0) ) {
		lwork = 2*n;
	} else {
		lwork = 2*n*n+2*n;
	}
	info = -100;
	
	// Allcoate Memory
	malloc_clapack_variables();
	eigen_solution = malloc_Eigen_Solution_t(n);
}

//----------------------------------------------------------------------------
void Eigen_Solver_t::malloc_clapack_variables()
{
	a =  malloc_1d_array_doublecomplex(lda*n);
	alpha = malloc_1d_array_doublecomplex(n);
	vl = malloc_1d_array_doublecomplex(ldvl*n);	// For left handed eigensolution
	b = malloc_1d_array_doublecomplex(ldb*n);
	beta = malloc_1d_array_doublecomplex(n);
	vr = malloc_1d_array_doublecomplex(ldvr*n);
	lscale = malloc_1d_array_doublereal(n);
	rscale = malloc_1d_array_doublereal(n);
	rconde = malloc_1d_array_doublereal(n);		// For the system Ax=b, the condition number is the effect that a small change in b has on x
	rcondv = malloc_1d_array_doublereal(n);
	work = malloc_1d_array_doublecomplex(lwork);
	rwork = malloc_1d_array_doublereal(6*n);
	iwork = malloc_1d_array_integer(n+2);
	bwork = malloc_1d_array_logical(n);
}

//----------------------------------------------------------------------------
void Eigen_Solver_t::free_clapack_variables()
{
	free_1d_array_doublecomplex(a);
	free_1d_array_doublecomplex(alpha);
	free_1d_array_doublecomplex(vl);
	free_1d_array_doublecomplex(b);
	free_1d_array_doublecomplex(beta);
	free_1d_array_doublecomplex(vr);
	free_1d_array_doublereal(lscale);
	free_1d_array_doublereal(rscale);
	free_1d_array_doublereal(rconde);
	free_1d_array_doublereal(rcondv);
	free_1d_array_doublecomplex(work);
	free_1d_array_doublereal(rwork);
	free_1d_array_integer(iwork);
	free_1d_array_logical(bwork);
}

//----------------------------------------------------------------------------
void Eigen_Solver_t::convert_blitz_to_clapack(Array<complex<double>,2> aM, Array<complex<double>,2> bM)
{
	int pos = 0;
	for (int ii=0; ii<n; ii++) {
		for (int jj=0; jj<n; jj++) {
			a[pos].r = real(aM(jj,ii));					// Swap c++ indexing fortran indexing
			a[pos].i = imag(aM(jj,ii));
			b[pos].r = real(bM(jj,ii)); 
			b[pos].i = imag(bM(jj,ii));
			pos++;
		}
	}
}

//----------------------------------------------------------------------------
void Eigen_Solver_t::calculate_eigensolution(Array<complex<double>,2> aM, Array<complex<double>,2> bM) 
{
	convert_blitz_to_clapack(aM, bM);
	zggevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alpha, beta, vl,
			&ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, rconde, rcondv, work,
			&lwork, rwork, iwork, bwork, &info);
	convert_clapack_to_blitz(aM, bM);
	calculate_eigenvalues();
}

//----------------------------------------------------------------------------
void Eigen_Solver_t::calculate_eigenvalues()
{
	double beta_real, beta_imag;
	for (int ii=0; ii<n; ii++) {
		eigen_solution[ii].eigenvalue = eigen_solution[ii].alpha * conj(eigen_solution[ii].beta) / pow(abs(eigen_solution[ii].beta),2.0);
		beta_real = real(eigen_solution[ii].beta);
		beta_imag = imag(eigen_solution[ii].beta);
		if ( (ABS(beta_real) < 1.0e-9) && (ABS(beta_imag) < 1.0e-9) ) {
			real(eigen_solution[ii].eigenvalue) = 1.0e15;
			imag(eigen_solution[ii].eigenvalue) = -1.0e15;
			eigen_solution[ii].rel_error_bound = 1.0e15;
		}
	}
}

//----------------------------------------------------------------------------
void Eigen_Solver_t::convert_clapack_to_blitz(Array<complex<double>,2> aM, Array<complex<double>,2> bM)
{
	double epsilon = 1.0e-15;		// Indication of numerical round off error
	int pos = 0;
	double u_mode_sum, w_mode_sum;
	
	aM = 0.0;	bM = 0.0;
	
	for (int ii=0; ii<n; ii++) {
		eigen_solution[ii].unsorted_mode_num = ii;
		real(eigen_solution[ii].alpha) = alpha[ii].r;
		imag(eigen_solution[ii].alpha) = alpha[ii].i;		
		real(eigen_solution[ii].beta) = beta[ii].r;
		imag(eigen_solution[ii].beta) = beta[ii].i;
		eigen_solution[ii].rel_error_bound = sqrt(pow(abnrm,2.0)+pow(bbnrm,2.0))/rconde[ii]*epsilon;
		eigen_solution[ii].cond = 1.0/rconde[ii];	// Note rconde - is the reciprical condition number
		
		for (int jj=0; jj<n; jj++) {
			real(aM(jj,ii)) = vr[pos].r;
			imag(aM(jj,ii)) = vr[pos].i;
			real(bM(jj,ii)) = vl[pos].r;
			imag(bM(jj,ii)) = vl[pos].i;
			pos++;
		}
		
		// Classify if it is a 'w' mode and if it is symmetric
		u_mode_sum = 0.0;	w_mode_sum = 0.0;
		for (int jj=0; jj<n_block; jj++) {
			u_mode_sum = u_mode_sum + real(aM(jj,ii)) * real(aM(jj,ii));
			w_mode_sum = w_mode_sum + real(aM(jj+2*n_block,ii)) * real(aM(jj+2*n_block,ii));
		}
		
		if ( u_mode_sum > w_mode_sum ) {
			eigen_solution[ii].w_mode = 0.0;
			if ( real(aM(1,ii)) * real(aM(n_block-2,ii)) > 0.0 ) {
				eigen_solution[ii].symmetric = 1.0;
			} else{ 
				eigen_solution[ii].symmetric = 0.0;
			}
		} else {	
			eigen_solution[ii].w_mode = 1.0;
			if ( real(aM(2*n_block+1,ii)) * real(aM(3*n_block-2,ii)) > 0.0 ) {
				eigen_solution[ii].symmetric = 1.0;
			} else{ 
				eigen_solution[ii].symmetric = 0.0;
			}
		}
		
	}
	
}

//----------------------------------------------------------------------------


