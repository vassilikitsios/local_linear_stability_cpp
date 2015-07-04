
//============================================================================
class Eigen_Solver_t
{
private:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	char balanc, jobvl, jobvr, sense;
	integer n, lda, ldvl, ldb, ldvr, ilo, ihi, lwork;
	integer *iwork;
	doublereal abnrm, bbnrm;
	doublereal *lscale, *rscale, *rconde, *rcondv, *rwork;
	doublecomplex *a, *alpha, *vl, *b, *beta, *vr, *work;
	logical *bwork;
	integer info;
	int n_block;
	//-------------------------------------------------------------------------
	// helper functions
	//-------------------------------------------------------------------------	
	void malloc_clapack_variables();
	void convert_blitz_to_clapack(Array<complex<double>,2> aM, Array<complex<double>,2> bM);
	void calculate_eigenvalues();
	void convert_clapack_to_blitz(Array<complex<double>,2> aM, Array<complex<double>,2> bM);
	void free_clapack_variables();	
	//-------------------------------------------------------------------------
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	Eigen_Solution_t *eigen_solution;
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Eigen_Solver_t();
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------	
	void allocate_memory(int n, int n_block);
	void calculate_eigensolution(Array<complex<double>,2> aM, Array<complex<double>,2> bM);
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Eigen_Solver_t();
	//-------------------------------------------------------------------------
};

//============================================================================

