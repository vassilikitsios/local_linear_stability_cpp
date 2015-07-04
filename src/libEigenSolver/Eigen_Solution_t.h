
//============================================================================
class Eigen_Solution_t
	{
	public:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		complex<double> alpha, beta, eigenvalue;
		double rel_error_bound;
		double cond;
		int discrete_continuous;
		int w_mode;
		int symmetric;
		int unsorted_mode_num;
		//-------------------------------------------------------------------------	
	};

//============================================================================
Eigen_Solution_t* malloc_Eigen_Solution_t(long int Nx);
void free_Eigen_Solution_t(Eigen_Solution_t *eigen_solution, long int Nx);

//============================================================================

