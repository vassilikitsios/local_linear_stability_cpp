
#include "EigenSolver.h"

#define artifical_eigenvalue_shift		200.0
#define rel_error_bound_tol				1.0e-6

//============================================================================
class Linear_Navier_Stokes_t
{
private:
	//-------------------------------------------------------------------------
	// problem specificiation
	//-------------------------------------------------------------------------
	char *equations;
	bool with_non_parallel_terms;
	bool with_w_modes;
	double Re;
	double dissipation_ratio_guess;
	complex<double> kx, kz;
	double k, k2;
	int N_size;									// Size of matrix
	char *top_boundary_condition;				// wall or freestream
	char *bottom_boundary_condition;
	void dimensionalise_solution();
	
	//-------------------------------------------------------------------------
	// linear operator
	//-------------------------------------------------------------------------
	Chebyshev_t cheb;
	complex<double> im, one, unit;									// helper variables
	Array<complex<double>,2> A, B;
	void build_laminar_matrix_B();
	void build_laminar_matrix_A();
	void build_turbulent_matrix_B();
	void build_turbulent_matrix_A_order1();
	void build_turbulent_matrix_A_order2();
	void build_turbulent_matrix_A_order3();
	
	//-------------------------------------------------------------------------
	// boundary conditions
	//-------------------------------------------------------------------------	
	void apply_boundary_conditions();
	void apply_boundary_conditions_for_top_wall();
	void apply_boundary_conditions_for_bottom_wall();
	void apply_boundary_conditions_for_top_freestream();
	void apply_boundary_conditions_for_bottom_freestream();
	
	//-------------------------------------------------------------------------
	// outputing
	//-------------------------------------------------------------------------
	int number_of_eigenvectors_to_write;
	void write_eigenvalues(int kx_real_num, int kx_imag_num, int iteration);
	void write_eigenvectors(int kx_real_num, int kx_imag_num, int iteration);
	void write_eigenvector_header(string filename_out, int eigenvector_num);
	void write_eigenvector(string filename_out, int eigenvector_num);
	
	string get_eigenvalue_filename(int kx_real_num, int kx_imag_num, int iteration);
	string get_eigenvector_filename(int eigenvector_num, int kx_real_num, int kx_imag_num, int iteration); 
	string get_padded_number_string(int num);
	
	//-------------------------------------------------------------------------
public:	
	//-------------------------------------------------------------------------
	Base_Flow_t base_flow;
	Eigen_Solver_t eigen_solver;
	int N;										// Number of colocation points
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Linear_Navier_Stokes_t();
	//-------------------------------------------------------------------------
	// preliminaries
	//-------------------------------------------------------------------------
	void allocate_memory(int N);
	void set_problem_parameters(double Re, complex<double> kx, complex<double> kz, double dissipation_ratio_guess, double L_0);
	void set_number_of_eigenvectors_to_write(int number_of_eigenvectors_to_write);
	void set_equations(char* equations);
	void set_boundary_conditions(char *top_boundary_condition, char *bottom_boundary_condition);
	void set_with_non_parallel_terms(bool with_non_parallel_terms);
	void generate_base_flow(double Re);
	void generate_linear_operator();
	void add_non_parallel_terms();
	
	//-------------------------------------------------------------------------
	// Local EVP
	//-------------------------------------------------------------------------
	void calculate_eigensolution();
	complex<double> get_eigensolution_closest_to_this_omega(complex<double> this_omega, int *ptr_eigenvalue_number);
	
	//-------------------------------------------------------------------------
	// I/O
	//-------------------------------------------------------------------------
	void write_output(int kx_real_num, int kx_imag_num);
	void write_output(int kx_real_num, int kx_imag_num, int iteration);
	
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Linear_Navier_Stokes_t();
	//-------------------------------------------------------------------------
};

//============================================================================
int qsort_compare_function(Eigen_Solution_t *es1, Eigen_Solution_t *es2);

//============================================================================
