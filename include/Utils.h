//==========================
//binary_read_or_write.h
//==========================

using namespace std;

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

enum ReadOrWrite {Read, Write};

//****************************************************************************

extern void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, int *val);
extern void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, double *val);
extern void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, char *val);
extern void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, bool *val);

//==========================
//file_opening_closing.h
//==========================

//****************************************************************************
// files
//****************************************************************************
/** 
Opens a file in the given mode.
**/
extern FILE *open_file(char *filename, char *mode);
/** 
Opens a file in the given mode. 
**/
extern FILE *open_file(char *filename, char *mode, bool verbose);
/** 
Closes a file.
**/
extern void close_file(FILE *fp, char *filename);
/** 
Closes a file.
**/
extern void close_file(FILE *fp, char *filename, bool verbose);


//========================
//memory_allocations.h
//========================

//****************************************************************************
// memory allocations
//****************************************************************************
// int
extern int    *malloc_1d_array_int(long int Nx);
extern int   **malloc_2d_array_int(long int Nx, long int Ny);
extern int  ***malloc_3d_array_int(long int Nx, long int Ny, long int Nz);
extern int ****malloc_4d_array_int(int num_dimensions, long int Nx, long int Ny, long int Nz);
// int - fortran indices
extern int    *malloc_1d_array_int_fortran_indices(long int i0, long int i1);
extern int   **malloc_2d_array_int_fortran_indices(long int i0, long int i1, long int j0, long int j1);
extern int  ***malloc_3d_array_int_fortran_indices(long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
extern int ****malloc_4d_array_int_fortran_indices(long int I0, long int I1, long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
// long int
extern long int    *malloc_1d_array_long_int(long int Nx);
extern long int   **malloc_2d_array_long_int(long int Nx, long int Ny);
// bool
extern bool *malloc_1d_array_bool(long int Nx);
// double
extern double    *malloc_1d_array_double(long int Nx);
extern double   **malloc_2d_array_double(long int Nx, long int Ny);
extern double  ***malloc_3d_array_double(long int Nx, long int Ny, long int Nz);
extern double ****malloc_4d_array_double(int num_dimensions, long int Nx, long int Ny, long int Nz);
// double - fortran indices
extern double    *malloc_1d_array_double_fortran_indices(long int i0, long int i1);
extern double   **malloc_2d_array_double_fortran_indices(long int i0, long int i1, long int j0, long int j1);
extern double  ***malloc_3d_array_double_fortran_indices(long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
extern double ****malloc_4d_array_double_fortran_indices(long int I0, long int I1, long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
// char
extern char *malloc_1d_array_char(long int Nx);
extern char **malloc_2d_array_char(long int Nx, long int Ny);

//****************************************************************************
// memory re-allocations
//****************************************************************************
// double
extern double *realloc_1d_array_double(double *f_old, long int Nx);

//****************************************************************************
// memory de-allocations
//****************************************************************************
// double
extern void free_1d_array_double(double *f);
extern void free_2d_array_double(double **f, long int Nx);
extern void free_3d_array_double(double ***f, long int Nx, long int Ny);
extern void free_4d_array_double(double ****f, int num_dimensions, long int Nx, long int Ny);
// double - fortran indices
extern void free_1d_array_double_fortran_indices(double *f,long int i0, long int i1);
extern void free_2d_array_double_fortran_indices(double **f,long int i0, long int i1, long int j0, long int j1);
extern void free_3d_array_double_fortran_indices(double ***f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
extern void free_4d_array_double_fortran_indices(double ****f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1, long int l0, long int l1);
// int
extern void free_1d_array_int(int *f);
extern void free_2d_array_int(int **f, long int Nx);
extern void free_3d_array_int(int ***f, long int Nx, long int Ny);
extern void free_4d_array_int(int ****f, int num_dimensions, long int Nx, long int Ny);
// int - fortran indices
extern void free_1d_array_int_fortran_indices(int *f,long int i0, long int i1);
extern void free_2d_array_int_fortran_indices(int **f,long int i0, long int i1, long int j0, long int j1);
extern void free_3d_array_int_fortran_indices(int ***f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
extern void free_4d_array_int_fortran_indices(int ****f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1, long int l0, long int l1);
// long int
extern void free_1d_array_long_int(long int *f);
extern void free_2d_array_long_int(long int **f, long int Nx);
// bool
extern void free_1d_array_bool(bool *f);
//char
extern void free_1d_array_char(char *f);
extern void free_2d_array_char(char **f, long int Nx);

//==========
//misc.h
//==========

#include <cstring>

//****************************************************************************
// defines
//****************************************************************************
//#define PI 3.14159265358979
#define MAX_VAL(X,Y) ((X)>=(Y)?(X):(Y))
#define MIN_VAL(X,Y) ((X)<=(Y)?(X):(Y))
#define MAX(X,Y) (((X) > (Y) ? (X) : (Y)))
#define MIN(X,Y) (((X) < (Y) ? (X) : (Y)))

#define true	1
#define false	0

//#define PI 3.1415926535897932384626433832795
#define PI			M_PI
//#define E 2.7182818284590452353602874713527
#define E			M_E

#define BIG_NUMBER		HUGE_VAL
#define ZERO				1.0e-14
#define SMALL_NUMBER		-HUGE_VAL

#define SQ(X) ( (X) * (X) )
#define CU(X) ( (X) * (X) * (X) )
#define ABS(X) sqrt( SQ(X) )
#define SIGN(X) ABS(X)/X

// The maximum allowable line length whenever fgets or whatever is used to read
// from a file. 
#define MAX_LINE_LENGTH 1000

// The maximum length for a string which is expected to be fairly short
#define MAX_STRING_LENGTH 1000

// Convert a string to upper case
#define ucase(STR) {long int __ucase_i;for(__ucase_i=0; __ucase_i<strlen(STR); __ucase_i++) STR[__ucase_i]=toupper(STR[__ucase_i]);}


//****************************************************************************
// general
//****************************************************************************

// Terminates the program with an error message. The arguments are in the same format as printf().
extern void quit_error(char *fmt, ...);
void quit_error_mpi(int myrank, int nprocs, char *fmt, ...);

// Prints a warning to stderr. The arguments are in the same format as printf().
extern void print_warning(char *fmt, ...);
//param.h

#define MAX_PARAM_LINE_LENGTH 256

//typedef class Param_t;
//typedef class ParamListHead_t;

//****************************************************************************
//Param
//****************************************************************************
class Param_t
{
private:
	
public:
	//-------------------------------------------------------------------------
	//variables
	//-------------------------------------------------------------------------
	char *name;
	char *line;
	int numberOfAccesses;
	Param_t *next;
	Param_t *prev;
	//-------------------------------------------------------------------------
	//functions
	//-------------------------------------------------------------------------
	//functions
	Param_t(char *name, char *line); // constructor
	//-------------------------------------------------------------------------
	int get_int_param(char *name);
	int get_int_param(char *name, int defaultValue);
	//-------------------------------------------------------------------------
	int get_long_int_param(char *name);
	//-------------------------------------------------------------------------
	double get_double_param(char *name);
	//-------------------------------------------------------------------------
	char * get_string_param(char *name);
	char * get_param_line(char *name);
	//-------------------------------------------------------------------------
	double get_double_param(char *name, double defaultValue);
	//-------------------------------------------------------------------------
	Param_t * get_param(char *name);
	Param_t * get_next_param(char *name);
};

//****************************************************************************
//Param
//****************************************************************************
class ParamListHead_t
{
private:
	
public:
	//-------------------------------------------------------------------------
	//variables
	//-------------------------------------------------------------------------
	Param_t *param_list;
	int num_params;
	//-------------------------------------------------------------------------
	//functions
	//-------------------------------------------------------------------------
	ParamListHead_t(); //contructor
	//-------------------------------------------------------------------------
	void addParam(char *param_line);
	//-------------------------------------------------------------------------
	void print_params();
	void print_param_accesses();
	//-------------------------------------------------------------------------
	void get_params_from_file(char *filename);
	void get_params_from_file(char *filename, bool verbose);
	void get_params_from_command_line(int argc, char *argv[]);
	//-------------------------------------------------------------------------
	bool check_param(char *name);
	//-------------------------------------------------------------------------
	int get_int_param(char *name);
	int get_int_param(char *name, int defaultValue);
	//-------------------------------------------------------------------------
	int get_long_int_param(char *name);
	//-------------------------------------------------------------------------
	double get_double_param(char *name);
	double get_double_param(char *name, double defaultValue);
	//-------------------------------------------------------------------------
	Param_t * get_param(char *name);
	Param_t * get_next_param(char *name);
	char * get_string_param(char *name);
	char * get_param_line(char *name);
	//-------------------------------------------------------------------------
	Param_t * get_nth_param(int n, char *name);
	int get_nth_int_param(int n, char *name);
	long int get_nth_long_int_param(int n, char *name);
	double get_nth_double_param(int n, char *name);
};



#include <sstream>  // Required for stringstreams
#include <string>

extern int strcmp_case_insensitive(char *s1, char *s2);

extern char *copy_string(char *str1);

extern char *concatenate_strings(char *str1, char *str2);

extern string convert_int_to_string(int number);
extern string convert_int_to_zero_padded_string(int number);

extern double char_to_double(char *input);
// time_stats.h
// ****************************************************************************
#include <ctime>
// ****************************************************************************
//typedef class TimeStats;
// ****************************************************************************
// TimeStats
// ****************************************************************************
class TimeStats
{
private:
	time_t time_start;
	long int computationStartStep;
	long int step_start;
	long int step_finish;
	long int total_steps;
	inline void convert_seconds(double total_seconds, int &days, int &hours, int &minutes, int &seconds);
public:
	char time_spent_string[40];
	char time_remaining_string[40];
	double percentage_complete;
	double timesteps_per_minute;
		
	TimeStats(long int step_start, long int step_finish);
	void reset(long int step_start, long int step_finish);
	void produce_stats(long int step_current);
};
// ****************************************************************************




// tridiagonal_solver.h

extern void tridiagonal_solver_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd);
extern void tridiagonal_solver_nondestructive_static(int N, double *a, double *b, double *c, double *d, double *x);
extern void tridiagonal_solver_destructive(int N, double *a, double *b, double *c, double *d, double *x);
extern void tridiagonal_solver_periodic_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd, double *C);
extern void tridiagonal_solver_periodic_nondestructive_static(int N, double *a, double *b, double *c, double *d, double *x);
extern void tridiagonal_solver_periodic_destructive(int N, double *a, double *b, double *c, double *d, double *x);

extern void tridiagonal_solver_test(void);
extern void tridiagonal_solver_periodic_test(void);


extern int matrix_invert(double *A, double *Ainv, long int N);
extern void matrix_display(double *A, long int R, long int C);
extern void matrix_multiply(double *A, long int Arows, long int Acols, double *B, long int Brows, long int Bcols, double *product);
extern double matrix_determinant(double *A, long int N);
extern void display_augmented(double *A, double *B, long int N);
extern void matrix_transpose(double *A, long int R, long int C);
extern int matrix_cholesky(double *A, double *L, long int N);
extern void matrix_copy(double *A, long int R, long int C, double *Cp);
extern int matrix_cholesky_invert(double *A, double *Ainv, long int N);
extern void matrix_constant_multiply(double *A, long int Arows, long int Acols, double constant, double *product);
extern void matrix_add(double *A, long int Arows, long int Acols, double *B, long int Brows, long int Bcols, double *sum);
extern long int find_null_space_vectors(double *A, long int rows, long int cols, double **P, double absolute_tolerance);
extern void svd(double *A, long int rows, long int cols, double *U, double *D, double *V);
extern int matrix_row_reduce(double *A, long int rows, long int cols, double tolerance);
extern long int matrix_nullity(double *A, long int rows, long int cols, double absolute_tolerance);
extern long int find_null_space_decomposition(double *A, long int rows, long int cols, double **D, double **V, double tolerance);
extern long int find_null_space_decomposition(double *A, long int rows, long int cols, double tolerance);


//****************************************************************************
// Cubic_Spline_t
//****************************************************************************
class Cubic_Spline_t
{
private:	
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	double *a, *b, *c, *d, *points;
	//-------------------------------------------------------------------------	
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	long int num_points;			// Numer of points in the spline
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Cubic_Spline_t(long int num_points);
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------	
	void calculate_interpolation_function(double *points, double *data);
	double evaluate_function(double interp_point);
	double evaluate_first_derivative(double interp_point);
	double evaluate_second_derivative(double interp_point);
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Cubic_Spline_t();
	//-------------------------------------------------------------------------
};
double calculate_second_order_centred_finite_difference(double delta_x, double f_i_minus_1, double f_i_plus_1);
double calculate_forth_order_centred_finite_difference(double delta_x, double f_i_minus_2, double f_i_minus_1, double f_i_plus_1, double f_i_plus_2);
double calculate_semi_centred_forward_finite_difference(double delta_x, double f_i_minus_1, double f_i, double f_i_plus_1, double f_i_plus_2);
double calculate_semi_centred_backward_finite_difference(double delta_x, double f_i_minus_2, double f_i_minus_1, double f_i, double f_i_plus_1);
double calculate_second_order_forward_finite_difference(double delta_x, double f_i, double f_i_plus_1, double f_i_plus_2);
double calculate_second_order_backward_finite_difference(double delta_x, double f_i_minus_3, double f_i_minus_2, double f_i_minus_1, double f_i);
double calculate_third_order_forward_finite_difference(double delta_x, double f_i, double f_i_plus_1, double f_i_plus_2, double f_i_plus_3);
double calculate_third_order_backward_finite_difference(double delta_x, double f_i_minus_3, double f_i_minus_2, double f_i_minus_1, double f_i);

//============================================================================
// Dist_t
//============================================================================
class Dist_t
	{
	public:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		double delta_x;
		double delta_y;
		double delta_z;
		double delta_mag;
		long int node_num;
		//-------------------------------------------------------------------------	
	};

//============================================================================
// Interpolation_t
//============================================================================
class Interpolation_t
{
private:	
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	int order;						// Polynomial Order
	int p;							// Number of terms in polynomial
	int max_redundancy;				// Maximum number of additional points required above that needed to calculate a polynomial surface of order "order"
	int m;							// Number of points used in stencil
	double epsilon;					// Constant used for Gaussian weights
	double singularity_tolerance;	// Tolerance for determining if a stencil is singular or not
	bool two_PI_periodic;			// Flag if mesh is 2PI periodic
	
	long int num_data_points;		// Number of points in the source mesh
	double *D;						// Multiply this matrix by the data points in the stencil to get the results of the interpolation
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------		
	void order_closest_nodes_on_extruded_2D_data(double **x_data, double x_interp_point, double y_interp_point);
	void order_closest_nodes_on_2D_data(double **x_data, double x_interp_point, double y_interp_point);
	
	void generate_stencil_and_test_stability(double x_interp_point, double y_interp_point);
	void generate_stencil(double x_interp_point, double y_interp_point);

	void generate_w_matrix(double *w);
	void generate_w_matrix_for_continuity(double *w);
	
	void generate_B_matrix(double *B);
	void generate_B_matrix_for_continuity(double *B);
	void generate_b_matrix(double *b);
	void generate_b_x_matrix(double *b_x);
	void generate_b_y_matrix(double *b_y);

	void generate_C_matrix(double *interpolation_data, double *C);
	void generate_C_matrix_for_continuity(double *u, double *v, double *C);	
	//-------------------------------------------------------------------------	
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	Dist_t *distance_to_points;		// Sorted list of distance to points used to determine the stencil of closest points to the desired interpolation point
	double f, df_dx, df_dy, d2f_dx2, d2f_dxdy, d2f_dy2;			// Results of interpolation
	double u, du_dx, du_dy, d2u_dx2, d2u_dxdy, d2u_dy2;			// Results of interpolation for continuity conservation
	double v, dv_dx, dv_dy, d2v_dx2, d2v_dxdy, d2v_dy2;
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Interpolation_t(long int num_data_points, int order, int max_redundancy, double epsilon, double singularity_tolerance, bool two_PI_periodic);
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------	
	void form_interpolation_function_on_extruded_2D_data(double **x_data, double x_interp_point, double y_interp_point);
	void form_interpolation_function_on_2D_data(double **x_data, double x_interp_point, double y_interp_point);

	void evaluate_interpolation_function(double *interpolation_data);
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Interpolation_t();
	//-------------------------------------------------------------------------
};

//============================================================================
// Non class function
//============================================================================
Dist_t* malloc_Dist_t(long int N);
void free_Dist_t(Dist_t *dist);
int qsort_distance_compare_function(Dist_t *dist1, Dist_t *dist2);

//============================================================================
// Interpolation_1D_t
//============================================================================
class Interpolation_1D_t
{
private:	
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	int order;						// Polynomial Order
	int p;							// Number of terms in polynomial
	int max_redundancy;				// Maximum number of additional points required above that needed to calculate a polynomial surface of order "order"
	int m;							// Number of points used in stencil
	double epsilon;					// Constant used for Gaussian weights
	double singularity_tolerance;	// Tolerance for determining if a stencil is singular or not
	
	long int num_data_points;		// Number of points in the source mesh
	double *D;						// Multiply this matrix by the data points in the stencil to get the results of the interpolation
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------
	void generate_stencil(double x_interp_point);
	void order_closest_points(double *x_data, double x_interp_point);
	void generate_w_matrix(double *w);
	void generate_B_matrix(double *B);
	void generate_C_matrix(double *interpolation_data, double *C);
	//-------------------------------------------------------------------------	
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	Dist_t *distance_to_points;			// Sorted list of distance to points used to determine the stencil of closest points to the desired interpolation point
	double f, df_dx, d2f_dx2;			// Results of interpolation
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Interpolation_1D_t(long int num_data_points, int order, int max_redundancy, double epsilon, double singularity_tolerance);
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------	
	void form_interpolation_function(double *x_data, double x_interp_point);
	void evaluate_interpolation_function(double *interpolation_data);
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Interpolation_1D_t();
	//-------------------------------------------------------------------------
};
