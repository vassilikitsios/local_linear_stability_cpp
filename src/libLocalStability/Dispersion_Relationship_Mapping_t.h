
//============================================================================
class Eigen_Values_t
	{
	public:
		double kx_r, kx_i, kz_r, kz_i;
		double Omega_r, Omega_i, c_r, c_i;
		double c_gr_r; // d Omega_r / d kx_r
		double c_gr_i; // d Omega_i / d kx_r
		double rel_error_bounds, cond;
		int discrete_continuous, w_mode, symmetric;
		bool sorted;
	};

//============================================================================
class Dispersion_Relationship_Mapping_t
	{
	private:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		int N_kx_real, N_kx_imag, N_kz_real;								// Number of eigensolutions
		int N_jj_loop;
		int N_groups;
		int N_eigenvalues;													// Number of eigenvalues per solution
		double Re;															// To report the Re in the final outputing
		double x;															// Streamwise position
		bool search_for_same_mode_type;
		
		Array<complex<double>,2> omega_read, kx, kz;						// Data from most unstable mode
		
		//-------------------------------------------------------------------------
		// helper mode grouping functions
		//-------------------------------------------------------------------------
		string get_eigenvalue_filename(int kx_real_num, int kx_imag_num);
		string get_complex_map_filename(int number);
		string get_real_omega_filename(int number);
		string get_real_omega_raw_filename(int number);
		string get_real_omega_max_spatial_growth_filename(int number);
		string get_pinch_point_filename(int number);
		string get_padded_number_string(int number);
		string get_max_omega_imag_for_kx_real_filename(int number);
		
		void read_data(Eigen_Values_t **eigenvalue_data);
		int read_file(string filename, int file, Eigen_Values_t **eigenvalue_data);

		void sort_modes_old(Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes);
		void sort_modes(Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes);
		
		bool determine_if_this_a_valid_mode(int mode, Eigen_Values_t **sorted_modes);
		void search_forward(int first_group, int mode, Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes);
		void search_backward(int first_group, int mode, Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes);
		void find_closet_mode_in_group(int mode, int group, int group_check, Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes);
		
		void calculate_group_velocity(int mode, Eigen_Values_t **sorted_modes);
		void find_pinch_point(int mode, double x_pos, Eigen_Values_t **sorted_modes);
		
		void write_complex_map(int mode, double nu, Eigen_Values_t **sorted_modes);
		
		void write_zero_line();
		
		void write_real_omega(int mode, double nu, Eigen_Values_t **sorted_modes, double x_pos);
		int find_where_omega_is_completely_real(int N_interp_points, int mode, Eigen_Values_t **sorted_modes,
		double *omega_real_interp, double *kx_real_interp, double *kx_imag_interp, double *c_gr_real_interp, double *c_gr_imag_interp);
		void interpolate_along_spline(int N_spline_points, int num_zeros, int mode,
		double *omega_real_interp, double *kx_real_interp, double *kx_imag_interp, double *c_gr_real_interp, double *c_gr_imag_interp,
		double *omega_real_spline, double *kx_real_spline, double *kx_imag_spline, double *c_gr_real_spline, double *c_gr_imag_spline, double x_pos);
		
		void write_max_omega_imag_for_kx_real(int mode, Eigen_Values_t **sorted_modes, double x_pos);
		
		//-------------------------------------------------------------------------
	public:
		//-------------------------------------------------------------------------
		// constructor
		//-------------------------------------------------------------------------
		Dispersion_Relationship_Mapping_t(int N_eigenvalues, double Re, double x, bool search_for_same_mode_type,
				int N_kx_real, double kx_real_max, double kx_real_min,
				int N_kx_imag, double kx_imag_max, double kx_imag_min,
				int N_kz_real, double kz_real_max, double kz_real_min, double kz_imag);
		//-------------------------------------------------------------------------
		// functions
		//-------------------------------------------------------------------------
		complex<double> get_kx(int ii, int jj);
		complex<double> get_kz(int ii, int jj);
		double L_0;						// effective length of first wave
		
		int get_N_jj_loop();
		
		void group_modes();
		
		void read_dispersion_relatioship(string filename);
		complex<double> estimate_kx_closest_to_this_omega(complex<double> this_omega, double *ptr_min_omega_error, complex<double> *ptr_omega_orig);
		//-------------------------------------------------------------------------
		// destructor
		//-------------------------------------------------------------------------
		~Dispersion_Relationship_Mapping_t();
		//-------------------------------------------------------------------------
	};

//============================================================================
// Non-class functions
//============================================================================
Eigen_Values_t* malloc_1d_array_Eigen_Values_t(int Nx);
Eigen_Values_t** malloc_2d_array_Eigen_Values_t(int Nx, int Ny);
void free_1d_array_Eigen_Values_t(Eigen_Values_t *eigen_values);
void free_2d_array_Eigen_Values_t(Eigen_Values_t **eigen_values, int Nx);
void copy_2d_array_Eigen_Values_t_from_A_to_B(Eigen_Values_t **eigen_values_A, int A1, int A2, Eigen_Values_t **eigen_values_B, int B1, int B2);

//============================================================================
