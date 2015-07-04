
//============================================================================
class Base_Flow_t
{
private:	
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	bool unsteady;								// use triple decomposition corrections
	double lambda;								// for shear layer profile
	double u_deficit, wake_half_width;			// for wake profile
	double epsilon;								// cut off parameter
	int profile_num;							// for data read from file
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------
	void generate_poiseuille_profile(Chebyshev_t *cheb);
	void generate_couette_profile(Chebyshev_t *cheb);
	void generate_turbulent_channel_profile(Chebyshev_t *cheb);
	void generate_blasius_profile(Chebyshev_t *cheb);
	void generate_shear_layer_profile(Chebyshev_t *cheb);
	void generate_wake_profile(Chebyshev_t *cheb);
	void read_profile_from_file(Chebyshev_t *cheb);
	void read_profile_stresses_from_file(Chebyshev_t *cheb);
	void read_profile_eddy_viscosity_from_file(Chebyshev_t *cheb);
	void check_continuity();
	//-------------------------------------------------------------------------
public:
	//-------------------------------------------------------------------------
	// Variables
	//-------------------------------------------------------------------------
	char *base_flow_type;
	
	double tol;													// tolerence
	int max_iter, min_iter;										// max and min iterations
	
	int N;
	double u_ref, l_ref;
	double domain_size;
	
	double x;													// Streamwise position
	double *n;													// Wall normal position scaled by reference length
	double *u, *u_x, *u_y, *u_xx, *u_xy, *u_yy;					// Wall tangential velocity and derivatives
	double *v, *v_x, *v_y, *v_xx, *v_xy, *v_yy;					// Wall normal velocity and derivatives
	double *u_rms, *v_rms, *w_rms;								// rms profiles

	double nu;																	// kinematic viscosity
	double *Ek, *Ek_y, *Ek_yy;													// turbulent kinetic energy
	double *nuE1, *nuE1_y, *nuE1_original, *nuE1_y_original;					// representation of 1st order eddy viscosity and 'y' derivative
	double *constE2, *constE2_y, *constE2_original, *constE2_y_original;		// representation of 2nd order eddy viscosity and 'y' derivative
	double *constE3, *constE3_y, *constE3_original, *constE3_y_original;		// representation of 3rd order eddy viscosity and 'y' derivative
	
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Base_Flow_t();
	//-------------------------------------------------------------------------
	// Functions
	//-------------------------------------------------------------------------	
	void allocate_memory(int N);
	
	void set_base_flow_type(char* base_flow_type);
	void get_base_flow_type(char* base_flow_type_out);
	void set_profile_positions(double *n_input);
	void set_domain_size(double domain_size);
	void set_shear_layer_parameter(double lambda);
	void set_wake_parameter(double u_deficit, double wake_half_width);
	void set_iteration_parameters(double tol, int max_iter, int min_iter);
	void set_unsteady_bool(bool unsteady);
	void update_base_flow(double *u, double *u_y, double *v, double *v_y);

	void set_eddy_viscosity_parameters(double epsilon);
	void update_eddy_viscosity_fields(double dissipation_ratio);
	
	void generate_base_flow(Chebyshev_t *cheb, double nu);
	int get_number_of_colocation_points_in_input_data();
	void non_dimensionalise_base_flow();
	
	void write_profile();
	
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Base_Flow_t();
	//-------------------------------------------------------------------------
};

//============================================================================

