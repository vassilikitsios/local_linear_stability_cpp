
//----------------------------------------------------------------------------
class Input_Deck_t
{
public:

	char *equations;						// Either 'laminar' or 'turbulent'
	char *base_flow;						// Either name of analytical base_flow or name of file
	
	double epsilon;
	double dissipation_ratio;				// dissipation_ratio
	
	char *top_boundary_condition;
	char *bottom_boundary_condition;
	char *run_analysis;						// run 'DoE' or 'search' or 'Floquet
	char *dispersion_relationship;			// required for run_analysis = search
	
	double tol;								// Parameters for Blasius profile for testing, to search for specific eigenvalue
	int max_iter, min_iter;					// interations for blasius and to search for specific eigenvalue
	double domain_size;
	double lambda;							// shear layer parameter
	double u_deficit, wake_half_width;		// wake parameters
	
	int N;									// Number of collocation points N+1
    double Re;								// Reference Reynolds Number
	
	double kx_real_min;						// Streamwise wavenumber range
	double kx_real_max;
	double kx_imag_min;
	double kx_imag_max;
	int N_kx_real;
	int N_kx_imag;
	complex<double> kx;
	bool search_for_same_mode_type;			// mode grouping code
	
	double omega_real_search;				// Values to search for
	double omega_imag_search;
	
	double kz_real_min;						// Spanwise wavenumber range
	double kz_real_max;
	double kz_imag;
	int N_kz_real;
	complex<double> kz;
	
	int number_of_eigenvectors_to_write;	// Eigenvector outputing

	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Input_Deck_t();
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Input_Deck_t();
	//-------------------------------------------------------------------------
};

