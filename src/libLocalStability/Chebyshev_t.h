
//============================================================================
class Chebyshev_t
{
private:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	double l, s;																	// Required for grid clustering
	//-------------------------------------------------------------------------
	// helper functions
	//-------------------------------------------------------------------------
	void generate_chebyshev_D1();
	void generate_chebyshev_D2();
	void generate_chebyshev_D3();
	void generate_chebyshev_D4();
	void correct_diagonal_terms(Array<double,2>& D);
	void approximate_w();
	//-------------------------------------------------------------------------
	void scale_derivatives_for_finite_domain(double domain_size);
	void scale_derivatives_for_unclustered_finite_domain(double domain_size);
	void scale_derivatives_for_semi_finite_domain(double domain_size);
	void scale_derivatives_for_clustered_semi_finite_domain();
	void scale_derivatives_for_unclustered_semi_finite_domain(double domain_size);
	//-------------------------------------------------------------------------
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	int N;
	double *xi, *eta;
	Array<double,2> D1, D2, D3, D4, w;
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Chebyshev_t();
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------	
	void allocate_memory(int N);
	void generate_colocation_points();
	void generate_chebyshev_integral_weights();
	void generate_chebyshev_derivative_matrices();
	//-------------------------------------------------------------------------
	void scale_for_finite_domain(double domain_size);
	void scale_geometry_for_finite_domain(double domain_size);
	void scale_geometry_for_finite_domain(double y_min, double y_max);
	//-------------------------------------------------------------------------
	void scale_for_semi_finite_domain(double domain_size);
	void scale_geometry_for_semi_finite_domain(double domain_size);
	//-------------------------------------------------------------------------
	void scale_for_clustered_semi_finite_domain(double domain_size, double half_grid_point1);
	void scale_geometry_for_clustered_semi_finite_domain(double domain_size, double half_grid_point);
	//-------------------------------------------------------------------------
	void scale_for_unclustered_finite_domain(double domain_size);
	void scale_geometry_for_unclustered_finite_domain(double domain_size);
	void scale_geometry_for_unclustered_finite_domain(double y_min, double y_max);
	//-------------------------------------------------------------------------
	void scale_for_unclustered_semi_finite_domain(double domain_size);
	void scale_geometry_for_unclustered_semi_finite_domain(double domain_size);
	//-------------------------------------------------------------------------
	void test_chebyshev_matrices_with_exponential_function(double q);
	void test_chebyshev_matrices_with_eta_power_4();
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Chebyshev_t();
	//-------------------------------------------------------------------------
};

//============================================================================
