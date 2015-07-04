#include <blitz/array.h>
#include <fstream>

using namespace blitz;

#include "Utils.h"
#include "Chebyshev_t.h"
#include "Base_Flow_t.h"

//----------------------------------------------------------------------------
Base_Flow_t::Base_Flow_t()
{
}

//----------------------------------------------------------------------------
void Base_Flow_t::allocate_memory(int N)
{
	this->N = N;
	
	n = malloc_1d_array_double(N+1);
	
	u = malloc_1d_array_double(N+1);	
	u_x = malloc_1d_array_double(N+1);	u_y = malloc_1d_array_double(N+1);
	u_xx = malloc_1d_array_double(N+1);	u_xy = malloc_1d_array_double(N+1);	u_yy = malloc_1d_array_double(N+1);
	
	v = malloc_1d_array_double(N+1);	
	v_x = malloc_1d_array_double(N+1);	v_y = malloc_1d_array_double(N+1);
	v_xx = malloc_1d_array_double(N+1);	v_xy = malloc_1d_array_double(N+1);	v_yy = malloc_1d_array_double(N+1);
	
	u_rms = malloc_1d_array_double(N+1);		v_rms = malloc_1d_array_double(N+1);		w_rms = malloc_1d_array_double(N+1);
	
	Ek = malloc_1d_array_double(N+1);					Ek_y = malloc_1d_array_double(N+1);				Ek_yy = malloc_1d_array_double(N+1);
	
	nuE1_original = malloc_1d_array_double(N+1);		nuE1_y_original = malloc_1d_array_double(N+1);
	constE2_original = malloc_1d_array_double(N+1);		constE2_y_original = malloc_1d_array_double(N+1);
	constE3_original = malloc_1d_array_double(N+1);		constE3_y_original = malloc_1d_array_double(N+1);
	
	nuE1 = malloc_1d_array_double(N+1);		nuE1_y = malloc_1d_array_double(N+1);
	constE2 = malloc_1d_array_double(N+1);	constE2_y = malloc_1d_array_double(N+1);
	constE3 = malloc_1d_array_double(N+1);	constE3_y = malloc_1d_array_double(N+1);
}

//----------------------------------------------------------------------------
Base_Flow_t::~Base_Flow_t()
{
	free_1d_array_double(n);
	
	free_1d_array_double(u);	
	free_1d_array_double(u_x);	free_1d_array_double(u_y);
	free_1d_array_double(u_xx);	free_1d_array_double(u_xy);	free_1d_array_double(u_yy);
	
	free_1d_array_double(v);	
	free_1d_array_double(v_x);	free_1d_array_double(v_y);
	free_1d_array_double(v_xx);	free_1d_array_double(v_xy);	free_1d_array_double(v_yy);

	free_1d_array_double(u_rms);		free_1d_array_double(v_rms);			free_1d_array_double(w_rms);
	
	free_1d_array_double(Ek);	free_1d_array_double(Ek_y);		free_1d_array_double(Ek_yy);
	
	free_1d_array_double(nuE1_original);		free_1d_array_double(nuE1_y_original);
	free_1d_array_double(constE2_original);		free_1d_array_double(constE2_y_original);
	free_1d_array_double(constE3_original);		free_1d_array_double(constE3_y_original);
	
	free_1d_array_double(nuE1);		free_1d_array_double(nuE1_y);
	free_1d_array_double(constE2);	free_1d_array_double(constE2_y);
	free_1d_array_double(constE3);	free_1d_array_double(constE3_y);
}

//----------------------------------------------------------------------------
void Base_Flow_t::set_base_flow_type(char* base_flow_type)	{ this->base_flow_type = base_flow_type; }

//----------------------------------------------------------------------------
void Base_Flow_t::get_base_flow_type(char* base_flow_type)	{ base_flow_type = this->base_flow_type; }

//----------------------------------------------------------------------------
void Base_Flow_t::set_domain_size(double domain_size) { this->domain_size = domain_size; }

//----------------------------------------------------------------------------
void Base_Flow_t::set_shear_layer_parameter(double lambda) { this->lambda = lambda; }

//----------------------------------------------------------------------------
void Base_Flow_t::set_eddy_viscosity_parameters(double epsilon) { this->epsilon = epsilon; }

//----------------------------------------------------------------------------
void Base_Flow_t::set_wake_parameter(double u_deficit, double wake_half_width)
{
	this->u_deficit = u_deficit;
	this->wake_half_width = wake_half_width;
}

//----------------------------------------------------------------------------
void Base_Flow_t::set_iteration_parameters(double tol, int max_iter, int min_iter)
{
	this->tol = tol;
	this->max_iter = max_iter;
	this->min_iter = min_iter;
}

//----------------------------------------------------------------------------
void Base_Flow_t::set_profile_positions(double *n_input) { for (long int ii=0; ii<N+1; ii++) n[ii] = n_input[ii]; }

//----------------------------------------------------------------------------
void Base_Flow_t::set_unsteady_bool(bool unsteady) { this->unsteady = unsteady; }

//----------------------------------------------------------------------------
void Base_Flow_t::generate_base_flow(Chebyshev_t *cheb, double nu)
{
	this->nu = nu;
	
	if (strcmp(base_flow_type,"couette")==0) {
		cout << "generating couette profile..." << endl;
		generate_couette_profile(cheb);
	} else if (strcmp(base_flow_type,"channel_turb")==0) {
		cout << "generating turbulent channel profile..." << endl;
		generate_turbulent_channel_profile(cheb);
	} else if (strcmp(base_flow_type,"poiseuille")==0) {
		cout << "generating poiseuille profile..." << endl;
		generate_poiseuille_profile(cheb);
	} else if (strcmp(base_flow_type,"shear_layer")==0) {
		cout << "generating shear layer profile..." << endl;
		generate_shear_layer_profile(cheb);
	} else if (strcmp(base_flow_type,"wake")==0) {
		cout << "generating wake profile..." << endl;
		generate_wake_profile(cheb);
	} else if (strcmp(base_flow_type,"blasius")==0) {
		cout << "generating blasius profile..." << endl;
		generate_blasius_profile(cheb);
	} else {
		cout << "reading data from file..." << endl;
		read_profile_from_file(cheb);
		read_profile_stresses_from_file(cheb);
//		read_profile_eddy_viscosity_from_file(cheb);
	}
	
	check_continuity();
	
	cout << "u_ref = " << u_ref << endl;
	cout << "l_ref = " << l_ref << endl;
}

//----------------------------------------------------------------------------
void Base_Flow_t::generate_wake_profile(Chebyshev_t *cheb)
{
	double sech_n;
	u_ref = 1.0;
	l_ref = 1.0;
//	double wake_half_width = 1.0;// / log( 1.0 + sqrt(2.0) );
//	cheb->scale_for_finite_domain(domain_size/wake_half_width);
	cheb->scale_for_finite_domain(domain_size);
	for (long int ii=0; ii<=cheb->N; ii++) {
		n[ii] = cheb->eta[ii];
		sech_n = 1.0 / cosh(n[ii] / wake_half_width);
		u[ii] = 1.0 - u_deficit * SQ(sech_n);
		u_x[ii] = 0.0;
		u_y[ii] = 2.0 * u_deficit * SQ(sech_n) * tanh(n[ii]);
		v[ii] = 0.0;
		v_x[ii] = 0.0;
		v_y[ii] = 0.0;
	}
}

//----------------------------------------------------------------------------
void Base_Flow_t::generate_couette_profile(Chebyshev_t *cheb)
{
	for (long int ii=0; ii<=cheb->N; ii++) {
		cheb->eta[ii] = cheb->xi[ii];
		n[ii] = cheb->eta[ii];
		u[ii] = n[ii];
		u_x[ii] = 0.0;
		u_y[ii] = 1.0;
		v[ii] = 0.0;
		v_x[ii] = 0.0;
		v_y[ii] = 0.0;
	}
	u_ref = 1.0;	l_ref = 1.0;	domain_size = 2.0;
}

//----------------------------------------------------------------------------
void Base_Flow_t::generate_turbulent_channel_profile(Chebyshev_t *cheb)
{
	for (long int ii=0; ii<=cheb->N; ii++) {
		cheb->eta[ii] = cheb->xi[ii];
		n[ii] = cheb->eta[ii];
		u[ii] = 1.14*(2.01/PI/4.0*log(ABS(cos(n[ii]*PI/2.01))) / 7.7603807704e-01 + 1.0);
		u_x[ii] = 0.0;
		u_y[ii] = -tan(n[ii]*PI/2.01)/4.0 / 7.7603807704e-01 * 1.14;
		v[ii] = 0.0;
		v_x[ii] = 0.0;
		v_y[ii] = 0.0;
	}
	u[0] = 0.0;
	u[cheb->N] = 0.0;
	u_ref = 1.0;	l_ref = 1.0;	domain_size = 2.0;
}

//----------------------------------------------------------------------------
void Base_Flow_t::generate_poiseuille_profile(Chebyshev_t *cheb)
{
	u_ref = 1.0;	l_ref = 1.0;	domain_size = 2.0;
	//cheb->scale_for_unclustered_finite_domain(domain_size);
	for (long int ii=0; ii<=cheb->N; ii++) {
		cheb->eta[ii] = cheb->xi[ii];
		n[ii] = cheb->eta[ii];
		u[ii] = 1.0 - pow(n[ii],2.0);
		u_x[ii] = 0.0;
		u_y[ii] = -2.0 * n[ii];
		v[ii] = 0.0;
		v_x[ii] = 0.0;
		v_y[ii] = 0.0;
	}
}

//----------------------------------------------------------------------------
void Base_Flow_t::generate_shear_layer_profile(Chebyshev_t *cheb)
{
	quit_error((char*)"%s base flow not written yet", base_flow_type);
	cheb->scale_for_finite_domain(domain_size);
	for (int ii=0; ii<=N; ii++) {
		n[ii] = cheb->eta[ii];
		u[ii] = 1.0 + lambda * tanh(n[ii]);
		u_x[ii] = 0.0;
		u_y[ii] = lambda * ( 1.0 - pow(tanh(n[ii]), 2.0) );
		v[ii] = 0.0;
		v_x[ii] = 0.0;
		v_y[ii] = 0.0;
	}
}

//----------------------------------------------------------------------------
void Base_Flow_t::non_dimensionalise_base_flow() 
{	
	for(long int j=0 ; j<=N ; j++) {
		n[j] = n[j] / l_ref;
		u[j] = u[j] / u_ref;
		v[j] = v[j] / u_ref;
		u_y[j] = u_y[j] / u_ref * l_ref;
		u_x[j] = u_x[j] / u_ref * l_ref;
		v_y[j] = v_y[j] / u_ref * l_ref;
		v_x[j] = v_x[j] / u_ref * l_ref;
	}
}

//----------------------------------------------------------------------------
void Base_Flow_t::update_base_flow(double *u, double *u_y, double *v, double *v_y)
{
	for (long int ii=0; ii<=N; ii++) {
		this->u[ii] = u[ii];
		this->u_y[ii] = u_y[ii];
		this->v[ii] = v[ii];
		this->v_y[ii] = v_y[ii];
	}
}

//----------------------------------------------------------------------------
void Base_Flow_t::check_continuity()
{		
	double sum = 0.0;
	for(long int j=0 ; j<=N ; j++) sum = u_x[j] + v_y[j];
	cout << "continuity check (should be zero) = " << sum << endl;
}

//----------------------------------------------------------------------------
// blasius helper functions
//----------------------------------------------------------------------------
Array<double,2> evaluate_g(Array<double,2>& f, Chebyshev_t *cheb);
Array<double,2> evaluate_jacobian(Array<double,2>& f, Chebyshev_t *cheb);
Array<double,2> solve_linear_system(Array<double,2>& Bb, Array<double,2>& Ab, int num);
void test_solve_linear_system(Array<double,2>& b, Array<double,2>& A, int num);
void reconstruct_velocities(Chebyshev_t *cheb, Array<double,2>& f, 
							double *u, double *u_x, double *u_y, double *v, double *v_x, double *v_y);

//----------------------------------------------------------------------------
void Base_Flow_t::generate_blasius_profile(Chebyshev_t *cheb)
{
	//............................................................................	
	// Can only calculate the Blasius BL using a semi-finite coordinate transformation due to formulation
	// If want to use a semi-infinite domain for the OSSQ equations, must first calculate the Blasius profile on a finite domain, 
	// then interpolate it onto the semi-infinite domain
	
	//............................................................................		
	// Initialise variables
	Array<double,2> f(N+1,1), neg_g(N+1,1), df(N+1,1), dgdf(N+1,N+1);
	
	//............................................................................	
	// Initialise domain
	this->domain_size = domain_size;
	
//	cheb->scale_for_semi_finite_domain(domain_size);

//	double half_grid_point = domain_size / 4.0;
//	cheb->scale_for_clustered_semi_finite_domain(domain_size, half_grid_point);
	
	cheb->scale_for_unclustered_semi_finite_domain(domain_size);
	
//	cheb->test_chebyshev_matrices_with_eta_power_4();
	for (long int ii=0; ii<=cheb->N; ii++) n[ii] = cheb->eta[ii];

	//initialise domain
	for(int ii=0; ii<=cheb->N; ii++) f(ii) = n[ii] - 1.21 * (1.0 - exp(-1.0*n[ii]));
	
	//............................................................................	
	// Iterate solution
	int iter=0;	
	double c = BIG_NUMBER;
	while ( ( (c>tol) && (iter<max_iter)) || (iter<min_iter) ) {
		iter++;
		neg_g = -1.0*evaluate_g(f, cheb);
		dgdf = evaluate_jacobian(f, cheb);
		df = solve_linear_system(neg_g, dgdf, N+1);
		c = 0.5*sum(neg_g*neg_g);
		f = f+df;
		fprintf(stdout,"Iteration number = %d. Error = %16.10e.\n", iter, c);
	}
	fprintf(stdout,"Iterations completed.\n");
	
	//............................................................................	
	// Reconstruct velocities
	reconstruct_velocities(cheb, f, u, u_x, u_y, v, v_x, v_y);
	
	//............................................................................	
	// Actual displacement thickness of Blasius Profile for non-dimensionalisation
	u_ref = 1.0;
	l_ref = 1.72078764;	
	
	//............................................................................	
}

//----------------------------------------------------------------------------
Array<double,2> evaluate_g(Array<double,2>& f, Chebyshev_t *cheb) 
{
	firstIndex i;	secondIndex j;	thirdIndex k;
	
	int num = cheb->N;
	
	Array<double,2> g(num+1,1);
	
	//Blasius ODE
	g = sum(cheb->D3(i,k)*f(k,j),k) + 0.5*sum(cheb->D2(i,k)*f(k,j),k)*f;
	
	//Apply BC f(eta=0)=0
	g(0,0)=f(0,0);
	
	//Apply BC df(eta=0)/d_eta=0
	g(1,0)=0.0;
	for(int ii=0;ii<=num;ii++) g(1,0)=g(1,0)+cheb->D1(0,ii)*f(ii,0);
	
	//Apply BC df(eta=N)/d_eta=1
	g(num,0)=0.0;
	for(int ii=0;ii<=num;ii++) g(num,0)=g(num,0)+cheb->D1(num,ii)*f(ii,0);
	g(num,0)=g(num,0)-1.0;			
	
	return g;
}

//----------------------------------------------------------------------------
Array<double,2> evaluate_jacobian(Array<double,2>& f, Chebyshev_t *cheb) 
{
	int num = cheb->N;
	
	Array<double,2> dgdf(num+1,num+1);			dgdf=0.0;
	for(int ii=2; ii<=num-1; ii++) {
		for(int jj=0; jj<=num; jj++) {
			dgdf(ii,jj) = cheb->D3(ii,jj) + 0.5 * f(ii) * cheb->D2(ii,jj);
		}
	}
	
	double sumTerm;
	for(int ii=0; ii<=num; ii++) {
		sumTerm=0.0;
		for(int kk=0; kk<=num; kk++) {
			sumTerm = sumTerm + cheb->D2(ii,kk) * f(kk);
		}
		dgdf(ii,ii) = cheb->D3(ii,ii) + 0.5 * f(ii) * cheb->D2(ii,ii) + 0.5 * sumTerm;
	}
	
	//Apply BC f(eta=0)=0
	Range all = Range::all();
	dgdf(0,all) = 0.0;
	dgdf(0,0) = 1.0;
	
	//Apply BC df(eta=0)/d_eta=0
	dgdf(1,all) = cheb->D1(0,all);
	
	//Apply BC df(eta=N)/d_eta=2
	dgdf(num,all) = cheb->D1(num,all);
	
	return dgdf;
}

//----------------------------------------------------------------------------
Array<double,2> solve_linear_system(Array<double,2>& Bb, Array<double,2>& Ab, int num)
{
	Array<double,2> xb(num,1);
	
	double *B = malloc_1d_array_double(num);
	double *A = malloc_1d_array_double(num*num);
	double *A_inv = malloc_1d_array_double(num*num);
	double *x = malloc_1d_array_double(num);
	
	// Copy to standard arrays
	for (int ii = 0; ii < num; ii++) {
		B[ii] = Bb(ii,0);
		for (int jj = 0; jj < num; jj++) {
			A[jj+ii*num] = Ab(ii,jj); 
		}
	}
	
	// Solve system
	matrix_invert(A, A_inv, num);
	matrix_multiply(A_inv, num, num, B, num, 1, x);
	
	// Copy to blitxz array
	for(int ii = 0; ii < num; ii++) xb(ii) = x[ii];
	
	free_1d_array_double(B);
	free_1d_array_double(A);
	free_1d_array_double(A_inv);
	free_1d_array_double(x);
	
	return xb;
}

//----------------------------------------------------------------------------
void test_solve_linear_system(Array<double,2>& b, Array<double,2>& A, int num) 
{
	/*	This functions tests the LU decomposition solution of gsl.
	 It has been checked in Matlab and is correct.
	 The output should be:
	 df(0) = 1
	 df(1:N) = 0	*/
	
	for (int ii = 0; ii < num; ii++) b(ii)=1 + ii;
	
	for (int ii = 0; ii < num; ii++) {
		for (int jj = 0; jj < num; jj++) {
			A(ii,jj)=abs(ii - jj) + 1; 
		}
	}
	
	Array<double,2> x(num,1);	
	x = solve_linear_system(b, A, num);
	cout << "x " << x << endl;
}

//----------------------------------------------------------------------------
void reconstruct_velocities(Chebyshev_t *cheb, Array<double,2>& f, 
							double *u, double *u_x, double *u_y, double *v, double *v_x, double *v_y) 
{
	// Reconstruct velocities
	Array<double,2> u_temp(cheb->N+1,1), u_y_temp(cheb->N+1,1);
	Array<double,2> v_temp(cheb->N+1,1), v_y_temp(cheb->N+1,1);
	firstIndex i;	secondIndex j;	thirdIndex k;
	
	u_temp = sum(cheb->D1(i,k)*f(k,j),k);
	u_y_temp = sum(cheb->D2(i,k)*f(k,j),k);
	for(int ii=0; ii<=cheb->N; ii++) v_temp(ii,0) = (u_temp(ii) * cheb->eta[ii] - f(ii))/2.0;
	v_y_temp = sum(cheb->D1(i,k)*v_temp(k,j),k);
	
	// Copy to standard matrices
	for(int ii=0; ii<=cheb->N; ii++) {
		u[ii] = u_temp(ii,0);
		u_x[ii] = 0.0;
		u_y[ii] = u_y_temp(ii,0);
		v[ii] = v_temp(ii,0);
		v_x[ii] = 0.0;
		v_y[ii] = v_y_temp(ii,0);
	}
}

//----------------------------------------------------------------------------
int Base_Flow_t::get_number_of_colocation_points_in_input_data()
{
	double x_read, y_read, xCart_read, yCart_read;
	double u_read, v_read, du_dx_read, du_dy_read, dv_dx_read, dv_dy_read;
	double d2u_dx2_read, d2u_dxdy_read, d2u_dy2_read, d2v_dx2_read, d2v_dxdy_read, d2v_dy2_read;
	double Ek_read, Ek_y_read, Ek_yy_read, nu_read;
	N = 0;	
	
	// Open file
	string file_name = base_flow_type;
	file_name = file_name + ".dat";
	ifstream fin(file_name.c_str());
	if (fin.fail()) quit_error((char*)"file %s not found\n", file_name.c_str());
	
	// Remove initial comment lines
	char next_char;
	char comment_char;		comment_char = '#';	
	next_char = fin.peek();
	while(next_char == comment_char) {
		fin.ignore(INT_MAX, '\n');
		next_char = fin.peek();
	}
	
	// Get data
	cout << "Reading file " << base_flow_type << endl;
	while(fin){
		fin >> x_read >> y_read >> xCart_read >> yCart_read >> u_read >> v_read >> du_dx_read >> du_dy_read >> dv_dx_read >> dv_dy_read >> d2u_dx2_read >> d2u_dxdy_read >> d2u_dy2_read >> d2v_dx2_read >> d2v_dxdy_read >> d2v_dy2_read >> Ek_read >> Ek_y_read >> Ek_yy_read >> nu_read;
		N++;
	}
	
	fin.close();
	
	return N-2;
}

//----------------------------------------------------------------------------
void Base_Flow_t::read_profile_from_file(Chebyshev_t *cheb)
{		
	double x_read, y_read, xCart_read, yCart_read;
	double u_read, v_read, du_dx_read, du_dy_read, dv_dx_read, dv_dy_read;
	double d2u_dx2_read, d2u_dxdy_read, d2u_dy2_read, d2v_dx2_read, d2v_dxdy_read, d2v_dy2_read;
	double Ek_read, Ek_y_read, Ek_yy_read, nu_read;
	
	// Open file
	string file_name = base_flow_type;
	file_name = file_name + ".dat";
	ifstream fin(file_name.c_str());
	if (fin.fail()) quit_error((char*)"file %s not found\n", file_name.c_str());
	
	// Remove initial comment lines
	char next_char;
	char comment_char;		comment_char = '#';	
	next_char = fin.peek();
	while(next_char == comment_char) {
		fin.ignore(INT_MAX, '\n');
		next_char = fin.peek();
	}
	
	// Get data
	cout << "Reading file " << file_name << endl;
	for (int j=0; j<=N; j++) {
		fin >> x_read >> y_read >> xCart_read >> yCart_read >> u_read >> v_read >> du_dx_read >> du_dy_read >> dv_dx_read >> dv_dy_read >> d2u_dx2_read >> d2u_dxdy_read >> d2u_dy2_read >> d2v_dx2_read >> d2v_dxdy_read >> d2v_dy2_read >> Ek_read >> Ek_y_read >> Ek_yy_read >> nu_read;
		
		n[j] = y_read;
		u[j] = u_read;				v[j] = v_read;
		u_x[j] = du_dx_read;		u_y[j] = du_dy_read;		
		v_x[j] = dv_dx_read;		v_y[j] = dv_dy_read;
		u_xx[j] = d2u_dx2_read;		u_xy[j] = d2u_dxdy_read;		u_yy[j] = d2u_dy2_read;
		v_xx[j] = d2v_dx2_read;		v_xy[j] = d2v_dxdy_read;		v_yy[j] = d2v_dy2_read;
		Ek[j] = Ek_read;			Ek_y[j] = Ek_y_read;			Ek_yy[j] = Ek_yy_read;
	}
	fin.close();
	
	u_ref = 1.0;	l_ref = 1.0;
//	non_dimensionalise_base_flow();
	
	x = x_read;
	this->domain_size = ABS(n[N] - n[0]);
	if ( n[0] * n[N] < 0.0) {	// if they are the different signs, must be wake or channel
		if ( ( ABS(n[0]) == 1.0) && ( ABS(n[N]) == 1.0) ) {
			cheb->scale_for_finite_domain(domain_size);					// Assuming a channel
			cout << "Assuming channel flow" << endl;
		} else {
			cheb->scale_for_unclustered_finite_domain(domain_size);		// Assuming a wake
			cout << "Assuming wake flow" << endl;			
		}
	} else {					// if they are the same sign, must be wall bounded
		cheb->scale_for_semi_finite_domain(domain_size);
		cout << "Assuming semi-bounded flow" << endl;
	}
//	for (int i=0; i<=N; i++) n[i] = cheb->eta[i];
}

//----------------------------------------------------------------------------
void Base_Flow_t::read_profile_eddy_viscosity_from_file(Chebyshev_t *cheb)
{		
	double x_read, y_read, xCart_read, yCart_read;
	double nuE1_read, nuE1_x_read, nuE1_y_read;
	double constE2_read, constE2_x_read, constE2_y_read;
	double constE3_read, constE3_x_read, constE3_y_read;
	double nu_sgs_avg_read, nu_sgs_rms_read, nu_read;
	
	// Open file
	string file_name = base_flow_type;
	file_name = file_name + ".eddy_viscosity.dat";
	ifstream fin(file_name.c_str());
	if (fin.fail()) quit_error((char*)"file %s not found\n", file_name.c_str());
	
	// Remove initial comment lines
	char next_char;
	char comment_char;		comment_char = '#';	
	next_char = fin.peek();
	while(next_char == comment_char) {
		fin.ignore(INT_MAX, '\n');
		next_char = fin.peek();
	}
	
	// Get data
	cout << "Reading file " << file_name << endl;

	for (int j=0; j<=N; j++) {
		fin >> x_read >> y_read >> xCart_read >> yCart_read >> nuE1_read >> nuE1_x_read >> nuE1_y_read >> constE2_read >> constE2_x_read >> constE2_y_read >> constE3_read >> constE3_x_read >> constE3_y_read >> nu_sgs_avg_read >> nu_sgs_rms_read >> nu_read;
		
		nuE1_original[j] =			nuE1_read;
		nuE1_y_original[j] =		nuE1_y_read;
		constE2_original[j] =		constE2_read;
		constE2_y_original[j] =		constE2_y_read;
		constE3_original[j] =		constE3_read;
		constE3_y_original[j] =		constE3_y_read;
	}
	
	for (int j=0; j<=N; j++) {
		nuE1[j] =		nuE1_original[j];
		nuE1_y[j] =		nuE1_y_original[j];
		constE2[j] =	constE2_original[j];
		constE2_y[j] =	constE2_y_original[j];
		constE3[j] =	constE3_original[j];
		constE3_y[j] =	constE3_y_original[j];
	}
	
	fin.close();
}

//----------------------------------------------------------------------------
void Base_Flow_t::read_profile_stresses_from_file(Chebyshev_t *cheb)
{		
	// Open file
	string file_name = base_flow_type;
	file_name = file_name + ".stresses.dat";
	ifstream fin(file_name.c_str());
	if (fin.fail()) quit_error((char*)"file %s not found\n", file_name.c_str());
	
	// Remove initial comment lines
	char next_char;
	char comment_char;		comment_char = '#';	
	next_char = fin.peek();
	while(next_char == comment_char) {
		fin.ignore(INT_MAX, '\n');
		next_char = fin.peek();
	}
		
	// Get data
	double x_read, y_read, xCart_read, yCart_read;
	double u_rms_read, v_rms_read, w_rms_read, uv_avg_read, uw_avg_read, vw_avg_read;
	double duu_dx_read, duu_dy_read, duv_dx_read, duv_dy_read, dvv_dx_read, dvv_dy_read, dww_dx_read, dww_dy_read;
	double I2_read, I3_read, eta_read, zeta_read, nu_read;
	double G, dG_dy;
	double F, dF_dy;
	
	// Need cut-offs to avoid singularity at centre-line of channel
	double nu_cut_off = 0.0;
	double K2_cut_off = epsilon;
	double K3_cut_off = epsilon;
	
	cout << "Reading file " << base_flow_type << endl;
	for (int j=0; j<=N; j++) {
		fin >> x_read >> y_read >> xCart_read >> yCart_read >> u_rms_read >> v_rms_read >> w_rms_read >> uv_avg_read >> uw_avg_read >> vw_avg_read >> duu_dx_read >> duu_dy_read >> duv_dx_read >> duv_dy_read >> dvv_dx_read >> dvv_dy_read >> dww_dx_read >> dww_dy_read >> I2_read >> I3_read >> eta_read >> zeta_read >> nu_read;
		
		u_rms[j] = u_rms_read;
		v_rms[j] = v_rms_read;
		w_rms[j] = w_rms_read;
		
		G = 4.0 * SQ(u_x[j]) + SQ(v_x[j]) + 2.0 * v_x[j] * u_y[j] + SQ(u_y[j]);
		dG_dy = 8.0 * u_x[j] * u_xy[j] + 2.0 * v_x[j] * v_xy[j] + 2.0 * v_xy[j] * u_y[j] + 2.0 * v_x[j] * u_yy[j] + 2.0 * u_y[j] * u_yy[j];
        if (G>0.0) {
            G = G + nu_cut_off;
        } else {
            G = G - nu_cut_off;
        }
		
		nuE1_original[j] = -( u_x[j] * (SQ(u_rms_read) - SQ(v_rms_read)) + uv_avg_read * ( u_y[j] + v_x[j] ) ) / G;
		nuE1_y_original[j] = -( u_xy[j] * (SQ(u_rms_read) - SQ(v_rms_read)) 
					  + u_x[j] * (duu_dy_read - dvv_dy_read) 
					  + duv_dy_read * (u_y[j] + v_x[j]) 
					  + uv_avg_read * (u_yy[j] + v_xy[j]) ) / G - nuE1_original[j] * dG_dy / G;
		
		F = G - nu_cut_off + K2_cut_off;
		dF_dy = dG_dy;
		constE2_original[j] = -( (u_y[j] + v_x[j]) * ( SQ(v_rms_read) - SQ(u_rms_read) ) + 4.0 * uv_avg_read * u_x[j] ) / F;
		constE2_y_original[j] = -( (u_yy[j] + v_xy[j]) * ( SQ(v_rms_read) - SQ(u_rms_read) ) 
						 + (u_y[j] + v_x[j]) * ( dvv_dy_read - duu_dy_read ) 
						 + 4.0 * duv_dy_read * u_x[j]
						 + 4.0 * uv_avg_read * u_xy[j]  ) / F
						- constE2_original[j] * dF_dy / F;
		constE2_original[j] = constE2_original[j] / 2.0;
		constE2_y_original[j] = constE2_y_original[j] / 2.0;
		
		G = G - nu_cut_off + K3_cut_off;
		constE3_original[j] =	- 2.0 * ( SQ(u_rms_read) + SQ(v_rms_read) - 2.0 * SQ(w_rms_read) ) * u_y[j] / G;
		constE3_y_original[j] =	- 2.0 * ( duu_dy_read + dvv_dy_read - 2.0 * dww_dy_read ) * u_y[j] / G
						- 2.0 * ( SQ(u_rms_read) + SQ(v_rms_read) - 2.0 * SQ(w_rms_read) ) * u_yy[j] / G
						- constE3_original[j] * dG_dy / G;
		constE3_original[j] = constE3_original[j] / 2.0;
		constE3_y_original[j] = constE3_y_original[j] / 2.0;

	}
	fin.close();
	
	for (int j=0; j<=N; j++) {
		nuE1[j] =		nuE1_original[j];
		nuE1_y[j] =		nuE1_y_original[j];
		constE2[j] =	constE2_original[j];
		constE2_y[j] =	constE2_y_original[j];
		constE3[j] =	constE3_original[j];
		constE3_y[j] =	constE3_y_original[j];
	}
}

//----------------------------------------------------------------------------
void Base_Flow_t::update_eddy_viscosity_fields(double dissipation_ratio)
{	
	double phi = 1.0 / ( 1.0 + dissipation_ratio );
	for (int j=0; j<=N; j++) {		
		nuE1[j] = nuE1_original[j] * phi;
		nuE1_y[j] = nuE1_y_original[j] * phi;
		constE2[j] = constE2_original[j] * SQ(phi);
		constE2_y[j] = constE2_y_original[j] * SQ(phi);
		constE3[j] = constE3_original[j] * SQ(phi);
		constE3_y[j] = constE3_y_original[j] * SQ(phi);
	}
}

//----------------------------------------------------------------------------
void Base_Flow_t::write_profile() 
{
	FILE *fout;
    char filename_out[] = "./results/profile.dat";
    int ii;
	fout=fopen(filename_out,"w");
	fprintf(fout,"#\n");
	fprintf(fout,"# x, n, u, u_y, u_yy, nu, nuE1, nuE1_y, constE2, constE2_y, constE3, constE3_y, Ek, Ek_y, Ek_yy\n");
	fprintf(fout,"# l_ref = %16.10e\n", l_ref);
	fprintf(fout,"# u_ref = %16.10e\n", u_ref);
	for (ii=0; ii<=N; ii++)
		fprintf(fout,"%17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e\n", 
				x, n[ii], u[ii], u_y[ii], u_yy[ii], nu, 
				nuE1_original[ii], nuE1_y_original[ii], constE2_original[ii], constE2_y_original[ii], constE3_original[ii], constE3_y_original[ii], 
				Ek[ii], Ek_y[ii], Ek_yy[ii]);
	fclose(fout);
}

//----------------------------------------------------------------------------
