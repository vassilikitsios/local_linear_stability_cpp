
#include <blitz/array.h>

using namespace blitz;

#include "Utils.h"
#include "Chebyshev_t.h"

//----------------------------------------------------------------------------
Chebyshev_t::Chebyshev_t()
{
}

//----------------------------------------------------------------------------
Chebyshev_t::~Chebyshev_t()
{
	free_1d_array_double(xi);
	free_1d_array_double(eta);
}

//----------------------------------------------------------------------------
void Chebyshev_t::allocate_memory(int N)
{
	this->N = N;
	xi = malloc_1d_array_double(N+1);
	eta = malloc_1d_array_double(N+1);
	allocateArrays(shape(N+1,N+1), D1, D2, D3, D4, w);
}

//----------------------------------------------------------------------------
void Chebyshev_t::generate_colocation_points()
{
	for (long int ii=0; ii<=N; ii++) {
		xi[ii] = cos(ii * M_PI / N);
		eta[ii] = xi[ii];
	}
}

//----------------------------------------------------------------------------
void Chebyshev_t::generate_chebyshev_integral_weights() 
{
	double c;
	if (N%2==0) {		// EVEN
		w(0,0) = 1.0 / (1.0 * SQ(N) - 1.0);
		for (int ii=1; ii<=N-1; ii++) {
			for (int kk=0; kk<=N/2; kk++) {
				c = 1.0;
				if (kk==0) c = 2.0;
				w(ii,ii) = w(ii,ii) + 4.0 / N / c * cos(2.0*PI*ii*kk/N) / (1.0-4.0*SQ(kk));
			}
		}
		w(N,N) = w(0,0);
	} else {			// ODD
		w(0,0) = 1.0 / (1.0 * SQ(N));
		for (int ii=1; ii<=N-1; ii++) {
			for (int kk=0; kk<=(N-1)/2; kk++) {
				c = 1.0;
				if ( (kk==0) || (kk==(N-1)/2) ) c = 2.0;
				w(ii,ii) = w(ii,ii) + 4.0 / N / c * cos(2.0*PI*ii*kk/N) / (1.0-4.0*SQ(kk))
				+ 4.0 / N * pow(-1.0,ii) * 2.0 * cos(PI*ii/N) / (1.0*SQ(N)*(2.0-N*1.0));
			}
		}
		w(N,N) = w(0,0);
	}
}
//----------------------------------------------------------------------------
void Chebyshev_t::generate_chebyshev_derivative_matrices() 
{
	generate_chebyshev_D1();
	generate_chebyshev_D2();
	generate_chebyshev_D3();
	generate_chebyshev_D4();
	generate_chebyshev_integral_weights();
}

//----------------------------------------------------------------------------
void Chebyshev_t::generate_chebyshev_D1() 
{
	Range all=Range::all();
	D1 = 0.0;
	
	for (int ii=0; ii<=N; ii++)
		for (int jj=0; jj<=N; jj++)
			D1(ii,jj)=pow(-1.0,ii+jj)/(xi[ii]-xi[jj]);
	
	D1(0,all)=2.0*D1(0,all);
	D1(N,all)=2.0*D1(N,all);
	D1(all,0)=0.5*D1(all,0);
	D1(all,N)=0.5*D1(all,N);	
	
	correct_diagonal_terms(D1);	
}

//----------------------------------------------------------------------------
void Chebyshev_t::generate_chebyshev_D2() 
{
	Range all=Range::all();
	
	D2=0.0;
	for(int ii=1;ii<=N-1;ii++)
		for(int jj=0;jj<=N;jj++)
			D2(ii,jj)=pow(-1.0,ii+jj)*(pow(xi[ii],2.0)+xi[ii]*xi[jj]-2) / ( (1-pow(xi[ii],2.0))*pow((xi[ii]-xi[jj]),2.0));

	D2(all,0)=D2(all,0)/2.0;
	D2(all,N)=D2(all,N)/2.0;
	
	for(int jj=1;jj<=N;jj++) D2(0,jj)=2.0/3.0*pow(-1.0,jj)*((2*N*N+1)*(1-xi[jj])-6)/pow(1-xi[jj],2.0);
	D2(0,N)=D2(0,N)/2.0;
	
	for(int jj=0;jj<=N-1;jj++) D2(N,jj)=2.0/3.0*pow(-1.0,jj+N)*((2*N*N+1)*(1+xi[jj])-6)/pow(1+xi[jj],2.0);
	D2(N,0)=D2(N,0)/2.0;
	
	correct_diagonal_terms(D2);
}

//----------------------------------------------------------------------------
void Chebyshev_t::generate_chebyshev_D3() {
	// Not completed need to determine the analytical form of the matrix elements
	firstIndex i;
	secondIndex j;
	thirdIndex k;
	Range all=Range::all();
	D3=0.0;
	D3=sum(D1(i,k)*D2(k,j),k);
	correct_diagonal_terms(D3);
}

//----------------------------------------------------------------------------
void Chebyshev_t::generate_chebyshev_D4() {
	// Not completed need to determine the analytical form of the matrix elements
	firstIndex i;
	secondIndex j;
	thirdIndex k;
	Range all=Range::all();
	D4=0.0;
	D4=sum(D2(i,k)*D2(k,j),k);
	correct_diagonal_terms(D4);
}

//----------------------------------------------------------------------------
void Chebyshev_t::correct_diagonal_terms(Array<double,2>& D) 
{
	int ii=0;
	Range all=Range::all();
	for(ii=0;ii<=N;ii++) {
		D(ii,ii)=0.0;
		D(ii,ii)=-sum(D(ii,all));
	}
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_for_semi_finite_domain(double domain_size) 
{
	scale_geometry_for_semi_finite_domain(domain_size);
	scale_derivatives_for_semi_finite_domain(domain_size);
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_for_finite_domain(double domain_size) 
{
	scale_geometry_for_finite_domain(domain_size);
	scale_derivatives_for_finite_domain(domain_size);
}


//----------------------------------------------------------------------------
void Chebyshev_t::scale_geometry_for_semi_finite_domain(double domain_size) 
{
	for (int ii=0; ii<=N; ii++) eta[ii] = domain_size/2.0*(1.0-xi[ii]);
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_geometry_for_finite_domain(double domain_size) 
{
	for (int ii=0; ii<=N; ii++) eta[ii] = -domain_size/2.0*xi[ii];
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_geometry_for_finite_domain(double y_min, double y_max)
{
	for (int ii=0; ii<=N; ii++) eta[ii] = (y_min - y_max)/2.0*xi[ii] + (y_max + y_min)/2.0;
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_derivatives_for_finite_domain(double domain_size) 
{
	scale_derivatives_for_semi_finite_domain(domain_size);
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_derivatives_for_semi_finite_domain(double domain_size) 
{
	D1=-2.0/domain_size*D1;
	D2=4.0/pow(domain_size,2)*D2;
	D3=-8.0/pow(domain_size,3)*D3;
	D4=16.0/pow(domain_size,4)*D4;
	w = w * domain_size/2.0;
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_for_unclustered_finite_domain(double domain_size) 
{
	scale_geometry_for_unclustered_finite_domain(domain_size);
	scale_derivatives_for_unclustered_finite_domain(domain_size);
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_geometry_for_unclustered_finite_domain(double domain_size) 
{
	double a = cos(8*PI/N);
	for (int ii=0; ii<=N; ii++) eta[ii] = -domain_size * asin(a * xi[ii]) / asin(a) / 2.0;
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_geometry_for_unclustered_finite_domain(double y_min, double y_max)
{
	double a = cos(8*PI/N);
	for (int ii=0; ii<=N; ii++) eta[ii] = (y_min - y_max) * asin(a * xi[ii]) / asin(a) / 2.0  + (y_max + y_min)/2.0;
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_derivatives_for_unclustered_finite_domain(double domain_size) 
{
	double a = cos(8*PI/N);
	double b = asin(a);
	double d_xi_d_eta, d2_xi_d_eta2;
	
	// scale D2
	for (int ii=0; ii<=N; ii++) {
		d_xi_d_eta = -2.0 * b / domain_size / a * cos(-b * 2.0 * eta[ii] / domain_size);
		d2_xi_d_eta2 = -SQ(2.0 * b / domain_size) / a * sin(-2.0 * b * eta[ii] / domain_size);
		for (int jj=0; jj<=N; jj++) D2(ii,jj) = D2(ii,jj) * SQ(d_xi_d_eta) + D1(ii,jj) * d2_xi_d_eta2;
	}
	correct_diagonal_terms(D2);
	
	// Scale D1
	for (int ii=0; ii<=N; ii++) {
		d_xi_d_eta = -2.0 * b / domain_size / a * cos(b*(-2.0 * eta[ii]/domain_size));
		for (int jj=0; jj<=N; jj++) D1(ii,jj) = D1(ii,jj) * d_xi_d_eta;
	}
	correct_diagonal_terms(D1);
	
	generate_chebyshev_D3();
	generate_chebyshev_D4();
	
	approximate_w();
}

//----------------------------------------------------------------------------
void Chebyshev_t::approximate_w() 
{
	for (int ii=1; ii<=N; ii++) w(ii,ii) = ABS(eta[ii] - eta[ii-1]);
	w(0,0) = w(1,1);	
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_for_unclustered_semi_finite_domain(double domain_size) 
{
	scale_geometry_for_unclustered_semi_finite_domain(domain_size);
	scale_derivatives_for_unclustered_semi_finite_domain(domain_size);
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_geometry_for_unclustered_semi_finite_domain(double domain_size) 
{
	double a = cos(8*PI/N);
	for (int ii=0; ii<=N; ii++) eta[ii] = domain_size * (1.0 - asin(a * xi[ii]) / asin(a)) / 2.0;
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_derivatives_for_unclustered_semi_finite_domain(double domain_size) 
{
	double a = cos(8*PI/N);
	double b = asin(a);
	double d_xi_d_eta, d2_xi_d_eta2;
	
	// scale D2
	for (int ii=0; ii<=N; ii++) {
		d_xi_d_eta = -2.0 * b / domain_size / a * cos(b*(1.0 - 2.0 * eta[ii]/domain_size));
		d2_xi_d_eta2 = -SQ(2.0 * b / domain_size) / a * sin(b*(1.0 - 2.0 * eta[ii]/domain_size));
		for (int jj=0; jj<=N; jj++) D2(ii,jj) = D2(ii,jj) * SQ(d_xi_d_eta) + D1(ii,jj) * d2_xi_d_eta2;
	}
	correct_diagonal_terms(D2);
	
	// Scale D1
	for (int ii=0; ii<=N; ii++) {
		d_xi_d_eta = -2.0 * b / domain_size / a * cos(b*(1.0 - 2.0 * eta[ii]/domain_size));
		for (int jj=0; jj<=N; jj++) D1(ii,jj) = D1(ii,jj) * d_xi_d_eta;
	}
	correct_diagonal_terms(D1);
	
	generate_chebyshev_D3();
	generate_chebyshev_D4();
	
	approximate_w();
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_for_clustered_semi_finite_domain(double domain_size, double half_grid_point) 
{
	scale_geometry_for_clustered_semi_finite_domain(domain_size, half_grid_point);
	scale_derivatives_for_clustered_semi_finite_domain();
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_geometry_for_clustered_semi_finite_domain(double domain_size, double half_grid_point) 
{
	double denom = 1.0 - 2.0 * half_grid_point / domain_size;
	if (ABS(denom)<1.0e-6) quit_error((char*)"cannot set half_grid_point at the mid point of the domain");
	l = half_grid_point / denom;
	s = 2.0 * l / domain_size;
	cout << "Chebyshev scaling of clustered semi infinite domain: s = " << s << ", l = " << l << endl;
	for (int ii=0; ii<=N; ii++) eta[ii] = l * (1.0 - xi[ii]) / (1.0 + s + xi[ii]);
	
	approximate_w();
}

//----------------------------------------------------------------------------
void Chebyshev_t::scale_derivatives_for_clustered_semi_finite_domain() 
{
	double d_xi_d_eta, d2_xi_d_eta2;

	// scale D2
	for (int ii=0; ii<=N; ii++) {
		d_xi_d_eta = -SQ(1.0 + s + xi[ii]) / l / (2.0 + s);
		d2_xi_d_eta2 = 2.0 * CU(1.0 + s + xi[ii]) / SQ(l) / SQ(2.0 + s);
		for (int jj=0; jj<=N; jj++) D2(ii,jj) = D2(ii,jj) * SQ(d_xi_d_eta) + D1(ii,jj) * d2_xi_d_eta2;
	}
	correct_diagonal_terms(D2);
	
	// Scale D1
	for (int ii=0; ii<=N; ii++) {
		d_xi_d_eta = -SQ(1.0 + s + xi[ii]) / l / (2.0 + s);
		for (int jj=0; jj<=N; jj++) D1(ii,jj) = D1(ii,jj) * d_xi_d_eta;
	}
	correct_diagonal_terms(D1);
	
	generate_chebyshev_D3();
	generate_chebyshev_D4();
	
	approximate_w();
}

//----------------------------------------------------------------------------
void Chebyshev_t::test_chebyshev_matrices_with_eta_power_4() 
{
	// Polynomial test function eta^4, such that the forth derivataive is a constant
	Array<double,2>	z, test, check, diff, dz_deta;
	allocateArrays(shape(N+1,1), z, test, check, diff, dz_deta);
	double diff_sum;
	
	firstIndex i;	secondIndex j;	thirdIndex k;
	
	for (int ii=0; ii<=N; ii++) z(ii,0) = pow(eta[ii],4.0);
	
	test=0.0;
	test=sum(D1(i,k)*z(k,j),k);
	for (int ii=0; ii<=N; ii++) check(ii,0) = 4.0*pow(eta[ii],3.0);
	diff = test-check;
	diff_sum = 0.0;
	for (int ii=0; ii<=N; ii++) diff_sum = diff_sum + diff(ii)*diff(ii);
	diff_sum = diff_sum / N;
	std::cout << "Average Squared Error for D1 is " << diff_sum << endl;		
	
	test=0.0;
	test=sum(D2(i,k)*z(k,j),k);
	for (int ii=0; ii<=N; ii++) check(ii,0) = 12.0*pow(eta[ii],2.0);
	diff = test-check;
	diff_sum = 0.0;
	for (int ii=0; ii<=N; ii++) diff_sum = diff_sum + diff(ii)*diff(ii);
	diff_sum = diff_sum / N;
	std::cout << "Average Squared Error for D2 is " << diff_sum << endl;	
	
	test=0.0;
	test=sum(D3(i,k)*z(k,j),k);
	for (int ii=0; ii<=N; ii++) check(ii,0) = 24.0*eta[ii];
	diff = test-check;
	diff_sum = 0.0;
	for (int ii=0; ii<=N; ii++) diff_sum = diff_sum + diff(ii)*diff(ii);
	diff_sum = diff_sum / N;
	std::cout << "Average Squared Error for D3 is " << diff_sum << endl;	
	
	test=0.0;
	test=sum(D4(i,k)*z(k,j),k);
	for (int ii=0; ii<=N; ii++) check(ii,0) = 24.0;
	diff = test-check;
	diff_sum = 0.0;
	for (int ii=0; ii<=N; ii++) diff_sum = diff_sum + diff(ii)*diff(ii);
	diff_sum = diff_sum / N;
	std::cout << "Average Squared Error for D4 is " << diff_sum << endl;
	
	double integral = 0.0;
	for (int ii=0; ii<=N; ii++) integral = integral + w(ii,ii) * z(ii,0);
	cout << "Integral W, analytical : " <<integral << " " << 0.2*pow(eta[0],5.0) - 0.2*pow(eta[N],5.0) << endl;
}

//----------------------------------------------------------------------------
