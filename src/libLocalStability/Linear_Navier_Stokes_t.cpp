
#include <cstdio>
#include <cstdlib>

#include <blitz/array.h>

using namespace blitz;

#include "Chebyshev_t.h"
#include "Base_Flow_t.h"
#include "Linear_Navier_Stokes_t.h"

#include "Utils.h"

//----------------------------------------------------------------------------
Linear_Navier_Stokes_t::Linear_Navier_Stokes_t()
{
	// helper variables
	real(im) = 0.0;			imag(im) = 1.0;
	real(one) = 1.0;		imag(one) = 0.0;
	real(unit) = 1.0;		imag(unit) = 1.0;
}

//----------------------------------------------------------------------------
Linear_Navier_Stokes_t::~Linear_Navier_Stokes_t() {
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::allocate_memory(int N_input_deck)
{
	int N_data = 0;
	if ( (strcmp(base_flow.base_flow_type,"couette")!=0) && (strcmp(base_flow.base_flow_type,"poiseuille")!=0) && 
		(strcmp(base_flow.base_flow_type,"shear_layer")!=0) && (strcmp(base_flow.base_flow_type,"wake")!=0) && 
		(strcmp(base_flow.base_flow_type,"blasius")!=0) && (strcmp(base_flow.base_flow_type,"channel_turb")!=0) ) {
		N_data = base_flow.get_number_of_colocation_points_in_input_data();
		cout << "N_data = " << N_data << endl;
		cout << "N_input_deck = " << N_input_deck << endl;
	}
	
	if (N_data > 0) {
		this->N = N_data;
	} else {
		this->N = N_input_deck;
	}
	cout << "N = " << N << endl;
	
	if ( (strcmp(equations,"laminar")!=0) 
		&& (strcmp(equations,"turbulent_non_linear_order1")!=0) 
		&& (strcmp(equations,"turbulent_non_linear_order2")!=0) 
		&& (strcmp(equations,"turbulent_non_linear_order3")!=0) )
		quit_error((char*)"Set equation type before allocating memory\n");
	
	N_size = 4*(N+1);
	allocateArrays( shape(N_size, N_size), A, B);
	eigen_solver.allocate_memory(N_size, N+1);	
	cheb.allocate_memory(N);
	base_flow.allocate_memory(N);
	
	with_non_parallel_terms = false;
	with_w_modes = true;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::set_problem_parameters(double Re, complex<double> kx, complex<double> kz, double dissipation_ratio_guess, double L_0)
{
	this->Re = Re;
	this->kx = kx;
	this->kz = kz;
	k2 = real(kx * conj(kx) + kz * conj(kz));
	k = sqrt(k2);
	
	if ( (strcmp(equations,"turbulent_non_linear_order1")==0) 
		|| (strcmp(equations,"turbulent_non_linear_order2")==0) 
		|| (strcmp(equations,"turbulent_non_linear_order3")==0) )
		base_flow.update_eddy_viscosity_fields(dissipation_ratio_guess);
	
	this->dissipation_ratio_guess = dissipation_ratio_guess;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::set_number_of_eigenvectors_to_write(int number_of_eigenvectors_to_write)
{
	this->number_of_eigenvectors_to_write = number_of_eigenvectors_to_write;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::set_boundary_conditions(char *top_boundary_condition, char *bottom_boundary_condition)
{
	if ( (strcmp(top_boundary_condition,"wall")!=0) && (strcmp(top_boundary_condition,"freestream")!=0) ) {
		quit_error((char*)"<top_boundary_condition> = 'wall' or 'freestream''\n");
	} else {
		this->top_boundary_condition = top_boundary_condition;
	}
	
	if ( (strcmp(bottom_boundary_condition,"wall")!=0) && (strcmp(bottom_boundary_condition,"freestream")!=0) ) {
		quit_error((char*)"<bottom_boundary_condition> = 'wall' or 'freestream''\n");
	} else {
		this->bottom_boundary_condition = bottom_boundary_condition;
	}
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::set_equations(char* equations)
{
	if ( (strcmp(equations,"laminar")!=0)  
		&& (strcmp(equations,"turbulent_non_linear_order1")!=0)
		&& (strcmp(equations,"turbulent_non_linear_order2")!=0)
		&& (strcmp(equations,"turbulent_non_linear_order3")!=0) ) {
		quit_error((char*)"<equation> = 'laminar' or 'turbulent_non_linear_order1' or 'turbulent_non_linear_order2' or 'turbulent_non_linear_order3'\n");
	} else {
		this->equations = equations;
		if (strcmp(equations,"laminar")==0) {
			base_flow.set_unsteady_bool(false);
		} else {
			base_flow.set_unsteady_bool(true);
		}
	}
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::set_with_non_parallel_terms(bool with_non_parallel_terms)
{
	this->with_non_parallel_terms = with_non_parallel_terms;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::generate_base_flow(double Re)
{
	this->Re = Re;
	cheb.generate_colocation_points();
	cheb.generate_chebyshev_derivative_matrices();
	base_flow.generate_base_flow(&cheb, 1.0/Re);
//	cheb.test_chebyshev_matrices_with_eta_power_4();
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::generate_linear_operator()
{
	// Note all matrices have been multiplied through by hp^3 * hq^3 to avoid any division in the matrices terms
	if (strcmp(equations,"laminar")==0) {
		cout << "	Building laminar operator..." << endl;
		build_laminar_matrix_A();
		build_laminar_matrix_B();
	} else if (strcmp(equations,"turbulent_non_linear_order1")==0) {
		build_laminar_matrix_A();
		build_turbulent_matrix_A_order1();
		build_turbulent_matrix_B();
	} else if (strcmp(equations,"turbulent_non_linear_order2")==0) {
		build_laminar_matrix_A();
		build_turbulent_matrix_A_order1();
		build_turbulent_matrix_A_order2();
		build_turbulent_matrix_B();
	} else if (strcmp(equations,"turbulent_non_linear_order3")==0) {
		build_laminar_matrix_A();
		build_turbulent_matrix_A_order1();
		build_turbulent_matrix_A_order2();
		build_turbulent_matrix_A_order3();
		build_turbulent_matrix_B();
	}
	apply_boundary_conditions();	
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::build_laminar_matrix_B()
{
	Array<double,2> I(N+1,N+1);
	I = 0.0;			for (int ii=0; ii<=N; ii++) I(ii,ii) = 1.0;
	
	B = 0.0;
	B( Range(0,N), Range(0,N) ) = I * one;
	B( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) = I * one;
	if (with_w_modes) {
		B( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) = I * one;
	} else {
		B( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) = 0.0;
	}
	B( Range(3*(N+1),4*(N+1)-1), Range(3*(N+1),4*(N+1)-1) ) = 0.0;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::build_laminar_matrix_A()
{
	Array<double,2> I(N+1,N+1), diag_u(N+1,N+1), diag_u_y(N+1,N+1);
	I = 0.0;	diag_u = 0.0;	diag_u_y = 0.0;
	for (int ii=0; ii<=N; ii++) {
		I(ii,ii) = 1.0;
		diag_u(ii,ii) = base_flow.u[ii];
		diag_u_y(ii,ii) = base_flow.u_y[ii];
	}
	
	Array<complex<double>,2> Z(N+1,N+1);
	Z = diag_u * kx + im*1.0/Re * (cheb.D2 - k2 * I);
	
	A=0.0;
	A( Range(0,N), Range(0,N) ) = Z;
	A( Range(0,N), Range(N+1,2*(N+1)-1) ) = -im*diag_u_y;
	A( Range(0,N), Range(2*(N+1),3*(N+1)-1) ) = 0.0;
	A( Range(0,N), Range(3*(N+1),4*(N+1)-1) ) = kx * I;
	
	A( Range(N+1,2*(N+1)-1), Range(0,N) ) = 0.0;
	A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) = Z;
	A( Range(N+1,2*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) = 0.0;
	A( Range(N+1,2*(N+1)-1), Range(3*(N+1),4*(N+1)-1) ) = -im*cheb.D1;
	
	A( Range(2*(N+1),3*(N+1)-1), Range(0,N) ) = 0.0;
	A( Range(2*(N+1),3*(N+1)-1), Range(N+1,2*(N+1)-1) ) = 0.0;
	if (with_w_modes) {
		A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) = Z;
		A( Range(2*(N+1),3*(N+1)-1), Range(3*(N+1),4*(N+1)-1) ) = kz * I;
	} else {
		A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) = I * one;
		A( Range(2*(N+1),3*(N+1)-1), Range(3*(N+1),4*(N+1)-1) ) = I * one;
	}
	
	A( Range(3*(N+1),4*(N+1)-1), Range(0,N) ) = kx * I;
	A( Range(3*(N+1),4*(N+1)-1), Range(N+1,2*(N+1)-1) ) = -im * cheb.D1;
	A( Range(3*(N+1),4*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) = kz * I;
	A( Range(3*(N+1),4*(N+1)-1), Range(3*(N+1),4*(N+1)-1) ) = 0.0;
	
	if (with_non_parallel_terms) add_non_parallel_terms();
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::add_non_parallel_terms()
{
	Array<double,2> diag_v(N+1,N+1), diag_v_y(N+1,N+1), diag_v_D(N+1,N+1);
	diag_v = 0.0;	diag_v_y = 0.0;		diag_v_D = 0.0;
	for (int ii=0; ii<=N; ii++) {
		diag_v(ii,ii) = base_flow.v[ii];
		diag_v_y(ii,ii) = base_flow.v_y[ii];
		for (int jj=0; jj<=N; jj++) diag_v_D(ii,jj) = base_flow.v_y[ii] * cheb.D1(ii,jj);
	}
	
	A( Range(0,N), Range(0,N) ) = A( Range(0,N), Range(0,N) ) + diag_v_D * one;
	
	A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) = A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) + (diag_v_D + diag_v_y) * one;
	
	A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) = A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) + diag_v_D * one;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::build_turbulent_matrix_B()
{
	build_laminar_matrix_B();
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::build_turbulent_matrix_A_order1()
{
	// Additional turbulent terms
	Array<double,2> diag_nuE_y(N+1,N+1), diag_nuE_y_D(N+1,N+1);
	Array<complex<double>,2> Ze(N+1,N+1);
	
	diag_nuE_y = 0.0;	diag_nuE_y_D = 0.0;	Ze = 0.0;
	for (int ii=0; ii<=N; ii++) {
		diag_nuE_y(ii,ii) = base_flow.nuE1_y[ii];
		for (int jj=0; jj<=N; jj++) {
			diag_nuE_y_D(ii,jj) = base_flow.nuE1_y[ii] * cheb.D1(ii,jj);
			if (ii==jj) {
				Ze(ii,jj) = im * base_flow.nuE1[ii] * (cheb.D2(ii,jj) - k2) + diag_nuE_y_D(ii,jj) * im;
			} else {
				Ze(ii,jj) = im * base_flow.nuE1[ii] * cheb.D2(ii,jj) + diag_nuE_y_D(ii,jj) * im;
			}
		}
	}								
	
	A( Range(0,N), Range(0,N) ) = A( Range(0,N), Range(0,N) ) + Ze;
	A( Range(0,N), Range(N+1,2*(N+1)-1) ) = A( Range(0,N), Range(N+1,2*(N+1)-1) ) - kx * diag_nuE_y;
	
	A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) = A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) + Ze + im * diag_nuE_y_D;
	
	A( Range(2*(N+1),3*(N+1)-1), Range(N+1,2*(N+1)-1) ) =  A( Range(2*(N+1),3*(N+1)-1), Range(N+1,2*(N+1)-1) ) - kz * diag_nuE_y;
	A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) = A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) + Ze;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::build_turbulent_matrix_A_order2()
{
	double Z_nu, Z_U;
	Array<complex<double>,2> Z_xu(N+1,N+1), Z_xv(N+1,N+1), Z_yu(N+1,N+1), Z_yv_zw(N+1,N+1), Z_zu(N+1,N+1);
	Z_xu = 0.0;		Z_xv = 0.0;		Z_yu = 0.0;		Z_yv_zw = 0.0;		Z_zu = 0.0;

	for (int ii=0; ii<=N; ii++) {
		Z_nu = base_flow.constE2[ii];
		Z_U = base_flow.constE2_y[ii];
		for (int jj=0; jj<=N; jj++) {
			Z_xv(ii,jj) = -Z_U * cheb.D1(ii,jj);
			Z_yv_zw(ii,jj) = -Z_nu * cheb.D1(ii,jj) * kx * im;
			if (ii==jj) {
				Z_xu(ii,jj) = Z_U * kx * im;
				Z_yu(ii,jj) = 2.0 * Z_U * cheb.D1(ii,jj) + Z_nu * ( 2.0 * cheb.D2(ii,jj) - k2 );
				Z_zu(ii,jj) = ( Z_U + Z_nu * cheb.D1(ii,jj) ) * kz * im;
			} else {
				Z_xu(ii,jj) = 0.0;
				Z_yu(ii,jj) =  2.0 * Z_U * cheb.D1(ii,jj) + Z_nu * 2.0 * cheb.D2(ii,jj);
				Z_zu(ii,jj) = Z_nu * cheb.D1(ii,jj) * kz * im;
			}
		}
	}								
	
	A( Range(0,N), Range(0,N) ) =								A( Range(0,N), Range(0,N) ) + Z_xu * im;
	A( Range(0,N), Range(N+1,2*(N+1)-1) ) =						A( Range(0,N), Range(N+1,2*(N+1)-1) ) + Z_xv * im;
	
	A( Range(N+1,2*(N+1)-1), Range(0,N) ) =						A( Range(N+1,2*(N+1)-1), Range(0,N) ) + Z_yu * im;
	A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) =			A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) + Z_yv_zw * im;
	
	A( Range(2*(N+1),3*(N+1)-1), Range(0,N) ) =					A( Range(2*(N+1),3*(N+1)-1), Range(0,N) ) + Z_zu * im;
	A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) =	A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) + Z_yv_zw * im;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::build_turbulent_matrix_A_order3()
{
	double Z_nu, Z_U;
	Array<complex<double>,2> Z_xu(N+1,N+1), Z_xv(N+1,N+1), Z_yu(N+1,N+1), Z_yv(N+1,N+1), Z_zu(N+1,N+1), Z_zv(N+1,N+1), Z_zw(N+1,N+1);
	Z_xu = 0.0; Z_xv = 0.0; Z_yu = 0.0; Z_yv = 0.0; Z_zu = 0.0; Z_zv = 0.0; Z_zw = 0.0;

	for (int ii=0; ii<=N; ii++) {
		Z_nu = base_flow.constE3[ii];
		Z_U = base_flow.constE3_y[ii];
		for (int jj=0; jj<=N; jj++) {
			if (ii==jj) {
				Z_xu(ii,jj) = ( 6.0 * Z_U + 5.0 * Z_nu * cheb.D1(ii,jj) ) * kx * im;
				Z_xv(ii,jj) = 6.0 * Z_U * cheb.D1(ii,jj) + Z_nu * ( 3.0 * cheb.D2(ii,jj) - 3.0 * k2 + kx * kx );
				Z_yu(ii,jj) = 2.0 * Z_U * cheb.D1(ii,jj) + Z_nu * ( 2.0 * cheb.D2(ii,jj) - 3.0 * k2 );
				Z_yv(ii,jj) = ( 2.0 * Z_U + 5.0 * Z_nu * cheb.D1(ii,jj) ) * kx * im;
				Z_zu(ii,jj) = ( 3.0 * Z_U - Z_nu * cheb.D1(ii,jj) ) * kz * im;
				Z_zv(ii,jj) = Z_nu * kx * kz;
				Z_zw(ii,jj) = ( 3.0 * Z_U + 6.0 * Z_nu * cheb.D1(ii,jj) ) * kx * im;
			} else {
				Z_xu(ii,jj) = 5.0 * Z_nu * cheb.D1(ii,jj) * kx * im;
				Z_xv(ii,jj) = 6.0 * Z_U * cheb.D1(ii,jj) + Z_nu * 3.0 * cheb.D2(ii,jj);
				Z_yu(ii,jj) = 2.0 * Z_U * cheb.D1(ii,jj) + Z_nu * 2.0 * cheb.D2(ii,jj);
				Z_yv(ii,jj) = 5.0 * Z_nu * cheb.D1(ii,jj) * kx * im;
				Z_zu(ii,jj) = -Z_nu * cheb.D1(ii,jj) * kz * im;
				Z_zv(ii,jj) = 0.0;
				Z_zw(ii,jj) = 6.0 * Z_nu * cheb.D1(ii,jj) * kx * im;
			}
		}
	}								
	
	A( Range(0,N), Range(0,N) ) =								A( Range(0,N), Range(0,N) ) + Z_xu * im / 6.0;
	A( Range(0,N), Range(N+1,2*(N+1)-1) ) =						A( Range(0,N), Range(N+1,2*(N+1)-1) ) + Z_xv * im / 6.0;

	A( Range(N+1,2*(N+1)-1), Range(0,N) ) =						A( Range(N+1,2*(N+1)-1), Range(0,N) ) + Z_yu * im / 6.0;
	A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) =			A( Range(N+1,2*(N+1)-1), Range(N+1,2*(N+1)-1) ) + Z_yv * im / 6.0;
	
	A( Range(2*(N+1),3*(N+1)-1), Range(0,N) ) =					A( Range(2*(N+1),3*(N+1)-1), Range(0,N) ) + Z_zu * im / 6.0;
	A( Range(2*(N+1),3*(N+1)-1), Range(N+1,2*(N+1)-1) ) =		A( Range(2*(N+1),3*(N+1)-1), Range(N+1,2*(N+1)-1) ) + Z_zv * im / 6.0;
	A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) =	A( Range(2*(N+1),3*(N+1)-1), Range(2*(N+1),3*(N+1)-1) ) + Z_zw * im / 6.0;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::apply_boundary_conditions()
{	
	if (strcmp(top_boundary_condition,"wall")==0) {
		cout << "			apply_boundary_conditions_for_top_wall" << endl; 
		apply_boundary_conditions_for_top_wall();
	} else if (strcmp(top_boundary_condition,"freestream")==0) {
		cout << "			apply_boundary_conditions_for_top_freestream" << endl; 
		apply_boundary_conditions_for_top_freestream();
	}
	
	if (strcmp(bottom_boundary_condition,"wall")==0) {
		cout << "			apply_boundary_conditions_for_bottom_wall" << endl; 
		apply_boundary_conditions_for_bottom_wall();
	} else if (strcmp(bottom_boundary_condition,"freestream")==0) {
		cout << "			apply_boundary_conditions_for_bottom_freestream" << endl; 
		apply_boundary_conditions_for_bottom_freestream();
	}
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::apply_boundary_conditions_for_bottom_wall()
{
	// BC: real[u(xi=1)=0] + imag[u(xi=1)=0]						wall
	A(0, Range(0,4*(N+1)-1)) = 0.0;						A(0, 0) = one;
	B(0, Range(0,4*(N+1)-1)) = 0.0;		
	
	// BC: real[v(xi=1)=0] + imag[v(xi=1)=0]						wall
	A(N+1, Range(0,4*(N+1)-1)) = 0.0;					A(N+1, N+1) = one;
	B(N+1, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[w(xi=1)=0] + imag[w(xi=1)=0]						wall
	A(2*(N+1), Range(0,4*(N+1)-1)) = 0.0;				A(2*(N+1), 2*(N+1)) = one;
	B(2*(N+1), Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[dpi_dxj(xi=1)] + imag [dpi_dxj(xi=1)] = 0			wall
	A(3*(N+1), Range(0,4*(N+1)-1)) = 0.0;				A(3*(N+1), Range(3*(N+1),4*(N+1)-1) ) = cheb.D1(0,Range(0,N)) * one;
	//	A(3*(N+1), 3*(N+1)) = A(3*(N+1), 3*(N+1)) + (im*kx + im*kz) * one;
	B(3*(N+1), Range(0,4*(N+1)-1)) = 0.0;
	
	// Shift eigenvalues
	B(0, Range(0,4*(N+1)-1)) = A(0, Range(0,4*(N+1)-1)) * im / artifical_eigenvalue_shift;
	B(N+1, Range(0,4*(N+1)-1)) = A(N+1, Range(0,4*(N+1)-1)) * im / artifical_eigenvalue_shift;
	B(2*(N+1), Range(0,4*(N+1)-1)) = A(2*(N+1), Range(0,4*(N+1)-1)) * im / artifical_eigenvalue_shift;
	B(3*(N+1), Range(0,4*(N+1)-1)) = A(2*(N+1), Range(0,4*(N+1)-1)) * im / artifical_eigenvalue_shift;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::apply_boundary_conditions_for_top_wall()
{
	// BC: real[u(xi=-1)=0] + imag[u(xi=-1)=0]						wall
	A(N, Range(0,4*(N+1)-1)) = 0.0;						A(N, N) = one;
	B(N, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[v(xi=-1)=0] + imag[v(xi=-1)=0]						wall
	A(2*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;				A(2*(N+1)-1, 2*(N+1)-1) = one;
	B(2*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[w(xi=-1)=0] + imag[w(xi=-1)=0]						wall
	A(3*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;				A(3*(N+1)-1, 3*(N+1)-1) = one;
	B(3*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[dpi_dxj(xi=-1)] + imag [dpi_dxj(xi=-1)] = 0			wall
	A(4*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;				A(4*(N+1)-1, Range(3*(N+1),4*(N+1)-1) ) = cheb.D1(N,Range(0,N)) * one;	
	//	A(4*(N+1)-1, 4*(N+1)-1) = A(4*(N+1)-1, 4*(N+1)-1) + (im*kx + im*kz) * one;
	B(4*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;
	
	// Shift eigenvalues
	B(N, Range(0,4*(N+1)-1)) = A(N, Range(0,4*(N+1)-1)) * im / artifical_eigenvalue_shift;
	B(2*(N+1)-1, Range(0,4*(N+1)-1)) = A(2*(N+1)-1, Range(0,4*(N+1)-1)) * im / artifical_eigenvalue_shift;
	B(3*(N+1)-1, Range(0,4*(N+1)-1)) = A(3*(N+1)-1, Range(0,4*(N+1)-1)) * im / artifical_eigenvalue_shift;
	B(4*(N+1)-1, Range(0,4*(N+1)-1)) = A(4*(N+1)-1, Range(0,4*(N+1)-1)) * im / artifical_eigenvalue_shift;		
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::apply_boundary_conditions_for_top_freestream()
{
	// Note sign of 'k' is changed to ensure the BC decay's to zero ?????
	
	// BC: real[du_dy(xi=-1) = -k u(xi=-1)] + imag[du_dy(xi=-1) = -k u(xi=-1)]				free stream
	A(N, Range(0,4*(N+1)-1)) = 0.0;					A(N, Range(0,N)) = cheb.D1(N,Range(0,N)) * one;
	A(N, N) = A(N, N) + k * one;
	B(N, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[dv_dy(xi=-1) = -k u(xi=-1)] + imag[dv_dy(xi=-1) = -k u(xi=-1)]				free stream
	A(2*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;				A(2*(N+1)-1, Range(N+1,2*(N+1)-1) ) = cheb.D1(N,Range(0,N)) * one;
	A(2*(N+1)-1, 2*(N+1)-1) = A(2*(N+1)-1, 2*(N+1)-1) + k * one;
	B(2*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[dw_dy(xi=-1) = -k u(xi=-1)] + imag[dw_dy(xi=-1) = -k u(xi=-1)]				free stream
	A(3*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;				A(3*(N+1)-1, Range(2*(N+1),3*(N+1)-1) ) = cheb.D1(N,Range(0,N)) * one;
	A(3*(N+1)-1, 3*(N+1)-1) = A(3*(N+1)-1, 3*(N+1)-1) + k * one;
	B(3*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[d2pi_dxj2(xi=-1)] + imag [d2pi_dxj2(xi=-1)] = 0								freestream	
	A(4*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;				A(4*(N+1)-1, Range(3*(N+1),4*(N+1)-1) ) = cheb.D2(N,Range(0,N)) * one;	
	A(4*(N+1)-1, 4*(N+1)-1) = A(4*(N+1)-1, 4*(N+1)-1) - k2 * one;
	B(4*(N+1)-1, Range(0,4*(N+1)-1)) = 0.0;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::apply_boundary_conditions_for_bottom_freestream()
{	
	// BC: real[du_dy(xi=1) = -k u(xi=1)]													free stream
	A(0, Range(0,4*(N+1)-1)) = 0.0;					A(0, Range(0,N)) = cheb.D1(0,Range(0,N)) * one;
	A(0, 0) = A(0, 0) + k * one;
	B(0, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[dv_dy(xi=1) = -k u(xi=1)]													free stream
	A(N+1, Range(0,4*(N+1)-1)) = 0.0;				A(N+1, Range(N+1,2*(N+1)-1) ) = cheb.D1(0,Range(0,N)) * one;
	A(N+1, N+1) = A(N+1, N+1) + k * one;
	B(N+1, Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[dw_dy(xi=1) = -k u(xi=1)]													free stream
	A(2*(N+1), Range(0,4*(N+1)-1)) = 0.0;				A(2*(N+1), Range(2*(N+1),3*(N+1)-1) ) = cheb.D1(0,Range(0,N)) * one;
	A(2*(N+1), 2*(N+1)) = A(2*(N+1), 2*(N+1)) + k * one;
	B(2*(N+1), Range(0,4*(N+1)-1)) = 0.0;
	
	// BC: real[d2pi_dxj2(xi=1)] = 0														freestream
	A(3*(N+1), Range(0,4*(N+1)-1)) = 0.0;				A(3*(N+1), Range(3*(N+1),4*(N+1)-1) ) = cheb.D2(0,Range(0,N)) * one;	
	A(3*(N+1), 3*(N+1)) = A(3*(N+1), 3*(N+1)) - k2 * one;
	B(3*(N+1), Range(0,4*(N+1)-1)) = 0.0;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::calculate_eigensolution() 
{	
	eigen_solver.calculate_eigensolution(A, B);
	qsort(eigen_solver.eigen_solution, N_size, sizeof(Eigen_Solution_t), (int (*)(const void*, const void*)) qsort_compare_function);
}

//----------------------------------------------------------------------------
int qsort_compare_function(Eigen_Solution_t *es1, Eigen_Solution_t *es2)
{
	if(imag(es1->eigenvalue) > imag(es2->eigenvalue)) return -1;
	if(imag(es1->eigenvalue) < imag(es2->eigenvalue)) return 1;
	return 0;
}


//----------------------------------------------------------------------------
// Specific search code
//----------------------------------------------------------------------------
complex<double> Linear_Navier_Stokes_t::get_eigensolution_closest_to_this_omega(complex<double> this_omega, int *ptr_eigenvalue_number)
{
	double error, min_error = BIG_NUMBER*1.0;
	int closest_mode_index = 0;
	for (int ii=0; ii<N_size; ii++) {
		error = SQ(real(eigen_solver.eigen_solution[ii].eigenvalue) - real(this_omega)) + SQ(imag(eigen_solver.eigen_solution[ii].eigenvalue) - imag(this_omega));
		if (error < min_error) {
			min_error = error;
			closest_mode_index = ii;
		}
	}
	*ptr_eigenvalue_number = closest_mode_index+1;
	return eigen_solver.eigen_solution[closest_mode_index].eigenvalue;
}


//----------------------------------------------------------------------------
// I/O helper functions
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::write_output(int kx_real_num, int kx_imag_num)
{
	write_output(kx_real_num, kx_imag_num, -1);
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::write_output(int kx_real_num, int kx_imag_num, int iteration)
{
	write_eigenvalues(kx_real_num, kx_imag_num, iteration);
	write_eigenvectors(kx_real_num, kx_imag_num, iteration);
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::write_eigenvalues(int kx_real_num, int kx_imag_num, int iteration) 
{
	double c_r, c_i, abs_kx;
	string filename_out;
	filename_out = get_eigenvalue_filename(kx_real_num, kx_imag_num, iteration);
	
	FILE *fout;				fout=fopen(filename_out.c_str(),"a");
	fprintf(fout,"#\n");
	fprintf(fout,"# equations = %s\n", equations);
	fprintf(fout,"# base flow = %s\n", base_flow.base_flow_type);
	fprintf(fout,"# Re = \"%8.6f\"\n", Re);
	fprintf(fout,"# domain_size = %16.10e\n", ABS(base_flow.n[N] - base_flow.n[0]));	
	fprintf(fout,"# N = \"%d\"\n", N);
	fprintf(fout,"# N_eigenvalues = \"%d\"\n", N_size);
	fprintf(fout,"#\n");
	fprintf(fout,"# num, kx_r, kx_i, kz_r, kz_i, a_r, a_i, b_r, b_i, Omega_r, Omega_i, c_r, c_i, discrete=0 continuous=1, u_mode=0 w_mode=1, asymmetic=0 symmetic=1,  rel_error_bounds, condition number\n");
	fprintf(fout,"#\n");
	
	for (int ii=0; ii<N_size; ii++) {
		abs_kx = SQ(real(kx)) + SQ(imag(kx));
		if ( abs_kx > 0 ) {
			c_r = real(eigen_solver.eigen_solution[ii].eigenvalue) * real(kx) + imag(eigen_solver.eigen_solution[ii].eigenvalue) * imag(kx);
			c_r = c_r / abs_kx;
			c_i = imag(eigen_solver.eigen_solution[ii].eigenvalue) * real(kx) - real(eigen_solver.eigen_solution[ii].eigenvalue) * imag(kx);
			c_i = c_i / abs_kx;
		} else {
			c_r = 0.0;
			c_i = 0.0;
		}
		fprintf(fout,"%d %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %d %d %d %16.10e %16.10e \n",
				ii+1,
				real(kx), imag(kx), real(kz), imag(kz),
				real(eigen_solver.eigen_solution[ii].alpha), imag(eigen_solver.eigen_solution[ii].alpha),
				real(eigen_solver.eigen_solution[ii].beta), imag(eigen_solver.eigen_solution[ii].beta),
				real(eigen_solver.eigen_solution[ii].eigenvalue), imag(eigen_solver.eigen_solution[ii].eigenvalue), c_r, c_i,
				eigen_solver.eigen_solution[ii].discrete_continuous,
				eigen_solver.eigen_solution[ii].w_mode, eigen_solver.eigen_solution[ii].symmetric,
				eigen_solver.eigen_solution[ii].rel_error_bound, eigen_solver.eigen_solution[ii].cond);
	}
	fclose(fout);
	
	// Write filename to file
	string eigenvalue_files_filename;
	if (iteration==-1) {
		eigenvalue_files_filename = "./results/eigenvalue_files.list";
	} else {
		eigenvalue_files_filename = "./results.iter/eigenvalue_files.list";
	}
	fout=fopen(eigenvalue_files_filename.c_str(),"a");
	fprintf(fout,"%s\n", filename_out.c_str());
	fclose(fout);
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::write_eigenvectors(int kx_real_num, int kx_imag_num, int iteration) 
{	
	string filename_out;
	if (number_of_eigenvectors_to_write>2*N) number_of_eigenvectors_to_write = 2*N;
	for (int eigen_vector=0; eigen_vector<number_of_eigenvectors_to_write; eigen_vector++) {
		filename_out = get_eigenvector_filename(eigen_vector, kx_real_num, kx_imag_num, iteration);
		write_eigenvector_header(filename_out, eigen_vector);
		write_eigenvector(filename_out, eigen_vector);
	}
}

//----------------------------------------------------------------------------
string Linear_Navier_Stokes_t::get_eigenvalue_filename(int kx_real_num, int kx_imag_num, int iteration) 
{	
	string eigenvalue_filename;
	string kx_real_number_string;
	string kx_imag_number_string;
	string iteration_number_string;
	kx_real_number_string = get_padded_number_string(kx_real_num+1);
	kx_imag_number_string = get_padded_number_string(kx_imag_num+1);
	iteration_number_string = get_padded_number_string(iteration+1);
	
	if (iteration==-1) {
		eigenvalue_filename = "./results/kx_" + kx_real_number_string + "_" + kx_imag_number_string + ".eigen_values.dat";
	} else {
		eigenvalue_filename = "./results.iter/iter_" + iteration_number_string + ".kx_" + kx_real_number_string + "_" + kx_imag_number_string + ".eigen_values.dat";
	}
	
	return eigenvalue_filename;
}

//----------------------------------------------------------------------------
string Linear_Navier_Stokes_t::get_eigenvector_filename(int eigenvector_num, int kx_real_num, int kx_imag_num, int iteration) 
{	
	string eigenvector_filename;
	string mode_number_string;
	string kx_real_number_string;
	string kx_imag_number_string;
	string iteration_number_string;
	
	mode_number_string = get_padded_number_string(eigenvector_num+1);
	kx_real_number_string = get_padded_number_string(kx_real_num+1);
	kx_imag_number_string = get_padded_number_string(kx_imag_num+1);
	iteration_number_string = get_padded_number_string(iteration+1);
	
	if (iteration==-1) {
		eigenvector_filename = "./results/kx_" + kx_real_number_string + "_" + kx_imag_number_string + ".eigen_vector.mode_" + mode_number_string + ".dat";
	} else {
		eigenvector_filename = "./results.iter/iter_" + iteration_number_string + ".kx_" + kx_real_number_string + "_" + kx_imag_number_string + ".eigen_vector.mode_" + mode_number_string + ".dat";
	}
	
	return eigenvector_filename;
}

//----------------------------------------------------------------------------
string Linear_Navier_Stokes_t::get_padded_number_string(int number)
{
	string number_string, pad_zeros_string, padded_number_string;
	number_string = convert_int_to_string(number);
	if (number+1>9999) {
		quit_error((char*)"Too many modes to output, increase digits in output file names");
	} else if ( (number<9999) && (number>=1000) ) {
		pad_zeros_string = "";
	} else if ( (number<1000) && (number>=100) ) {
		pad_zeros_string = "0";
	} else if ( (number<100) && (number>=10) ) {
		pad_zeros_string = "00";
	} else if (number<10) {
		pad_zeros_string = "000";
	}
	
	padded_number_string = pad_zeros_string + number_string;
	return padded_number_string;
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::write_eigenvector_header(string filename_out, int eigenvector_num) 
{	
	FILE *fout=fopen(filename_out.c_str(),"w");
	
	fprintf(fout,"#\n");
	fprintf(fout,"# equations = %s\n", equations);
	fprintf(fout,"# base flow = %s\n", base_flow.base_flow_type);
	fprintf(fout,"# Re = \"%8.6f\"\n", Re);
	fprintf(fout,"# kx = \"(%8.6f, %8.6f i)\"\n", real(kx), imag(kx));
	fprintf(fout,"# kz = \"(%8.6f, %8.6f i)\"\n", real(kz), imag(kz));
	fprintf(fout,"# domain size = %16.10e\n", ABS(base_flow.n[N] - base_flow.n[0]));	
	fprintf(fout,"# N = \"%d\"\n", N);
	fprintf(fout,"# eigenvalue: Omega = \"(%16.10e, %16.10e i)\"\n", 
			real(eigen_solver.eigen_solution[eigenvector_num].eigenvalue), 
			imag(eigen_solver.eigen_solution[eigenvector_num].eigenvalue));
	if ( abs(real(kx))>0 ) {
		fprintf(fout,"# wave speed: c = \"(%16.10e, %16.10e i)\"\n", 
				real(eigen_solver.eigen_solution[eigenvector_num].eigenvalue)/real(kx), 
				imag(eigen_solver.eigen_solution[eigenvector_num].eigenvalue)/real(kx));
	}
	if (eigen_solver.eigen_solution[eigenvector_num].discrete_continuous==0) {
		fprintf(fout,"# discrete mode\n");
	} else if (eigen_solver.eigen_solution[eigenvector_num].discrete_continuous==1) {
		fprintf(fout,"# continuous mode\n");
	} else {
		fprintf(fout,"# discrete / continuous mode classificaiton failed\n");
	}
	fprintf(fout,"# relative error bounds = %16.10e\n", eigen_solver.eigen_solution[eigenvector_num].rel_error_bound);
	fprintf(fout,"#\n");
	fprintf(fout,"# n, u_r, u_i, uAdj_r, uAdj_i, v_r, v_i, vAdj_r, vAdj_i, w_r, w_i, wAdj_r, wAdj_i, p_r, p_i, pAdj_r, pAdj_i\n");
	fprintf(fout,"#\n");
	
	fclose(fout);
}

//----------------------------------------------------------------------------
void Linear_Navier_Stokes_t::write_eigenvector(string filename_out, int sorted_eigenvector_num) 
{	
	FILE *fout;			fout=fopen(filename_out.c_str(),"a");
	
	int unsorted_mode_num = eigen_solver.eigen_solution[sorted_eigenvector_num].unsorted_mode_num;
	
	for (int jj=0; jj<=N; jj++) {
		fprintf(fout,"%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n",
				base_flow.n[jj],
				real(A(jj,unsorted_mode_num)),
				imag(A(jj,unsorted_mode_num)),
				real(B(jj,unsorted_mode_num)),
				imag(B(jj,unsorted_mode_num)),
				real(A(jj+N+1,unsorted_mode_num)),
				imag(A(jj+N+1,unsorted_mode_num)),
				real(B(jj+N+1,unsorted_mode_num)),
				imag(B(jj+N+1,unsorted_mode_num)),
				real(A(jj+2*(N+1),unsorted_mode_num)),
				imag(A(jj+2*(N+1),unsorted_mode_num)),
				real(B(jj+2*(N+1),unsorted_mode_num)),
				imag(B(jj+2*(N+1),unsorted_mode_num)),
				real(A(jj+3*(N+1),unsorted_mode_num)),
				imag(A(jj+3*(N+1),unsorted_mode_num)),
				real(B(jj+3*(N+1),unsorted_mode_num)),
				imag(B(jj+3*(N+1),unsorted_mode_num)) );
	}
	
	fclose(fout);
}

//----------------------------------------------------------------------------
