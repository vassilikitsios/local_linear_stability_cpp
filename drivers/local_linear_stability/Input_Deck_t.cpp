#include <blitz/array.h>

using namespace blitz;

#include "Utils.h"
#include "Input_Deck_t.h"

//----------------------------------------------------------------------------
Input_Deck_t::Input_Deck_t()
{
	cout << "Reading input deck..." << endl;
	
	ParamListHead_t param_list_head;
	
	param_list_head.get_params_from_file((char*)"local_linear_stability.in");
	param_list_head.print_params();

	run_analysis = param_list_head.get_string_param((char*)"run_analysis");
	if ( (strcmp(run_analysis,"EVP_DoE")!=0) 
		&& (strcmp(run_analysis,"group_modes")!=0)  
		&& (strcmp(run_analysis,"EVP_search")!=0) )
		quit_error((char*)"<run_analysis> = 'EVP_DoE' or 'group_modes' or 'EVP_search'\n");
	
	equations = param_list_head.get_string_param((char*)"equations");
	
	base_flow = param_list_head.get_string_param((char*)"base_flow");
	epsilon = param_list_head.get_double_param((char*)"epsilon");			// to handle singularity for eddy viscosity fields
	dissipation_ratio = param_list_head.get_double_param((char*)"dissipation_ratio");		// initial estimate of dissipation_ratio for unsteady stability analysis
	
	top_boundary_condition = param_list_head.get_string_param((char*)"top_boundary_condition");
	bottom_boundary_condition = param_list_head.get_string_param((char*)"bottom_boundary_condition");
	
	
	// Generating profile geometry
	domain_size = param_list_head.get_double_param((char*)"domain_size");	
	
	// Parameters for Blasius profile for testing
	tol = param_list_head.get_double_param((char*)"tol");						if(tol < 0.0) quit_error((char*)"tol must be > 0"); 
	max_iter = param_list_head.get_int_param((char*)"max_iter");				if(max_iter < 0) quit_error((char*)"max_iter must be > 0"); 
	min_iter = param_list_head.get_int_param((char*)"min_iter");				if(min_iter < 0) quit_error((char*)"min_iter must be > 0"); 
	
	// Parameters for shear layer profile
	lambda = param_list_head.get_double_param((char*)"lambda");			// shear later parameter lambda = (u2-u1)/(u2+u1)
	u_deficit = param_list_head.get_double_param((char*)"u_deficit");		// wake parameter
	wake_half_width = param_list_head.get_double_param((char*)"wake_half_width");
		
	// Stabilitiy Analysis Parameters
	N = param_list_head.get_int_param((char*)"N");								if(N < 0) quit_error((char*)"N must be > 0"); 
	Re = param_list_head.get_double_param((char*)"Re");						if(Re < 0) quit_error((char*)"Re must be > 0"); 
	
	
	char *search_for_same_mode_type_char;
	search_for_same_mode_type_char = param_list_head.get_string_param((char*)"search_for_same_mode_type");	
	if (strcmp(search_for_same_mode_type_char,"true")==0) {
		search_for_same_mode_type = true;
	} else if (strcmp(search_for_same_mode_type_char,"false")==0) {
		search_for_same_mode_type = false;
	} else {
		quit_error((char*)"<search_for_same_mode_type> = 'true' or 'false'\n");
	}
	
	N_kx_real = param_list_head.get_int_param((char*)"N_kx_real");
	kx_real_min = param_list_head.get_double_param((char*)"kx_real_min");
	kx_real_max = param_list_head.get_double_param((char*)"kx_real_max");
	N_kx_imag = param_list_head.get_int_param((char*)"N_kx_imag");
	kx_imag_min = param_list_head.get_double_param((char*)"kx_imag_min");
	kx_imag_max = param_list_head.get_double_param((char*)"kx_imag_max");
	real(kx) = kx_real_min;
	imag(kx) = kx_imag_min;
	
	N_kz_real = param_list_head.get_int_param((char*)"N_kz_real");
	kz_real_max = param_list_head.get_double_param((char*)"kz_real_max");
	kz_real_min = param_list_head.get_double_param((char*)"kz_real_min");
	kz_imag = param_list_head.get_double_param((char*)"kz_imag");
	real(kz) = kz_real_min;	
	imag(kz) = kz_imag;

	// Input data for eigenvalue search
	dispersion_relationship = param_list_head.get_string_param((char*)"dispersion_relationship");
	
	omega_real_search = param_list_head.get_double_param((char*)"omega_real_search");
	omega_imag_search = param_list_head.get_double_param((char*)"omega_imag_search");
	
	// Outputing
	number_of_eigenvectors_to_write = param_list_head.get_double_param((char*)"number_of_eigenvectors_to_write");
}

//----------------------------------------------------------------------------
Input_Deck_t::~Input_Deck_t()
{
}

//----------------------------------------------------------------------------
