/*=========================================================================
 local_linear_stability
 =========================================================================
 
 Solves the local linear stability eigenvalue problem for a series of
 analytical profiles and also profiles input from file.
 
 The program maps from the complex wave number plane to the complex frequency
 plane in either a Design of Experiment (DoE) manner, or it can search
 for a wave number that matches a particular complex frequency. It can also 
 group modes from previous DoE runs.
 
 =========================================================================*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <blitz/array.h>
using namespace blitz;

#include "Utils.h"
#include "Input_Deck_t.h"
#include "LocalStability.h"

//----------------------------------------------------------------------------
// helper functions
//----------------------------------------------------------------------------
void run_EVP_DoE(Input_Deck_t input_deck, Linear_Navier_Stokes_t lnse_stability, Dispersion_Relationship_Mapping_t mapping);
void run_EVP_search(Input_Deck_t input_deck, Linear_Navier_Stokes_t lnse_stability, Dispersion_Relationship_Mapping_t mapping);

//----------------------------------------------------------------------------
int main (int argc, char * const argv[]) 
{
	cout << "Reading input deck ... " << endl;
	Input_Deck_t input_deck;
	
	cout << "Set problem parameters ... " << endl;
	Linear_Navier_Stokes_t lnse_stability;
	lnse_stability.set_equations(input_deck.equations);
	lnse_stability.set_number_of_eigenvectors_to_write(input_deck.number_of_eigenvectors_to_write);
	lnse_stability.set_boundary_conditions(input_deck.top_boundary_condition, input_deck.bottom_boundary_condition);
	
	lnse_stability.base_flow.set_base_flow_type(input_deck.base_flow);
	lnse_stability.base_flow.set_domain_size(input_deck.domain_size);
	lnse_stability.base_flow.set_shear_layer_parameter(input_deck.lambda);
	lnse_stability.base_flow.set_wake_parameter(input_deck.u_deficit, input_deck.wake_half_width);
	lnse_stability.base_flow.set_iteration_parameters(input_deck.tol, input_deck.max_iter, input_deck.min_iter);
	lnse_stability.base_flow.set_eddy_viscosity_parameters(input_deck.epsilon);
	
	cout << "Allocating memory ... " << endl;
	lnse_stability.allocate_memory(input_deck.N);
	
	cout << "Generating base flow ... " << endl;
	lnse_stability.generate_base_flow(input_deck.Re);
	
	cout << endl << "Setting dispersion relationship parameters ... " << endl;
	Dispersion_Relationship_Mapping_t mapping((lnse_stability.N+1)*2,
		input_deck.Re, lnse_stability.base_flow.x, input_deck.search_for_same_mode_type,
		input_deck.N_kx_real, input_deck.kx_real_max, input_deck.kx_real_min, 
		input_deck.N_kx_imag, input_deck.kx_imag_max, input_deck.kx_imag_min,
		input_deck.N_kz_real, input_deck.kz_real_max, input_deck.kz_real_min, input_deck.kz_imag);

	if (strcmp(input_deck.run_analysis,"EVP_DoE")==0) {
		cout << endl << "Running Eigenvalue Problem Design of Experiment ... " << endl;
		run_EVP_DoE(input_deck, lnse_stability, mapping);

	} else if (strcmp(input_deck.run_analysis,"group_modes")==0) {
		cout << endl << "Grouping modes from previous EVP_DoE run ... " << endl;
		mapping.group_modes();
		
	} else if (strcmp(input_deck.run_analysis,"EVP_search")==0) {
		cout << endl << "Running search for specific eigenvalue ... " << endl;
		run_EVP_search(input_deck, lnse_stability, mapping);
	}
	
	cout << "Analysis complete. " << endl;

	return 0;
}

//----------------------------------------------------------------------------
void run_EVP_DoE(Input_Deck_t input_deck, Linear_Navier_Stokes_t lnse_stability, Dispersion_Relationship_Mapping_t mapping)
{
	rmdir("./results");				mkdir("./results", 0777);
	lnse_stability.base_flow.write_profile();
	
	for (int ii=0; ii < input_deck.N_kx_real; ii++) {
		for (int jj=0; jj < mapping.get_N_jj_loop(); jj++) {
			fprintf(stdout, "Eigenvalue problem (%d, %d) of (%d, %d) ...\n", ii+1, jj+1, input_deck.N_kx_real, mapping.get_N_jj_loop());
			cout << " kx = " << mapping.get_kx(ii,jj) << endl;
			cout << " kz = " << mapping.get_kz(ii,jj) << endl;
			
			cout << "	Generating linear operator ... " << endl;
			lnse_stability.set_problem_parameters(input_deck.Re, mapping.get_kx(ii,jj), mapping.get_kz(ii,jj), input_deck.dissipation_ratio, mapping.L_0);
			lnse_stability.generate_linear_operator();			
			
			cout << "	Calculating eigensolution ... " << endl;
			lnse_stability.calculate_eigensolution();
			
			cout << "	Writing eigensolution output ... " << endl;
			lnse_stability.write_output(ii, jj);
			cout << endl;
		}
	}

	if ( input_deck.N_kx_real * mapping.get_N_jj_loop() > 1 ) {
		cout << "Grouping eigenvalues into like modes ... " << endl;
		mapping.group_modes();
	}
}

//----------------------------------------------------------------------------
void run_EVP_search(Input_Deck_t input_deck, Linear_Navier_Stokes_t lnse_stability, Dispersion_Relationship_Mapping_t mapping)
{	
	rmdir("./results.iter");		mkdir("./results.iter", 0777);
	
	complex<double> omega_search;
	real(omega_search) = input_deck.omega_real_search;
	imag(omega_search) = input_deck.omega_imag_search;
	mapping.read_dispersion_relatioship(input_deck.dispersion_relationship);
	
	double min_omega_error, *ptr_min_omega_error;	ptr_min_omega_error = &min_omega_error;
	complex<double> omega_orig, *ptr_omega_orig;	ptr_omega_orig = &omega_orig;
	complex<double> kx = mapping.estimate_kx_closest_to_this_omega(omega_search, ptr_min_omega_error, ptr_omega_orig);
	complex<double> final_kx = kx;
	complex<double> closest_omega = omega_orig;
	
	double kx_real_delta = (input_deck.kx_real_max - input_deck.kx_real_min) / (input_deck.N_kx_real * 1.0 - 1.0) / 2.0;
	double kx_imag_delta = (input_deck.kx_imag_max - input_deck.kx_imag_min) / (input_deck.N_kx_imag * 1.0 - 1.0) / 2.0;
	
	complex<double> next_kx, closest_omega_in_this_spectra;
	double omega_error;
	int iter=0, best_real_index=-1, best_imag_index=-1;
	
	int closest_eigenvalue_number = 0, eigenvalue_number, *ptr_eigenvalue_number;	ptr_eigenvalue_number = &eigenvalue_number;
	
	string filename("./results.iter/search_history.dat");
	FILE *fout=fopen(filename.c_str(),"w");
	fprintf(fout,"# target omega = %g + %g i \n",real(omega_search), imag(omega_search));
	fprintf(fout,"# iter, eigenvalue_number, best_real_index, best_imag_index, kx_imag, kx_real, omega_real, omega_imag, omega_error \n");
	fprintf(fout,"%d %d %d %d %g %g %g %g %g\n", 0, 0, 0, 0, real(kx), imag(kx), real(closest_omega), imag(closest_omega), min_omega_error);
	
	lnse_stability.set_number_of_eigenvectors_to_write(0);
	
	while ( ( (min_omega_error>input_deck.tol) && (iter<input_deck.max_iter)) || (iter<input_deck.min_iter) ) {
		
		best_real_index = -1;		// If better solution not found the same centre used with refined delta
		best_imag_index = -1;
		
		for (int r=0; r<2; r++) {
			for (int i=0; i<2; i++) {
				real(next_kx) = real(kx) + kx_real_delta * (2*r-1);
				imag(next_kx) = imag(kx) + kx_imag_delta * (2*i-1);
				
				fprintf(stdout, "Eigenvalue problem (%d, %d) of (2, 2) at iteration level %d ...\n", r+1, i+1, iter+1);
				cout << " kx = " << next_kx << endl;
				cout << " kz = " << input_deck.kz << endl;
				
				cout << "	Generating linear operator ... " << endl;
				lnse_stability.set_problem_parameters(input_deck.Re, next_kx, input_deck.kz, input_deck.dissipation_ratio, mapping.L_0);
				lnse_stability.generate_linear_operator();			
				
				cout << "	Calculating eigensolution ... " << endl;
				lnse_stability.calculate_eigensolution();
				closest_omega_in_this_spectra = lnse_stability.get_eigensolution_closest_to_this_omega(omega_search, ptr_eigenvalue_number);
				omega_error = SQ(real(closest_omega_in_this_spectra) - real(omega_search)) + SQ(imag(closest_omega_in_this_spectra) - imag(omega_search));
				
				if (omega_error < min_omega_error) {
					min_omega_error = omega_error;
					best_real_index = r;
					best_imag_index = i;
					closest_omega = closest_omega_in_this_spectra;
					closest_eigenvalue_number = eigenvalue_number;
					final_kx = next_kx;
					cout << "		min_omega_error = " << min_omega_error << endl;
					cout << "		final_kx = " << final_kx << endl;
				}
				
				cout << "	Writing eigensolution output ... " << endl << endl;
				lnse_stability.write_output(r,i,iter);
			}
		}
		
		if(best_real_index==-1) {
			fprintf(fout,"%d %d %d %d %g %g %g %g %g\n", 
					iter+1, closest_eigenvalue_number, 0, 0, 
					real(final_kx), imag(final_kx), real(closest_omega), imag(closest_omega), min_omega_error);
		} else {	
			fprintf(fout,"%d %d %d %d %g %g %g %g %g\n", 
					iter+1, closest_eigenvalue_number, best_real_index+1, best_imag_index+1, 
					real(final_kx), imag(final_kx), real(closest_omega), imag(closest_omega), min_omega_error);
			real(kx) = real(kx) + kx_real_delta * (2*best_real_index-1);
			imag(kx) = imag(kx) + kx_imag_delta * (2*best_imag_index-1);
		}
		
		iter++;
		kx_real_delta = kx_real_delta * 0.5;
		kx_imag_delta = kx_imag_delta * 0.5;
	}
	
	cout << "Repeating analysis of final kx writing out all eigenvectors ... " << endl;
	cout << "		final_kx = " << final_kx << endl;
	
	lnse_stability.set_number_of_eigenvectors_to_write(input_deck.N);	
	cout << "	Generating linear operator ... " << endl;
	lnse_stability.set_problem_parameters(input_deck.Re, final_kx, input_deck.kz, input_deck.dissipation_ratio, mapping.L_0);
	lnse_stability.generate_linear_operator();
	cout << "	Calculating eigensolution ... " << endl;
	lnse_stability.calculate_eigensolution();
	cout << "	Writing eigensolution output ... " << endl;
	lnse_stability.write_output(best_real_index,best_imag_index,iter-1);
	
	fclose(fout);
}

//----------------------------------------------------------------------------
