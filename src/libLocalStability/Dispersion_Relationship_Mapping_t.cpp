
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <math.h>
#include <iostream>
#include<fstream>

#include <blitz/array.h>

using namespace blitz;

#include "Utils.h"
#include "EigenSolver.h"
#include "Dispersion_Relationship_Mapping_t.h"
#define C_GR_MAG_MAX			1.0
#define REL_ERROR_BOUNDS_TOL	1.0e-8
#define DIST_TOL				0.5

//----------------------------------------------------------------------------
Dispersion_Relationship_Mapping_t::Dispersion_Relationship_Mapping_t(int N_eigenvalues, double Re, double x, bool search_for_same_mode_type,
																	 int N_kx_real, double kx_real_max, double kx_real_min,
																	 int N_kx_imag, double kx_imag_max, double kx_imag_min,
																	 int N_kz_real, double kz_real_max, double kz_real_min,	double kz_imag)
{
	this->N_eigenvalues = N_eigenvalues;	
	this->Re = Re;
	this->x = x;
	this->search_for_same_mode_type = search_for_same_mode_type;
	
	if (N_kx_real < 1) quit_error((char*)"N_kx_real must be >= 1");	this->N_kx_real = N_kx_real;
	if (N_kx_imag < 1) quit_error((char*)"N_kx_imag must be >= 1");	this->N_kx_imag = N_kx_imag;
	if (N_kz_real < 1) quit_error((char*)"N_kz_real must be >= 1");	this->N_kz_real = N_kz_real;
	if ( (N_kz_real > 1) && (N_kx_imag > 1) ) quit_error((char*)"N_kz_real > 1 and N_kx_imag > 1, pick one.");
	if (N_kz_real==1) {
		N_jj_loop = N_kx_imag;
	} else {
		N_jj_loop = N_kz_real;
	}
	
	allocateArrays( shape(N_kx_real, N_jj_loop), omega_read, kx, kz);
	N_groups = N_kx_real * N_jj_loop;
	
	double denom_kx_real = 1.0*(N_kx_real-1);			if (N_kx_real == 1) denom_kx_real = 1.0;
	double denom_kx_imag = 1.0*(N_kx_imag-1);			if (N_kx_imag == 1) denom_kx_imag = 1.0;
	double denom_kz_real = 1.0*(N_kz_real-1);			if (N_kz_real == 1) denom_kz_real = 1.0;
	double kx_real_delta = (kx_real_max - kx_real_min) / denom_kx_real;
	double kx_imag_delta = (kx_imag_max - kx_imag_min) / denom_kx_imag;
	double kz_real_delta = (kz_real_max - kz_real_min) / denom_kz_real;

	cout << "ii jj kx(ii,jj) kz(ii,jj)" << endl;
	for (int ii = 0; ii < N_kx_real; ii++) {
		for (int jj = 0; jj < N_jj_loop; jj++) {
			
			real(kx(ii,jj)) = kx_real_min + kx_real_delta*ii;
			imag(kz(ii,jj)) = kz_imag;

			if (N_kz_real==1) {
				imag(kx(ii,jj)) = kx_imag_min + kx_imag_delta*jj;
				real(kz(ii,jj)) = kz_real_min;
			} else {
				imag(kx(ii,jj)) = kx_imag_min;
				real(kz(ii,jj)) = kz_real_min + kz_real_delta*jj;
			}
			cout << ii+1 << " " << jj+1 << " " << kx(ii,jj) << " " << kz(ii,jj) << endl;
		}
	}
	
	double kx2_0 = real( kx(0,0) * conj(kx(0,0)) );
	double kz2_0 = real( kz(0,0) * conj(kz(0,0)) );
	double Lx_0 = 0.0;
	if (kx2_0>0.0) Lx_0 = 2.0 * PI / sqrt(kx2_0);
	double Lz_0 = 0.0;
	if (kz2_0>0.0) Lz_0 = 2.0 * PI / sqrt(kz2_0);
	L_0 = sqrt( SQ(Lx_0) + SQ(Lz_0) );
}

//----------------------------------------------------------------------------
Dispersion_Relationship_Mapping_t::~Dispersion_Relationship_Mapping_t()
{
}

//----------------------------------------------------------------------------
int Dispersion_Relationship_Mapping_t::get_N_jj_loop()
{
	return N_jj_loop;
}

//----------------------------------------------------------------------------
complex<double> Dispersion_Relationship_Mapping_t::get_kx(int ii, int jj)
{
	complex<double> kx_out;
	kx_out = kx(ii,jj);
	return kx_out;
}

//----------------------------------------------------------------------------
complex<double> Dispersion_Relationship_Mapping_t::get_kz(int ii, int jj)
{
	complex<double> kz_out;
	kz_out = kz(ii,jj);
	return kz_out;
}

//----------------------------------------------------------------------------
complex<double> Dispersion_Relationship_Mapping_t::estimate_kx_closest_to_this_omega(complex<double> this_omega,
																					 double *ptr_min_omega_error, complex<double> *ptr_omega_orig)
{
	complex<double> closest_kx;
	
	double omega_error, min_omega_error = BIG_NUMBER * 1.0;
	int best_real_index = 0, best_imag_index = 0;
	
	for (int r=0; r < N_kx_real; r++) {
		for (int i=0; i < N_jj_loop; i++) {
			omega_error = SQ(real(this_omega) - real(omega_read(r,i))) + SQ(imag(this_omega) - imag(omega_read(r,i)));
			if (omega_error < min_omega_error) {
				best_real_index = r;
				best_imag_index = i;
				min_omega_error = omega_error;
			}
		}
	}
	*ptr_min_omega_error = min_omega_error;
	*ptr_omega_orig = omega_read(best_real_index, best_imag_index);
	
	return closest_kx = kx(best_real_index, best_imag_index);
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::read_dispersion_relatioship(string filename)
{
	char next_char;
	char comment_char;		comment_char = '#';	
	double kx_r_read, kx_i_read, kz_r_read, kz_i_read, Omega_r_read, Omega_i_read, c_r_read, c_i_read, c_gr_r_read, c_gr_i_read, rel_error_bounds_read, dissipation_ratio;
	
	ifstream fin(filename.c_str());
	if (fin.fail()) quit_error((char*)"file %s not found\n", filename.c_str());
	
	// Remove initial comment lines
	next_char = fin.peek();
	while(next_char == comment_char) {
		fin.ignore(INT_MAX, '\n');
		next_char = fin.peek();
	}
	
	// Get data
	for (int r=0; r < N_kx_real; r++ ) {
		for (int i=0; i < N_jj_loop; i++ ) {
			fin >> kx_r_read >> kx_i_read >> kz_r_read >> kz_i_read >> Omega_r_read >> Omega_i_read >> c_r_read >> c_i_read >> c_gr_r_read >> c_gr_i_read >> rel_error_bounds_read >> dissipation_ratio;
			real(kx(r,i)) = kx_r_read;
			imag(kx(r,i)) = kx_i_read;
			real(kz(r,i)) = kz_r_read;
			imag(kz(r,i)) = kz_i_read;
			real(omega_read(r,i)) = Omega_r_read;
			imag(omega_read(r,i)) = Omega_i_read;
		}
		//fin.ignore(INT_MAX, '\n');
	}
	fin.close();
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::group_modes()
{	
	Eigen_Values_t **eigenvalue_data = malloc_2d_array_Eigen_Values_t(N_groups, N_eigenvalues);
	read_data(eigenvalue_data);
	
	Eigen_Values_t **sorted_modes = malloc_2d_array_Eigen_Values_t(N_eigenvalues, N_groups);
//	sort_modes(eigenvalue_data, sorted_modes);
	sort_modes_old(eigenvalue_data, sorted_modes);
	
	for (int mode=0; mode < N_eigenvalues; mode++ ) {
		if (sorted_modes[mode][0].sorted==true) {
			cout << "Valid mode found : " << mode+1 << endl;
			if (N_kz_real==1) {
				calculate_group_velocity(mode, sorted_modes);
				find_pinch_point(mode, x, sorted_modes);
				write_real_omega(mode, 1.0/Re, sorted_modes, x);
				write_max_omega_imag_for_kx_real(mode, sorted_modes, x);
				write_zero_line();
			}
			write_complex_map(mode, 1.0/Re, sorted_modes);
		}
	}

	free_2d_array_Eigen_Values_t(eigenvalue_data, N_groups);
	free_2d_array_Eigen_Values_t(sorted_modes, N_eigenvalues);
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::read_data(Eigen_Values_t **eigenvalue_data)
{
	int file, N_modes;
	string filename;
	
	for (int ii=0; ii < N_kx_real; ii++ ) {
		for (int jj=0; jj < N_jj_loop; jj++ ) {
	
			file = jj + ii * N_jj_loop;
			//filename = get_eigenvalue_filename(ii, jj);
			filename = get_eigenvalue_filename(N_kx_real-1-ii, N_jj_loop-1-jj);
			N_modes = read_file(filename, file, eigenvalue_data);
			
			real(kx(ii,jj)) = eigenvalue_data[file][0].kx_r;
			imag(kx(ii,jj)) = eigenvalue_data[file][0].kx_i;
			real(kz(ii,jj)) = eigenvalue_data[file][0].kz_r;
			imag(kz(ii,jj)) = eigenvalue_data[file][0].kz_i;
		}
	}
}

//----------------------------------------------------------------------------
int Dispersion_Relationship_Mapping_t::read_file(string filename, int file, Eigen_Values_t **eigenvalue_data)
{
	char next_char;
	char comment_char;		comment_char = '#';	
	
	int num_read, discrete_read, w_mode_read, symmetic_read;
	double kx_r_read, kx_i_read, kz_r_read, kz_i_read, alpha_r_read, alpha_i_read, beta_r_read, beta_i_read, Omega_r_read, Omega_i_read, c_r_read, c_i_read, rel_error_bounds_read, cond_read;
	
	ifstream fin(filename.c_str());
	if (fin.fail()) quit_error((char*)"file %s not found\n", filename.c_str());
	
	// Remove initial comment lines
	next_char = fin.peek();
	while(next_char == comment_char) {
		fin.ignore(INT_MAX, '\n');
		next_char = fin.peek();
	}
	
	// Get data
	for (int i=0; i<N_eigenvalues; i++) {
		fin >> num_read >> kx_r_read >> kx_i_read >> kz_r_read >> kz_i_read >> alpha_r_read >> alpha_i_read >> beta_r_read >> beta_i_read >> Omega_r_read >> Omega_i_read >> c_r_read >> c_i_read >> discrete_read >> w_mode_read >> symmetic_read >> rel_error_bounds_read >> cond_read;
		
		eigenvalue_data[file][i].kx_r = kx_r_read;
		eigenvalue_data[file][i].kx_i = kx_i_read;
		eigenvalue_data[file][i].kz_r = kz_r_read;
		eigenvalue_data[file][i].kz_i = kz_i_read;
		eigenvalue_data[file][i].Omega_r = Omega_r_read;
		eigenvalue_data[file][i].Omega_i = Omega_i_read;
		eigenvalue_data[file][i].c_r = c_r_read;
		eigenvalue_data[file][i].c_i = c_i_read;
		eigenvalue_data[file][i].discrete_continuous = discrete_read;
		eigenvalue_data[file][i].w_mode = w_mode_read;
		eigenvalue_data[file][i].symmetric = symmetic_read;
		eigenvalue_data[file][i].rel_error_bounds = rel_error_bounds_read;
		eigenvalue_data[file][i].cond = cond_read;
	}
	
	fin.close();
	return num_read;
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::sort_modes_old(Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes)
{
	int group, group_check;
	
	for (int mode=0; mode < N_eigenvalues; mode++ ) {
		
		eigenvalue_data[0][mode].sorted = true;
		copy_2d_array_Eigen_Values_t_from_A_to_B(eigenvalue_data, 0, mode, sorted_modes, mode, 0);
		
		for (int ii=0; ii < N_kx_real; ii++ ) {
			for (int jj=0; jj < N_jj_loop; jj++ ) {
				group = jj + ii * N_jj_loop;
				if ( group > 0 ) {
					if (jj==0) {
						group_check = group - N_jj_loop;
					} else {
						group_check = group - 1;
					}
					find_closet_mode_in_group(mode, group, group_check, eigenvalue_data, sorted_modes);
				}
			}
		}
		
	}
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::sort_modes(Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes)
{
	int first_group, first_unsorted_mode;

	for (int mode=0; mode < N_eigenvalues; mode++ ) {
		
		first_group = mode % N_groups;
		first_unsorted_mode = 0;
		while( (first_unsorted_mode<N_eigenvalues) 
			  && ( (eigenvalue_data[first_group][first_unsorted_mode].sorted==true) 
				  || (eigenvalue_data[first_group][first_unsorted_mode].rel_error_bounds>REL_ERROR_BOUNDS_TOL)
				  || (eigenvalue_data[first_group][first_unsorted_mode].rel_error_bounds<0.0) ) )
			first_unsorted_mode++;
		
		eigenvalue_data[first_group][first_unsorted_mode].sorted = true;
		copy_2d_array_Eigen_Values_t_from_A_to_B(eigenvalue_data, first_group, first_unsorted_mode, sorted_modes, mode, first_group);

		if( ( (first_group+1) % N_jj_loop!=0) || (first_group==0) ) {
			search_backward(first_group, mode, eigenvalue_data, sorted_modes);
			search_forward(first_group, mode, eigenvalue_data, sorted_modes);
		} else {
			search_forward(first_group, mode, eigenvalue_data, sorted_modes);
			search_backward(first_group, mode, eigenvalue_data, sorted_modes);
		}
	}
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::search_forward(int first_group, int mode, Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes)
{
	int group_check;
	for (int group=first_group+1; group < N_groups; group++ ) {
		if (group % N_jj_loop==0) {
			group_check = group - N_jj_loop;
		} else {
			group_check = group - 1;
		}
		find_closet_mode_in_group(mode, group, group_check, eigenvalue_data, sorted_modes);
	}
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::search_backward(int first_group, int mode, Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes)
{
	int group_check;
	if (first_group>0) {
		for (int group=first_group-1; group>=0; group-- ) {
			if ( ((group+1) % N_jj_loop!=0) || (group==0) ){
				group_check = group + 1;
			} else {
				if (first_group < N_jj_loop) quit_error((char*)"Problem with logic");
				group_check = group + N_jj_loop;
			}
			find_closet_mode_in_group(mode, group, group_check, eigenvalue_data, sorted_modes);
		}
	}
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::find_closet_mode_in_group(int mode, int group, int group_check, 
																  Eigen_Values_t **eigenvalue_data, Eigen_Values_t **sorted_modes)
{
	double dist_to_mode, min_dist_to_mode = 1.0 * BIG_NUMBER;
	int closest_mode = -1;
	bool same_mode_type;
	double epsilon = 1.0e-6;
	
	for (int next_mode=0; next_mode < N_eigenvalues; next_mode++ ) {
		same_mode_type = false;
		if ( (eigenvalue_data[group][next_mode].rel_error_bounds<REL_ERROR_BOUNDS_TOL) &&
			(eigenvalue_data[group][next_mode].rel_error_bounds>0.0) ) {
			
			if (   (eigenvalue_data[group][next_mode].discrete_continuous	== sorted_modes[mode][group_check].discrete_continuous)
				&& (eigenvalue_data[group][next_mode].w_mode				== sorted_modes[mode][group_check].w_mode)
				&& (eigenvalue_data[group][next_mode].symmetric				== sorted_modes[mode][group_check].symmetric) )
				same_mode_type = true;
				
			if ( (same_mode_type==true) || (search_for_same_mode_type==false) ) {
				
				if ( (SQ(eigenvalue_data[group][next_mode].kx_r)+SQ(eigenvalue_data[group][next_mode].kx_i)<=epsilon) 
					|| (SQ(sorted_modes[mode][group_check].kx_r)+SQ(sorted_modes[mode][group_check].kx_i)<=epsilon) ) {
					dist_to_mode = 
					SQ(eigenvalue_data[group][next_mode].Omega_r / (eigenvalue_data[group][next_mode].kx_r + epsilon) 
						 - sorted_modes[mode][group_check].Omega_r / (sorted_modes[mode][group_check].kx_r + epsilon) ) 
					 + SQ(eigenvalue_data[group][next_mode].Omega_i / (eigenvalue_data[group][next_mode].kx_r + epsilon)
						  - sorted_modes[mode][group_check].Omega_i / (sorted_modes[mode][group_check].kx_r + epsilon) );
				} else {
					dist_to_mode = 
					( SQ(eigenvalue_data[group][next_mode].c_r - sorted_modes[mode][group_check].c_r) 
					 + SQ(eigenvalue_data[group][next_mode].c_i - sorted_modes[mode][group_check].c_i) );
				}
				
				if (dist_to_mode < min_dist_to_mode) {
					min_dist_to_mode = dist_to_mode;
					closest_mode = next_mode;
				}
			}
		}
	}
	
	eigenvalue_data[group][closest_mode].sorted = true;
	if (closest_mode==-1) {
		copy_2d_array_Eigen_Values_t_from_A_to_B(sorted_modes, mode, group_check, sorted_modes, mode, group);
		sorted_modes[mode][group].sorted = false;
		sorted_modes[mode][group].discrete_continuous = -1;
		sorted_modes[mode][group].w_mode = -1;
		sorted_modes[mode][group].symmetric = -1;
		sorted_modes[mode][group_check].discrete_continuous = -1;
		sorted_modes[mode][group_check].w_mode = -1;
		sorted_modes[mode][group_check].symmetric = -1;
	} else {
		eigenvalue_data[group][closest_mode].sorted = true;
		copy_2d_array_Eigen_Values_t_from_A_to_B(eigenvalue_data, group, closest_mode, sorted_modes, mode, group);
	}
}

//----------------------------------------------------------------------------
bool Dispersion_Relationship_Mapping_t::determine_if_this_a_valid_mode(int mode, Eigen_Values_t **sorted_modes)
{
	bool write_mode = true;
	int group;
	for (int ii=0; ii < N_kx_real; ii++ ) {
		for (int jj=0; jj < N_jj_loop; jj++ ) {
			group = jj + ii * N_jj_loop;
			if ( (sorted_modes[mode][group].rel_error_bounds > REL_ERROR_BOUNDS_TOL) ||
				(sorted_modes[mode][group].rel_error_bounds < 0.0) )
				write_mode = false;
		}
	}
	return write_mode;
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::calculate_group_velocity(int mode, Eigen_Values_t **sorted_modes)
{
	int group, other_group;
	
	double delta_Omega_r, delta_Omega_i, delta_kx_r, delta_kx_i, mag_delta_kx;
	for (int i=0; i < N_kx_imag; i++ ) {
		for (int r=0; r < N_kx_real; r++ ) {
			group = i + r * N_kx_imag;
			if (r>0) {
				other_group = i + (r-1) * N_kx_imag;
			} else {
				other_group = i + (r+1) * N_kx_imag;
			}
			if (i>0) {
				other_group = other_group - 1;
			} else {
				other_group = other_group + 1;
			}
			delta_Omega_r = sorted_modes[mode][group].Omega_r - sorted_modes[mode][other_group].Omega_r;
			delta_Omega_i = sorted_modes[mode][group].Omega_i - sorted_modes[mode][other_group].Omega_i;
			delta_kx_r = sorted_modes[mode][group].kx_r - sorted_modes[mode][other_group].kx_r;
			delta_kx_i = sorted_modes[mode][group].kx_i - sorted_modes[mode][other_group].kx_i;
			mag_delta_kx = SQ(delta_kx_r) + SQ(delta_kx_i);
			sorted_modes[mode][group].c_gr_r = (delta_Omega_r * delta_kx_r + delta_Omega_i * delta_kx_i) / mag_delta_kx;
			sorted_modes[mode][group].c_gr_i = (delta_Omega_i * delta_kx_r - delta_Omega_r * delta_kx_i) / mag_delta_kx;
		}
	}
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::find_pinch_point(int mode, double x_pos, Eigen_Values_t **sorted_modes)
{
	FILE *fout;
	string filename;
	int group, group_min_c_gr_mag=-1;
	double c_gr_mag = 0.0, min_c_gr_mag = BIG_NUMBER*1.0;
	
	for (int i=0; i < N_kx_imag; i++ ) {
		for (int r=0; r < N_kx_real; r++ ) {
			group = i + r * N_kx_imag;
			c_gr_mag = sqrt(SQ(sorted_modes[mode][group].c_gr_r) + SQ(sorted_modes[mode][group].c_gr_i));
			if (min_c_gr_mag > c_gr_mag) {
				min_c_gr_mag = c_gr_mag;
				group_min_c_gr_mag = group;
			}			
		}
		
		if ( (min_c_gr_mag<C_GR_MAG_MAX) && (group_min_c_gr_mag != -1) ) {
			filename = get_pinch_point_filename(mode+1);
			fout=fopen(filename.c_str(),"w");
			fprintf(fout,"# x_pos, mode, kx_real, kx_imag, Omega_real, Omega_imag, c_gr_real, c_gr_imag, kz_real, kz_imag\n");
			fprintf(fout,"%g %d %g %g %g %g %g %g %g %g\n", x_pos, mode, 
					sorted_modes[mode][group_min_c_gr_mag].kx_r, sorted_modes[mode][group_min_c_gr_mag].kx_i,
					sorted_modes[mode][group_min_c_gr_mag].Omega_r, sorted_modes[mode][group_min_c_gr_mag].Omega_i,
					sorted_modes[mode][group_min_c_gr_mag].c_gr_r, sorted_modes[mode][group_min_c_gr_mag].c_gr_i,
					sorted_modes[mode][group_min_c_gr_mag].kz_r, sorted_modes[mode][group_min_c_gr_mag].kz_i);
			fclose(fout);
		}
	}
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::write_complex_map(int mode, double nu, Eigen_Values_t **sorted_modes) 
{
	string filename_out = get_complex_map_filename(mode+1);
	FILE *fout=fopen(filename_out.c_str(),"w");
	
	fprintf(fout,"# Dispersion relationship\n");
	fprintf(fout,"# kx_r, kx_i, kz_r, kz_i, Omega_r, Omega_i, c_r, c_i, c_gr_r, c_gr_i, discrete=0 continuous=1, u_mode=0 w_mode=1, asymmetic=0 symmetic=1, rel_error_bound, condition number\n");
	fprintf(fout,"# Re = \"%8.6f\"\n", 1.0/nu);
	
	int group;	
	for (int ii=0; ii < N_kx_real; ii++ ) {
		for (int jj=0; jj < N_jj_loop; jj++ ) {
			group = jj + ii * N_jj_loop;
			
			fprintf(fout,"%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %d %d %d %16.10e %16.10e\n",
					sorted_modes[mode][group].kx_r, sorted_modes[mode][group].kx_i, sorted_modes[mode][group].kz_r, sorted_modes[mode][group].kz_i,
					sorted_modes[mode][group].Omega_r, sorted_modes[mode][group].Omega_i,
					sorted_modes[mode][group].c_r, sorted_modes[mode][group].c_i, 
					sorted_modes[mode][group].c_gr_r, sorted_modes[mode][group].c_gr_i,
					sorted_modes[mode][group].discrete_continuous, sorted_modes[mode][group].w_mode, sorted_modes[mode][group].symmetric,
					sorted_modes[mode][group].rel_error_bounds, sorted_modes[mode][group].cond); 
		}
		fprintf(fout,"\n");
	}
	fclose(fout);
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::write_zero_line()
{		
	FILE *fout;			
	string filename_out = "./images/zero_line.dat";
	fout=fopen(filename_out.c_str(),"w");
	fprintf(fout,"# zero line for complex mapping\n");
	fprintf(fout,"0.0     0.0\n");
	fprintf(fout,"0.0     0.0\n");
	fprintf(fout,"\n");
	fprintf(fout,"1.0     0.0\n");
	fprintf(fout,"1.0     0.0\n");
	fclose(fout);
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::write_max_omega_imag_for_kx_real(int mode, Eigen_Values_t **sorted_modes, double x_pos)
{
	double max_omega = -1.0*BIG_NUMBER;
	int group, max_group = -1;
	
	for (int i=0; i < N_kx_imag; i++ ) {
		for (int r=0; r < N_kx_real; r++ ) {
			group = i + r * N_kx_imag;
			if(ABS(sorted_modes[mode][group].kx_i)<1.0e-6) {
				if(max_omega < sorted_modes[mode][group].Omega_i) {
					max_omega = sorted_modes[mode][group].Omega_i;
					max_group = group;
				}
			}
		}
	}
	
	if (max_group>0) {
		string filename;
		filename = get_max_omega_imag_for_kx_real_filename(mode+1);
		FILE *fout=fopen(filename.c_str(),"w");
		fprintf(fout,"# x_pos, mode, kx_real, kx_imag=0, Omega_real, MAX Omega_imag, kz_real, kz_imag\n");
		fprintf(fout,"%g %d %g %g %g %g %g %g\n", 
				x_pos, mode, sorted_modes[mode][max_group].kx_r, sorted_modes[mode][max_group].kx_i,
				sorted_modes[mode][max_group].Omega_r, sorted_modes[mode][max_group].Omega_i, real(kz(0,0)), imag(kz(0,0)));
		fclose(fout);
	}
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::write_real_omega(int mode, double nu, Eigen_Values_t **sorted_modes, double x_pos)
{		
	FILE *fout;			
	string filename_out;
	int num_zeros;

	int N_interp_points = N_kx_real * N_kx_imag;
	// Data from most unstable mode interpolated onto line imag(omega)=0
	double *omega_real_interp = malloc_1d_array_double(N_interp_points);
	double *kx_real_interp = malloc_1d_array_double(N_interp_points);
	double *kx_imag_interp = malloc_1d_array_double(N_interp_points);
	double *c_gr_real_interp = malloc_1d_array_double(N_interp_points);
	double *c_gr_imag_interp = malloc_1d_array_double(N_interp_points);
	
	int N_spline_points = N_interp_points * 2;
	// Cubic spline used to add additional points for plotting
	double *omega_real_spline = malloc_1d_array_double(N_spline_points);
	double *kx_real_spline = malloc_1d_array_double(N_spline_points);
	double *kx_imag_spline = malloc_1d_array_double(N_spline_points);
	double *c_gr_real_spline = malloc_1d_array_double(N_spline_points);
	double *c_gr_imag_spline = malloc_1d_array_double(N_spline_points);
	
	num_zeros = find_where_omega_is_completely_real(N_interp_points, mode, sorted_modes,
													omega_real_interp, c_gr_real_interp, c_gr_imag_interp, kx_real_interp, kx_imag_interp);

	double c_r, c_i, c_gr_r, c_gr_i, mag_delta_kx;
	int pos;
	if (num_zeros>3) {

		filename_out = get_real_omega_raw_filename(mode+1);
		fout=fopen(filename_out.c_str(),"w");
		fprintf(fout,"# Dispersion relationship\n");
		fprintf(fout,"# kx_r, kx_i, -kx_i, Omega_r, Omega_i, c_r, c_i, c_gr_r, c_gr_i\n");
		fprintf(fout,"# Re = \"%8.6f\"\n", 1.0/nu);	
		fprintf(fout,"# k_z = \"(%8.6f, i%8.6f)\"\n", sorted_modes[mode][0].kz_r, sorted_modes[mode][0].kz_i);
		fprintf(fout,"# NOTE: data is repeated so that it can be plotted as a surface.\n");

		for (int ii=0; ii<num_zeros; ii++) {
			c_r = omega_real_interp[ii] * kx_real_interp[ii] / (SQ(kx_real_interp[ii]) + SQ(kx_imag_interp[ii]) );
			c_i = -omega_real_interp[ii] * kx_imag_interp[ii] / (SQ(kx_real_interp[ii]) + SQ(kx_imag_interp[ii]) );
			pos = 1;
			if (ii>0) pos = ii + 1;
			mag_delta_kx = SQ(kx_real_interp[ii] - kx_real_interp[pos]) + SQ(kx_imag_interp[ii] - kx_imag_interp[pos]);
			c_gr_r = ( omega_real_interp[ii] - omega_real_interp[pos] ) * (kx_real_interp[ii] - kx_real_interp[pos]) / mag_delta_kx;
			c_gr_i = -( omega_real_interp[ii] - omega_real_interp[pos] ) * (kx_imag_interp[ii] - kx_imag_interp[pos]) / mag_delta_kx;
			fprintf(fout,"%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n", 
					kx_real_interp[ii], kx_imag_interp[ii], -kx_imag_interp[ii], 
					omega_real_interp[ii], 0.0, c_r, c_i, c_gr_r, c_gr_i);
		}
		fprintf(fout,"\n");

		for (int ii=0; ii<num_zeros; ii++) {
			c_r = omega_real_interp[ii] * kx_real_interp[ii] / (SQ(kx_real_interp[ii]) + SQ(kx_imag_interp[ii]) );
			c_i = -omega_real_interp[ii] * kx_imag_interp[ii] / (SQ(kx_real_interp[ii]) + SQ(kx_imag_interp[ii]) );
			pos = 1;
			if (ii>0) pos = ii + 1;
			mag_delta_kx = SQ(kx_real_interp[ii] - kx_real_interp[pos]) + SQ(kx_imag_interp[ii] - kx_imag_interp[pos]);
			c_gr_r = ( omega_real_interp[ii] - omega_real_interp[pos] ) * (kx_real_interp[ii] - kx_real_interp[pos]) / mag_delta_kx;
			c_gr_i = -( omega_real_interp[ii] - omega_real_interp[pos] ) * (kx_imag_interp[ii] - kx_imag_interp[pos]) / mag_delta_kx;
			fprintf(fout,"%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n", 
					kx_real_interp[ii], kx_imag_interp[ii], -kx_imag_interp[ii], 
					omega_real_interp[ii], 0.0, c_r, c_i, c_gr_r, c_gr_i);
		}
		fclose(fout);
		
		interpolate_along_spline(N_spline_points, num_zeros, mode,
								 omega_real_interp, kx_real_interp, kx_imag_interp, c_gr_real_interp, c_gr_imag_interp,
								 omega_real_spline, kx_real_spline, kx_imag_spline, c_gr_real_spline, c_gr_imag_spline,
								 x_pos);

		filename_out = get_real_omega_filename(mode+1);
		fout=fopen(filename_out.c_str(),"w");
		fprintf(fout,"# Dispersion relationship\n");
		fprintf(fout,"# kx_r, kx_i, -kx_i, Omega_r, Omega_i, c_r, c_i, c_gr_r, c_gr_i\n");
		fprintf(fout,"# Re = \"%8.6f\"\n", 1.0/nu);	
		fprintf(fout,"# k_z = \"(%8.6f, i%8.6f)\"\n", sorted_modes[mode][0].kz_r, sorted_modes[mode][0].kz_i);
		fprintf(fout,"# NOTE: data is repeated so that it can be plotted as a surface.\n");

		for (int ii=0; ii<N_spline_points; ii++) {
			c_r = omega_real_spline[ii] * kx_real_interp[ii] / (SQ(kx_real_spline[ii]) + SQ(kx_imag_spline[ii]) );
			c_i = -omega_real_spline[ii] * kx_imag_spline[ii] / (SQ(kx_real_spline[ii]) + SQ(kx_imag_spline[ii]) );
			pos = 1;
			if (ii>0) pos = ii + 1;
			mag_delta_kx = SQ(kx_real_spline[ii] - kx_real_spline[pos]) + SQ(kx_imag_spline[ii] - kx_imag_spline[pos]);
			c_gr_r = ( omega_real_spline[ii] - omega_real_spline[pos] ) * (kx_real_spline[ii] - kx_real_spline[pos]) / mag_delta_kx;
			c_gr_i = -( omega_real_spline[ii] - omega_real_spline[pos] ) * (kx_imag_spline[ii] - kx_imag_spline[pos]) / mag_delta_kx;
			fprintf(fout,"%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n", 
					kx_real_spline[ii], kx_imag_spline[ii], -kx_imag_spline[ii], 
					omega_real_spline[ii], 0.0, c_r, c_i, c_gr_r, c_gr_i);
		}
		fprintf(fout,"\n");

		for (int ii=0; ii<N_spline_points; ii++) {
			c_r = omega_real_spline[ii] * kx_real_interp[ii] / (SQ(kx_real_spline[ii]) + SQ(kx_imag_spline[ii]) );
			c_i = -omega_real_spline[ii] * kx_imag_spline[ii] / (SQ(kx_real_spline[ii]) + SQ(kx_imag_spline[ii]) );
			pos = 1;
			if (ii>0) pos = ii + 1;
			mag_delta_kx = SQ(kx_real_spline[ii] - kx_real_spline[pos]) + SQ(kx_imag_spline[ii] - kx_imag_spline[pos]);
			c_gr_r = ( omega_real_spline[ii] - omega_real_spline[pos] ) * (kx_real_spline[ii] - kx_real_spline[pos]) / mag_delta_kx;
			c_gr_i = -( omega_real_spline[ii] - omega_real_spline[pos] ) * (kx_imag_spline[ii] - kx_imag_spline[pos]) / mag_delta_kx;
			fprintf(fout,"%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n", 
					kx_real_spline[ii], kx_imag_spline[ii], -kx_imag_spline[ii], 
					omega_real_spline[ii], 0.0, c_r, c_i, c_gr_r, c_gr_i);
		}
		fclose(fout);

	}
	
	free_1d_array_double(omega_real_interp);
	free_1d_array_double(kx_real_interp);
	free_1d_array_double(kx_imag_interp);
	free_1d_array_double(c_gr_real_interp);
	free_1d_array_double(c_gr_imag_interp);
	
	free_1d_array_double(omega_real_spline);
	free_1d_array_double(kx_real_spline);
	free_1d_array_double(kx_imag_spline);
	free_1d_array_double(c_gr_real_spline);
	free_1d_array_double(c_gr_imag_spline);
}

//----------------------------------------------------------------------------
int Dispersion_Relationship_Mapping_t::find_where_omega_is_completely_real(int N_interp_points, int mode, Eigen_Values_t **sorted_modes,
																		   double *omega_real_interp, 
																		   double *c_gr_real_interp, double *c_gr_imag_interp,
																		   double *kx_real_interp, double *kx_imag_interp)
{
	bool zero_found;
	int group, num_zeros=0, num_zeros_H=0, num_zeros_V=0;
	
	// zero arrays
	for (int point=0; point < N_interp_points; point++ ) {
		omega_real_interp[point] = 0.0;
		kx_real_interp[point] = 0.0;
		kx_imag_interp[point] = 0.0;
		c_gr_real_interp[point] = 0.0;
		c_gr_imag_interp[point] = 0.0;
	}
	
	// search_vertical
	for (int r = 0; r < N_kx_real; r++) {
		zero_found = false;
		
		
		if(sorted_modes[mode][r * N_kx_imag].Omega_i>0.0) {			
			for (int i = 1; i < N_kx_imag; i++) {
				group = i + r * N_kx_imag;
				
				if ( (zero_found==false) && (sorted_modes[mode][group].Omega_i<0.0) && (sorted_modes[mode][group-1].kx_r>0.0) ) {
					zero_found = true;
					
					kx_real_interp[num_zeros_V] = sorted_modes[mode][group].kx_r; 
					
					kx_imag_interp[num_zeros_V] = sorted_modes[mode][group-1].kx_i
					- sorted_modes[mode][group-1].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-1].Omega_i) 
					* ( sorted_modes[mode][group].kx_i - sorted_modes[mode][group-1].kx_i );
					
					c_gr_imag_interp[num_zeros_V] = sorted_modes[mode][group-1].c_gr_i
					- sorted_modes[mode][group-1].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-1].Omega_i) 
					* ( sorted_modes[mode][group].c_gr_i - sorted_modes[mode][group-1].c_gr_i );
					
					c_gr_real_interp[num_zeros_V] = sorted_modes[mode][group-1].c_gr_r
					- sorted_modes[mode][group-1].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-1].Omega_i) 
					* ( sorted_modes[mode][group].c_gr_r - sorted_modes[mode][group-1].c_gr_r );
					
					omega_real_interp[num_zeros_V] = sorted_modes[mode][group-1].Omega_r
					- sorted_modes[mode][group-1].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-1].Omega_i) 
					* ( sorted_modes[mode][group].Omega_r - sorted_modes[mode][group-1].Omega_r );
					
					num_zeros_V++;
				}
			}
		}	
		
		
		
		if(sorted_modes[mode][r * N_kx_imag].Omega_i<0.0) {			
			for (int i = 1; i < N_kx_imag; i++) {
				group = i + r * N_kx_imag;
				
				if ( (zero_found==false) && (sorted_modes[mode][group].Omega_i>0.0) && (sorted_modes[mode][group-1].kx_r>0.0) ) {
					zero_found = true;
					
					kx_real_interp[num_zeros_V] = sorted_modes[mode][group].kx_r; 
					
					c_gr_imag_interp[num_zeros_V] = sorted_modes[mode][group-1].c_gr_i
					- sorted_modes[mode][group-1].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-1].Omega_i) 
					* ( sorted_modes[mode][group].c_gr_i - sorted_modes[mode][group-1].c_gr_i );
					
					c_gr_real_interp[num_zeros_V] = sorted_modes[mode][group-1].c_gr_r
					- sorted_modes[mode][group-1].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-1].Omega_i) 
					* ( sorted_modes[mode][group].c_gr_r - sorted_modes[mode][group-1].c_gr_r );
					
					kx_imag_interp[num_zeros_V] = sorted_modes[mode][group-1].kx_i
					- sorted_modes[mode][group-1].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-1].Omega_i) 
					* ( sorted_modes[mode][group].kx_i - sorted_modes[mode][group-1].kx_i );
					
					omega_real_interp[num_zeros_V] = sorted_modes[mode][group-1].Omega_r
					- sorted_modes[mode][group-1].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-1].Omega_i) 
					* ( sorted_modes[mode][group].Omega_r - sorted_modes[mode][group-1].Omega_r );
					
					num_zeros_V++;
				}
			}
		}	
		
		
	}
	
	// search_horizontal
	// Data from most unstable mode interpolated onto line imag(omega)=0
	double *omega_real_interp_H = malloc_1d_array_double(N_interp_points);
	double *kx_real_interp_H = malloc_1d_array_double(N_interp_points);
	double *kx_imag_interp_H = malloc_1d_array_double(N_interp_points);
	double *c_gr_real_interp_H = malloc_1d_array_double(N_interp_points);
	double *c_gr_imag_interp_H = malloc_1d_array_double(N_interp_points);
	
	for (int i = 0; i < N_kx_imag; i++) {
		zero_found = false;
		
		if(sorted_modes[mode][i].Omega_i>0.0) {	
			for (int r = 1; r < N_kx_real; r++) {
				group = i + r * N_kx_imag;
				
				if ( (zero_found==false) && (sorted_modes[mode][group].Omega_i<0.0) && (sorted_modes[mode][group-1].kx_r>0.0) ) {
					zero_found = true;
					
					kx_imag_interp_H[num_zeros_H] = sorted_modes[mode][group].kx_i; 
					
					kx_real_interp_H[num_zeros_H] = sorted_modes[mode][group-1].kx_r
					- sorted_modes[mode][group-N_kx_imag].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-N_kx_imag].Omega_i) 
					* ( sorted_modes[mode][group].kx_r - sorted_modes[mode][group-N_kx_imag].kx_r );
					
					c_gr_real_interp_H[num_zeros_H] = sorted_modes[mode][group-1].c_gr_r
					- sorted_modes[mode][group-N_kx_imag].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-N_kx_imag].Omega_i) 
					* ( sorted_modes[mode][group].c_gr_r - sorted_modes[mode][group-N_kx_imag].c_gr_r );
					
					c_gr_imag_interp_H[num_zeros_H] = sorted_modes[mode][group-1].c_gr_i
					- sorted_modes[mode][group-N_kx_imag].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-N_kx_imag].Omega_i) 
					* ( sorted_modes[mode][group].c_gr_i - sorted_modes[mode][group-N_kx_imag].c_gr_i );
					
					omega_real_interp_H[num_zeros_H] = sorted_modes[mode][group-1].Omega_r
					- sorted_modes[mode][group-N_kx_imag].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-N_kx_imag].Omega_i) 
					* ( sorted_modes[mode][group].Omega_r - sorted_modes[mode][group-N_kx_imag].Omega_r );
					
					num_zeros_H++;
				}
			}
		}
		
		
		if(sorted_modes[mode][i].Omega_i<0.0) {	
			for (int r = 1; r < N_kx_real; r++) {
				group = i + r * N_kx_imag;
				
				if ( (zero_found==false) && (sorted_modes[mode][group].Omega_i>0.0) && (sorted_modes[mode][group-1].kx_r>0.0) ) {
					zero_found = true;
					
					kx_imag_interp_H[num_zeros_H] = sorted_modes[mode][group].kx_i; 
					
					kx_real_interp_H[num_zeros_H] = sorted_modes[mode][group-1].kx_r
					- sorted_modes[mode][group-N_kx_imag].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-N_kx_imag].Omega_i) 
					* ( sorted_modes[mode][group].kx_r - sorted_modes[mode][group-N_kx_imag].kx_r );
					
					c_gr_real_interp_H[num_zeros_H] = sorted_modes[mode][group-1].c_gr_r
					- sorted_modes[mode][group-N_kx_imag].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-N_kx_imag].Omega_i) 
					* ( sorted_modes[mode][group].c_gr_r - sorted_modes[mode][group-N_kx_imag].c_gr_r );
					
					c_gr_imag_interp_H[num_zeros_H] = sorted_modes[mode][group-1].c_gr_i
					- sorted_modes[mode][group-N_kx_imag].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-N_kx_imag].Omega_i) 
					* ( sorted_modes[mode][group].c_gr_i - sorted_modes[mode][group-N_kx_imag].c_gr_i );
					
					omega_real_interp_H[num_zeros_H] = sorted_modes[mode][group-1].Omega_r
					- sorted_modes[mode][group-N_kx_imag].Omega_i 
					/ (sorted_modes[mode][group].Omega_i - sorted_modes[mode][group-N_kx_imag].Omega_i) 
					* ( sorted_modes[mode][group].Omega_r - sorted_modes[mode][group-N_kx_imag].Omega_r );
					
					num_zeros_H++;
				}
			}
		}
		
	}
	
	if (num_zeros_H > num_zeros_V) {
		num_zeros = num_zeros_H;
		for (int k = 0; k < num_zeros; k++) {
			omega_real_interp[k] = omega_real_interp_H[k];
			kx_real_interp[k] = kx_real_interp_H[k];
			kx_imag_interp[k] = kx_imag_interp_H[k];
			c_gr_real_interp[k] = c_gr_real_interp_H[k];
			c_gr_imag_interp[k] = c_gr_imag_interp_H[k];
		}
	} else {
		num_zeros = num_zeros_V;
	}
	
	double temp_omega_real, temp_kx_real, temp_kx_imag, temp_c_gr_real, temp_c_gr_imag;
	for (int i = 0; i < num_zeros; i++) {
		for (int j = i; j < num_zeros; j++) {
			if (omega_real_interp[j] < omega_real_interp[i]) {
				temp_omega_real = omega_real_interp[j];
				temp_kx_real = kx_real_interp[j];
				temp_kx_imag = kx_imag_interp[j];
				temp_c_gr_real = c_gr_real_interp[j];
				temp_c_gr_imag = c_gr_imag_interp[j];
				for (int k = j; k > i; k--) {
					omega_real_interp[k] = omega_real_interp[k-1];
					kx_real_interp[k] = kx_real_interp[k-1];
					kx_imag_interp[k] = kx_imag_interp[k-1];
					c_gr_real_interp[k] = c_gr_real_interp[k-1];
					c_gr_imag_interp[k] = c_gr_imag_interp[k-1];
				}
				omega_real_interp[i] = temp_omega_real;
				kx_real_interp[i] = temp_kx_real;
				kx_imag_interp[i] = temp_kx_imag;
				c_gr_real_interp[i] = temp_c_gr_real;
				c_gr_imag_interp[i] = temp_c_gr_imag;
			}
		}
	}
	
	free_1d_array_double(omega_real_interp_H);
	free_1d_array_double(kx_real_interp_H);
	free_1d_array_double(kx_imag_interp_H);
	free_1d_array_double(c_gr_real_interp_H);
	free_1d_array_double(c_gr_imag_interp_H);
	
	return num_zeros;
}

//----------------------------------------------------------------------------
void Dispersion_Relationship_Mapping_t::interpolate_along_spline(int N_spline_points, int num_zeros, int mode, 
																 double *omega_real_interp, double *kx_real_interp, double *kx_imag_interp, double *c_gr_real_interp, double *c_gr_imag_interp,
																 double *omega_real_spline, double *kx_real_spline, double *kx_imag_spline, double *c_gr_real_spline, double *c_gr_imag_spline,
																 double x_pos)
{
	double min_kx_imag = BIG_NUMBER, associated_kx_real = 0.0, associated_omega_real = 0.0;

	for (int point=0; point < N_spline_points; point++ ) {
		omega_real_spline[point] = 0.0;
		c_gr_real_spline[point] = 0.0;
		c_gr_imag_spline[point] = 0.0;
		kx_real_spline[point] = 0.0;
		kx_imag_spline[point] = 0.0;
	}

	Cubic_Spline_t c_gr_real_cubic_spline(num_zeros);
	c_gr_real_cubic_spline.calculate_interpolation_function(omega_real_interp, c_gr_real_interp);
	
	Cubic_Spline_t c_gr_imag_cubic_spline(num_zeros);	
	c_gr_imag_cubic_spline.calculate_interpolation_function(omega_real_interp, c_gr_imag_interp);
	
	Cubic_Spline_t kx_real_cubic_spline(num_zeros);
	kx_real_cubic_spline.calculate_interpolation_function(omega_real_interp, kx_real_interp);
	
	Cubic_Spline_t kx_imag_cubic_spline(num_zeros);	
	kx_imag_cubic_spline.calculate_interpolation_function(omega_real_interp, kx_imag_interp);
	
	for (int ii = 0; ii < N_spline_points; ii++) {
		omega_real_spline[ii] = omega_real_interp[0] + ii * (omega_real_interp[num_zeros-1] - omega_real_interp[0]) / (N_spline_points-1);
		c_gr_real_spline[ii] = c_gr_real_cubic_spline.evaluate_function(omega_real_spline[ii]);
		c_gr_imag_spline[ii] = c_gr_imag_cubic_spline.evaluate_function(omega_real_spline[ii]);
		kx_real_spline[ii] = kx_real_cubic_spline.evaluate_function(omega_real_spline[ii]);
		kx_imag_spline[ii] = kx_imag_cubic_spline.evaluate_function(omega_real_spline[ii]);
		if (kx_imag_spline[ii] < min_kx_imag) {
			min_kx_imag = kx_imag_spline[ii];
			associated_kx_real = kx_real_spline[ii];
			associated_omega_real = omega_real_spline[ii];
		}
	}
	
	string filename;
	filename = get_real_omega_max_spatial_growth_filename(mode+1);
	FILE *fout=fopen(filename.c_str(),"w");
	fprintf(fout,"# x_pos, mode, kx_real, MIN kx_imag, Omega_real, Omega_imag=0, kz_real, kz_imag\n");
	fprintf(fout,"%g %d %g %g %g %g %g %g\n", 
			x_pos, mode, associated_kx_real, min_kx_imag, associated_omega_real, 0.0, real(kz(0,0)), imag(kz(0,0)));
	fclose(fout);
}

//----------------------------------------------------------------------------
string Dispersion_Relationship_Mapping_t::get_padded_number_string(int number)
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
string Dispersion_Relationship_Mapping_t::get_eigenvalue_filename(int kx_real_num, int kx_imag_num) 
{	
	string eigenvalue_filename;
	string output_location("./results/");
	string name("eigen_values");
	string kx_real_number_string;
	string kx_imag_number_string;
	string extension(".dat");
	
	kx_real_number_string = get_padded_number_string(kx_real_num+1);
	kx_imag_number_string = get_padded_number_string(kx_imag_num+1);	
	eigenvalue_filename = output_location + "kx_" + kx_real_number_string + "_" + kx_imag_number_string + "." + name + extension;
	
	return eigenvalue_filename;
}

//----------------------------------------------------------------------------
string Dispersion_Relationship_Mapping_t::get_complex_map_filename(int number) 
{	
	string filename;
	string name("./results/eigen_mapping_sorted.mode_");
	string padded_number_string;
	string extension(".dat");
	
	padded_number_string = get_padded_number_string(number);
	filename = name + padded_number_string + extension;
	
	return filename;
}

//----------------------------------------------------------------------------
string Dispersion_Relationship_Mapping_t::get_real_omega_raw_filename(int number) 
{	
	string filename;
	string name("./results/eigen_mapping_sorted.mode_");
	string padded_number_string;
	string extension(".real_omega_raw.dat");
	
	padded_number_string = get_padded_number_string(number);
	filename = name + padded_number_string + extension;
	
	return filename;
}

//----------------------------------------------------------------------------
string Dispersion_Relationship_Mapping_t::get_real_omega_filename(int number) 
{	
	string filename;
	string name("./results/eigen_mapping_sorted.mode_");
	string padded_number_string;
	string extension(".real_omega.dat");
	
	padded_number_string = get_padded_number_string(number);
	filename = name + padded_number_string + extension;
	
	return filename;
}

//----------------------------------------------------------------------------
string Dispersion_Relationship_Mapping_t::get_real_omega_max_spatial_growth_filename(int number) 
{	
	string filename;
	string name("./results/eigen_mapping_sorted.mode_");
	string padded_number_string;
	string extension(".real_omega.max_spatial_growth.dat");
	
	padded_number_string = get_padded_number_string(number);
	filename = name + padded_number_string + extension;
	
	return filename;
}

//----------------------------------------------------------------------------
string Dispersion_Relationship_Mapping_t::get_pinch_point_filename(int number) 
{	
	string filename;
	string name("./results/eigen_mapping_sorted.mode_");
	string padded_number_string;
	string extension(".pinch_point.dat");
	
	padded_number_string = get_padded_number_string(number);
	filename = name + padded_number_string + extension;
	
	return filename;
}

//----------------------------------------------------------------------------
string Dispersion_Relationship_Mapping_t::get_max_omega_imag_for_kx_real_filename(int number) 
{	
	string filename;
	string name("./results/eigen_mapping_sorted.mode_");
	string padded_number_string;
	string extension(".max_omega_imag_for_kx_real.dat");
	
	padded_number_string = get_padded_number_string(number);
	filename = name + padded_number_string + extension;
	
	return filename;
}


//============================================================================
// Non-class functions
//============================================================================

//----------------------------------------------------------------------------
Eigen_Values_t* malloc_1d_array_Eigen_Values_t(int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	Eigen_Values_t *eigen_values = (Eigen_Values_t *) malloc(sizeof(Eigen_Values_t) * Nx);
	if( !eigen_values ) quit_error((char*)"Malloc failed in malloc_Eigen_Values_t()");
	for(int i=0 ; i<Nx ; i++) eigen_values[i].sorted = false;
	return eigen_values;
}

//----------------------------------------------------------------------------
Eigen_Values_t** malloc_2d_array_Eigen_Values_t(int Nx, int Ny)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	Eigen_Values_t **eigen_values = (Eigen_Values_t **) malloc(sizeof(Eigen_Values_t*) * Nx);
	if( ! eigen_values ) quit_error((char*)"Malloc failed in malloc_2d_array_Eigen_Values_t()");
	for(int i=0 ; i<Nx ; i++) eigen_values[i] = malloc_1d_array_Eigen_Values_t(Ny);
	return eigen_values;
}

//----------------------------------------------------------------------------
void free_1d_array_Eigen_Values_t(Eigen_Values_t *eigen_values)
{
	if(!eigen_values)
		quit_error((char*)"Memory already deallocated in free_Eigen_Values_t()");
	free(eigen_values);
}

//----------------------------------------------------------------------------
void free_2d_array_Eigen_Values_t(Eigen_Values_t **eigen_values, int Nx)
{
	if(!eigen_values) quit_error((char*)"Memory already deallocated in free_2d_array_Eigen_Values_t()");
	for(int i=0;i<Nx;i++) free_1d_array_Eigen_Values_t(eigen_values[i]);
	free(eigen_values);
}

//----------------------------------------------------------------------------
void copy_2d_array_Eigen_Values_t_from_A_to_B(Eigen_Values_t **eigen_values_A, int A1, int A2, Eigen_Values_t **eigen_values_B, int B1, int B2)
{
	eigen_values_B[B1][B2].kx_r =						eigen_values_A[A1][A2].kx_r;
	eigen_values_B[B1][B2].kx_i =						eigen_values_A[A1][A2].kx_i;
	eigen_values_B[B1][B2].kz_r =						eigen_values_A[A1][A2].kz_r;
	eigen_values_B[B1][B2].kz_i =						eigen_values_A[A1][A2].kz_i;
	eigen_values_B[B1][B2].Omega_r =					eigen_values_A[A1][A2].Omega_r;
	eigen_values_B[B1][B2].Omega_i =					eigen_values_A[A1][A2].Omega_i;
	eigen_values_B[B1][B2].c_r =						eigen_values_A[A1][A2].c_r;
	eigen_values_B[B1][B2].c_i =						eigen_values_A[A1][A2].c_i;
	eigen_values_B[B1][B2].discrete_continuous =		eigen_values_A[A1][A2].discrete_continuous;
	eigen_values_B[B1][B2].w_mode =						eigen_values_A[A1][A2].w_mode;
	eigen_values_B[B1][B2].symmetric =					eigen_values_A[A1][A2].symmetric;
	eigen_values_B[B1][B2].rel_error_bounds =			eigen_values_A[A1][A2].rel_error_bounds;
	eigen_values_B[B1][B2].cond =						eigen_values_A[A1][A2].cond;
	eigen_values_B[B1][B2].sorted =						eigen_values_A[A1][A2].sorted;
}

//----------------------------------------------------------------------------
