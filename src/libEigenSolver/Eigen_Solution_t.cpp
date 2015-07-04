
#include <iostream>
#include <math.h>

using namespace std;

#include <stdio.h>
#include <stdarg.h>

#include <blitz/array.h>
using namespace blitz;

#include "Utils.h"
#include "lapack_mallocs.h"
#include "Eigen_Solution_t.h"

//----------------------------------------------------------------------------
Eigen_Solution_t* malloc_Eigen_Solution_t(long int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	Eigen_Solution_t *eigen_solution = (Eigen_Solution_t *) malloc(sizeof(Eigen_Solution_t) * Nx);
	
	if( !eigen_solution ) quit_error((char*)"Malloc failed in malloc_Eigen_Solution_t()");
	
	return eigen_solution;
}

//----------------------------------------------------------------------------
void free_Eigen_Solution_t(Eigen_Solution_t *eigen_solution, long int Nx)
{
	if(!eigen_solution)
		quit_error((char*)"Memory already deallocated in free_Eigen_Solution_t()");
	
	free(eigen_solution);
}

//----------------------------------------------------------------------------
