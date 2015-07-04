#include <iostream>
#include <math.h>

using namespace std;

#include <stdio.h>
#include <stdarg.h>

#include "lapack_mallocs.h"
#include "Utils.h"

//---------------------------------------------------------------------------------------------------
doublecomplex *malloc_1d_array_doublecomplex(integer Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	doublecomplex *f = (doublecomplex *) malloc(sizeof(doublecomplex) * Nx);

	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_array_doublecomplex()");

	return f;
}

//---------------------------------------------------------------------------------------------------
doublereal *malloc_1d_array_doublereal(integer Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	doublereal *f = (doublereal *) malloc(sizeof(doublereal) * Nx);

	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_array_doublereal()");

	return f;
}

//---------------------------------------------------------------------------------------------------
integer *malloc_1d_array_integer(integer Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	integer *f = (integer *) malloc(sizeof(integer) * Nx);

	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_array_integer()");

	return f;
}

//---------------------------------------------------------------------------------------------------
logical *malloc_1d_array_logical(integer Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	logical *f = (logical *) malloc(sizeof(logical) * Nx);

	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_array_logical()");

	return f;
}

//---------------------------------------------------------------------------------------------------
void free_1d_array_doublereal(doublereal *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_doublereal()");
	free(f);

	return;
}

//---------------------------------------------------------------------------------------------------
void free_1d_array_doublecomplex(doublecomplex *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_doublecomplex()");
	free(f);

	return;
}

//---------------------------------------------------------------------------------------------------
void free_1d_array_integer(integer *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_integer()");
	free(f);

	return;
}

//---------------------------------------------------------------------------------------------------
void free_1d_array_logical(logical *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_logical()");
	free(f);

	return;
}

//---------------------------------------------------------------------------------------------------
