
typedef long int integer;
typedef double doublereal;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;

doublecomplex *malloc_1d_array_doublecomplex(integer Nx);
doublereal *malloc_1d_array_doublereal(integer Nx);
integer *malloc_1d_array_integer(integer Nx);
logical *malloc_1d_array_logical(integer Nx);

void free_1d_array_doublereal(doublereal *f);
void free_1d_array_doublecomplex(doublecomplex *f);
void free_1d_array_integer(integer *f);
void free_1d_array_logical(logical *f);


