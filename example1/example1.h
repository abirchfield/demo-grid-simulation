#include<stdio.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>

#include <sundials/sundials_core.h>
#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_dense.h>

int example1();
int compress_matrix(SUNMatrix A, int n, int nz, int* Ti, int* Tj, double* Tx);