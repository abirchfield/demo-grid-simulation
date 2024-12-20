#include "example1.h"

int compress_matrix(SUNMatrix A, int n, int nz, int* Ti, int* Tj, double* Tx)
{
	/* These algorithms come from CSparse (Tim Davis's book) */

	/* Allocate variables, temporary matrix arrays, and workspace */
	int i,j,k,p,q,cnz,retval;
	int* Wc = (int*)malloc(n*sizeof(int));
	int* Wp = (int*)malloc((n + 1)*sizeof(int));
	int* Wi = (int*)malloc(nz*sizeof(int));
	double* Wx = (double*)malloc(nz*sizeof(double));
	if (Wp == NULL || Wi == NULL || Wx == NULL || Wc == NULL) return 200;

	/* Compress the transpose of the matrix */
	for (i = 0; i < n; ++i) Wc[i] = 0;
	for (k = 0; k < nz; ++k) Wc[Ti[k]]++;
	cnz = 0;
	for (i = 0; i < n; ++i)
	{
		Wp[i] = cnz;
		cnz += Wc[i];
		Wc[i] = Wp[i];
	}
	Wp[n] = cnz;
	for (k = 0; k < nz; ++k)
	{
		p = Wc[Ti[k]]++;
		Wi[p] = Tj[k];
		Wx[p] = Tx[k];
	}

	/* Sum up duplicate entries */
	cnz = 0;
	for (i = 0; i < n; ++i) Wc[i] = -1;
	for (j = 0; j < n; ++j)
	{
		q = cnz;
		for (p = Wp[j]; p < Wp[j + 1]; p++)
		{
			i = Wi[p];
			if (Wc[i] >= q)
			{
				Wx[Wc[i]] += Wx[p];
			}
			else
			{
				Wc[i] = cnz;
				Wi[cnz] = i;
				Wx[cnz] = Wx[p];
				cnz++;
			}
		}
		Wp[j] = q;
	}
	Wp[n] = cnz;

	/* Allocate resulting matrix */
    retval = SUNSparseMatrix_Reallocate(A, cnz);
    sunindextype *Ap = SUNSparseMatrix_IndexPointers(A);
    sunindextype *Ai = SUNSparseMatrix_IndexValues(A);
    double *Ax = SUNSparseMatrix_Data(A);
	if (retval < 0 || Ap == NULL || Ai == NULL || Ax == NULL) return -1;

	/* Transpose results, sorting the columns */
	for (i = 0; i < n; ++i) Wc[i] = 0;
	for (p = 0; p < Wp[n]; ++p) Wc[Wi[p]]++;
	cnz = 0;
	for (i = 0; i < n; ++i)
	{
		Ap[i] = cnz;
		cnz += Wc[i];
		Wc[i] = Ap[i];
	}
	Ap[n] = cnz;
	for (j = 0; j < n; ++j)
	{
		for (p = Wp[j]; p < Wp[j + 1]; ++p)
		{
			Ai[q = Wc[Wi[p]]++] = j;
			Ax[q] = Wx[p];
		}
	}

	/* Free workspace */
	free(Wc);
	free(Wp);
	free(Wi);
	free(Wx);
	return 0;
}