#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#define HOWMANYPROPERTIES 10


double Max_Eigenvalue_RealPart( gsl_matrix * MatrixFix ){ // returns the eigenvalue with the largest real part

	int S = MatrixFix->size1;

    gsl_matrix * Matrix = gsl_matrix_calloc(S, S);
    gsl_matrix_memcpy(Matrix, MatrixFix);
    
    double largest_eigenvalue;
    
	gsl_vector_complex * eval = gsl_vector_complex_calloc(S); // store eigenvalues
	gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc(S);
    gsl_eigen_nonsymm (Matrix, eval, w);
//	double * largest = malloc(3 * sizeof(double));
//	gsl_sort_largest(largest, 3, eval->data, 1, S);

    gsl_vector * real_eval = gsl_vector_calloc(S);
    int i;
    for (i = 0; i < S; i += 1)
    {
        gsl_vector_set ( real_eval , i, GSL_REAL( gsl_vector_complex_get (eval, i) ) );
    }

    largest_eigenvalue = gsl_vector_max( real_eval );

	gsl_vector_complex_free (eval);
	gsl_vector_free (real_eval);
	gsl_eigen_nonsymm_free (w);
	gsl_matrix_free(Matrix);

    return largest_eigenvalue;

}


double Min_Eigenvalue_RealPart( gsl_matrix * MatrixFix ){ // returns the eigenvalue with the smallest real part

	int S = MatrixFix->size1;

    gsl_matrix * Matrix = gsl_matrix_calloc(S, S);
    gsl_matrix_memcpy(Matrix, MatrixFix);
    
    double smallest_eigenvalue;
    
	gsl_vector_complex * eval = gsl_vector_complex_calloc(S); // store eigenvalues
	gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc(S);
    gsl_eigen_nonsymm (Matrix, eval, w);
//	double * largest = malloc(3 * sizeof(double));
//	gsl_sort_largest(largest, 3, eval->data, 1, S);

    gsl_vector * real_eval = gsl_vector_calloc(S);
    int i;
    for (i = 0; i < S; i += 1)
    {
        gsl_vector_set ( real_eval , i, GSL_REAL( gsl_vector_complex_get (eval, i) ) );
    }

    smallest_eigenvalue = gsl_vector_min( real_eval );

	gsl_vector_complex_free (eval);
	gsl_vector_free (real_eval);
	gsl_eigen_nonsymm_free (w);
	gsl_matrix_free(Matrix);

    return smallest_eigenvalue;

}


double Min_Eigenvalue_Reactivity( gsl_matrix * MatrixFix ){ // returns the min eigenvalue of the symmetric part of a matrix

	int S = MatrixFix->size1;

    gsl_matrix * Matrix = gsl_matrix_calloc(S, S);
    gsl_matrix_memcpy(Matrix, MatrixFix);
    
    gsl_matrix * MatrixTranspose = gsl_matrix_calloc(S, S);
    gsl_matrix_transpose_memcpy(MatrixTranspose, MatrixFix);

    gsl_matrix_add( Matrix, MatrixTranspose );

	gsl_matrix_free(MatrixTranspose);


    double smallest_eigenvalue;
    
	gsl_vector * eval = gsl_vector_calloc(S); // store eigenvalues
	gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc(S);
    gsl_eigen_symm(Matrix, eval, w);


    smallest_eigenvalue = gsl_vector_min( eval );

	gsl_vector_free (eval);
	gsl_eigen_symm_free (w);
	gsl_matrix_free(Matrix);

    return smallest_eigenvalue * 0.5;

}

void Sym_matrix_eigenvalues( const gsl_matrix * MatrixFix, gsl_vector * eigv ){

	int S = MatrixFix->size1;

    gsl_matrix * Matrix = gsl_matrix_calloc(S, S);
    gsl_matrix_memcpy(Matrix, MatrixFix);
    
    
	gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc(S);
    gsl_eigen_symm(Matrix, eigv, w);



	gsl_eigen_symm_free (w);
	gsl_matrix_free(Matrix);

    return;

}

double get_det(gsl_matrix *A) { // return the determinant of a matrix

    double det;
    int signum;
    gsl_permutation *p = gsl_permutation_alloc(A->size1);

    gsl_matrix *tmpA = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(tmpA , A);


    gsl_linalg_LU_decomp(tmpA , p , &signum);
    det = gsl_linalg_LU_det(tmpA , signum);
    gsl_permutation_free(p);
    gsl_matrix_free(tmpA);

    return det;

}



