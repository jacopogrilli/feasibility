#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include "SpectralStats.h"

#define ERROR_THRESHOLD 0.05
// error threshold sets the max ratio between the error on the integration and the value of the integration

void PositiveSphericalVector( int n, gsl_rng * r, gsl_vector * v ){ // generate a random vector with norm 1 in the positive orthant
    double * x = malloc(n * sizeof(double));
    gsl_ran_dir_nd( r, n, x); // generate a random vector with unit norm
    int i;
    for (i = 0; i < n; i += 1) gsl_vector_set( v, i, fabs(x[i]) );
    free( x );
    
    return;
}

void PolytopeGeneratorMatrix( gsl_matrix * GeneratorMat, gsl_matrix * Matrix ){
    // generates the matrix of polytope generators

    int S = GeneratorMat->size1;

    gsl_matrix_memcpy( GeneratorMat, Matrix);
    int i, j;
    double rowsum = 0.; 
    for (i = 0; i < S; i += 1){
        rowsum = 0.;
        for (j = 0; j < S; j += 1){
            rowsum += pow( gsl_matrix_get( Matrix, j,i ), 2.);
        }
        rowsum = sqrt( 1. / rowsum );
        for (j = 0; j < S; j += 1){
            gsl_matrix_set( GeneratorMat, j,i, gsl_matrix_get( Matrix, j,i ) * rowsum );
        }
    }

    return;

}


void StructuralStabilityPolytope(gsl_rng * r , gsl_matrix * Matrix,  double * StructStab, double * StructStabErr ){
    // numerical integration: returns the solid angle (size of the feasibility domain) and an associated error

    int S = Matrix->size1;

    // polytope generator matrix
    gsl_matrix * GeneratorMat = gsl_matrix_calloc(S, S);
    PolytopeGeneratorMatrix( GeneratorMat, Matrix ); // generate the matrix of polytope generators


    double detG = get_det( GeneratorMat ); // calculate the determinant



    double stab = 0, staberr = 0;
    double rs;
    gsl_vector * y = gsl_vector_alloc(S);  
    gsl_vector * u = gsl_vector_alloc(S);  

    int sample = 0;
    double relative_error = 1.;
    while( relative_error > ERROR_THRESHOLD ){ // it stops when reached a given precision
   
        PositiveSphericalVector( S,  r, u ); // generate a random vector on the sphere in the positive cone
        
        gsl_blas_dgemv( CblasNoTrans, 1., GeneratorMat, u, 0., y); // y = Matrix * u
        rs = gsl_blas_dnrm2(y); // calculate the norm of y

        rs = pow(1./rs,S);
        
        stab += rs ;
        staberr += rs * rs ;

        sample++;
        if( sample % 10*S == 0 ){ // every 10*S vectors calculate the error and check the precision

            relative_error = ( ( sample * staberr ) / ( stab * stab ) - 1 ) / ( sample - 1 );
        
            relative_error = sqrt( relative_error ); // relative error

            fprintf( stderr, "%f\n", relative_error); 

        }
                
    }

    stab = exp( log( detG ) + log( stab ) - log( sample ) ); // solid angle


    staberr = exp( 2* log(detG ) + log( staberr ) - log( sample ) ); // error on the solid angle

    *StructStab = stab;
    *StructStabErr = sqrt( staberr - pow(stab,2) ) / sqrt( sample - 1 );

    gsl_vector_free( y );
    gsl_vector_free( u );


    return ;
}



void StructuralStability3D( gsl_matrix * Matrix,  double * StructStab ){
    // if the matrix is 3x3 this function returns the solid angle computed using the analytic formula

    int S = Matrix->size1;

    // polytope generator matrix
    gsl_matrix * GeneratorMat = gsl_matrix_calloc(S, S);
    PolytopeGeneratorMatrix( GeneratorMat, Matrix );


    double detG = get_det( GeneratorMat );

    double scalarprod = 0.; // check this!
    double normg = 0.;
    int i;
    for (i = 0; i < 3; i += 1){
        normg += gsl_matrix_get( GeneratorMat,i, 0) * gsl_matrix_get( GeneratorMat,i, 0);
        scalarprod += gsl_matrix_get( GeneratorMat,i,0) * gsl_matrix_get( GeneratorMat,i,1);
        scalarprod += gsl_matrix_get( GeneratorMat, i,1) * gsl_matrix_get( GeneratorMat,i,2);
        scalarprod += gsl_matrix_get( GeneratorMat, i,2) * gsl_matrix_get( GeneratorMat, i,0);
    }


    if( 1. + scalarprod > 0 ) *StructStab = 4. * atan( fabs(detG) / ( 1. + scalarprod ) ) / M_PI;
    else *StructStab = 4. * ( atan( fabs(detG) / ( 1. + scalarprod ) )  / M_PI + 1.);

    return ;
}







