///////////////////////////////////////////
// StructStab ///////////////////////////////////
// Jacopo Grilli jgrilli@uchicago.edu
// Aug/Set 2014
///////////////////////////////////////////

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_inline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include "StructStability.h"
#include "SpectralStats.h"


#define GSL_MAX(a, b)   ((a) > (b) ? (a) : (b))

int main(int argc, char *argv[]){

    int S = atoi(argv[1]); // number of species
    char * FileInput = argv[2]; // file input
    char * FileOut = argv[3]; // file output
    
    int RndSeed = atoi(argv[4]); // rnd seed (-1 use time machine)
    
    // Allocate input matrix
    gsl_matrix * IntMat = gsl_matrix_calloc(S, S);


    // read matrix
    FILE * F;
    F = fopen(FileInput, "rb");
    gsl_matrix_fscanf(F, IntMat );
    fclose(F);


    // output file
    FILE * Fout;
    Fout = fopen(FileOut, "a"); // append results

    
    // Random number generator
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    if( RndSeed == -1 ) gsl_rng_set(r, time(NULL) );
    else gsl_rng_set(r, RndSeed );


    // calculate eigenvalue
    double eval_emp = Min_Eigenvalue_RealPart( IntMat ); // compute eigenvalue
    double reactivity_eval = Min_Eigenvalue_Reactivity( IntMat ); // compute reactivity


    double StructStab, StructStabErr;


    if( reactivity_eval > 0. ){ // Check

        StructuralStabilityPolytope( r, IntMat,  &StructStab, &StructStabErr ); // compute size of feasibility domain
        fprintf( Fout, "%f %f %.6e %.6e\n", eval_emp, reactivity_eval, StructStab, StructStabErr ); // print the output
        // output: dominant eigenvalue, largest eigenvalue of the symmetric part, estimated size of feasibility domain, error on the estimation
    }

    fclose(Fout);


    gsl_matrix_free(IntMat);
    
    return 0;
    
}



