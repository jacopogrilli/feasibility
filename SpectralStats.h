extern double Max_Eigenvalue_RealPart( gsl_matrix * MatrixFix );

extern double Min_Eigenvalue_RealPart( gsl_matrix * MatrixFix );					 

extern double Min_Eigenvalue_Reactivity( gsl_matrix * MatrixFix );

extern void Sym_matrix_eigenvalues( const gsl_matrix * MatrixFix, gsl_vector * eigv );

extern double get_det( gsl_matrix * Matrix );
