extern void PositiveSphericalVector( int n, gsl_rng * r, gsl_vector * v );

extern void PolytopeGeneratorMatrix( gsl_matrix * GeneratorMat, gsl_matrix * Matrix );

extern void StructuralStabilityIntegration(gsl_rng * r , gsl_matrix * Matrix,  double * StructStab, double * StructStabErr );

extern void StructuralStabilityPolytope(gsl_rng * r , gsl_matrix * Matrix,  double * StructStab, double * StructStabErr );

extern void StructuralStability3D( gsl_matrix * Matrix,  double * StructStab );
