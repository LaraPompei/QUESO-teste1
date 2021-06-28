#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <queso/GslMatrix.h>

using namespace QUESO;

struct likelihoodRoutine_DataType{
    const GslVector* meanVector;
    const GslMatrix* covMatrix;
};

double likelihoodRountine(const GslVector& paramValues, const GslVector* paramDirection, const void* functionDataPtr,  GslVector* gradVector, GslMatrix* hessianMatrix, GslVector* hessianEffect);

#endif
