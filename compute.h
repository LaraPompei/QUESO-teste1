#ifndef COMPUTE_H
#define COMPUTE_H

#include <queso/Environment.h>
#include <cmath>
#include <queso/GslMatrix.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>

#include <random>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <cmath>

using namespace QUESO;
using namespace std;

void compute(const FullEnvironment& env);
void filling_matrix(int& number_samples, double* t, double* values_a, double** data);

#endif
