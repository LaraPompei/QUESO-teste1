#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <cmath>

using namespace std;
using namespace QUESO;

template<class V = QUESO::GslVector, class M = QUESO::GslMatrix>
class Likelihood : public QUESO::BaseScalarFunction<V,M>{
    private:
        vector<double> data_mean; //exponencial model results
        vector<double> t; //values of a
        vector<double> stdDevs; //account for uncertainties
        double size;
        int poi;
        const QUESO::BaseEnvironment* m_env;
   
   public:
        Likelihood(const char* prefix, const VectorSet<V, M>& domain, const double* m_data_mean, const double* m_t, const double* m_srdDevs, const double& poi);
        virtual ~Likelihood();
        virtual double lnValue(const V& domainVector) const;
        virtual double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
};

#endif
