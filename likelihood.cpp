#include "likelihood.h"

//Construtor
template<class V, class M>
Likelihood<V,M>::Likelihood(const char* prefix, const VectorSet<V, M>& domain, const double* m_data_mean, const double* m_t, const double* m_srdDevs, const double& poi):BaseScalarFunction<V,M>(prefix,domain), data_mean(0), values_a(0), stdDevs(0){
    size_t const size = sizeof(m_data_mean)/sizeof(*m_data_mean);
    m_data_mean.assign(data_mean,data_mean+size);
    m_t.assign(t,t+n);
    m_stdDevs.assign(stdDevs,stdDevs+size);
    this->poi = poi;
}

template<class V, class M>
double Likelihood<V, M>::lnValue(const V& domainVector) const{
    double a = domainVector[0];

    double misfitValue = 0.0;    
    for (unsigned int i = 0; i < m_data_mean.size(); ++i) {
        double model = exp(-a*t[poi]);
        double error_A = abs((model - m_data_mean[i])/m_stdDevs[i]);
        misfitValue += error_A * error_A;
    }

    return -0.5*misfitValue;
}


template<class V, class M>
Likelihood<V,M>::~Likelihood(){}

template<class V, class M>
double Likelihood<V,M>::actualValue(const V & domainVector, const V * domainDirection, V * gradVector, M * hessianMatrix, V * hessianEffect) const{
    return exp(thiss->lnValue(domainVector));
}

template class Likelihood<QUESO::GslVector, QUESO::GslMatrix>;
