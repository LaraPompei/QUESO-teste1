#include "compute.h"
#include "likelihood.h"

#define tam 12

void filling_matrix(int& number_samples, double* t, double* values_a, double** data){
    mt19937 engine; // uniform random bit engine

    // seed the URBG
    random_device dev{};
    engine.seed(dev());

    // setup a distribution:
    double mu    = 0.5;
    double sigma = 2.0;
    lognormal_distribution<double> dist(mu, sigma);
    for (int i =0; i<number_samples; i++) {
        values_a[i] = dist(engine);
    }
    for(int i=0; i<number_samples; i++){
        double* lane = new double[tam];
        for(int j=0;j<tam; j++){
            lane[j] = exp(-values_a[i]*t[j]);
            data[i][j] = lane[j];
        }
    }
}

void compute(const FullEnvironment& env){
     //SIP:
  
    //instantiating the parameter space
    
    VectorSpace<> paramSpace(env, "param_", 1, NULL);
    
    //instantiating the parameter domain
    
    GslVector paramMinValues(paramSpace.zeroVector());
    GslVector paramMaxValues(paramSpace.zeroVector());
    paramMinValues[0] = 0.0;
    paramMaxValues[0] = 6.0;
    BoxSubset<> paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);
    
    //instantiating the likelihood function object

    double spacing_points = 0.5;
    int number_samples    = 10000;
    double* t             = new double[tam];
    double* values_a      = new double[number_samples];
    double** data         = new double*[number_samples];
    double* data_mean     = new double[number_samples]{0};
    double* data_std      = new double[number_samples]{0};

    //use poi = 1 if t = 0.5, 2 if t = 1.0; 3 if t = 2.5, 4 if t=2.0,...
    int poi = 3; //point of interest

    //memory allocation for the matrix data [number_samplesxtam]
    for(int i = 0; i<number_samples; i++){
        data[i] = new double[tam];
    }

    //vector t is filled with [0.0,5.5] interval catching each number after 0.5
    for(int i =0; i<tam;i++){
        t[i] = (i)*0.5;
    }

    //Generating and fillling the matrix, which the lanes represent the samples(yvalues) and the columns represent the date in each point of time(t)
    filling_matrix(number_samples, t, values_a, data);
    //defining t value in which the function will be evaluated
    double t_ = t[poi];

    //mean of data
    for(int j=0; j < tam; j++){
        for(int i = 0; i < tam; i++){
            data_mean[j] += data[i][j];
        }
    }
    for(int i = 0; i < tam; i++){
        data_mean[i] = data_mean[i]/number_samples;
    }

    //data standard deviation
    for(int j = 0; j < tam; j++){
        for(int i = 0; i<number_samples; i++){
            data_std[j] += pow((data[i][j] - data_mean[j]),2);
        }
    }
    for(int i = 0; i <tam; i++){
        data_std[i] = sqrt(data_std[i]/number_samples);
        cout<<"Para t = "<<i<<endl<<"Data mean: "<<data_mean[i]<<endl<<"Data standard deviation: "<<data_std[i]<<endl;
    }

    Likelihood<> lhood("like_", paramDomain, data[poi], values_a, data_std);

    //define the prior RV
    UniformVectorRV<> priorRv("prior_", paramDomain);

    //instantiate the inverse problem
    GenericVectorRV<> postRv("post_", paramSpace);
    StatisticalInverseProblem<> ip("", NULL, priorRv, lhood, postRv);

    //Solve the inverse problem
    //set the 'pdf' and the 'realizer' of the posterior RV
    GslVector paramInitials(paramSpace.zeroVector());
    priorRv.realizer().realization(paramInitials);

    GslMatrix proposalCovMatrix(paramSpace.zeroVector());
    proposalCovMatrix(0,0) = std::pow(std::abs(paramInitials[0])/20.0, 2.0);

    ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);
    
}
