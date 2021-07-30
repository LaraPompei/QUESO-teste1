#include "compute.h"
#include "likelihood.h"

#define tam 12
#define poi 3
#define number_samples 10000

void filling_matrix(double* t, double* values_a, double** data){
    mt19937 engine; // uniform random bit engine

    // seed the URBG
    random_device dev{};
    engine.seed(dev());

    // setup a distribution:
    double mu    = 0.5; //Mean value of the normal distribution
    double sigma = 0.2; //Standard deviation of the normal distribution
    lognormal_distribution<double> dist(mu, sigma);
    for (int i =0; i<number_samples; i++) {
        values_a[i] = dist(engine);
        cerr<<values_a[i]<<" "<<i<<endl;
    }
    for(int i=0; i<number_samples; i++){
        for(int j=0;j<tam; j++){
            data[i][j] = exp(-values_a[i]*t[j]);
            cerr<<data[i][j]<<" ";
        }
        cerr<<"\t"<<i<<endl<<endl;
    }
}


void compute(const FullEnvironment& env){
    struct timeval timevalNow;

    gettimeofday(&timevalNow, NULL);
    if (env.fullRank() == 0){
        cout<<"\nBeginning run of 'Example 1: Log-Normal Distribution Function (0.5,.2)' example at "<<ctime(&timevalNow.tv_sec)<<endl;
    }  
    env.fullComm().Barrier();
    env.subComm().Barrier();

    //instantiating the parameter space
    cerr<<"instantiating the parameter space.."<<endl;    
    VectorSpace<> paramSpace(env, "param_", 1, NULL);
    
    //defining the domain of a
    cerr<<"instantiating the parameter domain"<<endl;
    GslVector paramMinValues(paramSpace.zeroVector());
    GslVector paramMaxValues(paramSpace.zeroVector());
    paramMinValues[0] = 0.0;
    paramMaxValues[0] = 6.0;
    BoxSubset<> paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);
    
    //instantiating the likelihood function object
    cerr<<"instantiating the likelihood function object and generating the samples"<<endl;
    double spacing_points = 0.5;
    double* t             = new double[tam];
    double* values_a      = new double[number_samples];
    double** data         = new double*[number_samples];
    double data_mean      = 0;
    double data_std       = 0;
    
    //memory allocation for the matrix data [number_samples;tam]
    cerr<<"allocating memory for the data matrix"<<endl;
    for(int i = 0; i<number_samples; i++){
        data[i] = new double[tam];
    }
    
    //vector t is filled with [0.0,5.5] interval catching each number after 0.5
    cerr<<"Filling vector t"<<endl<<"t[ ";
    for(int i =0; i<tam;i++){
        t[i] = (i)*spacing_points;
        cerr<<t[i]<<" ";
    }
    cerr<<"]"<<endl;

    //Generating and fillling the matrix, each element of data is a sample(y(t)=e^(-a*t), each row is a differente t and each column is a different a
    cerr<<"Generating and filling the matrix"<<endl;
    filling_matrix(t, values_a, data);

    //mean of the data
    cerr<<"Calculating the mean of the data"<<endl;
    cerr<<"data[:][poi] = [ ";
    for(int i = 0; i < number_samples; i++){
        data_mean += data[i][poi];
        cerr<<data[i][poi]<<" ";
    }
    cerr<<" ]"<<endl;
    data_mean = data_mean/number_samples;
    cerr<<"Data_mean = "<<data_mean<<endl;
    
    //data standard deviation
    cerr<<"Calculating the standard deviation"<<endl;
    for(int i = 0; i<number_samples; i++){
        data_std += pow((data[i][poi] - data_mean),2);
    }
    data_std = sqrt(data_std/number_samples);

    cerr<<"Para t = "<<poi<<endl<<"Data mean: "<<data_mean<<endl<<"Data standard deviation: "<<data_std<<endl;

    cerr<<"Creating the likelihood object"<<endl;
    Likelihood<> lhood("like_", paramDomain, &data_mean, t, &data_std, poi, tam);
    //define the prior RV
    cerr<<"Defining the prior RV"<<endl;
    UniformVectorRV<> priorRv("prior_", paramDomain);

    //instantiate the inverse problem
    cerr<<"Instantiating the inverse problem"<<endl;
    GenericVectorRV<> postRv("post_", paramSpace);
    StatisticalInverseProblem<> ip("", NULL, priorRv, lhood, postRv);

    //Solve the inverse problem
    cerr<<"Solving the inverse problem"<<endl;
    GslVector paramInitials(paramSpace.zeroVector());
    priorRv.realizer().realization(paramInitials);
    
    GslMatrix proposalCovMatrix(paramSpace.zeroVector());
    
    proposalCovMatrix(0,0) = pow(abs(paramInitials[0])/20.0, 2.0);
    cerr<<"a = [";
    ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);
    cerr<<"]"<<endl;
    if (env.fullRank() == 0) {
        cout << "Ending run of 'Example 1: Log-Normal Distribution Function (0.5,.2)' example at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
    }
    cerr<<"values of a = [";
    for(int i =0 ;i<number_samples ; i++){
        cerr<<values_a[i]<<" ";
    }
    cerr<<"]"<<endl;
    save_data(NULL,values_a,NULL);
}

void save_data(double* model, double* data, double* values_of_a){
    fstream a_data, a_mcmc;
    char fileName[100];

    sprintf(fileName, "a_data.m");
    a_data.open(fileName, ios_base::out);
    a_data<<"a_data = [";
    for(int i=0 ; i<number_samples; i++){
        a_data<<data[i]<<endl;
    }
    a_data<<"];"<<endl;
    a_data.close();
}
