#include <random>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <cmath>

#define tam 12

using namespace std;

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

double model(int& a, double& t, double& data_mean, double& data_std){
    
    //Function
    long double y = exp(-a*t);
    
    //Evaluated for a point, the error was calculated as: value of the function minuus the mean, devided by the standart deviation
    double error_A = abs((y-data_mean)/(data_std)); 
    
    if(isnan(error_A) || isinf(error_A){
        error_A = 1.0e12;
    }
    return error_A;

}

int main() {
    double spacing_points = 0.5;
    int number_samples    = 10000;
    double* t             = new double[tam];
    double* values_a      = new double[number_samples];
    double* values_a_aux  = new double[number_samples];
    double** data         = new double*[number_samples];
    double** data_aux     = new double*[number_samples];
    double data_mean      = 0;
    double data_std       = 0;
    //use poi = 1 if t = 0.5, 2 if t = 1.0; 3 if t = 2.5, 4 if t=2.0,...
    int poi = 3; //point of interest

    //memory allocation for the matrix data [number_samplesxtam]
    for(int i = 0; i<number_samples; i++){
        data[i] = new double[tam];
        data_aux[i] = new double[tam];
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
    for(int i = 0; i < number_samples; i++){
        data_mean += data[i][poi];
    }
    data_mean = data_mean/number_samples;

    //data standard deviation
    for(int i = 0; i<number_samples; i++){
        data_std += pow((data[i][poi] - data_mean),2);
    }
    data_std = sqrt(data_std/number_samples);
    cout<<"Data average: "<<data_mean<<endl<<"Data standard deviation: "<<data_std<<endl;
} 
