/**
 * File owner : Said Guouihaj
 * Description : 
 * * This is my C implementation of EM-GMM algorithm
 * * The program generate random test data return The following parameters :
 * * * Mixture coefficent : a list of K scalars
 * * * Means vectors : a list of dimension D vectors
 * * * Covariance matrices : a list K of diagonal matrices represented with their diagonal vectors of dimension D
 * * The output reprented as a standard output and stored to a file model.txt
 * * The data generated stored in points.csv
 * * This is a sequencial program
**/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define M_PI		3.14159265358979323846	/* pi */

#define K   5   /* The number of clusters */
#define D   3   /* The dimension of data */
#define N   50  /* The number of test vectors generated */

//global vectors
double w[K];    /* The mixture weights of the K clusters */
double mu[K][D];    /* The mean vectors of the K cluster with D dimension*/
double theta[K][D]; /* The diagonal of Covariance matrices */
double x[N][D]; /* Input test vectors */

//intermediate vectors
double R[N][K]; /* The list of responsibility values (Calculated at E-step) */

/* function prototypes */
double drand ( double low, double high);
void initParams(void);
double G(int n, int k);
double r(int n, int k);
void calcR(void);
double calcM(int k);
void calcParams(void);
void printParams(void);
double calcL(void);
void storePoints(void);
void storeParams(void);

void main(){
    initParams();               
    bool isConverges = false;
    double L, prevL;
    prevL = calcL();
    while (!isConverges){
        calcR();
        calcParams();
        L = calcL();
        if ((double)-0.000001<L-prevL && L-prevL<(double)0.000001) /* convergence precision */
            isConverges = true;
        prevL = L;
    }
    printParams();
    storePoints();
    storeParams();
}

/// @brief Random double generator 
/// @param low : Minimum value
/// @param high : Maximum value
/// @return 
double drand ( double low, double high )
{
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

/// @brief initalize parameters w, mu, theta
void initParams(){
    for (int k=0;k<K;k++){
        for (int d=0;d<D;d++){
            mu[k][d] = (double)drand(0.1,5.0);  /* Mean random range*/
            theta[k][d] = (double)drand(5.0, 20.0); /* Variance random range*/
        }
        
        w[k] = 1.0/K;
    }
    for (int n=0;n<N;n++){
        for (int d=0;d<D;d++)
            x[n][d] = (double)drand(0.1, 0.4);
    }
}
/// @brief Gaussian distribution formula
/// @param n : input vector index 
/// @param k : cluster index
/// @return Gaussian desity result
double G(int n, int k){
    double result;
    double thetaInv[D];
    double xn_k[D];
    double exponent = 0.0;
    double detTheta = 1.0;
    for (int d=0;d<D;d++){
        thetaInv[d] = 1/theta[k][d];
        xn_k[d] = x[n][d] - mu[k][d];
        exponent += xn_k[d]*xn_k[d]*thetaInv[d];
        detTheta *= theta[k][d];
    }
    result = exp((-0.5*exponent)/(sqrt(2*M_PI*detTheta)));
    return result;
}

/// @brief : Calculate the responsibility value
/// @param n : input vector index
/// @param k : index cluster index
/// @return : The responsibility value
double r(int n, int k){
    double result = 0.0;
    for (int kt=0;kt<K;kt++)
        result += w[kt]*G(n, kt);
    result = w[k]*G(n,k)/result;
    return result;
}


/// @brief : Calculate the responsibilities of all values to all clusters (E-step)
void calcR(){
    for (int n=0;n<N;n++){
        for (int k=0;k<K;k++){
            R[n][k] = r(n,k);
        }
    }
}


/// @brief : Calculate m value (intermdiste value) of the cluster k
/// @param k : the cluster index
/// @return : the m value of the cluster k
double calcM(int k){
    double result = 0.0;
    for (int n=0;n<N;n++){
        result += R[n][k];
    }
    return result;
}

/// @brief Evaluate new parmeters (M-Step)
void calcParams(){
    for (int k=0;k<K;k++){
        double m = calcM(k);
        w[k] = m/N;
        for (int d=0;d<D;d++){
            mu[k][d] = 0.0;
            theta[k][d] = 0.0;
            for (int n=0;n<N;n++){
                mu[k][d] += R[n][k]*x[n][d];
                theta[k][d] += R[n][k]*(x[n][d]-mu[k][d])*(x[n][d]-mu[k][d]);
            }
            mu[k][d] /= m;
            theta[k][d] /= m;
        }
    }
}

/// @brief Display parameters
void printParams(){
    double sum = 0.0;
    for (int k=0;k<K;k++){
        printf("w[%d] = %lf\n", k, w[k]);
        sum += w[k];
    }
    printf("\n");
    for (int k=0;k<1;k++){
        for (int d=0;d<D;d++){
            printf("mu[%d][%d] = %lf\n",k, d, mu[k][d]);
        }
        printf("\n");
    }
    printf("\n");
    for (int k=0;k<1;k++){
        for (int d=0;d<D;d++){
            printf("theta[%d][%d] = %lf\n", k, d, theta[k][d]);
        }
        printf("\n");
    }
}

/// @brief Calculate the log-likelihood (The convergence criteria)
/// @return the log-likelihood
double calcL(){
    double L = 0.0;
    for (int n=0;n<N;n++){
        double subSum = 0.0;
        for (int k=0;k<K;k++){
            subSum += G(n,k)*w[k];
        }
        L += log(subSum);
    }
    return L;
}

/// @brief Store points in a file (points.csv)
void storePoints(){
    FILE * points = fopen("points.csv", "w");
    fprintf(points, "x,y\n");
    for (int n = 0; n < N; n++){
       for (int d=0;d<D;d++){
           fprintf(points, "%lf", x[n][d]);
            if (d != D-1)
                fprintf(points, ",");
       }    
        fprintf(points, "\n");
    }
    fclose(points);
}

/// @brief Store model parameters in  a file
void storeParams(){
    // w[k], mu[k], sigma2[k]
    FILE * gnuplot = fopen("model.txt", "w");
    fprintf(gnuplot, "%d %d\n", K, D);
    for (int k=0;k<K;k++){
        fprintf(gnuplot, "%lf ", w[k]);
        for (int d=0;d<D;d++){
           fprintf(gnuplot, "%lf ", mu[k][d]);
        } 
        for (int d=0;d<D;d++){
           fprintf(gnuplot, "%lf ", theta[k][d]);
        } 
        fprintf(gnuplot, "\n");
    }
    fclose(gnuplot);
}
