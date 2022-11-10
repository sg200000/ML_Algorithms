
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

# define M_PI		3.14159265358979323846	/* pi */

#define K   5
#define D   3
#define N   50

//global vectors
double w[K];
double mu[K][D];
double theta[K][D][D];
double x[N][D];

//intermediate vectors
double R[N][K];

double drand ( double low, double high )
{
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

void testDrand(){
    for (int i=0;i<100;i++){
        printf("%lf\n", drand(5.0,25.0));
    }
}

void initParams(){
    for (int k=0;k<K;k++){
        for (int d=0;d<D;d++){
            mu[k][d] = (double)drand(0.1, 5.0);
            theta[k][d][d] = (double)drand(0.1, 2.0);
        }
        
        w[k] = 1.0/K;
    }
    for (int n=0;n<N;n++){
        for (int d=0;d<D;d++)
            x[n][d] = (double)drand(0.1, 0.4);
        //printf("%lf\n", x[n][0]);
    }
}
double G(int n, int k){
    double result;
    double thetaInv[D];
    double xn_k[D];
    double exponent = 0.0;
    double detTheta = 1.0;
    for (int d=0;d<D;d++){
        thetaInv[d] = 1/theta[k][d][d];
        xn_k[d] = x[n][d] - mu[k][d];
        exponent += xn_k[d]*xn_k[d]*thetaInv[d];
        detTheta *= theta[k][d][d];
        //printf("%lf %lf %lf %lf %lf | ", theta[k][d][d], thetaInv[d], xn_k[d], exponent, detTheta);
        //printf("%lf | ", detTheta);
    }
   // printf("%lf %lf |", sqrt(2*M_PI*fabs(detTheta)), fabs(detTheta));
    result = exp((-0.5*exponent)/(sqrt(2*M_PI*detTheta)));
    //printf("%lf ", result);
    return result;
}

double r(int n, int k){
    double result = 0.0;
    for (int kt=0;kt<K;kt++)
        result += w[kt]*G(n, kt);
    result = w[k]*G(n,k)/result;
    //printf("r[%d][%d] = %lf\n", n, k, result);
    return result;
}

void calcR(){
    for (int n=0;n<N;n++){
        for (int k=0;k<K;k++){
            R[n][k] = r(n,k);
            //printf("r[%d][%d] = %lf | ",n,k,R[n][k]);
        }
        //printf(",\n");
    }
}
double calcM(int k){
    double result = 0.0;
    for (int n=0;n<N;n++){
        result += R[n][k];
    }
    return result;
}

void calcParams(){
    for (int k=0;k<K;k++){
        double m = calcM(k);
       // printf("m[%d] = %lf\n", k,m);
        w[k] = m/N;
        for (int d=0;d<D;d++){
            mu[k][d] = 0.0;
            theta[k][d][d] = 0.0;
            for (int n=0;n<N;n++){
                mu[k][d] += R[n][k]*x[n][d];
                theta[k][d][d] += R[n][k]*(x[n][d]-mu[k][d])*(x[n][d]-mu[k][d]);
            }
            mu[k][d] /= m;
            theta[k][d][d] /= m;
        }
    }
}

void printParams(){
    double sum = 0.0;
    for (int k=0;k<K;k++){
        printf("w[%d] = %lf\n", k, w[k]);
        sum += w[k];
    }
    printf("\n");
    //printf("\nsum = %lf\n\n", sum);
    for (int k=0;k<1;k++){
        //printf("mu[%d] = ", k);
        for (int d=0;d<D;d++){
            printf("mu[%d][%d] = %lf ",k, d, mu[k][d]);
        }
        //printf("\n");
    }
    //printf("\n");
    for (int k=0;k<1;k++){
       // printf("theta[%d] = ", k);
        for (int d=0;d<D;d++){
        //    printf("%lf ", theta[k][d][d]);
        }
        //printf("\n");
    }
}

double calcL(){
    double L = 0.0;
    for (int n=0;n<N;n++){
        double subSum = 0.0;
        for (int k=0;k<K;k++){
            subSum += G(n,k)*w[k];
            //printf("%lf %lf %lf | ", subSum, w[k], G(n,k));
        }
        L += log(subSum);
    }
    return L;
}

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
       // printf("%lf", g);
    }
    fclose(points);
}

void storeParams(){
    // w[k], mu[k], sigma2[k]
    FILE * gnuplot = fopen("model.csv", "w");
    fprintf(gnuplot, "%d %d\n", K, D);
    for (int k=0;k<K;k++){
        fprintf(gnuplot, "%lf ");
        for (int d=0;d<D;d++){
           fprintf(gnuplot, "%lf ", mu[k][d]);
        } 
        for (int d=0;d<D;d++){
           fprintf(gnuplot, "%lf ", theta[k][d][k]);
        } 
        fprintf(gnuplot, "\n");
    }
    fclose(gnuplot);
}
void main(){
    initParams();
    //testDrand();
    //double g = G(0,0); 
   // printf("%lf\n", g);
    bool isConverges = false;
    double L, prevL;
    prevL = calcL();
    int i = 0;
    printParams();
    while (!isConverges){
        calcR();
        calcParams();
        printParams();
        L = calcL();
        if ((double)-0.000001<L-prevL && L-prevL<(double)0.000001)
            isConverges = true;
        prevL = L;
        printf("%lf\n", L);
      //  if (isnan(L))
        //    break;
        i++;
        //break;
    }
    storePoints();
    storeParams();
}