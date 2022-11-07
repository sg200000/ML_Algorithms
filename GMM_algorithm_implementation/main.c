#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>


#define K   6
#define D   1
#define N   50

//global vectors
double w[K];
double mu[K][D];
double theta[K][D][D];
double x[N][D];

//intermediate vectors
double R[N][K];

void initParams(){
    mu[0][0] = 15.0;
    theta[0][0][0] = 5.0;
    w[0] = 1.0/K;
    for (int k=1;k<K;k++){
        mu[k-1][0] = mu[k][0]+3.0;
        theta[k][0][0] = theta[k-1][0][0]+2.0;
        w[k] = 1.0/K;
    }
    x[0][0] = 0.0;
    for (int n=1;n<N;n++){
        x[n][0] = x[n-1][0]+1.0;
        //printf("%lf\n", x[n][0]);
    }
}
double G(int n, int k){
    double result;
    double thetaInv[D];
    double xn_k[D];
    double exponent = 0;
    double detTheta = 1;
    for (int d=0;d<D;d++){
        thetaInv[d] = 1/theta[k][d][d];
        xn_k[d] = x[n][d] - mu[k][d];
        exponent += xn_k[d]*xn_k[d]*thetaInv[d];
        detTheta *= theta[k][d][d];
    }
    result = exp(-0.5*exponent)/(sqrt(2*M_PI*detTheta));
    return result;
}

double r(int n, int k){
    double result = 0.0;
    for (int kt=0;kt<K;kt++)
        result += w[kt]*G(n, kt);
    //printf("%lf\n", result);
    result = w[k]*G(n,k)/result;
    return result;
}

void calcR(){
    for (int n=0;n<N;n++){
        for (int k=0;k<K;k++){
            R[n][k] = r(n,k);
            //printf("%lf ", R[n][k]);
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
        //    printf("%lf ", mu[k][d]);
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
void main(){
    initParams();
    double g;
    FILE * gnuplot = fopen("gnuplot.txt", "w");
    for (int n = 0; n < N; n++){
       // g = G(n, 0);
        fprintf(gnuplot, "%g %lf\n", x[n][0], g);
       // printf("%lf", g);
    }
    for (int i=0;i<100;i++){
        calcR();
        calcParams();
        printParams();
    }
    
}