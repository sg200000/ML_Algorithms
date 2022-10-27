#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define N   100
#define K   5
#define D   3

void generateRandomData(void);
void serialEM_GMM(void);
float G(float * x, float * mu, float theta[][D]);
void diagMatrixInverse(float matrix[][D], float matrixInversed[][D]);
float diagMatrixBy2VectorsMultiply(float matrix[][D], float * vector1, float * vector2);
float det(float matrix[][D]);
void getComponent(int k);

float x[N][D], w[K], mu[K][D] = {0}, theta[K][D][D] = {0};

int main(){
    for (int k=0;k<K;k++){
        w[k] = 0;
    }
    printf("%f\n", w[0]);
    generateRandomData();
    serialEM_GMM();
    getComponent(0);
    return 0;
}

void generateRandomData(){
    for (int n=0;n<N;n++){
        for (int d=0;d<D;d++){
            x[n][d] = rand()%30;
        }
    }
}

void serialEM_GMM(){
    int s,n,k,d;
    float g[K], r[N][K], Nk[K] = {0}, sigma2[K][D] = {0};
    bool stopCondition = false;
    while (!stopCondition){
        // Expectation step (E-step)
        for (n=0;n<N;n++){
            s = 0;
            for (k=0;k<K;k++){
                theta[0][0][0] = 1;
                printf("%f\n", theta[0][0][0]);
                g[k] = w[k]*G(x[n],mu[k], (theta[k]));
                s = s + g[k];
            }
            for (k=0;k<K;k++){
                r[n][k] = g[k]/s;
            }
        }
        // Maximization step (M-step)
        for (n=0;n<N;n++){
            for (k=0;k<K;k++)
                Nk[k] = Nk[k] + r[n][k];
        }
        for (k=0;k<K;k++){
            w[k] = Nk[k]/N;
        }
        for (n=0;n<N;n++){
            for (k=0;k<K;k++){
                for (d=0;d<D;d++)
                    mu[k][d] = mu[k][d] + (r[n][k]*x[n][d])/Nk[k];
            }
        }
        for (n=0;n<N;n++){
            for (k=0;k<K;k++){
                for (d=0;d<D;d++)
                    sigma2[k][d] = sigma2[k][d] + r[n][k]*(x[n][d]-mu[k][d])*(x[n][d]-mu[k][d])/Nk[k];
            }
        }
        for (k=0;k<K;k++){
            for (d=0;d<D;d++)
                theta[k][d][d] = sigma2[k][d];
        }
        break;

    }
}

float G(float * x, float * mu, float theta[][D]){
    printf("%f\n", theta[0][0]);
    float x_mu[D];
    float thetaInverse[D][D] = {0};
    float matrixMult = 0;
    float expEval;
    for (int d=0;d<D;d++){
        x_mu[d] = x[d] - mu[d];
    }
    diagMatrixInverse(theta, thetaInverse);
    matrixMult = diagMatrixBy2VectorsMultiply(thetaInverse, x_mu, x_mu);
    expEval = (float)exp(-0.5*matrixMult);
    return expEval/sqrt(pow(2*M_PI, D)*det(theta));
}

void diagMatrixInverse(float matrix[][D], float matrixInversed[][D]){
    for (int d=0;d<D;d++){
        matrixInversed[d][d] = 1/matrix[d][d];
    }
}

float diagMatrixBy2VectorsMultiply(float matrix[][D], float * vector1, float * vector2){
    float result = 0;
    for (int d=0;d<D;d++){
            result = matrix[d][d]*vector1[d]*vector2[d];
    }
    return result;
}

float det(float matrix[][D]){
    float result = 1;
    for (int d=0;d<D;d++){
        result *= matrix[d][d];
    }
    return result;
}
void getComponent(int k){
    printf("w[%d] = %f\n",k, w[k]);
}