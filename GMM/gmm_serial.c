#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define N   100
#define K   10
#define D   3

void generateRandomData(void);
void initParameters(void);
void serialEM_GMM(void);
float G(int n, int k);
void diagMatrixInverse(int k, float matrixInversed[][D]);
float diagMatrixBy2VectorsMultiply(float matrix[][D], float * vector1, float * vector2);
float det(int k);
void getComponent(int k);

float x[N][D], w[K], mu[K][D], theta[K][D][D];

int main(){
    int k;
    for (k=0;k<K;k++){
        w[k] = 0;
    }
    initParameters();
    generateRandomData();
    serialEM_GMM();
    /*
    for (k=0;k<K;k++){
        getComponent(k);
        printf("\n");
    }
    */
    return 0;
}

void generateRandomData(){
    for (int n=0;n<N;n++){
        for (int d=0;d<D;d++){
            x[n][d] = ((float)rand()/(float)(RAND_MAX))*20.0;
        }
    }
}

void initParameters(){
    int k,d;
    // init w
    float acc_w_sum = 0;
    for (k=0;k<K;k++){
        w[k] = ((float)rand()/(float)(RAND_MAX)) * (1-acc_w_sum);
        acc_w_sum += w[k];
    }
    // init mu, theta
    for (k=0;k<K;k++){
        for (d=0;d<D;d++){
            mu[k][d] = ((float)rand()/(float)(RAND_MAX))*20.0;
            theta[k][d][d] = ((float)rand()/(float)(RAND_MAX))*20.0;
        }
    }
}

void serialEM_GMM(){
    int n,k,d;
    float g[K], r[N][K], Nk[K] = {0}, sigma2[K][D],s;
    bool stopCondition = false;
    int counter = 0;
    while (!stopCondition){
        // Expectation step (E-step)
        for (n=0;n<N;n++){
            s = 0;
            for (k=0;k<K;k++){
                g[k] = w[k]*G(n,k); // n : vector index ; k : guassian component index
                s = s + g[k];
            }
            for (k=0;k<K;k++){
                r[n][k] = g[k]/s;
            }
        }
        // Maximization step (M-step)
        for (k=0;k<K;k++){
            for (d=0;d<D;d++){
                sigma2[k][d] = 0.0;
                mu[k][d] = 0.0;
            }
            w[k] = 0.0;
            Nk[k] = 0.0;
        }
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
        counter ++;
        getComponent(0);
        if (counter == 1000)    
            break;

    }
}

float G(int n, int k){
    float x_mu[D];
    float thetaInverse[D][D] = {0};
    float matrixMult = 0;
    float expEval;
    for (int d=0;d<D;d++){
        x_mu[d] = x[n][d] - mu[k][d];
    }
    diagMatrixInverse(k, thetaInverse);
    matrixMult = diagMatrixBy2VectorsMultiply(thetaInverse, x_mu, x_mu);
    expEval = (float)exp(-0.5*matrixMult);
    return expEval/sqrt(pow(2*M_PI, D)*det(k));
}

void diagMatrixInverse(int k, float matrixInversed[][D]){
    for (int d=0;d<D;d++){
        matrixInversed[d][d] = 1/theta[k][d][d];
    }
}

float diagMatrixBy2VectorsMultiply(float matrix[][D], float * vector1, float * vector2){
    float result = 0;
    for (int d=0;d<D;d++){
            result += matrix[d][d]*vector1[d]*vector2[d];
    }
    return result;
}

float det(int k){
    float result = 1;
    for (int d=0;d<D;d++){
        result *= theta[k][d][d];
    }
    return result;
}
void getComponent(int k){
    int d;
    printf("w[%d] = %f\n",k, w[k]);
    
    printf("mu[%d] = [", k);
    for (d=0;d<D;d++){
            printf("%f, ", mu[k][d]);
    }
    printf("]\n");
    printf("theta[%d] = [", k);
    for (d=0;d<D;d++){
            printf("%f, ", theta[k][d][d]);
    }
    printf("]\n");
    
}