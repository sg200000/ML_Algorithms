#include <stdio.h>
#include <stdlib.h>
#include <math.h>

# define M_PI		3.14159265358979323846	/* pi */

double * w;
double ** mu;
double ** theta;

double ** x;

int dim,  numComps;

int main(){
    /* begin parsing model*/
    FILE * model = fopen("model.csv", "r");

    /* Parse number of components and he dimension*/
    fscanf(model, "%d", &numComps);
    fscanf(model, "%d", &dim);
    printf("%d %d\n", dim, numComps);

    /* Allocate memory for data */
    w = (double *) malloc(numComps*sizeof(double));
    mu = (double **) malloc(numComps*sizeof(double *));
    for (int i=0;i<numComps;i++)
        mu[i] = (double *)malloc(dim*sizeof(double));
    theta = (double **) malloc(numComps*sizeof(double *));
    for (int i=0;i<numComps;i++)
        theta[i] = (double *)malloc(dim*sizeof(double));

    double tempBuffer;
    for (int k=0;k<numComps;k++){
        fscanf(model,"%lf", &tempBuffer);
        w[k] = tempBuffer;
        fgetc(model);
        printf("w[%d] = %lf\n", k, w[k]);
        for (int d=0;d<dim;d++){
            fscanf(model,"%lf", &tempBuffer);
            mu[k][d] = tempBuffer;
            printf("mu[%d][%d] = %lf\n",k,d, mu[k][d]);
        }
        for (int d=0;d<dim;d++){
            fscanf(model,"%lf", &tempBuffer);
            theta[k][d] = tempBuffer;
            printf("theta[%d][%d] = %lf\n",k, d, theta[k][d]);
        }
    }
    /* end parsing model */
    fclose(model);
    for (int i=0;i<numComps;i++){
        free(mu[i]);
        free(theta[i]);
    }
    free(w);
    free(mu);
    free(theta);
    return 0;
}