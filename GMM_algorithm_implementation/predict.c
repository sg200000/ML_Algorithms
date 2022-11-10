#include <stdio.h>
#include <math.h>

int main(){
    FILE * model = fopen("model.csv", "r");
    int dim,  numComps;
    fscanf(model, "%d", &numComps);
    fscanf(model, "%d", &dim);
    printf("%d %d\n", dim, numComps);
    double w[numComps];
    double mu[numComps][dim];
    double theta[numComps][dim];
    double tempBuffer;
    //fgetc(model);
    for (int k=0;k<numComps;k++){
        fscanf(model,"%lf", &tempBuffer);
        w[k] = tempBuffer;
        fgetc(model);
        printf("w[%d] = %lf\n", k, w[k]);
        for (int d=0;d<dim;d++){
            fscanf(model,"%lf", &tempBuffer);
            mu[numComps][dim] = tempBuffer;
            printf("mu[%d][%d] = %lf\n",k,d, mu[numComps][dim]);
        }
        for (int d=0;d<dim;d++){
            fscanf(model,"%lf", &tempBuffer);
            theta[numComps][dim] = tempBuffer;
            printf("theta[%d][%d] = %lf\n",k, d, theta[numComps][dim]);
        }
    }
    /*while(fgets(chunk, sizeof(chunk), fp) != NULL) {
        fputs(chunk, stdout);
        fputs("|*\n", stdout); 
16  }*/
    return 0;
}