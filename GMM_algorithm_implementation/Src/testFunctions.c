void testG(double mu, double theta){
    mu[0][0] = mu;
    theta[0][0][0] = theta;
    x[0][0] = 0.0;
    double g;
    for (int n=1;n<N;n++){
        x[n][0] = x[n-1][0]+1.0;
        printf("%lf\n", x[n][0]);
    }
    FILE * gnuplot = fopen("gnuplot.txt", "w");
    for (int n = 0; n < N; n++){
        g = G(n, 0);
        fprintf(gnuplot, "%g %lf\n", x[n][0], g);
        printf("%lf", g);
    }
    fflush(gnuplot);
}