void randomseed();
void fixedseed();
void noparallel();
void initialize(RngStream* myRng);
double runif(RngStream* myRng);
double runifab(RngStream* myRng, double a, double b);
double rnorm(RngStream* myRng, double mu, double sigma);
double rexp(RngStream* myRng, double theta);
double rgamma(RngStream* myRng, double alpha, double beta);
double rchisq(RngStream* myRng, int df);
double rbeta(RngStream* myRng, double alpha, double beta);
double gamma(double x);

