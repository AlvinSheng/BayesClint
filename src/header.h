#include <RcppArmadillo.h>
#include <RcppGSL.h>
// [[Rcpp::depends(RcppArmadillo, RcppGSL)]]

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

using namespace arma;
using namespace Rcpp;

#include <setjmp.h>
#define TRY do{ jmp_buf ex_buf__; if( !setjmp(ex_buf__) ){
#define CATCH } else {
#define ETRY } }while(0)
#define THROW longjmp(ex_buf__, 1)

void SampleUU(gsl_rng * rr,int r, int n,int Np,int * P,double *** A, double ** U,double *** X,double** s2);
double mean(int n,double * x);
void SampleIntercept(gsl_rng * rr,int n, int r, double * intercept, double* sigma2, double sigma20,double ** U, double ** A, double **y);
bool compar(int n,bool *u,bool *v);
void logPostGam(double *logpo,int * IndVar,int Np,int r, int n,int* P,double *** Tau, double ** U,double *** X,double **s2,bool ** rho,bool *** Gam,double** qv,double* q);
bool **  UniqueModel(int nbrsample, int p, bool ** rhosample,int * modelidx,int *countd1);
void SamplerhoGamma(gsl_rng * rr,int r, int n,int IndVar,int p, bool * rho,double ** Tau, double ** U,double ** X, double* q1,double q2,double* s2,double* quadForm,bool** Gam,double *loggauss);
void save1d(char *filename,int n,double *data);
void sampleGam(gsl_rng * rr,int r, int n,int p, bool * rho,double ** Tau, double ** U,double ** X,double *s2,double* quadForm,bool ** Gam,double * loggauss,double q);
void logGausQuadForm(int j,int r, int n,int p, double ** Tau, double ** U,double ** X,double s2,double* quadForm,bool * Gam,double * loggauss);
void save2d(char *filename,int n,int p,double ** data);
void SampleU(gsl_rng * rr,int r, int n,int p0,int p1,int p2,double *** A, double ** U,double *** X,double** s2);
void EffectZero(int l,gsl_rng * rr,int K,int p,bool rho, bool * R,double* B,double *B0,bool ** Path,double alphab0, double betab0,double alpha, double * lambda2,double *AcceptB0,bool ** Gam);
void TauLambda(int l,gsl_rng * rr,int K,int p,double *A,double* B,double B0,double * Tau,bool ** Path,double alpha, double * lambda2,double *s2,bool ** Gam);
void GroupEffect(int l,gsl_rng * rr,bool rho,int K,int p, bool * R,double *A,double* B,double B0,double * Tau,bool ** Path,double alphab, double betab,double w,double alpha, double * lambda2,double *s2,double * AcceptR,bool **Gam);
void sigma2(gsl_rng * rr, int p,int n,double a0, double b0, double* quadForm,double * s2);
void LoadAOther(gsl_rng * rr,int r, int n,int p, bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,bool ** Gam);
void EstimateLoad(gsl_rng * rr,int r, int n,int p, bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,bool ** Gam);
int multivariate_gaussian (const gsl_rng * r,
                               const gsl_vector * mu,
                               const gsl_matrix * L,
                               gsl_vector * result);

// void mcmcfn(int *Method1,int * n1,int *P,int *r1,int *Np1,double *datasets,int * IndVar,int *K,int * Paths,int* maxmodel1, int *nbrsample1,int *burninsample1,double *CompoSelMean,double *VarSelMean,double * VarSelMeanGlobal,double * GrpSelMean,double *GrpEffectMean,double * IntGrpMean,double * EstU, double * EstSig2,double *InterceptMean,double *EstLoadMod,double *EstLoad,int *nbrmodel1,double *postgam,double* priorcompsel,double * priorcompselo,double* priorb0,double* priorb,double *priorgrpsel,double *probvarsel);
void SampleGammaRho(gsl_rng * rr,int r, int n,int IndVar,int p, bool * rho,double ** Tau, double ** U,double ** X, double* q1,double q2,double* s2,double* quadForm,bool** Gam,double *loggauss);



int logmultigaussian(const gsl_vector * x,
                                       const gsl_vector * mu,
                                       const gsl_matrix * L,
                                       double * result,double *quadF,
                                       gsl_vector * work);

int logmultigaussianT(const gsl_vector * x, const gsl_vector * y,
                                       const gsl_matrix * L,
                                       double * result,double *quadF,
                                       gsl_vector * work);
double  logpost(int r, int n,int p, bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,bool ** Gam,double * loggauss);
void rho1(gsl_rng * rr,int r, int n,int p, bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,bool** Gam,double *loggauss);
void proposal(int n,bool *R,bool *R1,float phi, gsl_rng * r);
void findc(int n,bool *R,int a,int * IDX, int *nx);
double inverseGaussian(gsl_rng * r, double mu, double lambda);
void readBoolVector(char *filename, int nRows,  bool * data );
void sort(int n,double *x,int *idx);
// double auc(int n, double * esti,bool class[n]);
void readDoubleArray(char *filename, int nRows, int nCols, double ** data );
void readBoolArray(char *filename, int nRows, int nCols, bool ** data );
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
// void nrerror(char error_text[]);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
bool **bmatrix(int nrl, int nrh, int ncl, int nch);
void free_bmatrix(bool **m, int nrl, int nrh, int ncl, int nch);

// new functions within utils.h
arma::mat convertToArmaMat(double** A, int rows, int cols);
double** convertToPtrArray(const arma::mat& A);
double **rbind_array(double ** array, int Np, arma::vec n, int nbrft);
arma::mat concatenateRowWise(const arma::field<arma::mat>& fieldMatrices);
arma::field<arma::mat> unbind_array(const arma::mat& arma_catmat, const arma::vec& n);
