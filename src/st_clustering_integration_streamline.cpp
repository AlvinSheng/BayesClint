#include <RcppArmadillo.h>
#include <RcppGSL.h>
// [[Rcpp::depends(RcppArmadillo, RcppGSL)]]

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <time.h>
#include <random>

#include "header.h"
#include "utils.h"

using namespace arma;
using namespace Rcpp;



void my_set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}



// Function to sample from a normal distribution
vec rnorm_cpp(int n, double mean, double sd) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(mean, sd);

  vec samples(n);
  for (int i = 0; i < n; ++i) {
    samples[i] = d(gen);
  }
  return samples;
}



// Function to sample from an inverse gamma distribution
double rinvgamma_cpp(double shape, double scale) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::gamma_distribution<> d(shape, 1.0 / scale);
  return 1.0 / d(gen);
}



///// myfunction.c subroutines

void logPostGam(double *logpo, arma::vec IndVar, int Np, int r, int n, arma::vec P,
                double *** Tau, double ** U,double *** X, double **s2, bool ** rho,
                bool *** Gam, double** qv, double* q) {
  double logpost=0;
  int m,j,l;
  for (m=0;m<Np;m++){
    for (j=0;j<P[m];j++){
      double logp;double quad;
      logGausQuadForm(j,r, n,P[m], Tau[m], U,X[m], s2[m][j],&quad,Gam[m][j],&logp);
      logpost+=logp;
      for (l=0;l<r;l++){
        if (IndVar[m]!=2)
          logpost+=Gam[m][j][l]*log(qv[m][l])+(1-Gam[m][j][l])*log(1-qv[m][l]);
      }
    }

    if (IndVar[m]==0){
      for (l=0;l<r;l++){
        logpost+=rho[m][l]*log(q[m])+(1-rho[m][l])*log(1-q[m]);
      }
    }
  }
  *logpo=logpost;
}

void logGausQuadForm(int j,int r, int n,int p,double ** Tau, double ** U,double ** X,double s2,double* quadForm,bool * Gam,double * loggauss){
  int i,s,s1;
  int NZ1[r];
  int nx1=0;
  findc(r,Gam,0,NZ1, &nx1);

  double result=0;double quadF=0;
  if (nx1>0){
    gsl_vector *work1 =gsl_vector_alloc (nx1);
    double * Sigma1= static_cast<double*>(malloc(nx1*nx1*sizeof(double)));
    for (s=0;s<nx1;s++){
      for (s1=0;s1<=s;s1++){
        double a=0;
        for (i=0;i<n;i++){
          //if (Gam[NZ1[s]]*Gam[NZ1[s1]]==1)
          a+=U[i][NZ1[s]]*U[i][NZ1[s1]];
        }
        Sigma1[s*nx1+s1]=Sigma1[s1*nx1+s]=a;

      }
      Sigma1[s*nx1+s]+=(1.0/Tau[NZ1[s]][j]);
      //printf("q1=%lf\n",Tau[NZ1[s]][j]);
    }
    gsl_matrix_view m1  = gsl_matrix_view_array (Sigma1, nx1,nx1);
    //printf("Error on Gam\n");

    gsl_linalg_cholesky_decomp (&m1.matrix);
    gsl_vector *xi =gsl_vector_alloc (n);
    for (i = 0; i < n; i++) {
      gsl_vector_set (xi, i, X[i][j]/sqrt(s2));
    }

    gsl_vector *y =gsl_vector_alloc (nx1);
    double sumT=0;
    for (s=0;s<nx1;s++){
      double a=0;
      //if (Gam[NZ1[s]]==1){
      sumT+=log(Tau[NZ1[s]][j]);
      for (i=0;i<n;i++){
        a+=U[i][NZ1[s]]*X[i][j]/sqrt(s2);
      }
      //}
      gsl_vector_set (y, s, a);
    }
    logmultigaussianT(xi, y,  &m1.matrix,&result,&quadF, work1);
    result=result-0.5*sumT-(n/2.0)*log(s2);
    gsl_vector_free (y);gsl_vector_free (work1);gsl_vector_free (xi);
    free(Sigma1);
  } else {
    for (i = 0; i < n; i++){
      quadF+=pow(X[i][j],2);
    }
    quadF=quadF/s2;
    result=-(n/2.0)*log(s2)- 0.5*n*log(2.0*M_PI)-0.5*quadF;
  }

  *quadForm=quadF;
  *loggauss=result;

}



void SampleGammaRho(gsl_rng * rr,int r, int n,int IndVar,int p, bool * rho,double ** Tau, double ** U,double ** X, double* q1,double q2,double* s2,double* quadForm,bool** Gam,double *loggauss, arma::mat ups2_scaling){
  int l,j;
  int i,ll;
  double ** Ustar = dmatrix(0, n-1, 0, r-1);
  bool *rhonew= static_cast<bool*>(malloc(r*sizeof(bool)));
  bool Gamnew[p][r];
  for (l=0;l<r;l++){
    for (j=0;j<p;j++){
      Gamnew[j][l]=Gam[j][l];
    }
    rhonew[l]=rho[l];
  }

  if (IndVar==1){ // you can re-purpose the IndVar variable to indicate
    // whether you want to do simple column-wise selection or nested selection.
    // But it's not an exact corollary, since I'm still using a two-dimensional matrix A
    for (l=0;l<r;l++){
      double quad1=0;double quad2=0;
      double loggaussold=0; double loggaussnew=0;
      Gam[0][l]=0;
      logGausQuadForm(0,r, n,p, Tau,  U,X,s2[0],&quad1,Gam[0],&loggaussold);
      Gamnew[0][l]=1;
      logGausQuadForm(0,r, n,p,Tau,  U,X,s2[0],&quad2,Gamnew[0],&loggaussnew);
      double rat=loggaussnew+log(q2)-loggaussold-log(1-q2);
      double uni=gsl_ran_flat (rr, 0, 1);
      if (log(uni/(1-uni))<rat){
        Gamnew[0][l]=Gam[0][l]=rho[l]=1;
        loggauss[0]=loggaussnew;
        quadForm[0]=quad2;
      } else {
        Gamnew[0][l]=Gam[0][l]=rho[l]=0;
        loggauss[0]=loggaussold;
        quadForm[0]=quad1;
      }
    }

  } else {
    for (l=0;l<r;l++){
      double logq=0;
      if (rho[l]==0){
        // Propose to turn the component on
        rhonew[l]=1;
        double logpostnew=0;
        double logpostold=0;
        double * quadForm1= static_cast<double*>(malloc(p*sizeof(double)));
        double * loggausnew1= static_cast<double*>(malloc(p*sizeof(double)));
        double quad1=0; double quad2=0;
        double loggaussold=0; double loggaussnew=0;
        // Propose to turn on/ off features
        for (j=0;j<p;j++){

          // Transform U
          for (i=0;i<n;i++) {
            for (ll=0;ll<r;ll++) {
              Ustar[i][ll] = U[i][ll] * ups2_scaling(i, j);
            }
          }

          double logqj=0;double logmqj=0;
          logpostold+=loggauss[j];
          Gam[j][l]=0;
          logGausQuadForm(j,r, n,p, Tau,  Ustar,X,s2[j],&quad1,Gam[j],&loggaussold);
          Gamnew[j][l]=1;
          logGausQuadForm(j,r, n,p,Tau,  Ustar,X,s2[j],&quad2,Gamnew[j],&loggaussnew);
          double rat=-loggaussold+loggaussnew+log(q1[l])-log(1-q1[l]); // log(P_lj / Q_lj)

          double x1=loggaussnew+log(q1[l]); // log G_1
          double x2=loggaussold+log(1-q1[l]); // log G_0
          double maxq=std::max(x1,x2);
          double den=maxq+log(exp(x1-maxq)+exp(x2-maxq)); // log(G_0 + G_1)
          logqj=x1-den; // log P_lj
          logmqj=x2-den; // log Q_lj = 1 - log P_lj

          double uni=gsl_ran_flat (rr, 0, 1);
          if ((log(uni/(1-uni))<rat)){
            // Turn the feature on
            Gam[j][l]=1;Gamnew[j][l]=1;
            logpostnew+=loggaussnew+log(q1[l]);
            loggausnew1[j]=loggaussnew;
            // log proposal difference- add the probability that a feature is on given a component is on
            logq+=logqj;
            quadForm1[j]=quad2;
          } else {
            Gam[j][l]=0;Gamnew[j][l]=0;
            logq+=logmqj;
            // log proposal difference- add the probability that a feature is off given a component is on
            logpostnew+=loggaussold+log(1-q1[l]);
            quadForm1[j]=quad1;
            loggausnew1[j]=loggaussold;
          }
        }
        // Add the log probability that the component is on
        logpostnew+=log(q2);
        logpostold+=log(1-q2);
        double un=gsl_ran_flat (rr, 0, 1);
        double rat1=logpostnew-logpostold-logq; // log acceptance ratio
        if (log(un)<rat1){
          // accept having the component on
          rho[l]=1;
          // store the new log gauss and quadForms
          for (j=0;j<p;j++){
            quadForm[j]=quadForm1[j];
            loggauss[j]=loggausnew1[j];
          }
        } else {
          // stay off
          rho[l]=0;
          for (j=0;j<p;j++) Gam[j][l]=Gamnew[j][l]=0;
        }
        // initialize gamma new for the next iteration & remove quadform and loggaussnew since unneccessary for the next step?
        rhonew[l]=rho[l];
        free(quadForm1);free(loggausnew1);
      } else {
        // gamma is on and we are proposing to turn it off
        rhonew[l]=0;
        // initialize proposal data structures
        double logpostnew=0;
        double logpostold=0;
        double * quadForm1= static_cast<double*>(malloc(p*sizeof(double)));
        double * loggausnew1= static_cast<double*>(malloc(p*sizeof(double)));
        double quad1;
        double loggaussnew=0; double loggaussn=0;
        double logpq=0;
        for (j=0;j<p;j++){

          // Transform U
          for (i=0;i<n;i++) {
            for (ll=0;ll<r;ll++) {
              Ustar[i][ll] = U[i][ll] * ups2_scaling(i, j);
            }
          }

          logpostold+=loggauss[j];
          Gamnew[j][l]=1;
          logGausQuadForm(j,r, n,p,Tau,  Ustar,X,s2[j],&quad1,Gamnew[j],&loggaussn);
          Gamnew[j][l]=0;
          logGausQuadForm(j,r, n,p,Tau,  Ustar,X,s2[j],&quad1,Gamnew[j],&loggaussnew);
          if (p!=1){
            double x1=loggaussn+log(q1[l]);
            double x2=loggaussnew+log(1-q1[l]);
            double maxq=std::max(x1,x2);
            double den=maxq+log(exp(x1-maxq)+exp(x2-maxq));
            double logqj=x1-den;
            double logmqj=x2-den;
            logq+=Gam[j][l]*log(q1[l])+(1-Gam[j][l])*log(1-q1[l]);
            logpq+=Gam[j][l]*logqj+(1-Gam[j][l])*logmqj;
          }
          logpostnew+=loggaussnew;
          loggausnew1[j]=loggaussnew;
          quadForm1[j]=quad1;
          Gamnew[j][l]=Gam[j][l];
        }
        logpostnew+=log(1-q2)+logpq;
        logpostold+=log(q2)+logq;
        double un=gsl_ran_flat (rr, 0, 1);
        double rat1=logpostnew-logpostold;
        if (log(un)<rat1){
          rho[l]=0;
          for (j=0;j<p;j++){
            Gamnew[j][l]=Gam[j][l]=0;
            quadForm[j]=quadForm1[j];
            loggauss[j]=loggausnew1[j];
          }
        } else {
          for (j=0;j<p;j++){
            Gamnew[j][l]=Gam[j][l];
          }
        }
        free(quadForm1);free(loggausnew1);
        rhonew[l]=rho[l];
      }
      /* Within move*/
      if (rho[l]==1){
        double  logpostold=0;
        double * quadForm1= static_cast<double*>(malloc(p*sizeof(double)));
        double * loggausnew1= static_cast<double*>(malloc(p*sizeof(double)));
        double quad1,quad2;
        double loggaussold,loggaussnew;
        for (j=0;j<p;j++){

          // Transform U
          for (i=0;i<n;i++) {
            for (ll=0;ll<r;ll++) {
              Ustar[i][ll] = U[i][ll] * ups2_scaling(i, j);
            }
          }

          logpostold+=loggauss[j];
          logGausQuadForm(j,r, n,p,Tau,  Ustar,X,s2[j],&quad1,Gam[j],&loggaussold);
          Gamnew[j][l]=1-Gam[j][l];
          logGausQuadForm(j,r, n,p,Tau,  Ustar,X,s2[j],&quad2,Gamnew[j],&loggaussnew);
          double rat=-loggaussold+loggaussnew;
          double uni=gsl_ran_flat (rr, 0, 1);
          if (Gam[j][l]==0){
            quad1=quad2;
            loggaussold=loggaussnew;
            rat+=log(q1[l])-log(1-q1[l]);
          } else {
            quad2=quad1;
            loggaussnew=loggaussold;
            rat=-rat+log(q1[l])-log(1-q1[l]);
          }
          if (log(uni/(1-uni))<rat){
            Gam[j][l]=1;
            quadForm[j]=quad1;
            loggauss[j]=loggaussold;
          } else {
            Gam[j][l]=0;
            quadForm[j]=quad2;
            loggauss[j]=loggaussnew;
          }
          Gamnew[j][l]=Gam[j][l];
        }
        free(quadForm1);free(loggausnew1);
      }
    }
  }
  free(rhonew);
  free_dmatrix(Ustar,0, n-1, 0, r-1);
}



arma::ivec zetalabelupdate(const int nsites, arma::mat u_mat, arma::mat mu, arma::mat inv_Sigma,
                           int M, arma::ivec zeta_m, arma::ivec kappa_m, arma::mat theta)
{

  // Update the zeta labels
  // Create new objects
  int rowstart=0, rowend=0;
  double numnbr;
  double log_prioru_matkernel;
  IntegerVector one2M;
  NumericVector log_probmasses(M), log_probmasses_diff(M), probmasses_ratio(M);
  double max_log_pm;
  arma::ivec clusternew(nsites);

  one2M = seq_len(M);



  //  Update each zeta label in turn
  clusternew = zeta_m;
  for(int j = 0; j < nsites; j++)
  {

    // Find the probability masses for each possible zeta label m = 1, ..., M
    for(int m = 1; m <= M; m++) {

      // Evaluating kernel of u_mat likelihood
      arma::colvec diff = u_mat.row(j).t() - mu.col(m - 1);
      log_prioru_matkernel =
        -0.5 * as_scalar( diff.t() * inv_Sigma * diff );

        // Multiplying the kernels together
        log_probmasses[m - 1] = log_prioru_matkernel + log(theta(m - 1, kappa_m(j) - 1));

    }

    // Use the log_probmasses to generate new zeta label

    max_log_pm = *std::max_element(log_probmasses.begin(), log_probmasses.end());

    log_probmasses_diff = log_probmasses;
    for(auto& element : log_probmasses_diff)
      element -= max_log_pm;

    probmasses_ratio = log_probmasses_diff;
    for(auto& element : probmasses_ratio)
      element = exp(element);

    // propose a value
    clusternew[j] = sample(one2M, 1, false, probmasses_ratio)[0];

  }

  return clusternew;
}



arma::ivec kappalabelupdate(arma::imat Wtriplet, arma::imat Wbegfin, const int nsites,
                            int M, double beta_m, arma::ivec zeta_m, arma::ivec kappa_m, arma::mat theta)
{

  // Update the kappa labels
  // Create new objects
  int rowstart=0, rowend=0;
  double numnbr;
  double log_priorpottskernel;
  IntegerVector one2M;
  NumericVector log_probmasses(M), log_probmasses_diff(M), probmasses_ratio(M);
  double max_log_pm;
  arma::ivec clusternew(nsites);

  one2M = seq_len(M);



  //  Update each kappa label in turn
  clusternew = kappa_m;
  for(int j = 0; j < nsites; j++)
  {

    // Find the probability masses for each possible kappa label m = 1, ..., M
    for(int m = 1; m <= M; m++) {

      // Evaluating kernel of Spatial Potts Model prior
      // needs count of neighbors with kappa label m
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      numnbr = 0;
      for(int l = rowstart; l < rowend; l++) numnbr += (clusternew[(Wtriplet(l,1) - 1)] == m);
      log_priorpottskernel = beta_m * numnbr;

      // Multiplying the kernels together
      log_probmasses[m - 1] = log_priorpottskernel + log(theta(zeta_m(j) - 1, m - 1));

    }

    // Use the log_probmasses to generate new kappa label

    max_log_pm = *std::max_element(log_probmasses.begin(), log_probmasses.end());

    log_probmasses_diff = log_probmasses;
    for(auto& element : log_probmasses_diff)
      element -= max_log_pm;

    probmasses_ratio = log_probmasses_diff;
    for(auto& element : probmasses_ratio)
      element = exp(element);

    // propose a value
    clusternew[j] = sample(one2M, 1, false, probmasses_ratio)[0];

  }

  return clusternew;
}



// [[Rcpp::export]]
Rcpp::List mcmcfn(arma::vec n, int P, int r, int Np, arma::mat datasets, int IndVar, int nbrsample, // nbrsample=n_iter-n_burnin
                  int burninsample, arma::vec CompoSelMean, arma::vec VarSelMean, arma::vec VarSelMeanGlobal, arma::vec priorcompsel, double probvarsel,
                  int C, arma::mat EstMat1, arma::mat EstMat2, arma::field<arma::ivec> EstIMat1,
                  int K, arma::field<arma::ivec> EstIMat2, arma::mat EstMat3, arma::field<arma::imat> Wtriplet_list, arma::field<arma::imat> Wbegfin_list, arma::vec beta, arma::mat interp_coef_mat,
                  int chainNbr, bool streamline = true) {

  int i,l,j,c,k;
  int m;
  clock_t t1 = clock();

  printf("Number of MCMC samples after burn-in is %d\n",nbrsample);

  printf("Number of burn-in is %d\n",burninsample);

  for (m=0;m<Np;m++){
    int n_samples = n[m];
    printf("Number of spots in section %d is %d\n", m+1, n_samples);
  }

  printf("Number of components is %d\n", r);

  arma::field<arma::mat> X = unbind_array(datasets, n);



  // set seed for reproducible gsl_rng results
  long seed = chainNbr;
  gsl_rng * rr = gsl_rng_alloc (gsl_rng_rand48);
  gsl_rng_set (rr, seed);
  // set seed for armadillo random generation too, according to https://blog.thecoatlessprofessor.com/programming/cpp/set-r-s-seed-in-rcpp-sequential-case/
  my_set_seed(chainNbr);

  arma::mat A = arma::mat(r, P);

  double ***U= static_cast<double***>(malloc(Np*sizeof(double **)));
  double ***meanU= static_cast<double***>(malloc(Np*sizeof(double **)));
  for (m=0;m<Np;m++){
    U[m]=dmatrix(0,n[m]-1,0,r-1);
    meanU[m]=dmatrix(0,n[m]-1,0,r-1);
    for (i=0;i<n[m];i++){
      for (l=0;l<r;l++){
        U[m][i][l]=gsl_ran_ugaussian(rr);
        meanU[m][i][l]=0;
      }
    }
  }



  /* Hyperparameter*/
  double q=0.5;
  double * qv= static_cast<double*>(malloc(r*sizeof(double)));
  for (l=0;l<r;l++){
    qv[l]= probvarsel;
  }
  // Define hyperparameters
  double al=priorcompsel[0]; double bl=priorcompsel[1]; // Hyper for q. No longer used
  double a0=0.01; double b0=0.01;
  arma::vec dmean(r, fill::zeros);
  arma::mat Sigma0 = arma::mat(r, r, fill::eye);
  double df0 = r + 1;
  arma::mat Omega0(r, r, fill::eye);
  arma::vec c0(C, fill::ones);

  // Initialize parameters to estimate
  bool * rhoest= static_cast<bool*>(malloc(r*sizeof(bool)));
  double * rhomean= static_cast<double*>(malloc(r*sizeof(double)));
  double * Gamvs= static_cast<double*>(malloc(P*sizeof(double)));
  bool ** Gam=bmatrix(0,P-1,0,r-1);
  double ** Gammean=dmatrix(0,P-1,0,r-1);
  double * quadForm= static_cast<double*>(malloc(P*sizeof(double)));
  double * loggauss= static_cast<double*>(malloc(P*sizeof(double)));

  arma::mat ups2(Np, P);
  arma::mat ups2Mean(Np, P);
  arma::mat mu = EstMat1;
  arma::mat Sigma = EstMat2;
  arma::field<arma::ivec> zeta = EstIMat1;
  arma::field<arma::ivec> kappa = EstIMat2;
  arma::mat theta = EstMat3;

  double ** Tau= dmatrix(0,r-1,0,P-1); // bypass parameter
  double * s2= static_cast<double*>(malloc(P*sizeof(double))); // bypass parameter

  for (l=0;l<r;l++){
    rhomean[l]=0;
    for (j=0;j<P;j++) Gammean[j][l]=0;

    double  uni=gsl_ran_flat (rr, 0, 1);

    if (uni<0.5) {
      rhoest[l]=1;
      for (j=0;j<P;j++){
        if (gsl_ran_flat (rr, 0, 1)<qv[l]) {
          Gam[j][l]=1;A(l, j)=0.01;
        } else {
          Gam[j][l]=0;A(l, j)=0;
        }
      }
    } else {
      rhoest[l]=0;
      for (j=0;j<P;j++) {
        Gam[j][l]=0;A(l, j)=0;
      }
    }

    for (j=0;j<P;j++){
      Tau[l][j]=1; // Setting the bypass parameters Tau (i.e., parameters left over from BIP code, which are kept in order to use the SampleGammaRho function)
    }

  }

  for (j=0;j<P;j++) {
    loggauss[j]=0;
    quadForm[j]=Gamvs[j]=0;
    s2[j]=1; // Setting the bypass parameters s2 (i.e., parameters left over from BIP code, which are kept in order to use the SampleGammaRho function)
    for (m=0;m<Np;m++) {
      ups2(m,j)=0.1; ups2Mean(m,j)=0;
    }
  }

  // Dim of a model (i.e., number of gamma and eta elements)
  int dim=0;
  dim=P*r + r;

  int t;
  int n_iter=nbrsample+burninsample;

  bool ** rhomodel= static_cast<bool**>(malloc(nbrsample*sizeof(bool*)));
  for (t=0;t<nbrsample;t++){
    rhomodel[t]= static_cast<bool*>(malloc(dim*sizeof(bool)));
  }



  int nsum = arma::sum(n);
  double ** Xstar=dmatrix(0, nsum-1, 0, P-1);
  arma::mat Xstar_arma; // only used for the loadings matrix A update
  double ** Ucat = rbind_array(U, Np, n, r); // From this point on, U is no longer updated
  arma::mat Ucat_arma = convertToArmaMat(Ucat, nsum, r);
  arma::field<arma::mat> U_arma = unbind_array(Ucat_arma, n);
  double ** Ustar = dmatrix(0, nsum-1, 0, r-1); // only used for the initialized logGausQuadForm
  arma::mat Ustar_subset; // only used for the initialized logGausQuadForm
  arma::mat Ustar_arma; // only used for the loadings matrix A update
  arma::mat ups2_scaling = arma::mat(nsum, P);
  // do pre-multiplication of datasets to get Xstar
  // also initializing ups2_scaling with initial ups2
  k = 0;
  for (m=0;m<Np;m++) {
    for (i=0;i<n[m];i++) {
      for (j=0;j<P;j++) {
        Xstar[i+k][j] = datasets(i+k, j) / pow(ups2(m,j), 0.5);
        ups2_scaling(i+k, j) = 1 / pow(ups2(m,j), 0.5);
      }
    }
    k = k + n[m];
  }



  // Other initializations for MCMC computations
  arma::vec temp;
  // for loadings update step
  int ngamj;
  std::vector<arma::uword> gam_idxj;
  arma::mat fc_var;
  arma::vec fc_offset;
  // for mu
  int tmc;
  arma::vec tumic;
  // for Sigma
  arma::mat Sigma_update_outer_prod(r, r);
  // for theta
  arma::mat tmc_mat;
  arma::vec dirichlet_vec;
  arma::vec dirichlet_output(C);
  // for beta
  double min_beta = 0;
  double max_beta = 4;
  int p_order = 10;
  arma::vec logd_current(Np);
  arma::vec X2(p_order + 1);
  for (m=0;m<Np;m++) { // getting current log normalizing constant log(d(beta))
    for(int k = 1; k <= p_order+1; k++) {
      X2(k-1) = pow(beta[m], k) / double(k);
    }
    logd_current[m] = -arma::dot(X2, interp_coef_mat.col(m)); // -log(d(beta))
  }
  double temperature = 1 / burninsample; // for simulated annealing
  double prob; // for Metropolis sampling with symmetric normal candidate distributions
  double pairwise_count;
  double beta_proposal;
  double logprob_beta_current;
  double logd_proposal;
  double logprob_beta_proposal;

  // for several parameters
  arma::mat inv_Sigma = inv(Sigma);
  arma::mat inv_Sigma0 = inv(Sigma0);



  // Initializing the objects that will contain all of the MCMC draws, to be returned in the output
  arma::field<arma::imat> zeta2mcmc(Np);
  arma::field<arma::imat> kappa2mcmc(Np);
  for (m=0;m<Np;m++) {
    zeta2mcmc(m) = arma::imat(n[m], nbrsample);
    kappa2mcmc(m) = arma::imat(n[m], nbrsample);
  }
  arma::cube mu2mcmc(r, C, nbrsample);
  arma::cube theta2mcmc(C, K, nbrsample);
  arma::mat beta2mcmc(Np, nbrsample);

  // Creating VarSelMeanGlobal_powerset, which contains versions of VarSelMeanGlobal for every combo of components
  arma::vec Pzeros(P, arma::fill::zeros);
  std::vector<arma::Col<int>> rpowerset; // Use std::vector to store the collection of Armadillo subsets
  arma::mat VarSelMeanGlobal_powerset(P, pow(2, r));
  // Creating numdeg_powerset, which will contain the number of selected genes across every combo of components, for every iteration
  arma::vec PzerosN(nbrsample, arma::fill::zeros);
  arma::mat numdeg_powerset(pow(2, r), nbrsample);
  arma::vec PzerosP(pow(2, r), arma::fill::zeros);
  arma::vec numdeg_pow2r(pow(2,r));
  if (!streamline) {

    // For the selection of differentially expressed genes: retrieve a vector of vectors that represents the power set of r elements
    arma::vec one2r = arma::linspace<arma::vec>(1, r, r);
    // Outer Loop: Iterate through all 2^r possible subsets (0 to 2^r - 1)
    for (int i = 0; i < (1 << r); i++) {
      // We'll use a temporary std::vector to build the subset elements first,
      // then convert it to an arma::Col<int> for efficiency and ease of use.
      std::vector<int> temp_subset_elements;
      // Inner Loop: Check each element of the original vector 'one2r'
      for (int j = 0; j < r; j++) {
        // Bitwise Check: If the j-th bit in 'i' is set
        if (i & (1 << j)) {
          // Get the integer value at index j and add it to the temporary list
          temp_subset_elements.push_back(one2r(j)); // Armadillo uses () for element access
        }
      }
      // Convert the temporary std::vector<int> to an Armadillo vector (arma::Col<int>)
      if (temp_subset_elements.empty()) {
        // Handle the empty set (index i=0)
        rpowerset.push_back(arma::Col<int>());
      } else {
        // Create the Armadillo vector from the standard vector data
        arma::Col<int> subset(temp_subset_elements);
        rpowerset.push_back(subset);
      }
    }

    // pre-evaluating the first element of VarSelMeanGlobal_powerset, so you just need to evaluate the slots after the first
    VarSelMeanGlobal_powerset.col(0) = Pzeros;

    // // pre-evaluating the first element of numdeg_powerset, so you just need to evaluate the columns after the first
    // numdeg_powerset.row(0) = PzerosN;
    // // not necessary
  }



  //////// Begin MCMC

  for (t=0;t<n_iter;t++){

    double sumrho=0;

    for (l=0;l<r;l++){
      sumrho+=rhoest[l];
    }

    // Sampling q, the Bernoulli probability for rho
    q=1;//gsl_ran_beta (rr, al+sumrho, bl+r-sumrho);



    // Sampling rho_l and gamma_{jl}, the component and feature selection indicators
    for (j=0;j<P;j++){
      // Transform U according to ups2_scaling
      for (i=0;i<nsum;i++) {
        for (l=0;l<r;l++) {
          Ustar[i][l] = Ucat[i][l] * ups2_scaling(i, j);
        }
      }
      logGausQuadForm(j, r, nsum, P, Tau, Ustar, Xstar, s2[j], &quadForm[j],Gam[j],&loggauss[j]); // updates quadForm and loggauss
      if (t==0){
        loggauss[j]=-DBL_MAX;
      }
    }

    SampleGammaRho(rr, r, nsum, IndVar, P, rhoest,Tau, Ucat, Xstar, qv, q, s2, quadForm, Gam, loggauss, ups2_scaling); // Ucat will be transformed into Ustar for each j inside the SampleGammaRho function



    // Sampling upsilon^2_{mj}, the error variance
    for (m=0;m<Np;m++) {
      for (j=0;j<P;j++) {
        temp = X(m).col(j) - U_arma(m) * A.col(j);
        double inv=1/(b0+0.5*arma::dot(temp, temp));
        ups2(m,j) = 1 / gsl_ran_gamma(rr, a0+n(m)/2.0, inv);
      }
    }



    // Pre-multiplying datasets with Upsilon2^{-1/2} to get Xstar for standardizing away the ups2 error variances
    // also updating ups2_scaling with the updated ups2
    k = 0;
    for (m=0;m<Np;m++) {
      for (i=0;i<n[m];i++) {
        for (j=0;j<P;j++) {
          Xstar[i+k][j] = datasets(i+k, j) / pow(ups2(m,j), 0.5);
          ups2_scaling(i+k, j) = 1 / pow(ups2(m,j), 0.5);
        }
      }
      k = k + n[m];
    }
    Xstar_arma = convertToArmaMat(Xstar, nsum, P);

    // Sampling A, factor loadings matrix
    A = arma::mat(r, P, fill::zeros); // resetting A entries to zero, the default value
    for (j=0;j<P;j++) { // go through each column of A, and update the ones that are nonzero according to gamma_{j.}
      ngamj = 0;
      gam_idxj.clear();
      for (l=0;l<r;l++) {
        if (Gam[j][l]) {
          ngamj += 1;
          gam_idxj.push_back(static_cast<arma::uword>(l));
        }
      }
      Ustar_arma = Ucat_arma.each_col() % ups2_scaling.col(j);
      Ustar_subset = Ustar_arma.cols(arma::uvec(gam_idxj));
      fc_var = inv(Ustar_subset.t() * Ustar_subset + arma::mat(ngamj, ngamj, fill::eye));
      fc_offset = Ustar_subset.t() * Xstar_arma.col(j);

      A.submat(arma::uvec(gam_idxj), arma::uvec(1).fill(j)) = mvnrnd(fc_var * fc_offset, fc_var);
    }



    // Sampling the factors U
    for (m=0;m<Np;m++){
      arma::mat cov_matrix = (A.each_row() / ups2.row(m)) * A.t() + inv_Sigma;
      fc_var = inv(cov_matrix + 1e-6 * arma::eye(r, r));
      for (i=0;i<n[m];i++) {
        fc_offset = (A.each_row() / ups2.row(m)) * X(m).row(i).t() + inv_Sigma * mu.col(zeta(m)(i) - 1);
        U_arma(m).row(i) = mvnrnd(fc_var * fc_offset, fc_var).t();
      }
    }
    Ucat_arma = concatenateRowWise(U_arma);
    Ucat = convertToPtrArray(Ucat_arma);



    // Sampling the cell type labels zeta
    for (m=0;m<Np;m++) {
      zeta(m) = zetalabelupdate(n[m], U_arma(m), mu, inv_Sigma, C, zeta(m), kappa(m), theta);
    }



    // Sampling the cell type means mu
    for (c=1;c<=C;c++) {
      tmc = 0;
      tumic = arma::vec(r, fill::zeros);
      for (m=0;m<Np;m++) {
        for (i=0;i<n[m];i++) {
          if (zeta(m)(i) == c) {
            tmc += 1;
            tumic = tumic + U_arma(m).row(i).t();
          }
        }
      }

      arma::mat cov_matrix = inv_Sigma * tmc + inv_Sigma0;
      fc_var = inv(cov_matrix + 1e-6 * arma::eye(r, r));

      fc_offset = inv_Sigma * tumic + inv_Sigma0 * dmean;

      mu.col(c - 1) = mvnrnd(fc_var * fc_offset, fc_var);
    }



    // // Sampling the cell type hypermeans dmean
    // dmean = mvnrnd(arma::sum(mu, 1) / C, Sigma0 / C);



    // // Sampling the cell type mean hypercovariance matrix Sigma0
    // Sigma_update_outer_prod = arma::mat(r, r, fill::zeros);
    // for (c=1;c<=C;c++) {
    //   Sigma_update_outer_prod += (mu.col(c - 1) - dmean) * (mu.col(c - 1) - dmean).t();
    // }
    // Sigma0 = iwishrnd(Sigma_update_outer_prod + Omega0, C + df0);
    // inv_Sigma0 = inv(Sigma0);



    // Sampling the factor covariance matrix Sigma
    Sigma_update_outer_prod = arma::mat(r, r, fill::zeros);
    for (m=0;m<Np;m++) {
      for (i=0;i<n[m];i++) {
        Sigma_update_outer_prod += (U_arma(m).row(i).t() - mu.col(zeta(m)(i) - 1)) * (U_arma(m).row(i).t() - mu.col(zeta(m)(i) - 1)).t();
      }
    }
    Sigma = iwishrnd(Sigma_update_outer_prod + Omega0 + 1e-6 * arma::eye(r, r), nsum + df0);
    inv_Sigma = inv(Sigma);



    // Sampling the cell type proportions theta
    tmc_mat = arma::mat(C, K, fill::zeros);
    for (m=0;m<Np;m++) {
      for (i=0;i<n[m];i++) {
        tmc_mat(zeta(m)(i) - 1, kappa(m)(i) - 1) += 1;
      }
    }
    for (k=1;k<=K;k++) {
      dirichlet_vec = tmc_mat.col(k - 1) + c0;
      gsl_ran_dirichlet(rr, static_cast<std::size_t>(C), dirichlet_vec.memptr(), dirichlet_output.memptr());
      theta.col(k - 1) = dirichlet_output;
    }



    // Sampling the spatial domain labels kappa
    for (m=0;m<Np;m++) {
      kappa(m) = kappalabelupdate(Wtriplet_list(m), Wbegfin_list(m), n[m], K, beta[m], zeta(m), kappa(m), theta);
    }



    // Sampling the Potts smoothing parameter beta
    for (m=0;m<Np;m++) {

      arma::uvec idx1 = arma::conv_to<arma::uvec>::from( Wtriplet_list(m).col(0) - 1 );
      arma::uvec idx2 = arma::conv_to<arma::uvec>::from( Wtriplet_list(m).col(1) - 1 );
      pairwise_count = arma::accu( kappa(m).elem(idx1) == kappa(m).elem(idx2) ) / 2.0;

      logprob_beta_current = logd_current(m) + beta(m) * pairwise_count;
      // beta_proposal = arma::randu( distr_param(beta(m) - 0.1, beta(m) + 0.1) );
      double u = gsl_ran_flat(rr, 0, 1);
      beta_proposal = (beta(m) - 0.1) + 0.2 * u;
      for(int k = 1; k <= p_order+1; k++) {
        X2(k-1) = pow(beta_proposal, k) / double(k);
      }
      logd_proposal = -arma::dot(X2, interp_coef_mat.col(m)); // -log(d(beta))
      logprob_beta_proposal = logd_proposal + beta_proposal * pairwise_count;
      prob = exp(logprob_beta_proposal - logprob_beta_current);
      prob = pow(prob, temperature); // simulated annealing

      // Accept or reject the proposal
      if(prob > gsl_ran_flat(rr, 0, 1) && beta_proposal >= min_beta && beta_proposal <= max_beta)
      {
        beta[m] = beta_proposal;
        logd_current[m] = logd_proposal;
      }

    }



    if (t>=burninsample){
      // Update posterior mean estimates

      int rm=0;
      for (j=0;j<P;j++){
        for (l=0;l<r;l++){
          rhomodel[t-burninsample][rm]=Gam[j][l];
          rm++;
        }
      }

      for (l=0;l<r;l++){
        rhomodel[t-burninsample][rm]=rhoest[l];
        rm++;
      }

      for (m=0;m<Np;m++){
        for (l=0;l<r;l++){
          for (i=0;i<n[m];i++){
            meanU[m][i][l]+=U_arma(m)(i,l)/nbrsample;
          }
        }
      }

      numdeg_pow2r = PzerosP;
      for (j=0;j<P;j++){
        for (m=0;m<Np;m++){
          ups2Mean(m,j)+=ups2(m,j)/nbrsample;
        }

        if (!streamline) {
          for (int power_idx=1;power_idx<(rpowerset.size() - 1);power_idx++){ // skipping the first and last elements of rpowerset
            int xx=1;
            for (int comp_idx : rpowerset[power_idx]) {
              xx*=1-Gam[j][comp_idx];
            }
            if (xx==0) {
              VarSelMeanGlobal_powerset.col(power_idx)(j) += 1.0/nbrsample;
              numdeg_pow2r(power_idx) += 1.0;
            }
          }
        }

        int xx=1;
        for (l=0;l<r;l++){
          Gammean[j][l]+=(double) Gam[j][l]/nbrsample;
          xx*=1-Gam[j][l];
        }
        if (xx==0) {
          Gamvs[j]+=1.0/nbrsample;
          if (!streamline) {
            VarSelMeanGlobal_powerset.col(rpowerset.size() - 1)(j) += 1.0/nbrsample; // also assigning it into VarSelMeanGlobal_powerset's last vector too
            numdeg_pow2r(rpowerset.size() - 1) += 1.0;
          }
        }

      }
      if (!streamline) {
        numdeg_powerset.col(t-burninsample) = numdeg_pow2r;
      }

      for (l=0;l<r;l++){
        rhomean[l]+=(double) rhoest[l]/nbrsample;
      }

      // Including the MCMC draw in the respective storage objects
      for (m=0;m<Np;m++) {
        zeta2mcmc(m).col(t-burninsample) = zeta(m);
        kappa2mcmc(m).col(t-burninsample) = kappa(m);
      }
      mu2mcmc.slice(t-burninsample) = mu;
      theta2mcmc.slice(t-burninsample) = theta;
      beta2mcmc.col(t-burninsample) = beta;

      if (t % (n_iter/10) == 1) {
        printf("\n");
        printf("Current iteration is %d\n", t);
      }

      // if (t % (n_iter/100) == 1) { // beware: if you use this if-stmt, you need to increase iterations for the test R function call
      //   printf("Current beta_proposal is %f ", beta_proposal);
      //   printf("Current logd_proposal is %f\n", logd_proposal);
      // }

    } else { // t < burninsample

      // simulated annealing
      temperature = (t + 1) / burninsample;

    }

  }
  printf("\n");



  //////// Post MCMC loop things

  for (l=0;l<r;l++){
    CompoSelMean[l]=rhomean[l];
  }

  for (j=0;j<P;j++){
    for (l=0;l<r;l++){
      VarSelMean[j*r+l]=Gammean[j][l];
    }
  }

  for (j=0;j<P;j++){
    VarSelMeanGlobal[j]=Gamvs[j];
  }



  // Return all of the intermediary values

  if (streamline) {

    return Rcpp::List::create(Rcpp::Named("nsum") = nsum,
                              Rcpp::Named("U") = U_arma,
                              Rcpp::Named("A") = A,
                              Rcpp::Named("CompoSelMean") = CompoSelMean,
                              Rcpp::Named("VarSelMean") = VarSelMean,
                              Rcpp::Named("VarSelMeanGlobal") = VarSelMeanGlobal,
                              Rcpp::Named("ups2Mean") = ups2Mean,
                              Rcpp::Named("misc1") = Tau[1][1],
                                                           Rcpp::Named("misc2") = s2[0],
                                                                                    Rcpp::Named("samples") = List::create(
                                                                                      Rcpp::Named("zeta2mcmc") = zeta2mcmc,
                                                                                      Rcpp::Named("kappa2mcmc") = kappa2mcmc,
                                                                                      Rcpp::Named("mu2mcmc") = mu2mcmc,
                                                                                      Rcpp::Named("theta2mcmc") = theta2mcmc,
                                                                                      Rcpp::Named("beta2mcmc") = beta2mcmc
                                                                                    )
    );

  } else {

    return Rcpp::List::create(Rcpp::Named("nsum") = nsum,
                              Rcpp::Named("U") = U_arma,
                              Rcpp::Named("A") = A,
                              Rcpp::Named("CompoSelMean") = CompoSelMean,
                              Rcpp::Named("VarSelMean") = VarSelMean,
                              Rcpp::Named("VarSelMeanGlobal") = VarSelMeanGlobal,
                              Rcpp::Named("VarSelMeanGlobal_powerset") = VarSelMeanGlobal_powerset, // unstreamlined output
                              Rcpp::Named("numdeg_powerset") = numdeg_powerset, // unstreamlined output
                              Rcpp::Named("ups2Mean") = ups2Mean,
                              Rcpp::Named("misc1") = Tau[1][1],
                                                           Rcpp::Named("misc2") = s2[0],
                                                                                    Rcpp::Named("samples") = List::create(
                                                                                      Rcpp::Named("zeta2mcmc") = zeta2mcmc,
                                                                                      Rcpp::Named("kappa2mcmc") = kappa2mcmc,
                                                                                      Rcpp::Named("mu2mcmc") = mu2mcmc,
                                                                                      Rcpp::Named("theta2mcmc") = theta2mcmc,
                                                                                      Rcpp::Named("beta2mcmc") = beta2mcmc
                                                                                    )
    );

  }

}





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
// /*** R
// library(BIPnet)
// library(dplyr)
// library(here)
// source(here("scripts/real_data_helper_fns.R"))
// */
//
// /*** R
// set.seed(700)
//
// # portion prior to BIP() call
//
// BIPdat=Simulate(scenario=1,setting=1,seed=2) # even though there's a seed parameter for Simulate,
// # still need to set.seed before running it
// dataListP <- list(BIPdat$X1[1:5, 1:4], BIPdat$X2[1:6, 1:4])
//
// dataList=dataListP; r=3
// */
//
// /*** R
// # Function Parameters
// IndicVar = 0; nbrsample=50; burnin=10 # nbrsample=5000; burnin=1000
// # BIP priors
// priorcompselv=c(1,1)
// probvarsel=0.05 # probvarsel=0.5
// */
//
// /*** R
// # portion within BIP() call, but prior to mcmcfunction() call
//
// P <- sapply(dataList, ncol) %>% unique
// if (length(P) != 1) {
//   stop("Each element of dataList must have a number of columns equal to the number of features.")
// }
//
//
//
// Np=length(dataList)
// n=NULL
// MeanData=list()
// SD=list()
// for (i in 1:Np){
//   dataList[[i]]=as.matrix(dataList[[i]])
//   n[i]=nrow(dataList[[i]])
//
//   MeanData[[i]]=apply(dataList[[i]],2,mean)
//   SD[[i]]=apply(dataList[[i]],2,sd)
//   # dataList[[i]]=scale(dataList[[i]],T,T)
//
// }
// datasets=do.call("rbind", dataList)
//
//
//
// if (nbrsample<=20){
//   stop("Please specify a larger number of MCMC iterations")
// }
//
// if (is.null(r)){ # could rbind the datasets and do one PCA truncation instead
//   mysvd=lapply(1:Np, function(i)  svd(dataList[[i]]))
//   mysumsvd=lapply(1:Np, function(i) cumsum(mysvd[[i]]$d)/max(cumsum(mysvd[[i]]$d))*100)
//   KMax=max(unlist(lapply(1:Np, function(i) min(which(mysumsvd[[i]]>=80, arr.ind = TRUE))))) #chooses maximum from Np cumulative proportions
//   r=min(KMax+1,10)
// }
//
//
//
// W1 <- matrix(c(0, 1, 0, 0, 0,
//                1, 0, 1, 0, 0,
//                0, 1, 0, 1, 0,
//                0, 0, 1, 0, 1,
//                0, 0, 0, 1, 0), nrow = 5, ncol = 5, byrow = T)
//
// W2 <- matrix(c(0, 1, 0, 0, 0, 0,
//                1, 0, 1, 0, 0, 0,
//                0, 1, 0, 1, 0, 0,
//                0, 0, 1, 0, 1, 0,
//                0, 0, 0, 1, 0, 1,
//                0, 0, 0, 0, 1, 0), nrow = 6, ncol = 6, byrow = T)
//
// W1.quants <- common.Wcheckformat(W1)
// W2.quants <- common.Wcheckformat(W2)
//
// beta <- rep(2, Np)
//
// C <- 3
// K <- 3
//
// # corresponding to the mcmcfunction() call
//
// # n = as.integer(n); P = as.integer(P); r = as.integer(r); Np = as.integer(Np); datasets = datasets; IndVar = as.integer(IndicVar); nbrsample = as.integer(nbrsample);
// # burninsample = as.integer(burnin); CompoSelMean = as.double(rep(0,r)); VarSelMean = as.double(rep(0,r*P)); VarSelMeanGlobal = as.double(rep(0,P));
// # EstU1 = as.double(rep(0,sum(n)*r)); EstU2 = as.double(rep(0,sum(n)*r)); EstUps2 = matrix(c(rep(0.25, P), rep(1.21, P)), nrow = Np, ncol = P, byrow = T);
// # EstLoad=as.double(rep(0,r*P)); priorcompsel = priorcompselv;
// # probvarsel = as.double(probvarsel)
//
// seed <- 2
// results <- mcmcfn(n = as.integer(n), P = as.integer(P), r = as.integer(r), Np = as.integer(Np), datasets = datasets, IndVar = as.integer(IndicVar), nbrsample = as.integer(nbrsample),
//                   burninsample = as.integer(burnin), CompoSelMean = as.double(rep(0,r)), VarSelMean = as.double(rep(0,r*P)), VarSelMeanGlobal = as.double(rep(0,P)),
//                   priorcompsel = priorcompselv, probvarsel = as.double(probvarsel),
//                   C = C, EstMat1 = matrix(rnorm(r*C), nrow = r, ncol = C), EstMat2 = diag(r), EstIMat1 = list(sample(1:C, size = 5, replace = T), sample(1:C, size = 6, replace = T)),
//                   K = K, EstIMat2 = list(sample(1:K, size = 5, replace = T), sample(1:K, size = 6, replace = T)), EstMat3 = matrix(runif(C*K), nrow = C, ncol = K),
//                   list(W1.quants$W.triplet, W2.quants$W.triplet), list(W1.quants$W.begfin, W2.quants$W.begfin), beta, interp_coef_mat = cbind(1:11, 2:12),
//                   chainNbr = seed)
// */
//
// /*** R
// results$nsum
// results$misc1
// results$misc2
//
// results$CompoSelMean
// results$VarSelMeanGlobal
//
// VarSelMean <- matrix(results$VarSelMean,P,byrow=T)
// t(round(VarSelMean, digits = 3))
//
// results$ups2Mean
// */


