// MCMC for 1d ordinal IRT model
//
// Kevin M. Quinn
// UC Berkeley School of Law
// kquinn@law.berkeley.edu
// 
// Daniel E. Ho
// Stanford Law School
// dho@law.stanford.edu
//
//
// Copyright (C) 2007 Kevin M. Quinn and Daniel E. Ho
// 
//

#ifndef ORDRATING_CC
#define ORDRATING_CC

#include <iostream>

#include "matrix.h"
#include "algorithm.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;

// used to access Ydata like a 2d array 
#define M(ROW,COL,NROWS) (COL*NROWS+ROW)




template <typename RNGTYPE>
void ordrating_impl(rng<RNGTYPE>& stream,
		    double *sampledata, 
		    const int& samplerow, const int& samplecol, 
		    const int* Y, const int& Yrow, const int& Ycol,
		    const int& I, const int& J, const int& K,
		    const int& beta_constraint,
		    const int& theta_neg_index, 
		    const int& theta_pos_index,
		    const double& vinv1, const double& vinv2,
		    const double& m1, const double& m2,
		    const int& burnin, 
		    const int& mcmc, const int& thin, const double& tune,
		    const int& verbose, 
		    double *alpha,
		    double *beta,
		    double *theta,
		    double *gamma
		    ){
  

    // define constants and from cross-product matrices
    const int tot_iter = burnin + mcmc;  // total number of mcmc iterations
    //const int nstore = mcmc / thin;      // number of draws to store
    const int N = Yrow;
    const double vinv1m1 = vinv1 * m1;
    const double vinv2m2 = vinv2 * m2;

    
    // set up data structures 
    // R gives, for each i, the set of indices of things rated by i
    vector< vector<int> > R;
    // S gives, for each j, the set of indices of raters that rated j
    vector< vector<int> > S;

    //double Z[I][J];
    double **Z = new double*[I];
    for (int i=0; i<I; ++i){
      Z[i] = new double[J];
    }

    // fill in Z and R
    for (int i=0; i<I; ++i){
      vector<int> jholder;
      for (int ind=0; ind<N; ++ind){
	if (Y[M(ind,0,Yrow)] == i){
	  jholder.push_back(Y[M(ind,1, Yrow)]);
	  //Z[i][Y[M(ind, 1, *Yrow)]] = (Y[M(ind, 2, *Yrow)] - 1.5)*1;
	  Z[i][Y[M(ind, 1, Yrow)]] = 0.0;
	}
      }
      R.push_back(jholder);
    }
    // fill in S
    for (int j=0; j<J; ++j){
      vector<int> iholder;
      for (int ind=0; ind<N; ++ind){
	if (Y[M(ind,1, Yrow)] == j){
	  iholder.push_back(Y[M(ind,0, Yrow)]);
	}
      }
      S.push_back(iholder);
    }



    // gamma proposal 
    //double gamma_p[K + 1];
    double *gamma_p = new double[K+1];
    for (int k=0; k < (K+1); ++k){
      gamma_p[k] = gamma[k];
    }

    // Gibbs loop
    int count = 0;
    int accepts = 0;
    for (int iter = 0; iter < tot_iter; ++iter){
      
      
      
      if ( K > 2){
	// [gamma | Z, alpha, beta, theta]
	for (int i=2; i< K; ++i){
	  if (i==(K-1)){
	    gamma_p[i] = stream.rtbnorm_combo(gamma[i], std::pow(tune, 2.0), 
					       gamma_p[i-1]);
	  }
	  else {
	    gamma_p[i] = stream.rtnorm_combo(gamma[i], std::pow(tune, 2.0), 
					      gamma_p[i-1], 
					      gamma[i+1]);
	  }
	}
	double loglikerat = 0.0;
	double loggendenrat = 0.0;
	
	
	
	// loop over observations and construct the acceptance ratio
	for (int i=0; i<N; ++i){
	  const double mu = alpha[Y[M(i,0, Yrow)]] + 
	    beta[Y[M(i,0, Yrow)]] * theta[Y[M(i,1,Yrow)]];
	  const int Yi = Y[M(i,2, Yrow)];
	  if (Yi == K){
	    loglikerat = loglikerat 
	      + log(1.0  - 
		    pnorm(gamma_p[Yi-1] - mu, 0.0, 1.0 ) ) 
	      - log(1.0 - 
		    pnorm(gamma[Yi-1] - mu, 0.0, 1.0) );
	  }
	  else if (Yi == 1){
	    loglikerat = loglikerat 
	      + log(pnorm(gamma_p[Yi] - mu, 0.0, 1.0)  ) 
	      - log(pnorm(gamma[Yi] - mu, 0.0, 1.0) );
	  }
	  else{
	    loglikerat = loglikerat 
	      + log(pnorm(gamma_p[Yi] - mu, 0.0, 1.0) - 
		    pnorm(gamma_p[Yi-1] - mu, 0.0, 1.0) ) 
	      - log(pnorm(gamma[Yi] - mu, 0.0, 1.0) - 
		    pnorm(gamma[Yi-1] - mu, 0.0, 1.0) );
	  }
	}
	for (int j=2; j<K; ++j){	   
	  loggendenrat = loggendenrat
	    + log(pnorm(gamma[j+1], gamma[j], tune) - 
		  pnorm(gamma_p[j-1], gamma[j], tune) )  
	    - log(pnorm(gamma_p[j+1], gamma_p[j], tune) - 
		  pnorm(gamma[j-1], gamma_p[j], tune) );
	  
	}
	double logacceptrat = loglikerat + loggendenrat;
	if (stream.runif() <= exp(logacceptrat)){
	  for (int k=0; k < (K+1); ++k){
	    gamma[k] = gamma_p[k];
	  }
	  ++accepts;
	}
      } // end if *K > 2 


      //cout << "Z sampling starts now " << endl;
      
      // [Z| gamma, alpha, beta, theta, y] 
      //Matrix<double> Z_mean = X * beta;
      for (int i=0; i<N; ++i){
	const double mu = alpha[Y[M(i,0, Yrow)]] + 
	  beta[Y[M(i,0, Yrow)]] * theta[Y[M(i,1, Yrow)]];
	const int Yi = Y[M(i,2, Yrow)];
	Z[Y[M(i,0, Yrow)]][Y[M(i,1, Yrow)]] = 
	  stream.rtnorm_combo(mu, 1.0, gamma[Yi-1], gamma[Yi]);
      }
      
      
      // [alpha, beta|Z, theta, gamma]
      if (beta_constraint == 0){
	for (int i=0; i<I; ++i){
	  const int ni = R[i].size();
	  double sumthetajsq = 0.0;
	  double sumthetaj = 0.0;
	  double sumzj = 0.0;
	  double sumthetajzj = 0.0;

	  for (int jj=0; jj<ni; ++jj){
	    const int j = R[i][jj];
	    sumthetajsq += std::pow(theta[j], 2.0);
	    sumthetaj += theta[j];
	    sumzj += Z[i][j];
	    sumthetajzj += Z[i][j] * theta[j];
	  } // end jj loop
	  
	  const double sumthetajallsq = std::pow(sumthetaj, 2.0);
	  const double invdet = 1.0 / ((ni + vinv1) * (sumthetajsq + vinv2) - 
				       sumthetajallsq);
	  const double v11star = invdet * (sumthetajsq + vinv2);
	  const double v12star = invdet * (-1.0 * sumthetaj);
	  const double v22star = invdet * (ni + vinv2);
	  const double s1star = std::sqrt(v11star);
	  const double s2star = std::sqrt(v22star);
	  const double rho = v12star / (s1star * s2star);
	  const double holder1 = sumzj + vinv1m1;
	  const double holder2 = sumthetajzj + vinv2m2;
	  const double m1star = v11star * holder1 + v12star * holder2;
	  const double m2star = v12star * holder1 + v22star * holder2;
	  // (alpha[i], beta[i]) ~ 
	  //     N((m1star, m2star), c(v11star, v12star, v22star) )
	  alpha[i] = stream.rnorm(m1star, s1star);
	  const double cond_mean = m2star - m1star * (v12star / v11star) + 
	    alpha[i] * (v12star / v11star);
	  const double cond_sd = std::sqrt(v22star * ( 1 - std::pow(rho, 2.0)));
	  beta[i] = stream.rnorm(cond_mean, cond_sd);
	} // end i loop
      }// end if beta_constrain == 0 
      else if (beta_constraint > 0){
	for (int i=0; i<I; ++i){
	  const int ni = R[i].size();
	  double sumthetajsq = 0.0;
	  double sumthetaj = 0.0;
	  double sumzj = 0.0;
	  double sumthetajzj = 0.0;
	  
	  for (int jj=0; jj<ni; ++jj){
	    const int j = R[i][jj];
	    sumthetajsq += std::pow(theta[j], 2.0);
	    sumthetaj += theta[j];
	    sumzj += Z[i][j];
	    sumthetajzj += Z[i][j] * theta[j];
	  } // end jj loop
	  
	  const double sumthetajallsq = std::pow(sumthetaj, 2.0);
	  const double invdet = 1.0 / ((ni + vinv1) * (sumthetajsq + vinv2) - 
				       sumthetajallsq);
	  const double v11star = invdet * (sumthetajsq + vinv2);
	  const double v12star = invdet * (-1.0 * sumthetaj);
	  const double v22star = invdet * (ni + vinv2);
	  const double s1star = std::sqrt(v11star);
	  const double s2star = std::sqrt(v22star);
	  const double rho = v12star / (s1star * s2star);
	  const double holder1 = sumzj + vinv1m1;
	  const double holder2 = sumthetajzj + vinv2m2;
	  const double m1star = v11star * holder1 + v12star * holder2;
	  const double m2star = v12star * holder1 + v22star * holder2;
	  // (alpha[i], beta[i]) ~ 
	  //     N((m1star, m2star), c(v11star, v12star, v22star) )
	  alpha[i] = stream.rnorm(m1star, s1star);
	  const double cond_mean = m2star - m1star * (v12star / v11star) + 
	    alpha[i] * (v12star / v11star);
	  const double cond_var = v22star * ( 1 - std::pow(rho, 2.0));
	  beta[i] = stream.rtbnorm_combo(cond_mean, cond_var, 0);
	} // end i loop
      }// end if *beta_constraint > 0
      else if (beta_constraint < 0){
	for (int i=0; i<I; ++i){
	  const int ni = R[i].size();
	  double sumthetajsq = 0.0;
	  double sumthetaj = 0.0;
	  double sumzj = 0.0;
	  double sumthetajzj = 0.0;
	  
	  for (int jj=0; jj<ni; ++jj){
	    const int j = R[i][jj];
	    sumthetajsq += std::pow(theta[j], 2.0);
	    sumthetaj += theta[j];
	    sumzj += Z[i][j];
	    sumthetajzj += Z[i][j] * theta[j];
	  } // end jj loop
	  
	  const double sumthetajallsq = std::pow(sumthetaj, 2.0);
	  const double invdet = 1.0 / ((ni + vinv1) * (sumthetajsq + vinv2) - 
				       sumthetajallsq);
	  const double v11star = invdet * (sumthetajsq + vinv2);
	  const double v12star = invdet * (-1.0 * sumthetaj);
	  const double v22star = invdet * (ni + vinv2);
	  const double s1star = std::sqrt(v11star);
	  const double s2star = std::sqrt(v22star);
	  const double rho = v12star / (s1star * s2star);
	  const double holder1 = sumzj + vinv1m1;
	  const double holder2 = sumthetajzj + vinv2m2;
	  const double m1star = v11star * holder1 + v12star * holder2;
	  const double m2star = v12star * holder1 + v22star * holder2;
	  // (alpha[i], beta[i]) ~ 
	  //     N((m1star, m2star), c(v11star, v12star, v22star) )
	  alpha[i] = stream.rnorm(m1star, s1star);
	  const double cond_mean = m2star - m1star * (v12star / v11star) + 
	    alpha[i] * (v12star / v11star);
	  const double cond_var = v22star * ( 1 - std::pow(rho, 2.0));
	  beta[i] = stream.rtanorm_combo(cond_mean, cond_var, 0);
	} // end i loop
      }// end if *beta_constraint < 0
      



      // [theta | Z, alpha, beta, gamma]
      for (int j=0; j<J; ++j){
	const int nj = S[j].size();
	double sumbetaisq = 0.0;
	double sumbetaizi = 0.0;
	
	for (int ii=0; ii<nj; ++ii){
	  const int i = S[j][ii];
	  sumbetaisq += std::pow(beta[i], 2.0);
	  sumbetaizi += beta[i] * (Z[i][j] - alpha[i]);
	} // end ii loop

	const double vstar = 1.0 / (sumbetaisq + 1.0);
	const double sstar = std::sqrt(vstar);
	const double mstar = vstar * sumbetaizi;
	// theta[j] ~ N(mstar, vstar)
	theta[j] = stream.rnorm(mstar, sstar);

	if (j == theta_neg_index){
	  theta[j] = stream.rtanorm_combo(mstar, vstar, 0);
	}
	else if (j == theta_pos_index){
	  theta[j] = stream.rtbnorm_combo(mstar, vstar, 0);
	}
	
      } // end j loop

      
      
      
      // store values in matrices
      if (iter >= burnin && ((iter % thin)==0)){ 
	for (int i=0; i<I; ++i){
	  sampledata[M(count, i, samplerow)] = alpha[i];
	  sampledata[M(count, (i+ I), samplerow)] = beta[i]; 
	}
	for (int j=0; j<J; ++j){
	  sampledata[M(count, (j + 2* I), samplerow)] = theta[j];
	}
	for (int k=1; k< K; ++k){
	  sampledata[M(count, (k-1 +J + I*2), samplerow)] = gamma[k];
	}	
	++count;
      }
      
      // print output to stdout
      if(verbose > 0 && iter % verbose == 0){
	Rprintf("\n\niteration %i of %i \n", (iter+1), tot_iter);
	Rprintf("Metropolis acceptance rate for gamma = %3.5f\n\n", 
		static_cast<double>(accepts) / 
		static_cast<double>(iter+1));		
	/*
	cout << "theta = " << endl;
	for (int j=0; j<*J; ++j){
	  cout << theta[j] << "  " ;
	}
	cout << endl << "gamma[2] = " << gamma[2] << endl;
	cout << endl << endl;
	*/
      }
     
      R_CheckUserInterrupt(); // allow user interrupts           
    }

    for (int i=0; i<I; ++i){
      delete[] Z[i];
    }
    delete[] Z;
    delete[] gamma_p;
    
    Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    Rprintf("The Metropolis acceptance rate for gamma was %3.5f", 
	    static_cast<double>(accepts) / static_cast<double>(tot_iter));
    Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

  }












extern "C"{

  
  void ordratingpost(double *sampledata, const int *samplerow, 
		     const int *samplecol, 
		     const int *Y, const int *Yrow, const int *Ycol,
		     const int *I, const int *J, const int *K,
		     const int *beta_constraint,
		     const int *theta_neg_index, const int *theta_pos_index,
		     const double *vinv1, const double *vinv2,
		     const double *m1, const double *m2,
		     const int *burnin, 
		     const int *mcmc, const int *thin, const double* tune,
		     const int *verbose, 
		     double *alpha,
		     double *beta,
		     double *theta,
		     double *gamma,
		     const int *uselecuyer, 
		     const int *seedarray, 
		     const int *lecuyerstream) {  
    
    

    MCMCPACK_PASSRNG2MODEL(ordrating_impl, 
			   sampledata, *samplerow, *samplecol, 
			   Y, *Yrow, *Ycol,
			   *I, *J, *K,
			   *beta_constraint,
			   *theta_neg_index, 
			   *theta_pos_index,
			   *vinv1, *vinv2,
			   *m1, *m2,
			   *burnin, 
			   *mcmc, *thin, *tune,
			   *verbose, 
			   alpha,
			   beta,
			   theta,
			   gamma);
    
  }
  
}

#endif

