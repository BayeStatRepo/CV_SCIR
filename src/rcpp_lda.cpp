#include <RcppArmadillo.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <numeric>
#include <limits>


using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void print_it(int t) {
  Rcpp::Environment base_env("package:base");
  
  Rcpp::Function print_r = base_env["print"];
  Rcpp::Function paste_r = base_env["paste"];
  print_r(paste_r("Iteration:", t));
}

//[[Rcpp::export]]
Rcpp::NumericVector rcpp_words_idx(NumericVector doc){
  // Environment myEnv = Environment::global_env();
  Environment myEnv = Environment::namespace_env("LDARcppPackage");
  Function words_idx = myEnv["words_idx"];
  NumericVector idx = words_idx(doc);
  return idx;
}

//[[Rcpp::export]]
List rcpp_mat_init(int K, int W, int D){
  Environment myEnv = Environment::namespace_env("LDARcppPackage");
  Function mat_init = myEnv["mat_init"];
  List mats = mat_init(K, W, D);
  return mats;
}

//[[Rcpp::export]]
Rcpp::NumericVector rcpp_cumprob_ini(int K){
  // Environment myEnv = Environment::global_env();
  Environment myEnv = Environment::namespace_env("LDARcppPackage");
  Function cumprob_ini = myEnv["cumprob_ini"];
  NumericVector cumprob = cumprob_ini(K);
  return cumprob;
}

//[[Rcpp::export]]
Rcpp::NumericVector rcpp_runif(int n){
  // Environment myEnv = Environment::global_env();
  Environment myEnv = Environment::namespace_env("LDARcppPackage");
  Function runif_fun = myEnv["runif_fun"];
  NumericVector u = runif_fun(n);
  return u;
}

// [[Rcpp::export]]
List rcpp_sample_cnt(double alpha, NumericVector probs_init,
                     List docs_list, NumericMatrix theta,  
                     int gibbs_iters, int gibbs_burnin){
  
  // Obtaining namespace of base package
  Environment base = Environment::namespace_env("base");
  // Picking up sample() function from base package
  Function sample = base["sample"];
  
  NumericVector temp_doc = docs_list[0];
  int W = temp_doc.size();
  int D = docs_list.size();
  int K = probs_init.size();
  
  List mats=rcpp_mat_init(K, W, D);
  
  NumericMatrix Adk=mats[0];
  NumericMatrix Adk_avg=mats[1];
  NumericMatrix Bkw=mats[2]; 
  NumericMatrix Bkw_avg=mats[3]; 
 
  IntegerVector topic_vec = Rcpp::seq(0, K - 1);
  
  Rcpp::List z(D);
  
  // Initialise the z_id for each document in the batch
  for(int d = 0; d < D; ++d){
    NumericVector doc = docs_list[d];
    
    NumericVector idx_w_d=rcpp_words_idx(doc);
    
    Rcpp::List z_w(idx_w_d.size());
    
    for(int v=0; v < idx_w_d.size(); ++v){
      int w=idx_w_d[v];
      
      int word_cnt = doc[w];
      NumericVector zdw(word_cnt, 0.0);
      
      for(int i = 0; i < zdw.size(); ++i){
        int zInit = as<int>(sample(Named("x")=topic_vec,  Named("size")=1, 
                                   Named("replace", false), 
                                   Named("prob")=as<NumericVector>(probs_init)));
        Adk(d, zInit) = Adk(d, zInit) + 1;
        Bkw(zInit, w) = Bkw(zInit, w) + 1;
        zdw[i] = zInit;
      }
      z_w[v] = zdw;
    }
    z[d] = z_w;
  }  
  
  // Draw samples from the posterior on z_ids using Gibbs sampling
  for(int t=0; t < gibbs_iters; ++t){
    for(int d = 0; d < D; ++d){
      
      NumericVector doc = docs_list[d];
      
      NumericVector idx_w_d=rcpp_words_idx(doc);
      
      List z_w=z[d];
      
      for(int v=0; v < idx_w_d.size(); ++v){
        int w=idx_w_d[v];
        
        // int word_cnt = doc[w];
        
        NumericVector zdw = z_w[v];
        for(int i = 0; i < zdw.size(); ++i){
          int zOld = zdw[i];
          
          NumericVector un_probs(K);
          NumericVector probs(K);
          
          double sum_p=0;
          for(int k=0; k < K; ++k){
            un_probs[k] = (alpha + Adk(d,k) - (k == zOld)) * theta(k,w);
            sum_p = sum_p + un_probs[k];
          }
          for(int k=0; k < K; ++k){
            probs[k] = un_probs[k]/sum_p;
          } 
          
          int zNew = as<int>(sample(Named("x")=topic_vec,  Named("size")=1, 
                             Named("replace", false), 
                             Named("prob")=as<NumericVector>(probs)));
          zdw[i] = zNew;
          
          Adk(d,zOld)     = Adk(d,zOld) - 1;
          Adk(d,zNew)     = Adk(d,zNew) + 1;
          Bkw(zOld,w)     = Bkw(zOld,w) - 1;
          Bkw(zNew,w)     = Bkw(zNew,w) + 1;
          
        }
        z_w[v] = zdw;
      }
      z[d] = z_w;
    }
    
    if(t >= gibbs_burnin){
      
      for(int k=0; k < K; ++k){
        for(int d = 0; d < D; ++d){
          Adk_avg(d, k) = Adk_avg(d, k) + Adk(d, k);
        }
        for(int w=0; w < W; ++w){
          Bkw_avg(k, w) = Bkw_avg(k, w) + Bkw(k, w);
        }  
      }
      
    }
    print_it(t+1);
    
  }

  // Creating a List for output
  List out(2);
  
  out[0]=Adk_avg;
  out[1]=Bkw_avg;
  return out;
}







