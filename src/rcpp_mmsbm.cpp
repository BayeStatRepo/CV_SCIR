#include <RcppArmadillo.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <numeric>
#include <limits>
#include <math.h>
#include <Rmath.h>
#include <thread>
#include <random>
#include <mutex>  // For std::mutex

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
int sample_index(const arma::vec& P_marg_z_ij) {
  static thread_local std::mt19937 gen(std::random_device{}()); // Efficient and safe in RcppParallel
  
  std::discrete_distribution<> dist(P_marg_z_ij.begin(), P_marg_z_ij.end());
  return dist(gen);
}
//[[Rcpp::export]]
arma::uvec sample_index_from_joint(const arma::mat& P_joint) {
  arma::vec flat_P = arma::vectorise(P_joint);  // Flatten to 1D vector (column-major order)
  arma::vec cum_P = arma::cumsum(flat_P);       // Cumulative sum
  
  // double u = R::runif(0.0, 1.0);                // Uniform [0, 1)
  double u = arma::randu<double>();
  
  // arma::uword idx = arma::as_scalar(arma::find(cum_P >= u, 1));  // Find where u falls
  arma::uvec found_idx = arma::find(cum_P >= u, 1);
  arma::uword idx = found_idx.is_empty() ? cum_P.n_elem - 1 : found_idx(0);
  
  // Convert linear index to row and column (remember: column-major)
  arma::uword row = idx % P_joint.n_rows;
  arma::uword col = idx / P_joint.n_rows;
  
  return arma::uvec({row, col});  // Return as vector of 2 elements
}

// NEW FUNCTIONS CORRECTLY PARAMETRIZED
//[[Rcpp::export]]
Rcpp::List CIR_B_update_grad_ij(arma::mat theta_1, arma::mat theta_0, arma::vec pi_i, arma::vec pi_j, int Y_ij){
  
  int K = theta_1.n_rows;
  
  arma::mat P_joint_z_ij_zji_unnorm(K, K, arma::fill::zeros);  
  for(int k = 0; k < K; ++k){
    for(int h = 0; h < K; ++h){
      P_joint_z_ij_zji_unnorm(k,h) = ( Y_ij * theta_1(k,h) + (1-Y_ij) * theta_0(k,h) ) * pi_i[k] * pi_j[h];
    }
  }
  
  arma::mat P_joint_z_ij_zji = P_joint_z_ij_zji_unnorm / arma::accu(P_joint_z_ij_zji_unnorm);  // normalize by total sum
  
  
  arma::uvec zij_zji = sample_index_from_joint(P_joint_z_ij_zji);
  
  // Assign the row and column to two integers
  int row = zij_zji(0);  // First element of the vector
  int col = zij_zji(1);  // Second element of the vector
  
  arma::mat eta_par_mat_0(K, K, arma::fill::zeros);  
  arma::mat eta_par_mat_1(K, K, arma::fill::zeros);  
  
  eta_par_mat_0(row,col) = (1-Y_ij);
  eta_par_mat_1(row,col) = (Y_ij);
  
  List out(2);
  out[0]=eta_par_mat_0;
  out[1]=eta_par_mat_1;
  
  return out;
}
//[[Rcpp::export]]
Rcpp::List CIR_B_update_no_par(arma::mat theta_1, arma::mat theta_0, arma::mat pi, arma::mat Y, arma::mat subE){
  
  int K = theta_1.n_rows;
  int n = subE.n_rows;
  
  arma::mat eta_par_mat_0(K, K, arma::fill::zeros);  // Initialize the result matrices
  arma::mat eta_par_mat_1(K, K, arma::fill::zeros);
  
  for (int r = 0; r < n; ++r) {
    // Computation
    int i = subE(r,0);
    int j = subE(r,1);
    
    arma::vec pi_i = pi.row(i).t();
    arma::vec pi_j = pi.row(j).t();
    
    int Y_ij = Y(i,j);
    
    List grad_ij = CIR_B_update_grad_ij(theta_1, theta_0, pi_i, pi_j, Y_ij);
    
    eta_par_mat_0 += as<arma::mat>(grad_ij[0]);  // Accumulate the results across all threads
    eta_par_mat_1 += as<arma::mat>(grad_ij[1]);
  }
  

  List out(2);
  out[0] = eta_par_mat_0;
  out[1] = eta_par_mat_1;
  
  return out;
}
//[[Rcpp::export]]
Rcpp::List CIR_pi_update_grad_i(int i, arma::mat B, arma::mat phi, arma::vec Y_i_r, arma::vec Y_i_c, arma::vec subV_i){
  
  int K = B.n_rows;
  int n = subV_i.n_elem;
  
  arma::vec phi_i = phi.row(i).t();
  
  arma::vec sum_z_ij(K, arma::fill::zeros);
  arma::vec sum_z_ji(K, arma::fill::zeros);
  
  for(int r = 0; r < n; ++r){
    int j = subV_i[r];
    
    arma::vec phi_j = phi.row(j).t();
    
    int Y_ij = Y_i_r[j];
    int Y_ji = Y_i_c[j];
    
    arma::vec P_marg_z_ij_unnorm(K, arma::fill::zeros);
    arma::vec P_marg_z_ji_unnorm(K, arma::fill::zeros);
    
    for(int k = 0; k < K; ++k){
      double temp_ij = 0;
      double temp_ji = 0;
      for(int h = 0; h < K; ++h){
        temp_ij += ( Y_ij*B(k,h) + (1-Y_ij)*(1-B(k,h)) ) * phi_j[h];
        temp_ji += ( Y_ji*B(h,k) + (1-Y_ji)*(1-B(h,k)) ) * phi_j[h];
      }
      P_marg_z_ij_unnorm[k] = phi_i[k] * temp_ij;
      P_marg_z_ji_unnorm[k] = phi_i[k] * temp_ji;
    }
    
    arma::vec P_marg_z_ij = P_marg_z_ij_unnorm / arma::sum(P_marg_z_ij_unnorm);  // normalize by total sum
    arma::vec P_marg_z_ji = P_marg_z_ji_unnorm / arma::sum(P_marg_z_ji_unnorm);  // normalize by total sum
    
    int z_ij = sample_index(P_marg_z_ij);
    int z_ji = sample_index(P_marg_z_ji);
    
    // Rcpp::Rcout << "Sampled z_ij = " << z_ij << std::endl;
    // Rcpp::Rcout << "Sampled z_ji = " << z_ji << std::endl;
    
    // Rcpp::Rcout << "P_marg_z_ij = " << P_marg_z_ij_unnorm.t() << std::endl;
    // Rcpp::Rcout << "P_marg_z_ji = " << P_marg_z_ji_unnorm.t() << std::endl;
    
    arma::vec z_ij_vec(K, arma::fill::zeros);
    arma::vec z_ji_vec(K, arma::fill::zeros);
    
    z_ij_vec[z_ij] = 1;
    z_ji_vec[z_ji] = 1;
    
    sum_z_ij += z_ij_vec;
    sum_z_ji += z_ji_vec;
    
  }
  
  List out(2);
  out[0] = sum_z_ij;
  out[1] = sum_z_ji;
  
  return out;
}
//[[Rcpp::export]]
arma::mat CIR_pi_update_no_par(arma::mat B, arma::mat phi, arma::mat Y, arma::vec sub_V_xi, Rcpp::List sub_V_i_list){
  
  int K = B.n_rows;
  int N = Y.n_rows;
  
  int n = sub_V_xi.n_elem;
  
  arma::mat alpha_par_mat(N, K, arma::fill::zeros); 
  
  for (int r = 0; r < n; ++r) {
    // Computation
    int i = sub_V_xi[r];
    
    arma::vec subV_i = as<arma::vec>(sub_V_i_list[r]);
    
    arma::rowvec Y_i_r = Y.row(i); 
    arma::vec Y_i_c = Y.col(i);
    
    Rcpp::List grad_i = CIR_pi_update_grad_i(i, B, phi, Y_i_r.t(), Y_i_c, subV_i);
    
    arma::vec sum_z_ij = as<arma::vec>(grad_i[0]);
    arma::vec sum_z_ji = as<arma::vec>(grad_i[1]);
    
    arma::rowvec temp = sum_z_ij.t() + sum_z_ji.t();
    alpha_par_mat.row(i) = temp;
    
  }
  
  return alpha_par_mat;
  
}

//[[Rcpp::export]]
Rcpp::List grad_SG_global_update_ij( arma::mat theta_0, arma::mat theta_1, arma::mat B, arma::vec pi_i, arma::vec pi_j, int Y_ij){
  
  int K = B.n_rows;
  
  arma::mat P_joint_z_ij_zji_unnorm(K, K, arma::fill::zeros);
  for(int k = 0; k < K; ++k){
    for(int h = 0; h < K; ++h){
      P_joint_z_ij_zji_unnorm(k,h) = ( Y_ij * B(k,h) + (1-Y_ij) * (1-B(k,h) )) * pi_i[k] * pi_j[h];
    }
  }
  
  arma::mat P_joint_z_ij_zji = P_joint_z_ij_zji_unnorm / arma::accu(P_joint_z_ij_zji_unnorm);  // normalize by total sum
  
  
  arma::uvec zij_zji = sample_index_from_joint(P_joint_z_ij_zji);
  
  // Assign the row and column to two integers
  int row = zij_zji(0);  // First element of the vector
  int col = zij_zji(1);  // Second element of the vector
  
  arma::mat eta_par_mat_0(K, K, arma::fill::zeros);
  arma::mat eta_par_mat_1(K, K, arma::fill::zeros);
  
  eta_par_mat_0(row,col) = (1-Y_ij)/theta_0(row,col) - 1/(theta_0(row,col)+theta_1(row,col));
  eta_par_mat_1(row,col) = (  Y_ij)/theta_1(row,col) - 1/(theta_0(row,col)+theta_1(row,col));
  
  List out(2);
  out[0]=eta_par_mat_0;
  out[1]=eta_par_mat_1;      
        
  return(out);
}
//[[Rcpp::export]]
Rcpp::List SG_B_update_no_par(arma::mat theta_0, arma::mat theta_1, arma::mat B, arma::mat pi, arma::mat Y, arma::mat subE){
  
  int K = B.n_rows;
  int n = subE.n_rows;
  
  arma::mat eta_par_mat_0(K, K, arma::fill::zeros);  // Initialize the result matrices
  arma::mat eta_par_mat_1(K, K, arma::fill::zeros);
  
  for (int r = 0; r < n; ++r) {
    // Computation
    int i = subE(r,0);
    int j = subE(r,1);
    
    arma::vec pi_i = pi.row(i).t();
    arma::vec pi_j = pi.row(j).t();
    
    int Y_ij = Y(i,j);
    
    List grad_ij = grad_SG_global_update_ij(theta_0, theta_1, B, pi_i, pi_j, Y_ij);
    
    eta_par_mat_0 += as<arma::mat>(grad_ij[0]);  // Accumulate the results across all threads
    eta_par_mat_1 += as<arma::mat>(grad_ij[1]);
  }
  
  
  List out(2);
  out[0] = eta_par_mat_0;
  out[1] = eta_par_mat_1;
  
  return out;
}
//[[Rcpp::export]]
arma::mat grad_SG_local_update_i(int i, arma::vec phi_i, arma::mat B, arma::mat pi, arma::vec Y_i_p, arma::vec Y_p_i, arma::vec subV_i){
  
  int K = B.n_rows;
  int n = subV_i.n_elem;
  
  arma::vec pi_i = pi.row(i).t();
  
  arma::vec sum_z_ij(K);
  arma::vec sum_z_ji(K);
  sum_z_ij.fill(0);
  sum_z_ji.fill(0);
  sum_z_ij.fill(-n/arma::sum(phi_i));
  sum_z_ji.fill(-n/arma::sum(phi_i));

  for (int q = 0; q < n; ++q) {
    int j = subV_i[q];
    arma::vec pi_j = pi.row(j).t();
    
    int Y_ij = Y_i_p[j];
    int Y_ji = Y_p_i[j];
    
    arma::vec P_marg_z_ij_unnorm(K, arma::fill::zeros);
    arma::vec P_marg_z_ji_unnorm(K, arma::fill::zeros);
    
    for(int k = 0; k < K; ++k){
      double temp_ij = 0;
      double temp_ji = 0;
      for(int h = 0; h < K; ++h){
        temp_ij += ( Y_ij*B(k,h) + (1-Y_ij)*(1-B(k,h)) ) * pi_j[h];
        temp_ji += ( Y_ji*B(h,k) + (1-Y_ji)*(1-B(h,k)) ) * pi_j[h];
      }
      P_marg_z_ij_unnorm[k] = pi_i[k] * temp_ij;
      P_marg_z_ji_unnorm[k] = pi_i[k] * temp_ji;
    }
    
    arma::vec P_marg_z_ij = P_marg_z_ij_unnorm / arma::sum(P_marg_z_ij_unnorm);  // normalize by total sum
    arma::vec P_marg_z_ji = P_marg_z_ji_unnorm / arma::sum(P_marg_z_ji_unnorm);  // normalize by total sum
    
    int z_ij = sample_index(P_marg_z_ij);
    int z_ji = sample_index(P_marg_z_ji);
    
    arma::vec z_ij_vec(K, arma::fill::zeros);
    arma::vec z_ji_vec(K, arma::fill::zeros);
    
    z_ij_vec[z_ij] = 1/phi_i[z_ij];
    z_ji_vec[z_ji] = 1/phi_i[z_ji];
    
    sum_z_ij += z_ij_vec;
    sum_z_ji += z_ji_vec;
    
    // z_ij_vec[z_ij] = 1;
    // z_ji_vec[z_ji] = 1;
    // sum_z_ij += z_ij_vec/phi_i[z_ij];
    // sum_z_ji += z_ji_vec/phi_i[z_ji];
    
  }
  
  // arma::mat sum_z = arma::join_cols(sum_z_ij, sum_z_ji);
  arma::mat sum_z = arma::join_cols(sum_z_ij.t(), sum_z_ji.t());  // 2 x K
  
  return(sum_z);
      
}
//[[Rcpp::export]]
arma::mat SG_pi_update_no_par(arma::mat phi, arma::mat B, arma::mat pi, arma::mat Y, arma::vec V_all, Rcpp::List sub_V_xi){
  
  int N = Y.n_rows;
  int K = B.n_rows;
  
  int n = V_all.n_elem;
  
  arma::mat alpha_par_mat(N, K, arma::fill::zeros);

  for(int r = 0; r < n; ++r){
    
    int i = V_all[r];
    
    arma::vec subV_i = sub_V_xi[r];
    
    arma::vec phi_i = phi.row(i).t();

    arma::vec Y_i_p = Y.row(i).t();
    arma::vec Y_p_i = Y.col(i);
    
    arma::mat grad_i = grad_SG_local_update_i(i, phi_i, B, pi, Y_i_p, Y_p_i, subV_i);
    
    
    arma::rowvec sum_z_ij = grad_i.row(0);
    arma::rowvec sum_z_ji = grad_i.row(1);
    
    alpha_par_mat.row(i) = sum_z_ij + sum_z_ji;
    
  }
  return(alpha_par_mat);
}


double P_ij_fun(int Y_ij, arma::vec pi_i, arma::vec pi_j, arma::mat B) {
  
  double p_0_ij = arma::as_scalar(pi_i.t() * B * pi_j);
  
  double p_ij = Y_ij*( p_0_ij ) + (1-Y_ij)*( 1-p_0_ij );
  return(p_ij);
}

//[[Rcpp::export]]
double P_0_fun(arma::mat Y, arma::mat pi_0, arma::mat B_0){
  
  int N = Y.n_rows;
  
  int M = N*(N-1);
  
  arma::vec p_ij_0(M, arma::fill::zeros);
  
  
  int count = 0; 
  for (int i = 0; i < N; ++i){
    
    arma::vec pi_i_0 = pi_0.row(i).t(); // column vector
    
    for (int j = 0; j < N; ++j){
      
      if(i == j) continue;
      
      arma::vec p_j_0 = pi_0.row(j).t(); // column vector
      
      int Y_ij = Y(i,j);
      
      p_ij_0[count]       = P_ij_fun(Y_ij, pi_i_0, p_j_0, B_0);
      
      count = count + 1;
      
    }
  }
  double log_P_0_Y = arma::sum(arma::log(p_ij_0));
  
  return( log_P_0_Y );
}

arma::vec P_fun_old(arma::mat Y, arma::mat pi_scir, arma::mat B_scir, arma::mat pi_scir_cv, arma::mat B_scir_cv, arma::mat pi_sg, arma::mat B_sg){
  
  int N = Y.n_rows;
  
  int M = N*(N-1);
  
  arma::vec p_ij_scir(M, arma::fill::zeros);
  arma::vec p_ij_scir_cv(M, arma::fill::zeros);
  arma::vec p_ij_sg(M, arma::fill::zeros);
  
  int count = 0; 
  for (int i = 0; i < N; ++i){
    
    arma::vec pi_i_scir = pi_scir.row(i).t(); // column vector
    arma::vec pi_i_scir_cv = pi_scir_cv.row(i).t(); // column vector
    arma::vec pi_i_sg = pi_sg.row(i).t(); // column vector
    for (int j = 0; j < N; ++j){
      
      if(i == j) continue;
      
      arma::vec pi_j_scir = pi_scir.row(j).t(); // column vector
      arma::vec pi_j_scir_cv = pi_scir_cv.row(j).t(); // column vector
      arma::vec pi_j_sg = pi_sg.row(j).t(); // column vector
      
      
      int Y_ij = Y(i,j);
      
      p_ij_scir[count]    = P_ij_fun(Y_ij, pi_i_scir, pi_j_scir, B_scir);
      p_ij_scir_cv[count] = P_ij_fun(Y_ij, pi_i_scir_cv, pi_j_scir_cv, B_scir_cv);
      p_ij_sg[count]      = P_ij_fun(Y_ij, pi_i_sg, pi_j_sg, B_sg);
      
      count = count + 1;
      
    }
  }
  double log_P_scir_Y = arma::sum(arma::log(p_ij_scir));
  double log_P_scir_cv_Y = arma::sum(arma::log(p_ij_scir_cv));
  double log_P_sg_Y = arma::sum(arma::log(p_ij_sg));
  
  
  arma::vec p_out(3, arma::fill::zeros);
  p_out[0]= exp(log_P_scir_Y);
  p_out[1]= exp(log_P_scir_cv_Y);
  p_out[2]= exp(log_P_sg_Y);
  
  return( p_out );
}

//[[Rcpp::export]]
arma::vec P_fun(arma::mat Y, Rcpp::List pi_list, Rcpp::List B_list){
  
  int N = Y.n_rows;
  
  int M = N*(N-1);
  
  arma::mat B0 = B_list[0];
  int K = B0.n_rows;
  
  int n_methods = pi_list.size();
  
  arma::mat pij_mat(M, n_methods, arma::fill::zeros);
  
  int count = 0; 
  for (int i = 0; i < N; ++i){
    
    arma::mat p_i(K, n_methods);
    for(int m=0; m < n_methods; ++m ){
      arma::mat pi_m = pi_list[m];
      p_i.col(m)=pi_m.row(i).t(); // column vector
    }
    
    
    for (int j = 0; j < N; ++j){
      
      if(i == j) continue;
      
      int Y_ij = Y(i,j);
      arma::mat p_j(K, n_methods);
      
      for(int m=0; m < n_methods; ++m ){
        
        arma::mat B_m = B_list[m];
        
        arma::vec pi_m_i = p_i.col(m);
        
        arma::mat pj_m = pi_list[m];
        p_j.col(m)=pj_m.row(j).t(); // column vector
        
        arma::vec pi_m_j = p_j.col(m);
        
        pij_mat(count,m) = P_ij_fun(Y_ij, pi_m_i, pi_m_j, B_m);
    
      }
      
      count = count + 1;
      
    }
  }
  
  arma::vec log_P_m(n_methods);
  for(int m=0; m < n_methods; ++m ){
    
    arma::vec pij_m = pij_mat.col(m);
    
    log_P_m[m] = exp(arma::sum(arma::log(pij_m)));
    
  }

  return( log_P_m );
}
