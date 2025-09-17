library(CVSCIR)
library(tidyverse)
source(paste0(getwd(), "/Experiments/LDA/R/00_LDA_paths.R"))

load(paste0(env_dir, "/input.RData"))

train_list_all <- readRDS( paste0( data_batch, "/train_list.rds" ) )
valid_list_all <- readRDS( paste0( data_batch, "/valid_list.rds" ) )

D_partial<-round(length(train_list_all)/minibatch_size, 0)

tau=iterations
hh<-1
kappa<-log(0.1/hh)/log((1+(iterations)/tau)^(-1))# kappa s.t. for t = tau = T (numb of iterations) h_T=0.1
h_vec=hh*(1+(1:iterations)/tau)^(-kappa)
plot(h_vec, type="l" )

ho_train <- readRDS( paste0(data_batch, "/ho_train.rds") )
ho_test  <- readRDS( paste0(data_batch, "/ho_test.rds") )

z_up_t_idx<-seq(1, iterations, by=update_z_every )

time_arr<-array(NA, dim=c(iterations, 2, n_seeds))

for(s in 1:n_seeds){
  
  theta_SCIR<-theta0
  theta_SCIR_CV<-theta0
  theta_SGRLD<-theta0

  ho_count_SCIR<-ho_count_CV<-ho_count_SGRLD<-ho_count
  old_avg_SCIR<-old_avg_CV<-old_avg_SGRLD<-matrix(0, nrow = n_test_set, ncol = W)
  
  perp_valid<-matrix(NA, nrow=iterations, ncol = 3 )
  
  ho_perp_SCIR<-ho_perp_CV<-ho_perp_SGRLD<-c()
  
  st<-Sys.time()
  set.seed(s)
  for( t in 1:iterations ){
    
    st_t<-Sys.time()
    
    h = eps = h_vec[t]

    if( t %in% z_up_t_idx ){
      
      idx_t<-c(sample(1:length(train_list_all), size=D_partial, replace = F ))
      
      train_list <- train_list_all[idx_t]
      
      st_z<-Sys.time()
      sam_cnt_CV_tot <- sample_counts_Rcpp_parallel(n_cores_used, 
                                  theta_SCIR_CV, alpha, probs_init, train_list,
                                  gibbs_iters, gibbs_burnin)
      en_z<-Sys.time()
      
      print( paste0("Iteration ", t, "; Topic assignment latent updated; Time ", 
                    round(en_z-st_z,2), " ", attr(en_z-st_z, "units") ) )
      
      A_tot<-sam_cnt_CV_tot[[1]]
      B_tot<-sam_cnt_CV_tot[[2]]
      
    }
    
    idx_every_t<-sample(idx_t, size = minibatch_size, replace = F)
    
    train_t<-train_list_all[idx_every_t]
    valid_t<-valid_list_all[idx_every_t]
    
    sam_cnt   <-rcpp_sample_cnt(alpha, probs_init, train_t, theta_SCIR,  
                             gibbs_iters, gibbs_burnin)

    sam_cnt_CV<-rcpp_sample_cnt(alpha, probs_init, train_t, theta_SCIR_CV,
                                gibbs_iters, gibbs_burnin)
    
    sam_cnt_SGRLD<-rcpp_sample_cnt(alpha, probs_init, train_t, theta_SGRLD,  
                                   gibbs_iters, gibbs_burnin)
    
    # COMPUTE THE AVERAGE COUNTS
    A_avg<-sam_cnt[[1]]/(gibbs_iters - gibbs_burnin)
    B_avg<-sam_cnt[[2]]/(gibbs_iters - gibbs_burnin)
    
    A_avg_CV<-sam_cnt_CV[[1]]/(gibbs_iters - gibbs_burnin)
    B_avg_CV<-sam_cnt_CV[[2]]/(gibbs_iters - gibbs_burnin)
    
    A_avg_SGRLD<-sam_cnt_SGRLD[[1]]/(gibbs_iters - gibbs_burnin)
    B_avg_SGRLD<-sam_cnt_SGRLD[[2]]/(gibbs_iters - gibbs_burnin)
    
    omega_SGRLD<-theta_SGRLD/rowSums(theta_SGRLD)
    
    for(k in 1:K){
      # SCIR
      par1<-2*(beta + (D/minibatch_size)*B_avg[k,] )
      par2<-2*theta_SCIR[k,]*(exp(-h)/(1-exp(-h)))
      
      X<-rchisq(W, par1,  par2)
      theta_SCIR[k,] <- ((1-exp(-h))/2)*X
      
      # SCIR-CV
      a_hat_CV <- (beta + (D/minibatch_size)*B_avg_CV[k,] )
      a_CV <- beta + (D/D_partial)*B_tot[k, ]

      b_hat_CV <- (a_hat_CV-1)/(a_CV-1)
      par1_CV<-2*a_hat_CV
      par2_CV<-2*theta_SCIR_CV[k,]*( b_hat_CV*(exp(-h*b_hat_CV)/(1-exp(-h*b_hat_CV))) )

      X_CV<-rchisq(W, par1_CV,  par2_CV)
      theta_SCIR_CV[k,] <- ((1-exp(-h*b_hat_CV))/(2*b_hat_CV)) * X_CV

      # SGRLD
      z <- rnorm(W)
      # Update theta according to Equation 11 in paper;
      grad = beta - theta_SGRLD[k,] + (D/minibatch_size)*(B_avg_SGRLD[k,] - sum(B_avg_SGRLD[k,])*omega_SGRLD[k,])
      theta_SGRLD[k,] = abs(theta_SGRLD[k,] + 0.5*eps*grad + (theta_SGRLD[k,])^(0.5) * z*(eps)^(0.5) )
      
    }

    perp_valid[t, 1]<-exp(-log_pred(alpha, A_avg,  theta_SCIR, valid_t) )
    
    perp_valid[t, 2]<-exp(-log_pred(alpha, A_avg_CV,  theta_SCIR_CV, valid_t) )
    
    perp_valid[t, 3]<-exp(-log_pred(alpha, A_avg_SGRLD,  theta_SGRLD, valid_t) )
    
    print( paste0("Iteration ", t, "; Perplexity validation: ", round(perp_valid[t, 2], 2) ) )
    
    if( t %% perpl_every == 0 ){
      
      lp_SCIR<-hold_out_log_pred(alpha, probs_init, 
                                      ho_train, ho_test, 
                                      ho_count_SCIR, old_avg_SCIR, 
                                      theta_SCIR, 
                                      gibbs_iters, gibbs_burnin)
      
      ho_perp_SCIR<-c(ho_perp_SCIR, exp(-lp_SCIR[[1]]) )
      
      ho_count_SCIR<-lp_SCIR[[2]]
      old_avg_SCIR<-lp_SCIR[[3]]
      
      lp_CV<-hold_out_log_pred(alpha, probs_init, 
                                    ho_train, ho_test, 
                                    ho_count_CV, old_avg_CV, 
                                    theta_SCIR_CV, 
                                    gibbs_iters, gibbs_burnin)
      
      ho_perp_CV<-c(ho_perp_CV, exp(-lp_CV[[1]]) )
      
      ho_count_CV<-lp_CV[[2]]
      old_avg_CV<-lp_CV[[3]]
      
      
      lp_SGRLD<-hold_out_log_pred(alpha, probs_init, 
                                       ho_train, ho_test, 
                                       ho_count_SGRLD, old_avg_SGRLD, 
                                       theta_SGRLD, 
                                       gibbs_iters, gibbs_burnin)
      
      ho_perp_SGRLD<-c(ho_perp_SGRLD, exp(-lp_SGRLD[[1]]) )
      
      ho_count_SGRLD<-lp_SGRLD[[2]]
      old_avg_SGRLD<-lp_SGRLD[[3]]
      
      print( paste0("Iteration ", t, "; Perplexity TEST: ", 
                    round(ho_perp_CV[ho_count_CV], 2)) )
    }
   
    en_t<-Sys.time()
    
    time_arr[t,1,s]<-st_t
    time_arr[t,2,s]<-en_t
    
    time_diff_t<-en_t-st_t
    
    print(paste0( "Iteration ", t, "; Time: ", round(time_diff_t, 2), " ", 
                  attr(time_diff_t, "units") ))
    
    if( t %in% store_t_vec ){
      saveRDS(perp_valid, file = paste0(output_path, "/perp_valid_", s,".rds"  ) )
      
      ho_perp<-cbind(ho_perp_SCIR, ho_perp_CV, ho_perp_SGRLD)
      saveRDS(ho_perp,  file = paste0(output_path, "/ho_perp_", s,".rds"  ) )
      rm(ho_perp)
      
      saveRDS(time_arr, file = paste0(output_path, "/time_arr.rds"  ) )
    }
    
  }
  en<-Sys.time()
  time_diff_s<-en-st
  print(paste0( "Seed ", s, "; Time: ", round(time_diff_s, 2), " ", 
                attr(time_diff_s, "units") ))
  
  rm(batch_list)
  gc()
  

}











