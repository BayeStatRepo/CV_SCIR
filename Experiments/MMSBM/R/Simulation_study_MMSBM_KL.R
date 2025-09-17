library(gtools) #rdirichlet
options(scipen = 100) 

sim_path <- paste0(getwd(), "/Experiments/MMSBM/Simulation")
dir.create( sim_path )

# SIMULATE GRAPH
N <- 100
V_all <- 1:N

hold_out_perc <- 0.2
n_hold_out <- round(N*hold_out_perc, 0)

K <- 10
block_sizes <- c( rep(round(N/K,0), K-1), round(N-(K-1)*N/K,0) )
true_label <- rep(1:length(block_sizes), block_sizes)

sum_alpha <- 1000
perc_K2 <- 0.8 
perc_K3_on <- c(0.7, 0.2) 
pi_0 <- matrix(NA, nrow = N, ncol = K)
for(i in 1:length(true_label)){
  alpha_gen <- rep(NA, K)
  if( K==1 ){
    alpha_gen[true_label[i]] <- sum_alpha
  }
  else{
    if( K==2 ){
      alpha1 <- round(sum_alpha*perc_K2,0)
      alpha_gen[true_label[i]] <- alpha1
      alpha_gen[which(is.na(alpha_gen))] <- (sum_alpha-alpha1) 
    }
    if( K>2 ){
      alpha_first <- round(sum_alpha*perc_K3_on,0)
      alpha_on <- rep( (sum_alpha-sum(alpha_first))/(K-length(perc_K3_on)), K-length(perc_K3_on) )
      alpha_vec <- c(alpha_first, alpha_on)
      alpha_gen[true_label[i]] <- max(alpha_vec)
      alpha_gen[ which(is.na(alpha_gen)) ] <- alpha_vec[-which.max(alpha_vec)] 
    }
  }
  pi_0[i,] <- rdirichlet(1, alpha_gen)
}

set.seed(K)
offdiag_B0 <- 0.01
B_0<-matrix(offdiag_B0, nrow = K, ncol = K)
diag(B_0) <- round(rbeta(K, 10, 1), 2)

n_seeds <- 100
methods <- c("SCIR", "SCIR_CV", "SG_MMSBM")
n_methods <- length(methods)

n_sim <- 10000

n_iter <- 1000
burnin <- 1000
thin <- 1
iterations<-burnin+thin*n_iter

# HYPERPARAMETERS
eta_prior_0 <- diag(K) + 0.1
eta_prior_0[row(eta_prior_0)!=col(eta_prior_0)] = 10

eta_prior_1 <- diag(K) + 9
eta_prior_1[row(eta_prior_1)!=col(eta_prior_1)] = 1.1

alpha_prior<-matrix(1.1, N, K)
for( k in 1:K){
  alpha_prior[which(true_label==k), k] <- 10
}

tau<-iterations
h_1<-1
h_T<-0.1
kappa<-log(h_T/h_1)/log((1+(iterations)/tau)^(-1))# kappa s.t. for t = tau = T (numb of iterations) h_T=0.1
h_vec=h_1*(1+(1:iterations)/tau)^(-kappa)

ell<-5
t_start_every_ell <- seq(1, iterations, by=ell)
t_end_every_ell   <- seq(ell, iterations, by=ell)

update_every_ell <- seq(1, iterations, by=ell )

ell_perf<-5
performance_every_ell <- seq(ell, iterations, by=ell_perf )

rm(alpha_first, alpha_gen, alpha_on, alpha_vec, block_sizes, h_1, h_T, hold_out_perc,
   i, k, kappa, offdiag_B0, perc_K2, perc_K3_on, sum_alpha, tau)

save.image(paste0(sim_path,"/Simulation_env.RData"))

library(parallel)
library(doParallel)
library(foreach)
sim_path <- paste0(getwd(), "/Experiments/MMSBM/Simulation")
load(paste0(sim_path,"/Simulation_env.RData"))
n_cores <- parallel::detectCores() - 1

dir.create(paste0(sim_path , "/Results"))

myCluster <- makeCluster(n_cores, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)
res = foreach(s = 1:n_seeds,
        .packages = c("CVSCIR") 
        ) %dopar% {
          
          load(paste0(sim_path,"/Simulation_env.RData"))
          
          set.seed(s)
          
          V_hold_out <- as.numeric(sample(as.character(V_all), n_hold_out, replace = F))
          V_in <- setdiff(V_all, V_hold_out)
          
          # TRUE marginal pmf
          P_0 <- matrix(pi_0[V_hold_out,], ncol=K)%*%B_0%*%t(matrix(pi_0[V_hold_out,], ncol=K))
          diag(P_0) <- 0
          
          sim_Y <- array(NA, dim = c(length(V_hold_out),
                                     length(V_hold_out),
                                     n_sim))
          
          
          for(u in 1:length(V_hold_out)){
            sim_Y[u,u,] <- rep(0, n_sim)
            for(v in setdiff(1:length(V_hold_out), u) ){
              sim_Y[u,v,] <- rbinom(n_sim, 1, P_0[u,v])
            }
          }
          
          gc()
          
          Y_sim_list <- lapply(1:n_sim, function(s) sim_Y[,,s] )
          remove(sim_Y)
          gc()
          
          log_P0 <- sapply(1:n_sim, function(l){
            P_0_fun(Y_sim_list[[l]], matrix(pi_0[V_hold_out,], ncol=K), B_0)
          })
          
          log_P_0_mean <- mean(log_P0)
          
          Y<-matrix(NA, nrow = N, ncol = N)
          p<-matrix(NA, nrow = N, ncol = N)
          
          for (i in 1:N){
            for(j in 1:N){
              if(j==i){
                p[i,j] <- 0
                Y[i,j] <- 0
              }
              else{
                z_ij <- t(rmultinom(1, size = 1, pi_0[i,]))
                z_ji <- rmultinom(1, size = 1, pi_0[j,])
                
                p[i,j] <- z_ij %*% B_0 %*% z_ji
                Y[i,j] <- rbinom(1, 1, p[i,j])
              }
            }
          }
          
          Y_in <- Y[V_in, V_in]
          
          m <- ceiling( (2*choose(ncol(Y_in),2)-sum(Y_in))/sum(Y_in) )
        
          subE_lists <- minibatch_sampling(Y, V_in, m, iterations, ell, seed=s) # repeat for different seeds
          
          # MCMC INITIALIZATION
          theta_0 <- matrix(rgamma(K^2, 1, 1), K, K)
          theta_1 <- matrix(rgamma(K^2, 1, 1), K, K)
          
          phi_0 <- matrix(rgamma(N*K, 1, 1), N, K)
          
          # SCIR
          theta_SCIR <- array(NA, dim = c(K,K,2))
          theta_SCIR[,,1] <- theta_0
          theta_SCIR[,,2] <- theta_1
          
          B_SCIR<-theta_SCIR[,,2]/apply(theta_SCIR, c(1,2), sum)
          
          phi_SCIR <- phi_0
          pi_SCIR  <- phi_SCIR/rowSums(phi_SCIR)
          
          # SCIR-CV
          theta_SCIR_CV <- array(NA, dim = c(K,K,2))
          theta_SCIR_CV[,,1] <- theta_0
          theta_SCIR_CV[,,2] <- theta_1
          
          B_SCIR_CV<-theta_SCIR_CV[,,2]/apply(theta_SCIR_CV, c(1,2), sum)
          
          phi_SCIR_CV <- phi_0
          pi_SCIR_CV  <- phi_SCIR_CV/rowSums(phi_SCIR_CV)
          
          # SG_MMSBM
          theta_SG_MMSBM <- array(NA, dim = c(K,K,2))
          theta_SG_MMSBM[,,1] <- theta_0
          theta_SG_MMSBM[,,2] <- theta_1
          
          B_SG_MMSBM <- theta_SG_MMSBM[,,2]/apply(theta_SG_MMSBM, c(1,2), sum)
          
          phi_SG_MMSBM <- phi_0
          pi_SG_MMSBM  <- phi_SG_MMSBM/rowSums(phi_SG_MMSBM)
          
          perf_count <- 0
          update_count <- 0
          
          P_hat_Y_old_mat <- matrix(0, nrow = n_sim, ncol = length(methods) )
          
          cls_measures<-matrix(NA, nrow = length(performance_every_ell), ncol= n_methods)
          for(t in 1:iterations){
            
            st_t<-Sys.time()
            
            coint_t <- subE_lists[[3]][t]
            scaling_t <- ifelse(coint_t==0, length(V_in)*m, length(V_in))
            
            pos_t <- which( (t >= t_start_every_ell & t <= t_end_every_ell)==T )
            n_link_t <- sum(subE_lists[[3]][t_start_every_ell[pos_t]:t_end_every_ell[pos_t]])
            scaling_all_t <- ifelse(coint_t==0, (length(V_in)*m)-(ell-n_link_t), length(V_in)-n_link_t )
            
            h_par = eps = h_vec[t]
            
            if( t %in% update_every_ell ){
              
              update_count <- update_count + 1
              
              subE_all   <- subE_lists[[2]][[update_count]]
              subE_all_0 <- matrix(subE_all[ which(subE_all[,3]==0), -3 ], ncol = 2)
              subE_all_1 <- matrix(subE_all[ which(subE_all[,3]==1), -3 ], ncol = 2)
              
              subE_all_list <- list(subE_all_0, subE_all_1)
              
              eta_par_mat_all_0 <- CIR_B_update_no_par(theta_SCIR_CV[,,2], theta_SCIR_CV[,,1], 
                                                       pi_SCIR_CV, Y, 
                                                       subE_all_0-1)
              
              eta_par_mat_all_1 <- CIR_B_update_no_par(theta_SCIR_CV[,,2], theta_SCIR_CV[,,1], 
                                                       pi_SCIR_CV, Y, 
                                                       subE_all_1-1)
              
              eta_par_all_list <- list(eta_par_mat_all_0, eta_par_mat_all_1)
              
            }
            subE <- matrix(subE_lists[[1]][[t]][,-3], ncol=2)
            
            # UPDATE GLOBAL PARAMETER B
            eta_par_mat_SCIR   <- CIR_B_update_no_par(theta_SCIR[,,2], theta_SCIR[,,1], 
                                                      pi_SCIR, Y, subE-1)
            eta_par_mat_0_SCIR <- eta_par_mat_SCIR[[1]]
            eta_par_mat_1_SCIR <- eta_par_mat_SCIR[[2]]
            
            eta_par_mat_SCIR_CV  <- CIR_B_update_no_par(theta_SCIR_CV[,,2], theta_SCIR_CV[,,1], 
                                                        pi_SCIR_CV, Y, subE-1)
            eta_par_mat_0_SCIR_CV <- eta_par_mat_SCIR_CV[[1]]
            eta_par_mat_1_SCIR_CV <- eta_par_mat_SCIR_CV[[2]]
            
            eta_par_mat_all   <- eta_par_all_list[[coint_t+1]]
            eta_par_mat_0_all <- eta_par_mat_all[[1]]
            eta_par_mat_1_all <- eta_par_mat_all[[2]]
            
            eta_par_mat_SG_MMSBM   <- SG_B_update_no_par(theta_SG_MMSBM[,,1], theta_SG_MMSBM[,,2], 
                                                         B_SG_MMSBM, pi_SG_MMSBM, Y, subE-1)
            eta_par_mat_0_SG_MMSBM <- eta_par_mat_SG_MMSBM[[1]]
            eta_par_mat_1_SG_MMSBM <- eta_par_mat_SG_MMSBM[[2]]
            
            for(k in 1:K){
              for( h in 1:K){
                
                eta0 = eta_prior_0[k,h]
                eta1 = eta_prior_1[k,h]
                
                # SCIR
                # ell = 0
                eta_kh_0_SCIR <- eta0 + scaling_t *  eta_par_mat_0_SCIR[k,h]
                
                par_0_1_SCIR <- 2*eta_kh_0_SCIR
                par_0_2_SCIR <- 2*theta_SCIR[k,h,1] * ( exp(-h_par)/(1-exp(-h_par)) )
                
                W_kh_0_SCIR <- rchisq(1, par_0_1_SCIR,  par_0_2_SCIR) * (1-exp(-h_par))/(2)
                
                theta_SCIR[k,h,1] <- W_kh_0_SCIR 
                
                # ell = 1
                eta_kh_1_SCIR <- eta1 + scaling_t *  eta_par_mat_1_SCIR[k,h]
                
                par_1_1_SCIR <- 2*eta_kh_1_SCIR
                par_1_2_SCIR <- 2*theta_SCIR[k,h,2] * ( exp(-h_par)/(1-exp(-h_par)) )
                
                W_kh_1_SCIR <- rchisq(1, par_1_1_SCIR,  par_1_2_SCIR) * (1-exp(-h_par))/(2)
                
                theta_SCIR[k,h,2] <- W_kh_1_SCIR 
                
                B_SCIR[k,h] <- W_kh_1_SCIR/(W_kh_0_SCIR + W_kh_1_SCIR)
                
                # SCIR-CV
                # ell = 0
                eta_kh_0_SCIR_CV <- eta0 + scaling_t     *  eta_par_mat_0_SCIR_CV[k,h]
                eta_kh_0_all     <- eta0 + scaling_all_t *  eta_par_mat_0_all[k,h]
                b_eta_0 <- (eta_kh_0_SCIR_CV - 1)/(eta_kh_0_all - 1)
                
                par_0_1_SCIR_CV <- 2*eta_kh_0_SCIR_CV
                par_0_2_SCIR_CV <- 2*theta_SCIR_CV[k,h,1] * b_eta_0 * ( exp(-h_par*b_eta_0)/(1-exp(-h_par*b_eta_0)) )
                
                W_kh_0_SCIR_CV <- rchisq(1, par_0_1_SCIR_CV,  par_0_2_SCIR_CV) * (1-exp(-h_par*b_eta_0))/(2*b_eta_0)
                
                theta_SCIR_CV[k,h,1] <- W_kh_0_SCIR_CV
                
                # ell = 1
                eta_kh_1_SCIR_CV <- eta1 + scaling_t     * eta_par_mat_1_SCIR_CV[k,h]
                eta_kh_1_all     <- eta1 + scaling_all_t * eta_par_mat_1_all[k,h]
                b_eta_1 <- (eta_kh_1_SCIR_CV - 1)/(eta_kh_1_all - 1)
                
                par_1_1_SCIR_CV <- 2*eta_kh_1_SCIR_CV
                par_1_2_SCIR_CV <- 2*theta_SCIR_CV[k,h,2] * b_eta_1 * ( exp(-h_par*b_eta_1)/(1-exp(-h_par*b_eta_1)) )
                
                W_kh_1_SCIR_CV <- rchisq(1, par_1_1_SCIR_CV,  par_1_2_SCIR_CV) * (1-exp(-h_par*b_eta_1))/(2*b_eta_1)
                
                theta_SCIR_CV[k,h,2] <- W_kh_1_SCIR_CV
                
                B_SCIR_CV[k,h] <- W_kh_1_SCIR_CV/(W_kh_0_SCIR_CV + W_kh_1_SCIR_CV)
                
                # SG_MMSBM
                # ell = 0
                z <- rnorm(2)
                
                theta_SG_MMSBM[k,h,1] <- abs(theta_SG_MMSBM[k,h,1] +
                                               eps/2 * (eta0 - theta_SG_MMSBM[k,h,1] +
                                                          scaling_t * eta_par_mat_0_SG_MMSBM[k,h]) +
                                               (eps)^(0.5) * theta_SG_MMSBM[k,h,1]^(0.5) * z[1]
                )
                
                # ell = 1
                theta_SG_MMSBM[k,h,2] <- abs(theta_SG_MMSBM[k,h,2] +
                                               eps/2 * (eta1 - theta_SG_MMSBM[k,h,2] +
                                                          scaling_t * eta_par_mat_1_SG_MMSBM[k,h]) +
                                               (eps)^(0.5) * theta_SG_MMSBM[k,h,2]^(0.5) * z[2]
                )
                
                B_SG_MMSBM[k,h] <- theta_SG_MMSBM[k,h,2]/(theta_SG_MMSBM[k,h,1] + theta_SG_MMSBM[k,h,2])
                
              }
            }
            B_list <- list(B_SCIR, B_SCIR_CV, B_SG_MMSBM)
            
            # UPDATE LOCAL PARAMETER pi
            subV <- sort(unique(as.vector(subE)))
            
            subV_all <- sort(unique(as.vector(subE_all_list[[coint_t+1]])))
            
            sub_V_i <-lapply(1:N,
                             function(i){ setdiff( subV, i ) - 1 } )
            
            sub_V_i_all <-lapply(1:N,
                                 function(i){ setdiff( subV_all, i ) - 1 } )
            
            
            if( t %in% update_every_ell ){
              
              alpha_par_mat_all <- CIR_pi_update_no_par(B_SCIR_CV, phi_SCIR_CV, Y, V_all-1, sub_V_i_all)
              
            }
            
            # SCIR
            alpha_par_mat_SCIR <- CIR_pi_update_no_par(B_SCIR, phi_SCIR, Y, V_all-1, sub_V_i)
            
            # SCIR-CV
            alpha_par_mat_SCIR_CV <- CIR_pi_update_no_par(B_SCIR_CV, phi_SCIR_CV, Y, V_all-1, sub_V_i)
            
            # SG_MMSBM
            alpha_par_mat_SG_MMSBM <- SG_pi_update_no_par(phi_SG_MMSBM, B_SG_MMSBM, pi_SG_MMSBM, Y, V_all-1, sub_V_i)
            
            for( r in 1:length(V_all) ){
              
              i = V_all[r]
              
              n_V_i     <- length(sub_V_i[[r]])
              n_V_i_all <- length(sub_V_i_all[[r]])
              
              alpha <- alpha_prior[i,] 
              
              # SCIR
              alpha_i_SCIR <- alpha + ((N-1)/n_V_i) * alpha_par_mat_SCIR[i, ]
              
              par1_alpha_SCIR <- 2*alpha_i_SCIR
              par2_alpha_SCIR <- 2*phi_SCIR[i,] * exp(-h_par)/(1-exp(-h_par))
              
              phi_SCIR[i,] <- (1-exp(-h_par))/2 * rchisq(K, par1_alpha_SCIR,  par2_alpha_SCIR)
              pi_SCIR[i, ] <- phi_SCIR[i,]/sum(phi_SCIR[i,])
              
              # SCIR-CV
              alpha_i_SCIR_CV <- alpha + ((N-1)/n_V_i) * alpha_par_mat_SCIR_CV[i, ]
              alpha_i_all <- alpha + ((N-1)/n_V_i_all) * alpha_par_mat_all[i, ]
              
              b_alpha_i <- (alpha_i_SCIR_CV - 1)/(alpha_i_all - 1)
              
              par1_alpha_SCIR_CV <- 2*alpha_i_SCIR_CV
              par2_alpha_SCIR_CV <- 2*phi_SCIR_CV[i,] * b_alpha_i * ( exp(-h_par*b_alpha_i)/(1-exp(-h_par*b_alpha_i)) )
              
              phi_SCIR_CV[i,] <- (1-exp(-h_par*b_alpha_i))/(2*b_alpha_i )* rchisq(K, par1_alpha_SCIR_CV,  par2_alpha_SCIR_CV)
              pi_SCIR_CV[i, ] <- phi_SCIR_CV[i,]/sum(phi_SCIR_CV[i,])
              
              # SG_MMSBM
              z_i <- rnorm(K)
              
              phi_SG_MMSBM[i,] <- abs(phi_SG_MMSBM[i,] + eps/2 * (alpha - phi_SG_MMSBM[i,] +
                                                                    ((N-1)/n_V_i) * alpha_par_mat_SG_MMSBM[i,]) +
                                        eps^(0.5)*phi_SG_MMSBM[i,]^(0.5) * z_i)
              
              pi_SG_MMSBM[i,]  <- phi_SG_MMSBM[i,]/sum(phi_SG_MMSBM[i,])
              
            }
            
            pi_list <- list( matrix(    pi_SCIR[V_hold_out,], ncol=K),
                             matrix( pi_SCIR_CV[V_hold_out,], ncol=K),
                             matrix(pi_SG_MMSBM[V_hold_out,], ncol=K)
            )
            
            if( t %in% performance_every_ell ){
              
              P_hat_Y_old_mat   <- P_hat_t_Y_list_fun(perf_count, P_hat_Y_old_mat, 
                                                      Y_sim_list, pi_list, B_list)
              
              
              perf_count <- perf_count + 1
              
              KL_vec <- log_P_0_mean-colMeans(log( P_hat_Y_old_mat  ))
              
              cls_measures[perf_count,1] <- KL_vec[1]
              cls_measures[perf_count,2] <- KL_vec[2]
              cls_measures[perf_count,3] <- KL_vec[3]
              
            }
          }
          saveRDS(cls_measures, paste0(sim_path, "/Results/cls_measures_", s, ".rds"))

}
stopCluster(myCluster)


# PLOT RESULTS
library(tidyverse)
library(ggh4x) 
library(reporter)

sim_path <- paste0(getwd(), "/Experiments/MMSBM/Simulation")
load(paste0(sim_path,"/Simulation_env.RData"))

n_inter_eval <- 1:length(performance_every_ell)

burnin_perf <- 1:(burnin/ell_perf)
file_names <- list.files(paste0(sim_path, "/Results"))

cls_measures_df <- data.frame()
for(f in 1:length(file_names)){
  cls_s = readRDS( paste0(sim_path, "/Results/", file_names[f]) )
  
  cls_measures_df <- rbind(cls_measures_df,
                           data.frame(value_SCIR = mean(cls_s[-burnin_perf ,1]), 
                                      value_CV   = mean(cls_s[-burnin_perf ,2]), 
                                      value_SG   = mean(cls_s[-burnin_perf ,3]), Measure="KL") )
  
}

cls_measures_df_longer = cls_measures_df %>% pivot_longer(c(cls_measures_df %>% select(contains("value")) %>% names()),
                                                          names_to = c(".value", "Method"),
                                                          names_sep = "\\_")



pp<-ggplot(cls_measures_df_longer,  
           xlab=",",
           aes(x=Method, y=value, fill=Method)) +
  geom_violin()  + 
  theme_bw()  + 
  scale_fill_manual(values = c( rgb(0.9290, 0.6940, 0.1250), 4, rgb(0, 0, 0, 0.5) )
  ) +
  facet_grid2(~Measure,
              labeller = labeller(Measure = c("KL" = "KL divergence")),
              scales = "free_y" , independent="y" ) +
  theme( axis.text.x=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         legend.title=element_blank()) +
  theme( strip.text.x = element_text(size = 20),
         strip.text.y.right  = element_text(size = 20),
         axis.text=element_text(size=20)
  ) +
  theme(aspect.ratio = 1) 

pp

ggsave(paste0(sim_path,"/KL.pdf"), 
       plot = pp, width = 30, height = 15, 
       units = "cm", dpi = 300)
