data_path <- paste0(getwd(), "/Experiments/MMSBM/Data")

email.Eu.core <- read.table(paste0(data_path, "/email-Eu-core.txt"), 
                            quote="\"", comment.char="")

email.Eu.core.department.labels <- read.table(paste0(data_path, "/email-Eu-core-department-labels.txt"), 
                                              quote="\"", comment.char="")

true_label <- as.numeric( email.Eu.core.department.labels$V2 + 1 ) 

library(tidyverse)
library(gtools)
library(igraph)

K <- nrow( email.Eu.core.department.labels %>% group_by(V2) %>% 
             summarise(n=n()) )

N <- nrow(email.Eu.core.department.labels)
M <- 2*choose(N,2)

V_all <- 1:N

email.Eu.core_ori<-email.Eu.core
email.Eu.core<-data.frame(email.Eu.core) %>% filter(V1!=V2)

Y<-matrix(0, N, N)
for( r in 1:nrow(email.Eu.core) ){
  Y[email.Eu.core$V1[r]+1, email.Eu.core$V2[r]+1]<-1 
  if( r %in% seq(1000, floor(nrow(email.Eu.core)/1000)*1000, 1000) ){
    print(paste0("row: ", r, " out of ", nrow(email.Eu.core), "." ))
  }
}

m <- ceiling((M-sum(Y))/sum(Y))

methods <- c("SCIR", "SCIR_CV", "SG_MMSBM")
n_methods <- length(methods)
cls_measures <- c("aRI", "entropy", "NMI", "purity" )
n_cls_measures <- length(cls_measures)

n_iter <- 10000
burnin <- 0
thin <- 1
iterations<-burnin+thin*n_iter

# HYPERPARAMETERS
eta_prior_0 <- diag(K) + 0.1
eta_prior_0[row(eta_prior_0)!=col(eta_prior_0)] = 15

eta_prior_1 <- diag(K) + 14
eta_prior_1[row(eta_prior_1)!=col(eta_prior_1)] = 1.1

alpha_prior<-matrix(1.1, N, K)
for( k in 1:K){
  alpha_prior[which(true_label==k), k] <- 15
}

eps <- 0.01

tau<-iterations
h_1<-1
kappa=1.3
h_vec=h_1*(1+(1:iterations)/tau)^(-kappa)

ell<-5
t_start_every_ell <- seq(1, iterations, by=ell)
t_end_every_ell   <- seq(ell, iterations, by=ell)

update_every_ell <- seq(1, iterations, by=ell )

ell_perf<-5
performance_every_ell <- seq(ell, iterations, by=ell_perf )

rm(email.Eu.core, email.Eu.core_ori, email.Eu.core.department.labels, k)

email_eu_path <- paste0(getwd(), "/Experiments/MMSBM/Email_EU")
dir.create(email_eu_path)
save.image(paste0(email_eu_path, "/email_eu_pre_env.RData"))

library(CVSCIR)
library(fossil) # adj.rand.index
library(NMF) # purity and entropy of a clustering
library(aricode) # NMI
email_eu_path <- paste0(getwd(), "/Experiments/MMSBM/Email_EU")
load(paste0(email_eu_path, "/email_eu_pre_env.RData"))

subE_lists <- minibatch_sampling(Y, V_all, m, iterations, ell, seed=1) 

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

pi_cum_SCIR     <- matrix(0, nrow = N, ncol = K)
pi_cum_SCIR_CV  <- matrix(0, nrow = N, ncol = K)
pi_cum_SG_MMSBM <- matrix(0, nrow = N, ncol = K)

cls_measures<-array(NA, dim = c(length(performance_every_ell), 
                                n_cls_measures, 
                                n_methods))
for(t in 1:iterations){
  
  st_t<-Sys.time()
  
  h_par <- h_vec[t]
 
  coint_t <- subE_lists[[3]][t]
  scaling_t <- ifelse(coint_t==0, N*m, N)
  
  pos_t <- which( (t >= t_start_every_ell & t <= t_end_every_ell)==T )
  n_link_t <- sum(subE_lists[[3]][t_start_every_ell[pos_t]:t_end_every_ell[pos_t]])
  scaling_all_t <- ifelse(coint_t==0, (N*m)-(ell-n_link_t), N-n_link_t )
  
  if( t %in% update_every_ell ){
    update_count <- update_count + 1
    
    subE_all   <- subE_lists[[2]][[update_count]]
    subE_all_0 <- matrix(subE_all[ which(subE_all[,3]==0), -3 ], ncol = 2)
    subE_all_1 <- matrix(subE_all[ which(subE_all[,3]==1), -3 ], ncol = 2)
    
    subE_all_list <- list(subE_all_0, subE_all_1)
    
    st_ell<-Sys.time()
    eta_par_mat_all_0 <- CIR_B_update_no_par(theta_SCIR_CV[,,2], theta_SCIR_CV[,,1], 
                                             pi_SCIR_CV, Y, 
                                             subE_all_0-1)
    
    eta_par_mat_all_1 <- CIR_B_update_no_par(theta_SCIR_CV[,,2], theta_SCIR_CV[,,1], 
                                             pi_SCIR_CV, Y, 
                                             subE_all_1-1)
    
    eta_par_all_list <- list(eta_par_mat_all_0, eta_par_mat_all_1)
    en_ell<-Sys.time()
    
    print( paste0("Iteration ", t, "; global parameter all updated; Time ",
                  round(en_ell-st_ell,2), " ", attr(en_ell-st_ell, "units") ) )
    
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
  
  # UPDATE LOCAL PARAMETER pi
  subV <- sort(unique(as.vector(subE)))
  
  subV_all <- sort(unique(as.vector(subE_all_list[[coint_t+1]])))
  
  sub_V_i <-lapply(1:length(V_all),
                   function(i){ setdiff( subV, i ) - 1 } )
  
  sub_V_i_all <-lapply(1:length(V_all),
                       function(i){ setdiff( subV_all, i ) - 1 } )
  
  
  if( t %in% update_every_ell ){
    
    st_ell<-Sys.time()
    
    alpha_par_mat_all <- CIR_pi_update_no_par(B_SCIR_CV, phi_SCIR_CV, Y, V_all-1, sub_V_i_all)
  
    en_ell<-Sys.time()
    
    print( paste0("Iteration ", t, "; local parameter all updated; Time ",
                  round(en_ell-st_ell,2), " ", attr(en_ell-st_ell, "units") ) )
    
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
  
  
  if( t %in% performance_every_ell ){
    
    perf_count <- perf_count + 1
    
    # SCIR
    pi_cum_SCIR <- pi_cum_SCIR + pi_SCIR
    pi_mean_SCIR <-  pi_cum_SCIR/perf_count
    
    z_hat_SCIR <- apply(pi_mean_SCIR, 1, which.max)
    
    cls_measures[perf_count,1,1] <- adj.rand.index(true_label, z_hat_SCIR)
    cls_measures[perf_count,2,1] <- NMF::entropy(as.factor(z_hat_SCIR), true_label )
    cls_measures[perf_count,3,1] <- NMI(true_label, z_hat_SCIR)
    cls_measures[perf_count,4,1] <- purity(as.factor(z_hat_SCIR), true_label)
    
    # SCIR-CV
    pi_cum_SCIR_CV  <- pi_cum_SCIR_CV + pi_SCIR_CV
    pi_mean_SCIR_CV <-  pi_cum_SCIR_CV/perf_count
    
    z_hat_SCIR_CV <- apply(pi_mean_SCIR_CV, 1, which.max)
    
    cls_measures[perf_count,1,2] <- adj.rand.index(true_label, z_hat_SCIR_CV)
    cls_measures[perf_count,2,2] <- NMF::entropy(as.factor(z_hat_SCIR_CV), true_label)
    cls_measures[perf_count,3,2] <- NMI(true_label, z_hat_SCIR_CV)
    cls_measures[perf_count,4,2] <- purity(as.factor(z_hat_SCIR_CV), true_label)
    
    # SG_MMSBM
    pi_cum_SG_MMSBM  <- pi_cum_SG_MMSBM + pi_SG_MMSBM
    pi_mean_SG_MMSBM <- pi_cum_SG_MMSBM/perf_count
    
    z_hat_SG_MMSBM <- apply(pi_mean_SG_MMSBM, 1, which.max)
    
    cls_measures[perf_count,1,3] <- adj.rand.index(true_label, z_hat_SG_MMSBM)
    cls_measures[perf_count,2,3] <- NMF::entropy(as.factor(z_hat_SG_MMSBM), true_label)
    cls_measures[perf_count,3,3] <- NMI(true_label, z_hat_SG_MMSBM)
    cls_measures[perf_count,4,3] <- purity(as.factor(z_hat_SG_MMSBM), true_label)
    
  }
  
  en_t<-Sys.time()
  time_diff_t<-en_t-st_t
  print(paste0("Iteration ", t, "; Time: ", round(time_diff_t, 2), " ", 
                attr(time_diff_t, "units") ))
}
saveRDS(cls_measures, paste0(email_eu_path, "/email_eu_cls_measures.rds") )

# PLOT RESULTS
email_eu_path <- paste0(getwd(), "/Experiments/MMSBM/Email_EU")
load(paste0(email_eu_path, "/email_eu_pre_env.RData"))

cls_measures = readRDS(paste0(email_eu_path, "/email_eu_cls_measures.rds") )
  
library(tidyverse)
cls_measures_df <- data.frame()
cls_measures_df <- rbind(cls_measures_df, 
                         data.frame(Time_index = 1:dim(cls_measures)[1], 
                                    value_SCIR=cls_measures[,1,1], 
                                    value_CV=cls_measures[,1,2], 
                                    value_SG=cls_measures[,1,3], 
                                    Measure="ARI"),
                         data.frame(Time_index = 1:dim(cls_measures)[1],
                                    value_SCIR=cls_measures[,2,1], 
                                    value_CV=cls_measures[,2,2], 
                                    value_SG=cls_measures[,2,3], 
                                    Measure="Entropy"),
                         data.frame(Time_index = 1:dim(cls_measures)[1],
                                    value_SCIR=cls_measures[,3,1], 
                                    value_CV=cls_measures[,3,2], 
                                    value_SG=cls_measures[,3,3], 
                                    Measure="NMI"),
                         data.frame(Time_index = 1:dim(cls_measures)[1],
                                    value_SCIR=cls_measures[,4,1], 
                                    value_CV=cls_measures[,4,2], 
                                    value_SG=cls_measures[,4,3], 
                                    Measure="Purity")
)

cls_measures_df_longer <- cls_measures_df %>% pivot_longer(cols  = c(cls_measures_df %>% select( contains("value") ) %>% names()), 
                                                           names_to = c(".value", "Method"),
                                                           names_sep = "\\_")


library(ggh4x) 
df_running <- cls_measures_df_longer %>%
  group_by(Method, Measure) %>%
  arrange(Time_index) %>%
  mutate(running_value = cumsum(value) / seq_along(value)) %>%
  ungroup()

df_summary <- df_running %>%
  group_by(Method, Measure, Time_index) %>%
  summarise(
    mean = mean(running_value),
    .groups = "drop" )

# Plot with ggplot
(pp <- ggplot(df_summary, 
              aes(x = Time_index, y = mean, color = Method, fill = Method)) +
    geom_line(size = 0.2) +
    scale_color_manual(name="Method", values = c( rgb(0.9290, 0.6940, 0.1250), 4, rgb(0, 0, 0, 0.5) )            
    ) +
    scale_fill_manual(name="Method", values = c( rgb(0.9290, 0.6940, 0.1250), 4, rgb(0, 0, 0, 0.5) )
    ) +
    scale_x_continuous(breaks = seq(0, nrow(df_summary)/(n_methods*n_cls_measures), length.out=5),
                       labels = as.character( seq(0, nrow(df_summary)/(n_methods*n_cls_measures) * ell_perf, length.out=5) ) )  + # n_methods * n_cls_measures * ell_perf
    facet_grid2(~ Measure, 
      labeller = labeller(Measure = c(
        # "ARI" = "Adjusted Rand Index",
        "Entropy" = "Entropy",
        # "NMI" = "Normalized Mutual Information",
        "Purity" = "Purity"))) +
    theme_bw() +
    theme( 
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()  ) +
    theme( 
      strip.text.x = element_text(size = 15),
      strip.text.y.right  = element_text(size = 15),
      axis.text=element_text(size=10) ) + 
    theme(aspect.ratio = 1) 
)

ggsave(paste0(email_eu_path, "/plot_cls_measure_eu_email.pdf"),
       plot = pp, width = 30, height = 15, units = "cm", dpi = 300)


# # PLOT GRAPHS
library(igraph)
library(ggraph)
library(RColorBrewer)
email_eu_path <- paste0(getwd(), "/Experiments/MMSBM/Email_EU")
load(paste0(email_eu_path, "/email_eu_pre_env.RData"))

Y_graph <- graph_from_adjacency_matrix(Y)

V(Y_graph)$community <-true_label

# Plot the largest two communities only.
b_k <- as.numeric(attr(sort(table(true_label))[(K-1):K], "dimnames")$true_label)
v_keep = which(true_label %in% b_k )
sub_Y_graph <- induced_subgraph(Y_graph, vids = v_keep)

V(sub_Y_graph)$color <- ifelse(V(sub_Y_graph)$community == b_k[1], rgb(0.9290, 0.6940, 0.1250), 4)
V(sub_Y_graph)$name <- paste0(1:vcount(sub_Y_graph))

layout <- layout_with_fr(sub_Y_graph)

pdf(paste0(email_eu_path, "/plot_largest_communities_eu_email.pdf"), 
    width = 12, height = 10)

ggraph(sub_Y_graph, layout = "fr") +
  geom_edge_link(alpha =0.5, colour = "black",
                 arrow = arrow(length = unit(5, "mm"), 
                               type = "closed", ends = "last",
                               angle=20),  # add arrows
                 end_cap = circle(2, 'mm')  # spacing at target node
  ) +
  geom_node_point(aes(color = color), size = 5) +
  geom_node_text(aes(label = name), repel = F, size = 2, color="blue") + 
  scale_color_manual(values = c(rgb(0.9290, 0.6940, 0.1250), 4  ),
                     name = "Community") +
  theme_void() +
  theme(legend.position = "none")
dev.off()