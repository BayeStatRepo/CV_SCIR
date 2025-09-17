library(CVSCIR)
library(tidyverse)
source(paste0(getwd(), "/Experiments/LDA/R/00_LDA_paths.R"))

load(paste0(env_dir, "/minibatch_size.RData"))
# D=count_all_docs(data_in)
D = 3467794 # result of the count

batch_list<-readRDS(paste0(data_batch, "/batch_list_1.rds"))

test_list<-readRDS( paste0( data_batch, "/",
                 data.frame(File=list.files(data_batch)) %>% 
                   filter( str_detect(File, "test") &  str_detect(File, "parsed") ) %>% 
                   as.matrix() %>% as.character() ))

W <- length(batch_list[[1]]$train[[1]])
iterations<-length(batch_list)
n_test_set<-length(test_list)

store_every_perc <- 0.2
store_t_vec<-unique( c( seq( (store_every_perc)*100, iterations, 
                             by = (store_every_perc)*100 ),
                        iterations ) )

# HYPERPARAMETERS
K<-50 # latent topics
alpha=1.1
beta<-rep(1.1, W )

# STARTING VALUES
theta0<-matrix(NA, nrow = K, ncol = W)
phi0<-matrix(NA, nrow = K, ncol = W)
for(k in 1:K){
  theta0[k,]<-sapply(1:W, function(w){rgamma(1, shape = beta[w], rate = 1)})
  phi0[k,]<-theta0[k,]/sum(theta0[k,])
}

probs_init<-diff( seq(0, 1, length.out = K+1 ) )


n_iter<-100
gibbs_burnin<-100
thin<-1
gibbs_iters<-gibbs_burnin+thin*n_iter

avg_probs<-matrix(0, nrow = n_test_set, ncol = W)
ho_count<-0

n_cores <- parallel::detectCores()

n_cores_used <- n_cores - 1

perpl_every <- 5

rm(list=setdiff(ls(), c( "env_dir", "W", "K", "alpha", "beta", "probs_init", "theta0",
                         "gibbs_iters", "gibbs_burnin", "n_test_set",
                         "iterations", "minibatch_size", "minibatch_size_N", 
                         "update_z_every", "perpl_every", "store_t_vec", 
                         "avg_probs", "ho_count", "n_cores_used", "n_seeds",
                         "hh", "D", "tau", "kappa"
                          )))

save.image(paste0(env_dir, "/input.RData"))

