link_non_link_sets_fun <- function(Y){
  N <- dim(Y)[1]
  V_all <- 1:N
  
  Y_graph <- graph_from_adjacency_matrix(Y)
  
  link_non_link_sets <- lapply(V_all, function(i){ 
    
    # LINK SET
    neig_i_out <- as.numeric( neighbors(Y_graph, i, mode = "out") )
    if( length(neig_i_out)==0 ){
      neig_i_out <- -1
    }
    
    neig_i_in  <- as.numeric( neighbors(Y_graph, i, mode = "in" ) )
    if( length(neig_i_in)==0 ){
      neig_i_in <- -1
    }
    
    link_set_i <- rbind(
      cbind(i=i, j=neig_i_out),
      cbind(i=neig_i_in, j=i)
    )
    
    neig_to_rm <- union(which(link_set_i[,1]<0), which(link_set_i[,2]<0))
    if( length(neig_to_rm) > 0  ){
      link_set_i <- link_set_i[-neig_to_rm, ]
    }
    
    # NON-LINK SET
    non_neig_i_out <- as.numeric(neighbors(complementer(Y_graph, loops = FALSE), i, mode = "out" ))
    if( length(non_neig_i_out)==0 ){
      non_neig_i_out <- -1
    }
    
    non_neig_i_in  <- as.numeric(neighbors(complementer(Y_graph, loops = FALSE), i, mode = "in"  ))
    if( length(non_neig_i_in)==0 ){
      non_neig_i_in <- -1
    }
    
    non_link_set_i <- rbind(
      cbind(i=i, j=non_neig_i_out),
      cbind(i=non_neig_i_in, j=i)
    )
    
    non_neig_to_rm <- union(which(non_link_set_i[,1]<0), which(non_link_set_i[,2]<0))
    if( length(non_neig_to_rm) > 0  ){
      non_link_set_i <- non_link_set_i[-non_neig_to_rm, ]
    }
    
    out <- list( link=link_set_i,  non_link=non_link_set_i )
    return(out)  
  }  )
  
  return(link_non_link_sets)
}
minibatch_sampling <- function(Y, V_in, m, iterations, ell, seed){
  set.seed(seed)
  
  t_start_every_ell <- seq(1, iterations, by=ell)
  t_end_every_ell   <- seq(ell, iterations, by=ell)
  
  link_non_link_sets <- link_non_link_sets_fun(Y)
  
  V_all <- V_in
  
  coin_t <- rep(NA, iterations)
  node_t <- rep(NA, iterations)
  subE_list <- list()
  for(t in 1:iterations){
    
    feas_nodes <- V_all
    
    node_cand_t <- as.numeric(sample(as.character(feas_nodes), 1, replace = F ))
    node_cand_t_sets <- link_non_link_sets[[node_cand_t]]
    
    coin <- rbinom(1,1,0.5)
    coin_t[t] <- coin 
    
    if(coin==1){
      
      subE <- node_cand_t_sets$link 
      
      while( nrow(matrix(subE, ncol=2))==0 ){
        feas_nodes <- setdiff(feas_nodes, node_cand_t)
        
        node_cand_t <- as.numeric(sample(as.character(feas_nodes), 1, replace = F ))
        node_cand_t_sets <- link_non_link_sets[[node_cand_t]]
        subE <- node_cand_t_sets$link 
      }
      
      node_t[t] <- node_cand_t
    }
    else{
      idx <- 1:nrow(matrix(node_cand_t_sets$non_link, ncol = 2))
      
      non_link_sets_sizes <- c( rep(round(length(idx)/m,0), m-1), length(idx)-round((m-1)*length(idx)/m,0) )
      
      subset_pos <- as.numeric( sample( as.character(1:length(non_link_sets_sizes)), 1, replace = F  ) )
      
      subset_idx <- as.numeric( sample( as.character( idx ), non_link_sets_sizes[subset_pos], replace = F  ) )
      
      subE <- node_cand_t_sets$non_link[subset_idx,]
      
      node_t[t] <- node_cand_t
    }
    subE_list[[t]] <- cbind(matrix(subE, ncol=2), rep(coin, nrow(matrix(subE, ncol=2))))
    print(t)
  }
  
  subE_all_list <- list()
  for( r in 1:length(t_start_every_ell) ){
    subE_all_dup <- do.call(rbind, subE_list[ c( t_start_every_ell[r]:t_end_every_ell[r] ) ] )
    
    subE_all <- subE_all_dup[!duplicated(subE_all_dup), ]
    subE_all_list[[r]] <- subE_all
  }
  
  out_list <- list( subE_list,  subE_all_list, coin_t )
}
CIR_B_update_par <- function(n_cores, Y, B, pi, subE_all){
  
  n <- nrow(subE_all)
  
  n_edges_per_core<-c(rep(floor(n/n_cores), n_cores-1), 
                      n-sum(rep(floor(n/n_cores), n_cores-1)))
  
  st_idx<-c(1, cumsum(n_edges_per_core)[-length(cumsum(n_edges_per_core))]+1)
  en_idx<-cumsum(n_edges_per_core)
  
  core_edges_list<-list()
  for( c in 1:n_cores){
    core_edges_list[[c]]<-subE_all[c(st_idx[c]:en_idx[c]),]
  }
  
  myCluster <- makeCluster(n_cores, # number of cores to use
                           type = "PSOCK") # type of cluster
  
  registerDoParallel(myCluster)
  eta_par_list <- foreach(c = 1:n_cores, 
                          .combine = 'c', 
                          .packages = c("CVSCIR")
  ) %dopar% {
    
    subE_c <- matrix(core_edges_list[[c]], ncol = 2) 
    
    eta_par_mat <- CIR_B_update_no_par(B, pi, Y, subE_c-1)
    
  }
  stopCluster(myCluster)
  
  eta0_idx<-(1:length(eta_par_list))[((1:length(eta_par_list))%%2)==1]
  eta1_idx<-setdiff(1:length(eta_par_list), eta0_idx)
  
  c<-1
  eta0 <- eta_par_list[[eta0_idx[c]]]
  eta1 <- eta_par_list[[eta1_idx[c]]]
  for( c in 2:length(eta0_idx)){
    eta0 <- eta0 + eta_par_list[[eta0_idx[c]]]
    eta1 <- eta1 + eta_par_list[[eta1_idx[c]]]
  }
  
  return(list(eta0=eta0, eta1=eta1))
}
CIR_pi_update_par<- function(n_cores, Y, B, pi, sub_V_xi, sub_V_xi_j){
  
  n <- length(sub_V_xi)
  
  n_nodes_per_core<-c(rep(floor(n/n_cores), n_cores-1), 
                      n-sum(rep(floor(n/n_cores), n_cores-1)))
  
  st_idx<-c(1, cumsum(n_nodes_per_core)[-length(cumsum(n_nodes_per_core))]+1)
  en_idx<-cumsum(n_nodes_per_core)
  
  core_nodes_list<-list()
  core_pair_nodes_list<-list()
  for( c in 1:n_cores){
    core_nodes_list[[c]]<-sub_V_xi[c(st_idx[c]:en_idx[c])]
    core_pair_nodes_list[[c]] <- sub_V_xi_j[c(st_idx[c]:en_idx[c])]
  }
  
  myCluster <- makeCluster(n_cores, # number of cores to use
                           type = "PSOCK") # type of cluster
  
  registerDoParallel(myCluster)
  alpha_par_mat_sum <- foreach(c = 1:n_cores, 
                               .combine = '+', 
                               .packages = c("CVSCIR")
  ) %dopar% {
    
    sub_V_xi_c <- core_nodes_list[[c]]
    
    sub_V_xi_j_c <- core_pair_nodes_list[[c]] 
    
    alpha_par_mat <- CIR_pi_update_no_par(B, pi, Y, sub_V_xi_c-1, sub_V_xi_j_c)
    
  }
  stopCluster(myCluster)
  
  return(alpha_par_mat_sum)
}
P_hat_t_Y_fun <- function(t_old, P_hat_Y_old_vec, Y, pi_list, B_list){
  
  P_hat_Y_t <- P_fun(Y, pi_list, B_list)
  
  P_hat_Y_t_mean <- (t_old*P_hat_Y_old_vec + P_hat_Y_t)/(t_old+1)
  
  return(P_hat_Y_t_mean)
}
P_hat_t_Y_list_fun <- function(t_old, P_hat_Y_old_mat, Y_list, pi_list, B_list){
  
  L <- length(Y_list)
  
  P_t_Y_mat <- matrix(t(sapply(1:L, function(l){
    Y <- Y_list[[l]]
    P_hat_Y_old_vec <- P_hat_Y_old_mat[l,]
    P_hat_t_Y_fun(t_old, P_hat_Y_old_vec, Y, pi_list, B_list) })),
    ncol=length(pi_list))
  
  
  return(P_t_Y_mat)
  
}
P_hat_t_Y_list_fun_par  <- function(n_cores, t_old, P_hat_Y_old_mat, Y_list, pi_list, B_list) {
  
  L <- length(Y_list)
  
  n_Y_per_core<-c(rep(floor(L/n_cores), n_cores-1), 
                  L-sum(rep(floor(L/n_cores), n_cores-1)))
  
  
  st_idx<-c(1, cumsum(n_Y_per_core)[-length(cumsum(n_Y_per_core))]+1)
  en_idx<-cumsum(n_Y_per_core)
  
  core_Y_list<-list()
  core_P_hat_Y_old_list <- list()
  for( c in 1:n_cores){
    core_Y_list[[c]]<-Y_list[c(st_idx[c]:en_idx[c])]
    core_P_hat_Y_old_list[[c]] <- matrix( P_hat_Y_old_mat[c(st_idx[c]:en_idx[c]),], 
                                         nrow = length(st_idx[c]:en_idx[c]) )
  }
  
  myCluster <- makeCluster(n_cores, # number of cores to use
                           type = "PSOCK") # type of cluster
  
  registerDoParallel(myCluster)
  
  P_hat_list <- foreach(c = 1:n_cores, 
                          .combine = 'rbind', 
                          .packages = c("CVSCIR")
  ) %dopar% {
    
    Y_list_c <- core_Y_list[[c]] 
    P_hat_Y_old_mat_c <- core_P_hat_Y_old_list[[c]]
    
    P_hat_c <- P_hat_t_Y_list_fun(t_old, P_hat_Y_old_mat_c, Y_list_c, pi_list, B_list)
    
  }
  stopCluster(myCluster)
  
  return(P_hat_list)
  
}