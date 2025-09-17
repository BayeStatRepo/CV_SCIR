library(CVSCIR)
source(paste0(getwd(), "/Experiments/LDA/R/00_LDA_paths.R"))
load( paste0(env_dir, "/minibatch_size.RData") )

train_det<-data.frame(File=list.files(data_batch)) %>% 
  filter( str_detect(File, "train") &  str_detect(File, "det"))

train_l<-data.frame(File=list.files(data_batch)) %>% 
  filter( str_detect(File, "train") &  str_detect(File, "parsed"))
valid<-data.frame(File=list.files(data_batch)) %>% 
  filter( str_detect(File, "valid") )


train_det_all<-data.frame()
train_list<-list()
valid_list<-list()

for( r in 1:nrow(train_det) ){
  train_r <- readRDS( paste0( data_batch, "/",  
                              train_det %>% 
                                filter(row_number()==r ) %>% 
                                as.matrix() %>% as.character() ) ) 
  train_det_all<-rbind( train_det_all, train_r  )
  
  train_l_r<-readRDS( paste0( data_batch, "/",  
                              train_l %>% 
                                filter(row_number()==r ) %>% 
                                as.matrix() %>% as.character() ) )
  
  valid_r<- readRDS( paste0( data_batch, "/",  
                             valid %>% 
                               filter(row_number()==r ) %>% 
                               as.matrix() %>% as.character() ) )  
  
  train_list <- append(train_list, train_l_r[c(1:nrow(train_r))] )
  valid_list <- append(valid_list,  valid_r[c(1:nrow(train_r))] )
  print(r)
}

saveRDS(train_list, file = paste0( data_batch, "/train_list.rds" )  )
saveRDS(valid_list, file = paste0( data_batch, "/valid_list.rds" )  )

n_docs <-nrow( train_det_all )  
iterations <- floor(n_docs/minibatch_size)

n_seeds<-5
idx_arr<-array( NA, dim=c(iterations, minibatch_size, n_seeds) )
for(s in 1:n_seeds){
  set.seed(s)
  batch_list<-list()
  for(t in 1:iterations){
    idx_s_t<-sample(1:n_docs, minibatch_size, replace = F)
    idx_arr[t,,s]<-idx_s_t
    
    batch_list[[t]]<-list(train=train_list[c(idx_s_t)], 
                          valid=valid_list[c(idx_s_t)])
  }
  saveRDS(batch_list, paste0(data_batch, "/batch_list_", s, ".rds"))
}
saveRDS(idx_arr, paste0(data_batch, "/idx_array.rds"))

rm(list=setdiff(ls(), c("n_seeds", "env_dir") ))

load( paste0(env_dir, "/minibatch_size.RData") )
save.image(paste0(env_dir, "/minibatch_size.RData"))
