library(CVSCIR)
source(paste0(getwd(), "/Experiments/LDA/R/00_LDA_paths.R"))

minibatch_size<-50
n_test_set<-1000

# ell
update_z_every<-5

batch_output_name<-"raw_docs"
save.image(paste0(env_dir, "/minibatch_size.RData"))


burnin<-0
n_iter<-1000
thin<-1
iterations<-burnin+thin*n_iter

n_docs <- minibatch_size * iterations + n_test_set

seed=2

st<-Sys.time()
online_wiki(data_in, data_batch, batch_output_name, n_docs, seed=seed)
en<-Sys.time()
en-st

