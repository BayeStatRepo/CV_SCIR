lda_dir     <- paste0(getwd(), "/Experiments/LDA")
env_dir     <- paste0(lda_dir, "/Environments")
output_path <- paste0(lda_dir, "/Output")

data_dir   <- paste0(lda_dir, "/Data")
data_in    <- paste0(data_dir, "/wikidump")
data_batch <- paste0(data_dir, "/wiki_batch")

if( !( dir.exists(env_dir) ) ){
  dir.create(env_dir)
}
if( !( dir.exists(output_path) ) ){
  dir.create(output_path)
}
if( !( dir.exists(data_batch) ) ){
  dir.create(data_batch)
}
