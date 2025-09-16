lda_dir<-paste0(getwd(), "/Experiments/LDA")
env_dir<-paste0(lda_dir, "/Environments")
data_dir<-paste0(lda_dir, "/Data")
# data_in<-paste0(data_dir, "/wikidump")
data_in<-"/home/fb21093/Desktop/LDA/data/wikidump"
data_out<-paste0(data_dir, "/wikidump_df")
data_batch<-paste0(data_dir, "/wiki_batch")
output_path<-paste0(lda_dir, "/output")
if( !( dir.exists(data_batch) ) ){
  dir.create(data_batch)
}
if( !( dir.exists(output_path) ) ){
  dir.create(output_path)
}