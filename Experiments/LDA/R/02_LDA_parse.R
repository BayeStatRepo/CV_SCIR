library(CVSCIR)
library(tidyverse)
library(parallel)
source(paste0(getwd(), "/Experiments/LDA/R/00_LDA_paths.R"))

load( paste0(env_dir, "/minibatch_size.RData") )

vocab_raw <- read.table(paste0(data_dir, "/wiki.vocab"), quote="\"", comment.char="")
vocab_sort<-vocab_raw %>% rename(words=1) %>% arrange(words)
vocab<-vocab_sort %>% filter(words!="None")
rm(vocab_sort, vocab_raw)

# convert the vocabulary words into lower case and remove non-alphabetic characters
# in this case not necessary to run as the vocalubary is already clean but, in general, 
# it is worth to run as a check

vocab<-data.frame(words=sapply(1:nrow(vocab), function(w){
  text <- vocab %>% filter(row_number()==w) %>% as.matrix() %>% as.character()
  text_lower_case<-stringr::str_to_lower(text)
  # Remove non-alphabetic characters except space
  text_lower_case <- gsub('[^a-z ]', ' ', text_lower_case)
  print(paste0("word: ", w ) )
  return(text_lower_case)
}  ) )
saveRDS( vocab, paste0(data_dir, "/vocab.rds") )

vocab<-readRDS(paste0(data_dir, "/vocab.rds"))

test_perc=0.1

n_cores = parallel::detectCores() - 1

seed=2
parse_doc_parallel_train_test(data_in, data_batch, batch_output_name, test_perc, n_test_set, vocab, n_cores, seed)

parse_hold_out(data_batch, vocab, test_perc)
