# Function returning the position index (in C++) of positive (>0) words frequencies in a given document
words_idx<-function(doc){
  idx<-as.numeric(which(doc>0))-1
}
# Function returning matrices initialized to 0 for the sample counts C++ function
mat_init <- function(K, W, D) {
  Adk<-matrix(0, nrow = D, ncol = K)
  Adk_avg<-matrix(0, nrow = D, ncol = K)
  
  Bkw<-matrix(0, nrow = K, ncol = W)
  Bkw_avg<-matrix(0, nrow = K, ncol = W)
  
  return(list(Adk, Adk_avg, Bkw, Bkw_avg))
}
# Function returning the inizialization of the vector cumprobs
cumprob_ini<-function(K){
  cumprobs = seq(0, 1, length.out=(K+1) )[0:K]
  return(cumprobs)
}
# Function to generate Uniform rv
runif_fun<-function(n){
  u<-runif(n)
  return(u)
}
# Function to count all the documents in the dowloaded wikipedia corpora  
count_all_docs <- function(data_in){
  
  data_in_folder <- list.files(data_in)
  
  corpora_docs_count<-0
  
  for( f in 1:length(data_in_folder) ){
    
    f_list_files <- list.files( paste0(data_in, "/", data_in_folder[f]  ) )
    
    for(j in 1:length( f_list_files ) ){
      
      file_path <- paste0(data_in, "/", data_in_folder[f], "/", f_list_files[j]  )
      
      con <- file(file_path, open="r")
      
      all_docs_raw<-readLines(con, n=-1, warn = FALSE)
      close(con)
      
      corpora_docs_count <- corpora_docs_count + length(all_docs_raw)
      
      rm( all_docs_raw )
      
      print( paste0( "Folder ", data_in_folder[f], ", file ", f_list_files[j] ) )
      
    }
    gc()
    
  }
  return(corpora_docs_count)
}
# Function to check if the required number of documents have been sampled (used in the next function)
online_wiki_missing_docs_check <- function(data_in, data_out, batch_output_name, n_docs, seed) {
  
  batch_df <- readRDS(file = paste0(data_out, "/", batch_output_name, "_details.rds") )
  
  n_docs_new <- n_docs - nrow(batch_df)
  
  empty_file <- data.frame(Folder_name="", File_name="" )
  
  while( n_docs_new > 0 ){
    count_read<-0
    
    data_in_folder <- list.files(data_in)
    
    rep_flag <- ifelse(n_docs_new > length(data_in_folder), TRUE, FALSE )
    
    data_folder_sam <- sample(data_in_folder, size = n_docs_new, replace = rep_flag )
    sam_table<-data.frame(table(data_folder_sam))
    
    batch_df_new <- data.frame()
    
    for( f in 1:nrow(sam_table) ){
      f_det<-sam_table %>% filter( row_number()==f)
      folder_name <- f_det %>% select(data_folder_sam) %>% as.matrix() %>% as.character()
      
      f_list_files <- setdiff( list.files( paste0(data_in, "/", folder_name  ) ), 
                               empty_file %>% filter( Folder_name %in% folder_name ) %>% 
                                 select(File_name) %>% as.matrix() %>% as.character() )
      doc_freq <- f_det %>% select(Freq) %>% as.numeric()
      
      rep2_flag = ifelse( doc_freq > length(f_list_files), TRUE, FALSE )
      
      files_samp_long<-sample( f_list_files, doc_freq, replace = rep2_flag  )
      
      files_samp<-data.frame(table(files_samp_long))
      
      
      for(j in 1:nrow(files_samp)){
        
        file_det <- files_samp %>% filter(row_number()==j)
        file_name <- file_det %>% select(files_samp_long) %>% as.matrix() %>% as.character()
        
        file_path <- paste0(data_in, "/", f_det %>% 
                              select(data_folder_sam) %>% 
                              as.matrix() %>% 
                              as.character(), "/", file_name  )
        
        file_freq <- file_det %>% select(Freq) %>% as.numeric()
        
        con <- file(file_path, open="r")
        
        all_docs_raw<-readLines(con, n=-1, warn = FALSE)
        close(con)
        
        sampled_doc_set<-batch_df %>% filter(Folder_name %in% folder_name & File_name %in% file_name ) %>% 
          select(Doc_row) %>% as.matrix() %>% as.numeric() %>% as.character()
        
        sampled_doc_set_ori<-sampled_doc_set
        
        doc_set<-as.character( setdiff( 1:length(all_docs_raw), sampled_doc_set) )
        
        n_doc_ff<-min(file_freq, length(doc_set))
        
        doc_list <- list()
        
        i<-1
        while( (i %in% 1:n_doc_ff) & length(doc_set)>0 ){
          
          doc_idx<-sample(doc_set, 1, replace = F )
          doc_set<-setdiff(doc_set, doc_idx )
          doc_jas<-fromJSON(all_docs_raw[[ as.numeric( doc_idx ) ]])
          
          while(nchar(doc_jas$text)==0  & length(doc_set)>0 ){
            doc_idx<-sample(doc_set, 1)
            doc_set<-setdiff(doc_set, doc_idx )
            doc_jas<-fromJSON(all_docs_raw[[as.numeric(doc_idx)]])
          }
          if( nchar(doc_jas$text)>0 & !(doc_idx %in% sampled_doc_set) ){
            sampled_doc_set<-c(sampled_doc_set, doc_idx)
            doc_list <- append(doc_list, doc_jas)
            count_read <- count_read + 1
            print( paste0( "Processed ",  round(count_read/n_docs_new * 100, 2 ), "%" ) )
          }
          i <- i + 1
        }

        
        if( length(doc_list)==0 ){
          
          empty_file<-rbind(empty_file,
                            data.frame(Folder_name=folder_name, File_name=file_name ))
          empty_file<-empty_file %>% filter( !( Folder_name=="" & File_name==""  )  )
        }
        
        if( length(sampled_doc_set)>length(sampled_doc_set_ori) ){
          
          batch_df_new<-rbind(batch_df_new, data.frame( Folder_name=folder_name, 
                                                        File_name=file_name,
                                                        Doc_row=as.numeric(setdiff(sampled_doc_set,
                                                                                   sampled_doc_set_ori) )))
        }
      }  
    }
    
    batch_df_up <- rbind( batch_df_new, batch_df )
    
    batch_df_new <-  batch_df_up %>% arrange(Folder_name, File_name, Doc_row)
    rm( batch_df_up )
    
    saveRDS(batch_df_new, file = paste0(data_out, "/", batch_output_name, "_details.rds"))
    
    batch_df <- batch_df_new
    
    n_docs_new <- n_docs - nrow(batch_df)
    
    print(paste0( count_read, " missing documents have been sampled and added to the batch list. ", 
                  n_docs_new, " documents missing."  ))
  }
  
}
# Function to read (non-empty) "n_docs" documents at random in raw format from jason files. Output list and Output details data.frame stored in data_out path. 
online_wiki <- function(data_in, data_out, batch_output_name, n_docs, seed){
  data_in_folder <- list.files(data_in)
  
  rep_flag <- ifelse(n_docs > length(data_in_folder), TRUE, FALSE )
  
  data_folder_sam <- sample(data_in_folder, size = n_docs, replace = rep_flag )
  sam_table<-data.frame(table(data_folder_sam))
  
  count_read<-0
  batch_df <- data.frame()

  set.seed(seed)
  for( f in 1:nrow(sam_table) ){
    f_det<-sam_table %>% filter( row_number()==f)
    folder_name <- f_det %>% select(data_folder_sam) %>% as.matrix() %>% as.character()
    
    f_list_files <- list.files( paste0(data_in, "/", folder_name  ) )
    doc_freq <- f_det %>% select(Freq) %>% as.numeric()
    
    rep2_flag = ifelse( doc_freq > length(f_list_files), TRUE, FALSE )
    
    files_samp_long<-sample( f_list_files, doc_freq, replace = rep2_flag  )
    
    files_samp<-data.frame(table(files_samp_long))
    
    for(j in 1:nrow(files_samp)){
      
      file_det <- files_samp %>% filter(row_number()==j)
      file_name <- file_det %>% select(files_samp_long) %>% as.matrix() %>% as.character()
      
      file_path <- paste0(data_in, "/", f_det %>% 
                            select(data_folder_sam) %>% 
                            as.matrix() %>% 
                            as.character(), "/", file_name  )
      
      file_freq <- file_det %>% select(Freq) %>% as.numeric()
      
      con <- file(file_path, open="r")
      
      all_docs_raw<-readLines(con, n=-1, warn = FALSE)
      close(con)
      
      doc_set<- as.character( 1:length(all_docs_raw) )
      
      n_doc_ff<-min(file_freq, length(all_docs_raw))
      
      doc_list <- list()
      sampled_doc_set<-c()
      
      i<-1
      while( (i %in% 1:n_doc_ff) & length(doc_set)>0 ){
        
        doc_idx<-sample(doc_set, 1, replace = F )
        doc_set<-setdiff(doc_set, doc_idx )
        doc_jas<-fromJSON(all_docs_raw[[ as.numeric(doc_idx) ]])
        
        while(nchar(doc_jas$text)==0  & length(doc_set)>0 ){
          doc_idx<-sample(doc_set, 1, replace = F)
          doc_set<-setdiff(doc_set, doc_idx )
          doc_jas<-fromJSON(all_docs_raw[[as.numeric(doc_idx)]])
        }
        if( nchar(doc_jas$text)>0 & !(doc_idx %in% sampled_doc_set) ){
          sampled_doc_set<-c(sampled_doc_set, doc_idx)
          doc_list <- append(doc_list, doc_jas)
          count_read <- count_read + 1
          print( paste0( "Processed ",  round(count_read/n_docs * 100, 2 ), "%" ) )
        }
        i <- i + 1
      }
      
      # batch_list<-append( batch_list, doc_list)
      
      batch_df<-rbind(batch_df, data.frame( Folder_name=folder_name, 
                                            File_name=file_name,
                                            Doc_row=as.numeric(sampled_doc_set)
      ) )
    }
  }
  
  saveRDS(batch_df, file = paste0(data_out, "/", batch_output_name, "_details.rds"))
  
  online_wiki_missing_docs_check(data_in, data_out, batch_output_name, n_docs, seed)
}
# Function to parse the documents and store train/validation and test lists of vectors
parse_doc_parallel_train_test <- function(data_in, data_batch, batch_output_name, test_perc, n_test_set, vocab, n_cores, seed){
  parse_test_doc <- function(raw_doc, vocab){
    
    text<-raw_doc$text
    text_lower_case<-str_to_lower(text)
    # Replace '-' with space
    text_lower_case <- gsub('-', ' ', text_lower_case)
    # Remove non-alphabetic characters except space
    text_lower_case <- gsub('[^a-z ]', ' ', text_lower_case)
    # Replace multiple spaces with a single space
    text_clean <- gsub(' +', ' ', text_lower_case)
    
    doc_vec<-sapply(1:W, function(w){
      word <- paste0("\\b", vocab %>% filter(row_number()==w) %>% 
                       as.matrix() %>% as.character(), "\\b")
      str_count(text_clean, word)
    })
    
    return(doc_vec)
  }
  
  batch_df <- readRDS(file = paste0(data_batch, "/", batch_output_name, "_details.rds") )
  
  W<-nrow(vocab)
  test_every<-round(test_perc*100, 0)
  
  n_docs<-nrow(batch_df)
  n_docs_train<-n_docs-n_test_set
  
  set.seed(seed)
  test_set_idx<-sample(1:n_docs, n_test_set, replace = F)
  
  
  batch_test_df <- batch_df %>% filter( (row_number() %in% test_set_idx) )
  
  batch_test_df_short<-batch_test_df %>% 
    group_by(Folder_name, File_name) %>% 
    summarise(Rows = paste0(Doc_row, collapse = ", "))  %>% ungroup()
  
  batch_train_df <- batch_df %>% filter( !(row_number() %in% test_set_idx) ) %>%
    mutate(idx=sample(1:n_docs_train, size=n_docs_train, replace=F  )) %>%
    arrange(idx) %>% select(-idx)
  
  
  test_set_doc_det<-data.frame( )
  test_list<-list()
  count_parsed_test<-0
  
  for(d in 1:nrow(batch_test_df_short)){
    doc_file<-batch_test_df_short %>% filter(row_number() %in% d)
    
    folder_name <- doc_file %>% select(Folder_name) %>% as.matrix() %>% as.character()
    file_name <- doc_file %>% select(File_name) %>% as.matrix() %>% as.character()
    doc_rows <-  doc_file %>% select(Rows) %>% as.matrix() %>% as.character() %>% 
      str_split( pattern = ", ") %>% unlist() %>% as.numeric()
    
    file_path <- paste0(data_in, "/", folder_name, "/", file_name  )
    
    con <- file(file_path, open="r")
    all_docs_raw<-readLines(con, n=-1, warn = FALSE)
    close(con)  
    
    for(r in 1:length(doc_rows)){
      raw_doc<-fromJSON(all_docs_raw[[doc_rows[r]]])
      parsed_doc<-parse_test_doc(raw_doc, vocab)
      test_list<-append(test_list, list(parsed_doc))
      
      test_set_doc_det<-rbind( test_set_doc_det,
                               data.frame(Folder_name=folder_name, File_name=file_name, Row=doc_rows[r], 
                                          id=raw_doc$id, revid=raw_doc$revid, 
                                          url=raw_doc$url, title=raw_doc$title ) ) 
      saveRDS(test_list, paste0(data_batch, "/test_doc_parsed.rds" ))
      saveRDS(test_set_doc_det, paste0(data_batch, "/test_set_doc_det.rds" ))
      
      count_parsed_test <- count_parsed_test + 1
      print(paste0("Parsed ", round(count_parsed_test/n_test_set * 100, 2 ), "% of the test set documents"))
    }
  }
  
  print( "Starting parsing the train set documents in parallel" )
  
  n_docs_core<-c(rep(ceiling(n_docs_train/n_cores), n_cores-1),
                 n_docs_train-sum( rep(ceiling(n_docs_train/n_cores), n_cores-1) ))
  
  core_idx<-rbind( c(1, cumsum(n_docs_core)[-n_cores]+1), 
                   c( cumsum(n_docs_core)[-n_cores], n_docs_train ) )
  
  
  myCluster <- makeCluster(n_cores, # number of cores to use
                           type = "PSOCK") # type of cluster
  
  st_s<-Sys.time()
  registerDoParallel(myCluster)
  foreach(c = 1:n_cores, 
          .packages = c("stringr", 
                        "tidyverse", 
                        "jsonlite")
  ) %dopar% {
    
    parse_train_doc <- function(raw_doc, vocab, test_every){
      
      text<-raw_doc$text
      text_lower_case<-str_to_lower(text)
      # Replace '-' with space
      text_lower_case <- gsub('-', ' ', text_lower_case)
      # Remove non-alphabetic characters except space
      text_lower_case <- gsub('[^a-z ]', ' ', text_lower_case)
      # Replace multiple spaces with a single space
      text_clean <- gsub(' +', ' ', text_lower_case)
      
      space_pos<-str_locate_all(text_clean, " ")[[1]]
      
      if( nrow(space_pos)>=test_every){
        if(as.numeric(space_pos[nrow(space_pos),1])==nchar(text_clean)){
          text_clean <-substr( text_clean, start = 1, stop = (nchar(text_clean)-1) )
        }
        space_pos<-str_locate_all(text_clean, " ")[[1]][,1]
        
        word_pos_ini<-c(1, space_pos+1)
        word_pos_fin<-c(space_pos-1, nchar(text_clean))
        
        word_test_pos<-seq(test_every, length(word_pos_ini), test_every)
        
        ini_train_idx<-word_pos_ini[-word_test_pos]
        fin_train_idx<-word_pos_fin[-word_test_pos]
        
        ini_test_idx<-word_pos_ini[word_test_pos]
        fin_test_idx<-word_pos_fin[word_test_pos] 
        
        s<-1
        train_text<-substr( text_clean, start = ini_train_idx[s], stop = fin_train_idx[s] )
        if( length(ini_train_idx) > 1 ){
          for( s in 2:length(ini_train_idx) ){
            train_text<-paste(train_text, substr( text_clean, start = ini_train_idx[s], stop = fin_train_idx[s] ), sep = " " )
          }
        }
        
        s<-1
        test_text<-substr( text_clean, start = ini_test_idx[s], stop = fin_test_idx[s] )
        if( length(ini_test_idx) > 1 ){
          for( s in 2:length(ini_test_idx) ){
            test_text<-paste(test_text, substr( text_clean, start = ini_test_idx[s], stop = fin_test_idx[s] ), sep = " " )
          }
        }
        
        doc_vec<-sapply(1:W, function(w){
          word <- paste0("\\b", vocab %>% filter(row_number()==w) %>% 
                           as.matrix() %>% as.character(), "\\b")
          c(str_count(train_text, word), str_count(test_text, word))
        })
      } 
      else{
        if( nrow(space_pos)>0 ){
          if(as.numeric(space_pos[nrow(space_pos),1])==nchar(text_clean)){
            text_clean <-substr( text_clean, start = 1, stop = (nchar(text_clean)-1) )
          }
          train_text<-text_clean
          test_text<-""
          
          doc_vec<-sapply(1:W, function(w){
            word <- paste0("\\b", vocab %>% filter(row_number()==w) %>% 
                             as.matrix() %>% as.character(), "\\b")
            c(str_count(train_text, word), str_count(test_text, word))
          })
        }
        else{
          train_text<-text_clean
          test_text<-""
          
          doc_vec<-sapply(1:W, function(w){
            word <- paste0("\\b", vocab %>% filter(row_number()==w) %>% 
                             as.matrix() %>% as.character(), "\\b")
            c(str_count(train_text, word), str_count(test_text, word))
          })
        }
      }
      
      return(doc_vec)
    }
    
    batch_train_df_chunck <- batch_train_df %>% filter( row_number() %in% ( core_idx[1, c]:core_idx[2, c] )  )
    
    batch_train_df_short <- batch_train_df_chunck %>% 
      group_by(Folder_name, File_name) %>% 
      summarise(Rows = paste0(Doc_row, collapse = ", ")) %>% ungroup()
    
    
    train_set_doc_det<-data.frame( )
    train_list<-list()
    valid_list<-list()
    
    count_parsed_train<-0
    for(d in 1:nrow(batch_train_df_short)){
      
      doc_file<-batch_train_df_short %>% filter(row_number() %in% d)
      
      folder_name <- doc_file %>% select(Folder_name) %>% as.matrix() %>% as.character()
      file_name <- doc_file %>% select(File_name) %>% as.matrix() %>% as.character()
      doc_rows <-  doc_file %>% select(Rows) %>% as.matrix() %>% as.character() %>% 
        str_split( pattern = ", ") %>% unlist() %>% as.numeric()
      
      file_path <- paste0(data_in, "/", folder_name, "/", file_name  )
      
      con <- file(file_path, open="r")
      all_docs_raw<-readLines(con, n=-1, warn = FALSE)
      close(con)  
      
      for(r in 1:length(doc_rows)){
        raw_doc<-fromJSON(all_docs_raw[[doc_rows[r]]])
        
        parsed_doc<-parse_train_doc(raw_doc, vocab, test_every)
        
        train_list<-append(train_list, list(parsed_doc[1,]))
        valid_list<-append(valid_list, list(parsed_doc[2,]))
        
        
        train_set_doc_det<-rbind( train_set_doc_det,
                                  data.frame(Folder_name=folder_name, File_name=file_name, Row=doc_rows[r], 
                                             id=raw_doc$id, revid=raw_doc$revid, 
                                             url=raw_doc$url, title=raw_doc$title ) ) 
        
        saveRDS(train_list, paste0(data_batch, "/train_doc_parsed_", c, ".rds" ))
        saveRDS(valid_list, paste0(data_batch, "/valid_doc_parsed_", c, ".rds" ))
        
        saveRDS(train_set_doc_det, paste0(data_batch, "/train_set_doc_det_", c, ".rds" ))
        
        count_parsed_train <- count_parsed_train + 1
        print(paste0("Parsed ", round(count_parsed_train/nrow(batch_train_df_chunck) * 100, 2 ), 
                     "% of the train set documents by core ", c))
      }
    }
    
  }
  stopCluster(myCluster)
  en_s<-Sys.time()
  en_s-st_s
  print(paste0( "Total time: ", en_s-st_s ))
}
# Function to parse the hold_out documents and store in two lists of vectors
parse_hold_out<-function(data_batch, vocab, test_perc){
  parse_train_doc <- function(raw_doc, vocab, test_every){
    
    W<-nrow(vocab)
    
    text<-raw_doc$text
    text_lower_case<-str_to_lower(text)
    # Replace '-' with space
    text_lower_case <- gsub('-', ' ', text_lower_case)
    # Remove non-alphabetic characters except space
    text_lower_case <- gsub('[^a-z ]', ' ', text_lower_case)
    # Replace multiple spaces with a single space
    text_clean <- gsub(' +', ' ', text_lower_case)
    
    space_pos<-str_locate_all(text_clean, " ")[[1]]
    
    if( nrow(space_pos)>=test_every){
      if(as.numeric(space_pos[nrow(space_pos),1])==nchar(text_clean)){
        text_clean <-substr( text_clean, start = 1, stop = (nchar(text_clean)-1) )
      }
      space_pos<-str_locate_all(text_clean, " ")[[1]][,1]
      
      word_pos_ini<-c(1, space_pos+1)
      word_pos_fin<-c(space_pos-1, nchar(text_clean))
      
      word_test_pos<-seq(test_every, length(word_pos_ini), test_every)
      
      ini_train_idx<-word_pos_ini[-word_test_pos]
      fin_train_idx<-word_pos_fin[-word_test_pos]
      
      ini_test_idx<-word_pos_ini[word_test_pos]
      fin_test_idx<-word_pos_fin[word_test_pos] 
      
      s<-1
      train_text<-substr( text_clean, start = ini_train_idx[s], stop = fin_train_idx[s] )
      if( length(ini_train_idx) > 1 ){
        for( s in 2:length(ini_train_idx) ){
          train_text<-paste(train_text, substr( text_clean, start = ini_train_idx[s], stop = fin_train_idx[s] ), sep = " " )
        }
      }
      
      s<-1
      test_text<-substr( text_clean, start = ini_test_idx[s], stop = fin_test_idx[s] )
      if( length(ini_test_idx) > 1 ){
        for( s in 2:length(ini_test_idx) ){
          test_text<-paste(test_text, substr( text_clean, start = ini_test_idx[s], stop = fin_test_idx[s] ), sep = " " )
        }
      }
      
      doc_vec<-sapply(1:W, function(w){
        word <- paste0("\\b", vocab %>% filter(row_number()==w) %>% 
                         as.matrix() %>% as.character(), "\\b")
        c(str_count(train_text, word), str_count(test_text, word))
      })
    } 
    else{
      if( nrow(space_pos)>0 ){
        if(as.numeric(space_pos[nrow(space_pos),1])==nchar(text_clean)){
          text_clean <-substr( text_clean, start = 1, stop = (nchar(text_clean)-1) )
        }
        train_text<-text_clean
        test_text<-""
        
        doc_vec<-sapply(1:W, function(w){
          word <- paste0("\\b", vocab %>% filter(row_number()==w) %>% 
                           as.matrix() %>% as.character(), "\\b")
          c(str_count(train_text, word), str_count(test_text, word))
        })
      }
      else{
        train_text<-text_clean
        test_text<-""
        
        doc_vec<-sapply(1:W, function(w){
          word <- paste0("\\b", vocab %>% filter(row_number()==w) %>% 
                           as.matrix() %>% as.character(), "\\b")
          c(str_count(train_text, word), str_count(test_text, word))
        })
      }
    }
    
    return(doc_vec)
  }
  
  test_every<-round(test_perc*100, 0)
  
  test_set_doc_det <- readRDS(paste0(data_batch, "/test_set_doc_det.rds") )
  
  ho_train<-list()
  ho_test<-list()

  for(d in 1:nrow(test_set_doc_det)){
    doc_file<-test_set_doc_det  %>% filter(row_number() %in% d)
    
    folder_name <- doc_file %>% select(Folder_name) %>% as.matrix() %>% as.character()
    file_name <- doc_file %>% select(File_name) %>% as.matrix() %>% as.character()
    doc_rows <-  doc_file %>% select(Row) %>% as.matrix() %>% as.numeric()
    
    file_path <- paste0(data_in, "/", folder_name, "/", file_name  )
    
    con <- file(file_path, open="r")
    all_docs_raw<-readLines(con, n=-1, warn = FALSE)
    close(con)  
    
    raw_doc<-fromJSON(all_docs_raw[[doc_rows]])
    
    parsed_doc<-parse_train_doc(raw_doc, vocab, test_every)
    
    ho_train[[d]]<-parsed_doc[1,] 
    ho_test[[d]]<-parsed_doc[2,] 

    print(paste0("Document ", d))
  }  
  saveRDS( ho_train, paste0(data_batch, "/ho_train.rds") )
  saveRDS( ho_test, paste0(data_batch, "/ho_test.rds") )
}
# Function to estimate the sum of the topic assignment variables in the sample by parallelizing the task across cores
sample_counts_Rcpp_parallel<- function(n_cores_used, theta, alpha, 
                                       probs_init, docs_list, 
                                       gibbs_iters, gibbs_burnin) {
  
  N_doc<-length(docs_list)
  
  # K<-nrow(theta)
  # W<-ncol(theta)
  
  n_obs_per_core<-c(rep(ceiling(N_doc/n_cores_used), n_cores_used-1), 
                    N_doc-sum(rep(ceiling(N_doc/n_cores_used), n_cores_used-1)))
  
  st_idx<-c(1, cumsum(n_obs_per_core)[-length(cumsum(n_obs_per_core))]+1)
  en_idx<-cumsum(n_obs_per_core)
  
  core_docs_list_all<-list()
  for( m in 1:n_cores_used){
    core_docs_list_all[[m]]<-docs_list[c(st_idx[m]:en_idx[m])]
  }
  
  
  myCluster <- makeCluster(n_cores_used, # number of cores to use
                           type = "PSOCK") # type of cluster
  
  # st_s<-Sys.time()
  registerDoParallel(myCluster)
  out_list <- foreach(m = 1:n_cores_used, 
                      .combine = 'c', 
                      .packages = c("LDARcppPackage")
  ) %dopar% {
    
    # Bkw<-matrix(0, nrow = K, ncol = W)
    # Bkw_avg<-matrix(0, nrow = K, ncol = W)
    
    core_docs_list<-core_docs_list_all[[m]]
    
    # D<-length(core_docs_list)
    # Adk<-matrix(0, nrow = D, ncol = K)
    # Adk_avg<-matrix(0, nrow = D, ncol = K)
    
    set.seed(m)
    sam_cnt_Rcpp<-rcpp_sample_cnt(alpha, probs_init, core_docs_list, theta,  
                                  # Adk, Bkw, Adk_avg, Bkw_avg, 
                                  gibbs_iters, gibbs_burnin)
    return(sam_cnt_Rcpp)
  }
  stopCluster(myCluster)
  
  A_idx<-(1:length(out_list))[((1:length(out_list))%%2)==1]
  B_idx<-setdiff(1:length(out_list), A_idx)
  
  m<-1
  A_dk<-out_list[[A_idx[m]]]
  B_kw<-out_list[[B_idx[m]]]
  
  if(length(A_idx)){
    for(m in 2:length(A_idx)){
      A_dk<-rbind(A_dk, out_list[[A_idx[m]]])
      B_kw <- B_kw + out_list[[B_idx[m]]]
    }
  }
  
  A_dk <- A_dk/(gibbs_iters - gibbs_burnin)
  B_kw <- B_kw/(gibbs_iters - gibbs_burnin)
  
  return(list(Adk_avg=A_dk, Bkw_avg=B_kw))
  
}
# Function to compute perplexity on validation dataset
log_pred <- function(alpha, A_avg,  theta, valid_list){
  
  phi <- theta/rowSums(theta) 
  
  eta_hat <- alpha + A_avg
  eta_hat <- eta_hat/rowSums(eta_hat)
  
  lg_pr<-log(eta_hat %*% phi)
  
  valid_cnt<-sapply(1:length(valid_list), function(d){
    idx<-which(valid_list[[d]]>0)
    
    n_words<-valid_list[[d]][idx]
    
    log_p<-n_words*lg_pr[d,idx]
    # log_p<-valid_list[[d]][idx]*lg_pr[d,idx]
    
    return(c( sum(log_p), sum(n_words) ))
  })
  
  log_perp_vec<-rowSums(valid_cnt)
  log_perp_valid<-log_perp_vec[1]/log_perp_vec[2]
  return(log_perp_valid)
}
# Function to compute perplexity on test dataset
hold_out_log_pred <- function(alpha, probs_init, 
                              ho_train, ho_test, 
                              ho_count, old_avg, 
                              theta,  
                              gibbs_iters, gibbs_burnin){
  
  
  sam_cnt_Rcpp<-sample_counts_Rcpp_parallel(n_cores_used, 
                                            theta, 
                                            alpha, probs_init, 
                                            ho_train, 
                                            gibbs_iters, gibbs_burnin)
  
  
  A_avg<-sam_cnt_Rcpp[[1]]
  
  phi <- theta/rowSums(theta)
  
  eta_hat <- alpha + A_avg
  eta_hat <- eta_hat/rowSums(eta_hat)
  
  avg_probs<-(ho_count*old_avg + eta_hat %*% phi)/(ho_count+1)
  
  old_avg<-avg_probs
  ho_count <- ho_count + 1 
  
  test_cnt<-sapply(1:length(ho_test), function(d){
    idx<-which(ho_test[[d]]>0)
    
    n_words<-ho_test[[d]][idx]
    log_p<-n_words*log(avg_probs[d,idx])
    
    return(c( sum(log_p), sum(n_words) ))
  })
  
  log_perp_vec<-rowSums(test_cnt)
  log_perp_test<-log_perp_vec[1]/log_perp_vec[2]
  return(list( log_perp_test=log_perp_test, ho_count=ho_count, old_avg=old_avg))
}
