library(tidyverse)
source(paste0(getwd(), "/Experiments/LDA/R/00_LDA_paths.R"))

load(paste0(env_dir, "/input.RData"))
rm(list=setdiff(ls(), c("perpl_every", "minibatch_size")))
source(paste0(getwd(), "/Experiments/LDA/R/00_LDA_paths.R"))

ho_perp_names<-data.frame(File=list.files(output_path)) %>% filter( str_detect(File, "ho") & !str_detect(File, "sim") )

ho_perp_df<-data.frame()
for( r in 1:nrow(ho_perp_names) ){
  
  perp_df_r<-data.frame(readRDS( file = paste0( output_path, "/", ho_perp_names %>% 
                                                  filter( row_number()==r ) %>% 
                                                  as.matrix() %>% 
                                                  as.character() ) )) %>% 
    rename(perp_SCIR=1, perp_CV=2, perp_SGRLD=3 )
  
  ho_perp_df <- rbind(ho_perp_df,
  cbind(Time_index= 1:nrow(perp_df_r), perp_df_r, Run=as.character(r)) )
}

ho_perp_df_longer <- ho_perp_df %>% pivot_longer(cols  = c(ho_perp_df %>% select( contains("perp") ) %>% names()), 
                            names_to = c(".value", "Method"),
                            names_sep = "\\_")
  
ho_perp_df_longer <- ho_perp_df_longer %>% mutate(Method_lab=case_when( Method=="CV" ~ "SCIR-CV", 
                                                                        TRUE ~ Method))

ho_perp_df_longer <- ho_perp_df_longer %>% select(-Method) %>% rename(Method=Method_lab)


ho_perp_df_summary<-ho_perp_df_longer %>% group_by(Method, Run) %>%
  arrange(Time_index) %>%
  mutate(running_value = cumsum(perp) / seq_along(perp)) %>%
  ungroup()

# Now compute summaries across chains
ho_perp_df_summary <- ho_perp_df_summary %>%
  group_by(Method, Time_index) %>%
  summarise(
    perp_mean = mean(running_value),
    sd = sd(running_value),
    .groups = "drop" ) %>%
  mutate( perp_min = perp_mean-sd, perp_max = perp_mean+sd)



(  pp<-ggplot( ho_perp_df_summary, aes( Time_index, perp_mean, 
                                 fill = Method,
                                 color=Method 
                                 )  ) +
    xlab("Number of documents") + ylab("Perplexity") +
    theme_bw() +
    geom_line() + 
    geom_ribbon( aes( ymin=perp_min, ymax=perp_max, fill=Method ), alpha=0.5, linetype=1 ) + 
    scale_fill_manual( values = c( 4, rgb(0.9290, 0.6940, 0.1250),  rgb(0, 0, 0, 0.5) ) ) +
    scale_color_manual( values = c( 4, rgb(0.9290, 0.6940, 0.1250), rgb(0, 0, 0, 0.5) ) ) +
    scale_x_continuous(breaks = seq(0, nrow(ho_perp_df_summary)/3, length.out=5),
                     labels = as.character( seq(0, nrow(ho_perp_df_summary)/3*perpl_every*minibatch_size, 
                                                length.out=5) ) )  + #perpl_every*minibatch_size 
    theme( 
      strip.text.x = element_text(size = 20),
      strip.text.y  = element_text(size = 20),
      axis.text=element_text(size=20),
      axis.title.x = element_text(size = 20),  
      axis.title.y = element_text(size = 20)   
      ) +
  theme(legend.position = "none")
)


ggsave(paste0(output_path,"/lda_wiki_perp.pdf"), 
       plot = pp, width = 20, height = 15, units = "cm", dpi = 300)  
  
