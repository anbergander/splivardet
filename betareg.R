library(betareg)
library(rcompanion)
library(tidyverse)
library(car)
library(lmtest)
library(emmeans)

keys <- scan("isoform_keys", what="", sep="\n")
#keys <- commandArgs(trailingOnly = TRUE)
#keys <- c("00273bbc-2065-42be-b025-3a7a2408818c_chr13:41542000")
start_df <- read.csv(paste0(getwd(),"/",keys[1],".csv"))
i <- 1


for (key in keys) {
  i <- i+1
  #if(i%%10==0){
   # print(i)}
  #print(key)
  ex_df = read.csv(paste0(getwd(),"/",key,".csv"))
  
  if(sum(ex_df$value)/nrow(ex_df)!=1){
	
    model.beta <- try(betareg(transf_value ~ group, data=ex_df), silent=TRUE)
     if(class(model.beta)=="try-error") { 
	  print(key)
	  next 
	  }
	lrtest(model.beta)
	summary(model.beta)
  
	joint_tests(model.beta)
  
  resAn <- Anova(model.beta)
  if (resAn[1,3] < 1.05){
    Post_Hoc <- emmeans (model.beta, specs = pairwise ~ group, type = "response")
    resph <- summary(Post_Hoc)
    resph_df = resph[["contrasts"]]
    resph_df <- resph_df %>% add_column(isoform_key = key)
    if (exists("posthoc_df")){
      posthoc_df <- rbind(posthoc_df, resph_df)
    }
    else{
      posthoc_df <- resph_df
      
    }
  }
  }
}

write.csv(posthoc_df,paste0("posthoc",".csv"), row.names=TRUE)



