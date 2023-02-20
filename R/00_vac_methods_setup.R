###R library set up
library(dplyr)
library(data.table)
library(lubridate)
library(janitor)
library(haven)
library(cli)
library(arrow)
library(foreach)
library(ggplot2)
library(hash)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(ggalluvial)


###recurrent functions
#make sure that path for list.files is set correctly
#extract medcodes
extractor_med <- function(list.files, codelist){
  
  df <- as.data.table(read_parquet(list.files[[1]]))
  print(paste0(length(unique(df$patid)), " in original"))
  df <- df[medcodeid %chin% codelist[, medcodeid]]
  print(paste0(length(unique(df$patid)), " subjects with code extracted"))
  df <- merge(df, codelist, by ="medcodeid")
  
  for(i in 2:length(list.files)){
    df_temp <- as.data.table(read_parquet(list.files[[i]]))
    print(paste0(length(unique(df_temp$patid)), " in original"))
    df_temp <- df_temp[medcodeid %chin% codelist[, medcodeid]]
    print(paste0(length(unique(df_temp$patid)), " subjects with code extracted"))
    df_temp <- merge(df_temp, codelist, by = "medcodeid")
    
    df <- rbind(df, df_temp)
  }
  
  return(df)
}

#extract prodcodes
extractor_prod <- function(list.files, codelist){
  
  df <- as.data.table(read_parquet(list.files[[1]]))
  print(paste0(length(unique(df$patid)), " in original"))
  df <- df[prodcodeid %chin% codelist[, prodcodeid]]
  print(paste0(length(unique(df$patid)), " subjects with code extracted"))
  df <- merge(df, codelist, by ="prodcodeid")
  
  for(i in 2:length(list.files)){
    df_temp <- as.data.table(read_parquet(list.files[[i]]))
    print(paste0(length(unique(df_temp$patid)), " in original"))
    df_temp <- df_temp[prodcodeid %chin% codelist[, prodcodeid]]
    print(paste0(length(unique(df_temp$patid)), " subjects with code extracted"))
    df_temp <- merge(df_temp, codelist, by = "prodcodeid")
    
    df <- rbind(df, df_temp)
  }
  
  return(df)
}


#rounding the results
rounding <- function(X){
  n <- round(X, digits = 2)
  return(n)
}

#keeping the max exclugind NA values
max_NA <- function(x){
  return(sum(x, na.rm=T))
}

#---function for the data cleaning
#clean out date before minimum date
clean_min_vac_date <- function(df, dose, min_age){
  df[n_dose == dose, keep := age >=min_age]
  per <- round((nrow(df[keep ==F])/nrow(df)*100), digits = 2)
  print(paste0(nrow(df[keep ==F]), " vaccine records dropped (",
               per, "%)")) 
  df <- df[keep == T |is.na(keep)]
  df[, keep := NULL]
  print(paste0(nrow(df), " entries remaining in df")) #3334178
  #renumbering the vaccines
  df[DTP_vacc == "vaccinated", n_dose := seq_len(.N), by = "patid"]
  df[DTP_vacc == "vaccinated", n_dose := rowid(patid)]
  df[DTP_vacc == "vaccinated", by = .(patid, DTP_vacc),
     DTP_all_doses :=  .N]
  return(df)
}

#checking time diff between dose a and b
#depent on max_NA function from above
#checking now for the time difference between the two doses
check_min_diff <- function(df, dose_a, dose_b, min_diff){
  
  df[n_dose == dose_a, age1 := age]
  df[n_dose == dose_b, age2 := age]
  df[DTP_all_doses > dose_a, age_1 := lapply(.SD, max_NA), .SDcols = c("age1"), by="patid"]
  df[DTP_all_doses > dose_a, age_2 := lapply(.SD, max_NA), .SDcols = c("age2"), by="patid"]
  df[DTP_all_doses > dose_a, age_diff_12 := age_2-age_1]
  df[!is.na(age_diff_12), keep := age_diff_12 >=min_diff]
  per <- round((nrow(df[keep ==F])/nrow(df)*100), digits = 2)
  print(paste0(nrow(df[keep ==F]), " vaccine records dropped (",
               per, "%)")) 
  df <- df[keep == T |is.na(keep)]
  df[, keep := NULL]
  df[, age1:= NULL]
  df[, age2:= NULL]
  df[, age_1:= NULL]
  df[, age_2:= NULL]
  df[, age_diff_12:= NULL] 
  print(paste0(nrow(df), " entries remaining in df")) 
  
  #renumbering the vaccines
  df[DTP_vacc == "vaccinated", n_dose := seq_len(.N), by = "patid"]
  df[DTP_vacc == "vaccinated", n_dose := rowid(patid)]
  df[DTP_vacc == "vaccinated", by = .(patid, DTP_vacc),
     DTP_all_doses :=  .N]
  return(df)
}



#---vaccination coverage at ages 1, 2, and 5 plus 95%-CI
vacc_age <- function(data, vac, age_vac, start, end, vac_dose){
  
  n <- nrow(data[vaccine == vac &
                   age <= age_vac &
                   age_startfu <= start &
                   age_endfu >= end &
                   n_dose == vac_dose])
  return(n)
}

pop_vacc <- function(data_p, start, end){
  
  n <- nrow(data_p[age_startfu <= start & age_endfu >= end])
  return(n)
}


coverage_1y <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 365, start = 40, end = 365, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 365)
  
  return((a/b)*100)
}
coverage_1y_lb <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 365, start = 40, end = 365, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 365)
  p_hat <- a/b
  res <- p_hat -1.96*(sqrt(p_hat*(1-p_hat)/b))
  return(res*100)
}
coverage_1y_ub <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 365, start = 40, end = 365, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 365)
  p_hat <- a/b
  res <- p_hat +1.96*(sqrt(p_hat*(1-p_hat)/b))
  return(res*100)
}




coverage_2y <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 731, start = 40, end = 731, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 731)
  
  return((a/b)*100)
}
coverage_2y_lb <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 731, start = 40, end = 731, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 731)
  p_hat <- a/b
  res <- p_hat -1.96*(sqrt(p_hat*(1-p_hat)/b))
  return(res*100)
}
coverage_2y_ub <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 731, start = 40, end = 731, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 731)
  p_hat <- a/b
  res <- p_hat +1.96*(sqrt(p_hat*(1-p_hat)/b))
  return(res*100)
}


coverage_3y <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 1096, start = 40, end = 1096, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 1096)
  
  return((a/b)*100)
}

coverage_4y <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 1461, start = 40, end = 1461, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 1461)
  
  return((a/b)*100)
}


coverage_5y <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 1826, start = 40, end = 1826, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 1826)
  
  return((a/b)*100)
}
coverage_5y_lb <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 1826, start = 40, end = 1826, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 1826)
  p_hat <- a/b
  res <- p_hat -1.96*(sqrt(p_hat*(1-p_hat)/b))
  return(res*100)
}
coverage_5y_ub <- function(vacc, dose, df, df_p){
  
  a <- vacc_age(data = df, vac = vacc, age = 1826, start = 40, end = 1826, vac_dose = dose)
  b <-  pop_vacc(data_p = df_p, start = 40, end = 1826)
  p_hat <- a/b
  res <- p_hat +1.96*(sqrt(p_hat*(1-p_hat)/b))
  return(res*100)
}

