library(tidyverse)
library(stringr)
library(dplyr)
library(data.table)

#reading in the data dictionaries for CPRD Aurum release May 2022
med_aurum<- read.delim(".../CPRD Aurum/Code browsers/2022_05/CPRDAurumMedical.txt", 
                       sep= "\t",  colClasses="character")
prod_aurum <- read.delim("...CPRD Aurum/Code browsers/2022_05/CPRDAurumProduct.txt", 
                         sep= "\t",  colClasses="character")


#function to look up medical terms
term_lookup <- function(terms){
  #result list for all the searched terms
  list_terms <- list()
  
  for(i in 1:length(terms)){
    
    index <- grep(terms[i], med_aurum$Term, ignore.case=T)
    list_terms[[i]]<- med_aurum[index,]
  }
  
  #generate one data frame  after deduplicating entries
  df <- list_terms[[1]]
  for(i in 2:length(list_terms)){
    df <- unique(rbind(df, list_terms[[i]]))
  }
  
  df<-df%>%
    arrange(OriginalReadCode)
  return(df)
}


#function to look up product codes
prod_lookup <- function(terms){
  
  list_terms <- list()
  
  for(i in 1:length(terms)){
    
    index_1 <- grep(terms[i], prod_aurum$Term.from.EMIS, ignore.case=T)
    index_2 <- grep(terms[i], prod_aurum$DrugSubstanceName, ignore.case=T)
    index <- c(index_1, index_2)
    list_terms[[i]]<- prod_aurum[index,]
    
  }
  
  #generate one data frame  after deduplicating entries
  df <- list_terms[[1]]
  for(i in 2:length(list_terms)){
    df <- unique(rbind(df, list_terms[[i]]))
  }
  
  df<-df%>%
    arrange(Term.from.EMIS)
  return(df)
}

###---looking up the codes

terms1 <- c("vacc", "immunisation", "immunization", "booster", "course")
terms2 <- c("MMR", "DTP", "Pneumococcal", "measles", "rubella", "mumps",
            "dipht", "tetanus", "pert", "Pneu", "Rota", "hepa",
            "influenza", "MMR", "chickenpox", "papilloma", "rabies", "varicella",
            "polio", "tularemia", "bcg", "mengingo", "ncov", "corona", 
            "anthrax", "brucella", "chlamydia", "did not attend", "cholera", 
            "smallpox", "flu", "allergy", "adverse", "japanese", "lyme", 
            "meningo", "poisoning", "fever", "adult", "infanrix")


#vaccine terms
df1 <- term_lookup(terms = terms1)
df2 <- term_lookup(terms = terms2)
#terms only in df1 but not in df2
df_terms <- setdiff(df1, df2)


###saving the preliminary lists
write.csv2(df_terms, file = "results/General_terms_raw.csv")
