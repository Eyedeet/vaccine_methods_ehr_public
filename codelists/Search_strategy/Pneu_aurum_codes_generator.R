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

#function to look up snomed concept IDs
concept_lookup <- function(terms){
  #result list for all the searched terms
  list_terms <- list()
  
  for(i in 1:length(terms)){
    
    index <- grep(terms[i], med_aurum$SnomedCTConceptId, ignore.case=T)
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


###---looking up the initial codes

terms1 <- c("pneumococc", "pneumovax", "pnu-imune", "PPV", "PCV", "pneumococc",
            "pneumovax", "pnu-imune",  "PCV", "prevenar", "synflorix")
terms2 <- c("vac", "imm", "consent", "injec", "declin", "invit", "boost", 
            "not given", "call", "message", "dose", "syringe", "amp", "attenuated")


#vaccine terms
df1 <- term_lookup(terms = terms1)
df2 <- term_lookup(terms = terms2)
df_terms <- intersect(df1, df2)

#vaccine product codes
df1a <- prod_lookup(terms = terms1)
df2a <- prod_lookup(terms = terms2)

df_prods <- intersect(df1a, df2a)


#--- second step search, expanding for SNOMED concept IDs which weren't detected 
#related snomed concept IDs from pervious code list screening
snomed <- c("310578008", "12866006", "10100410000006110", "16660311000000100",
            "714821000000108","10245211000001108", "3439211000001108", "34783011000001102", "35111211000001108",
            "27396511000001105", "3439311000001100", "36461211000001108", "3005011000001103",
            "3017311000001106", "3018111000001105", "16649411000001104", "10231211000001106",
            "3021111000001108", "16660211000001102" )

concepts <- concept_lookup(terms = snomed)

#checking mismatch between already identified terms
missed <- setdiff(concepts, df_terms) # none are missing


###safing the preliminary lists
write.csv2(df_terms, file = "results/Pneu_terms_raw.csv")
write.csv2(df_prods, file = "results/Pneu_prod_raw.csv")