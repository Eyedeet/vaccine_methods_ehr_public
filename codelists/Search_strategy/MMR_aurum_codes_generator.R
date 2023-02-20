library(tidyverse)
library(stringr)
library(dplyr)
library(data.table)

#reading in the data dictionaries for CPRD Aurum release May 2022
med_aurum<- read.delim(".../CPRD Aurum/Code browsers/2022_05/CPRDAurumMedical.txt", 
                       sep= "\t",  colClasses="character")
prod_aurum <- read.delim("...CPRD Aurum/Code browsers/2022_05/CPRDAurumProduct.txt", 
                         sep= "\t",  colClasses="character")




###---filtering addintional terms in the aurum data set which might have been missed
#function for looking up the terms
term_lookup <- function(terms){
  #result list for all the searched terms
  list_terms <- list()
  
  for(i in 1:length(terms)){
    
    index <- grep(terms[i], med_aurum$Term, ignore.case=T)
    list_terms[[i]]<- med_aurum[index,]
  }
  
  #generate one data frame deduplicating entries
  df <- list_terms[[1]]
  for(i in 2:length(list_terms)){
    df <- unique(rbind(df, list_terms[[i]]))
  }
  
  df<-df%>%
    arrange(OriginalReadCode)
  return(df)
}


###---code to find the right products
prod_lookup <- function(terms){
  #result list for all the searched terms
  list_terms <- list()
  
  for(i in 1:length(terms)){
    
    index_1 <- grep(terms[i], prod_aurum$Term.from.EMIS, ignore.case=T)
    index_2 <- grep(terms[i], prod_aurum$DrugSubstanceName, ignore.case=T)
    index <- c(index_1, index_2)
    list_terms[[i]]<- prod_aurum[index,]
    
  }
  
  #generate one data frame deduplicating entries
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


###---looking up the codes

terms1 <- c("MMR", "measles", "rubella", "mumps", "Immrav", "priorix", "Mevilin",
            "M-M-R", "M-M-Rvaxpro" , "pluserix", "meruvax", "rubavax", "ervevax",
            "almevax", "attenuvax")
terms2 <- c("vac", "imm", "consent", "injec", "declin", "invit", "boost", 
            "not given", "call", "message", "dose", "attenuated")

#vaccine terms
df1 <- term_lookup(terms = terms1)
df2 <- term_lookup(terms = terms2)

df_terms <- intersect(df1, df2)

#vaccine product codes
df1a <- prod_lookup(terms = terms1)
df2a <- prod_lookup(terms = terms2)

df_prods <- intersect(df1a, df2a)

###related snomed concept IDs from pervious code list screening
snomed <- c("310578008", "12866006", "10100410000006110", "16660311000000100",
            "714821000000108", "432636005", "1128901000001103", 
            "13968211000001108",
             "14015211000001108", "2144301000001109",
             "3063401000001105", "3403401000001100",
             "34925111000001104",              "34938511000001103", 
             "1036101000001108",            "1036201000001101", 
             "1129001000001108",             "3403301000001108", 
             "4621611000001106",             "4639811000001109", 
             "4830211000001107",             "652601000001108", 
             "9927501000001103",             "9711000175104", 
             "505031000000103",            "150971000119104", 
             "308081000000105",            "432636005", 
             "571591000119106",             "38598009", 
             "505001000000109",              "572511000119105", 
             "871909005",              "150971000119104", 
             "170433008",              "433733003", 
             "505001000000109",          "572511000119105", 
             "170431005",            "170432003", 
             "38598009",           "432636005", 
             "571591000119106")
concepts_2 <- concept_lookup(terms = snomed_2)
missed_2 <- setdiff(concepts_2, df_terms)


concepts <- concept_lookup(terms = snomed)

#checking mismatch between already identified terms
missed <- setdiff(concepts, df_terms)  #no codes missed

###safing the preliminary lists
write.csv2(df_terms, file = "results/MMR_terms_raw.csv")
write.csv2(df_prods, file = "results/MMR_prod_raw.csv")
