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



###---looking up the codes

terms1 <- c("DT", "DTP", "diphteria", "tetanus", "pertussis",
            "diphth", "dipht", "trivax", "boostrix", "DTAP",
            "whooping", "infanrix", "tetavax")
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
snomed <- c("34771611000001100", "632481000119106", "770608009",
            "1082461000000110", "770616000", "412763007",
            "4689411000001100", "12218401000001100", "414004005",
            "10776411000001100",  "4695411000001100", "170429001",
            "170414003", "281040003",  "335051000000104",
            "924721000006106" ,"809002", "842801000000100")

concepts <- concept_lookup(terms = snomed)

#checking mismatch between already identified terms
missed <- setdiff(concepts, df_terms) # none are missing


###safing the preliminary lists
write.csv2(df_terms, file = "results/DTP_terms_raw.csv")
write.csv2(df_prods, file = "results/DTP_prod_raw.csv")
 