###Script to generate a variable for the Pneumococcal vaccine as an outcome

#--- reading in the the code lists for Pneu vaccine
csv_files <- list.files(paste0(codelists, "Pneu/"))
vacc_codes <- list(rep(data.table(), times = length(csv_files)))

for(i in 1:length(csv_files)){
  csv_files[[i]] <- paste(c(codelists, "Pneu/", csv_files[[i]]), collapse="")
  vacc_codes[[i]]<- as.data.table(read.csv(csv_files[[i]]))
}


#prepare the formatting of the code lists
#first using the medcode based ones
medcode_list <- list (vacc_codes[[1]], vacc_codes[[3]], vacc_codes[[4]])
for(i in 1:length(medcode_list)){
  medcode_list[[i]][, medcodeid := MedCodeId]
  medcode_list[[i]][, MedCodeId := NULL]
  medcode_list[[i]]$medcodeid <- as.character( medcode_list[[i]]$medcodeid)  
}
#the product code based one
vacc_codes[[2]] <- vacc_codes[[2]][,prodcodeid := ProdCodeId]
vacc_codes[[2]] <- vacc_codes[[2]][, ProdCodeId := NULL]
vacc_codes[[2]]$prodcodeid <- as.character(vacc_codes[[2]]$prodcodeid)

#attaching the right names to the codelists
Pneu_codes <- medcode_list
Pneu_codes[[4]] <- vacc_codes[[2]]
names(Pneu_codes) <- c("Pneu_history", "Pneu_terms", "Pneu_terms_declinded", 
                      "Pneu_products")

#---extracting the different codes
obs.files<- list.files(path = data, pattern = "\\Observation")
prod.files <- list.files(path = data, pattern = "\\DrugIssue")

setwd(data)
df_terms <- extractor_med(list.files = obs.files, codelist = Pneu_codes$Pneu_terms)
length(unique(df_terms$patid)) #1,663,291 babies with vaccine record
df_terms[given ==1, vac_status := "given"]
df_terms[is.na(given), vac_status := "neutral"]

df_declined <- extractor_med(list.files = obs.files, codelist = Pneu_codes$Pneu_terms_declinded)
length(unique(df_declined$patid)) # 5835 babies with declined vaccine record
df_declined[no_consent ==1, vac_status := "declined"]
df_declined[DNA ==1, vac_status := "declined"]
df_declined[contraindictaed ==1, vac_status := "declined"]
df_declined[not.immunised ==1, vac_status := "declined"]

df_products <- extractor_prod(list.files = prod.files, codelist = Pneu_codes$Pneu_products)
df_products <- df_products[, vac_status := "product"]
df_products <- df_products[, obsdate := issuedate]
df_products <- df_products[, issuedate := NULL] 
length(unique(df_products$patid))#21,574 product codes

df_history <- extractor_med(list.files = obs.files, codelist = Pneu_codes$Pneu_history)
length(unique(df_history$patid)) #99


#merging declined and term file
df <- merge(df_terms, df_declined, all = T, by = c("patid", "medcodeid", "dob", "consid",
                                                   "pracid", "obsid", "obsdate", "enterdate",
                                                   "staffid", "parentobsid", "value", "numunitid",
                                                   "obstypeid", "numrangelow", "numrangehigh",
                                                   "probobsid", "Term", 
                                                   "vac_status", "start_fu", "end_fu"))
#merging with product
df <- merge(df, df_products, all = T, by = c("patid", "dob", "pracid", "enterdate", "obsdate",
                                             "staffid", 
                                             "vac_status", "probobsid", "start_fu", "end_fu", "PPV"))

#merging with history codes
df <- merge(df, df_history, all = T, by = c("patid", "medcodeid", "dob", "consid",
                                            "pracid", "obsid", "obsdate", "enterdate",
                                            "staffid", "parentobsid", "value", "numunitid",
                                            "obstypeid", "numrangelow", "numrangehigh",
                                            "probobsid", "Term", 
                                            "start_fu", "end_fu", "PPV"))

write_parquet(df, paste0(data, "Pneu_outcome_uncleaned.parquet"))

rm(df_declined)
rm(df_history)
rm(df_products)
rm(df_terms)


#clean out PPV given only to older children
colnames(df)
length(which(df$PPV == 1)) #2,130 records of PPV vaccine to be dropped
df <- df[is.na(PPV)]
nrow(df) #4,735,843


###---Vaccine algorithm
#preparation
length(unique(df$patid)) #1,664,737 babies included in the cohort

#double checking the follow-up period of the study
df <- df[obsdate <= end_fu &
           obsdate >= start_fu]

length(unique(df$patid))  #1,663,639

################################################################################
###Step 1: finding out whether a vaccine was given
#flag vaccine as given if there is a product code on a day

tmp <- df[, .(patid, obsdate, enterdate, vac_status, Term)]
adverse_df <- tmp[is.na(vac_status)] # keep for memory list of adverse effects
length(unique(adverse_df$patid)) # 55
tmp <- tmp[!is.na(vac_status)]
janitor::tabyl(tmp$vac_status)
# tmp$vac_status       n     percent
# declined    6347 0.001343196
# given    8885 0.001880305
# neutral 4669056 0.988098100
# product   41008 0.008678398



length(unique(tmp$patid)) #1,663,636
length(unique(tmp$obsdate)) #4,281


tmp_wide <- dcast(tmp, patid + obsdate + enterdate+ Term  ~ vac_status, fun.aggregate = length,
                  value.var = "vac_status")
tmp_wide[, .(sum_declined = sum(declined),
             sum_given = sum(given), 
             sum_neutral = sum(neutral),
             sum_product = sum(product))]
# sum_declined sum_given sum_neutral sum_product
# 1:         6347      8885     4669056       41008

#summing up the same codes per day 
tmp_filtered <- tmp_wide[, by=c("patid", "obsdate"), lapply(.SD, sum),
                         .SDcols = c("declined", "given", "neutral", "product")]


#default - all codes count as vaccinated
tmp_filtered[, Pneu_vacc := "vaccinated"]
#declined when only declined codes or together with neutral
tmp_filtered[neutral > 0 & declined >0, Pneu_vacc:= "declined"]
tmp_filtered[declined >0 & product == 0 & given == 0 & neutral == 0,
             Pneu_vacc:= "declined"]
#conflict flag when declined and given or declined and product
tmp_filtered[declined >0 & given > 0, Pneu_vacc := "conflict"]
tmp_filtered[declined >0 & product > 0, Pneu_vacc := "conflict"]


janitor::tabyl(tmp_filtered$Pneu_vacc)
# tmp_filtered$Pneu_vacc       n      percent
# conflict       <10 6.532821e-07
# declined    5815 1.266278e-03
# vaccinated 4586379 9.987331e-01

tmp<- tmp_filtered[Pneu_vacc != "conflict"] # drop the conflicted ones 3
rm(tmp_filtered)

###Step 2: counting the number of doses
tmp[Pneu_vacc == "vaccinated", by = .(patid, Pneu_vacc),
    Pneu_all_doses :=  .N]
tmp_df <- unique(tmp[, .(patid, Pneu_all_doses)])
janitor::tabyl(tmp_df$Pneu_all_doses)
# tmp_df$Pneu_all_doses       n      percent valid_percent
# 1   82031 4.919713e-02  4.934771e-02
# 2  254629 1.527108e-01  1.531782e-01
# 3 1308522 7.847707e-01  7.871728e-01
# 4   16207 9.719958e-03  9.749709e-03
# 5     814 4.881870e-04  4.896812e-04
# 6      96 5.757487e-05  5.775110e-05
# 7     <10 NA  NA
# 8     <10 NA NA
# NA    5088 3.051468e-03            NA

#and just by events
janitor::tabyl(tmp$Pneu_all_doses)
# tmp$Pneu_all_doses       n      percent valid_percent
# 1   82031 1.786314e-02  1.788579e-02
# 2  509258 1.108964e-01  1.110371e-01
# 3 3925566 8.548345e-01  8.559184e-01
# 4   64828 1.411700e-02  1.413490e-02
# 5    4070 8.862866e-04  8.874103e-04
# 6     576 1.254302e-04  1.255893e-04
# 7      42 9.145955e-06  9.157551e-06
# 8       <10 NA  NA
# NA    5815 1.266279e-03            NA

#checking the distribution over years
hist(year(tmp$obsdate)) #increase from 2007 on, stable and declined from 2015, 2016
tmp <- tmp[order(patid, obsdate)]


#counting the number of doses 
tmp[Pneu_vacc == "vaccinated", n_dose := seq_len(.N), by = "patid"]
tmp[Pneu_vacc == "vaccinated", n_dose := rowid(patid)]

#now dob has be added as information again
dob <- unique(df[, list(patid, dob)])
tmp<- merge(tmp, dob, by = "patid")
tmp[, age := obsdate - dob]
tmp[, age := as.integer(age)] # make age an absolute value in days

###---quality checks
length(unique(tmp$patid)) #1,663,636
### Step 1: Excluding all the children who got their first dose of vaccine
###before they were born
q1 <- tmp[n_dose == 1 & age < 16]
patids_q1 <- unique(q1$patid)
length(unique(q1$patid))#929


#backdating as a cause?
q1vac <- tmp_wide[patid %chin% patids_q1]
q1vac[, n_dose := seq_len(.N), by = "patid"] #approx for qual check 
q1vac[, n_dose := rowid(patid)]
q1vac[n_dose ==1,late_entry := enterdate >= (obsdate+ years(1))]
tabyl(q1vac$late_entry)
# q1vac$late_entry   n    percent valid_percent
# FALSE 100 0.09182736      0.295858
# TRUE 238 0.21854913      0.704142
# NA 751 0.68962351            NA


#checking their date of birth
population <- read_parquet(paste0(data, "study_population.parquet"))
pop_dub <- population[patid %chin% patids_q1]
length(unique(pop_dub$patid))/length(unique(tmp$patid)) #0.01 % dropped
length(unique(pop_dub$patid)) #290 dropped
hist(pop_dub$yob) #less dropped born later 2013
write_parquet(pop_dub, paste0(data, "excluded_Pneu_clean.parquet"))

rm(population)
rm(pop_dub)
rm(q1)


#290 babies had their first vaccine before birth, 
#those will be dropped after the ones from DTP already dropped
tmp1 <- tmp_e[!patid %chin% patids_q1]
tabyl(tmp1$Pneu_all_doses)
length(unique(tmp1$patid)) #1,660,417




#at what age was the first Pneu dose given
summary(tmp1[n_dose ==1, list(age)])
# age        
# Min.   :  16.0  
# 1st Qu.:  56.0  
# Median :  65.0  
# Mean   :  76.9  
# 3rd Qu.:  74.0  
# Max.   :1825.0  


#defining the cut-off age for the first dose
plot <- tmp1[n_dose ==1]
hist_d1 <- plot %>%
  ggplot(aes(x=age))+
  geom_histogram(binwidth =1, color = "darkred", fill = "red", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(0, 500), 
                     breaks = seq(from= 0, to=500, by= 25))+
  theme_classic()+
  ggtitle("A. Frequency of first vaccine dose recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))

#at what age was the second Pneumo dose given
summary(tmp1[n_dose ==2, list(age)])
# age        
# Min.   :  23.0  
# 1st Qu.: 119.0  
# Median : 131.0  
# Mean   : 152.6  
# 3rd Qu.: 151.0  
# Max.   :1825.0  

#plotting only the second dose
plot <- tmp1[n_dose ==2]
hist_d2 <- plot %>%
  ggplot(aes(x=age))+
  geom_histogram(binwidth =1, color = "orange", fill = "yellow", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(0,500), 
                     breaks = seq(from= 0, to=500, by= 25))+
  theme_classic()+
  ggtitle("B. Frequency of second vaccine dose recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))


#at what age was the second Pneumo dose given
summary(tmp1[n_dose ==3, list(age)])
# age        
# Min.   :  34.0  
# 1st Qu.: 386.0  
# Median : 406.0  
# Mean   : 427.8  
# 3rd Qu.: 438.0  
# Max.   :1826.0 

#plotting the third dose
plot <- tmp1[n_dose ==3]
hist_d3 <-  plot %>%
  ggplot(aes(x=age))+
  geom_histogram(binwidth =1, color = "darkgreen", fill = "green", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(0,500), 
                     breaks = seq(from= 0, to=500, by= 25))+
  theme_classic()+
  ggtitle("C. Frequency of third vaccine dose recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))


###make the cumulative histograms for the cut-off ages for each of the three
#first doses 
#dose 1
plot <- tmp1[n_dose ==1]
histcum_d1 <- plot %>%
  ggplot(aes(x=age))+
  geom_histogram(aes(y = cumsum(..count..)),
                 binwidth =1, color = "darkred", fill = "red", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(25, 45), 
                     breaks = seq(from= 25, to=45, by= 1))+
  theme_classic()+
  ggtitle("A. First dose : Number of events excluded by age minium")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))

#dose 2
plot <- tmp1[n_dose ==2]
histcum_d2 <- plot %>%
  ggplot(aes(x=age))+
  geom_histogram(aes(y = cumsum(..count..)),
                 binwidth =1, color = "orange", fill = "yellow", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(85,110), 
                     breaks = seq(from= 85, to=110, by= 1))+
  theme_classic()+
  ggtitle("B. Second dose : Number of events excluded by age minium")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))

#dose 3
plot <- tmp1[n_dose ==3]
histcum_d3 <- plot %>%
  ggplot(aes(x=age))+
  geom_histogram(aes(y = cumsum(..count..)),
                 binwidth =1, color = "darkgreen", fill = "green", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(335, 360), 
                     breaks = seq(from= 335, to=360, by= 1))+
  theme_classic()+
  ggtitle("C. Third dose : Number of events excluded by age minium")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))



####now looking at the minimum  age difference of 26 days between each dose
###exploring the age difference between first and second dose
tmp3 <- tmp1[Pneu_all_doses >= 2 & (n_dose == 1 | n_dose == 2)]
tmp3[n_dose == 1, age2 := age]
tmp3[n_dose == 2, age3 := age]
tmp3[, age_2 := lapply(.SD, max_NA), .SDcols = c("age2"), by="patid"]
tmp3[, age_3 := lapply(.SD, max_NA), .SDcols = c("age3"), by="patid"]
tmp3[, error := age_2 > age_3]
tmp3[, age_diff := age_3-age_2]
tmp3 <- unique(tmp3[,list(patid, age_2, age_3, age_diff)])
summary(tmp3$age_diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   56.00   63.00   83.45   79.00 1753.00 

##graph of age difference between dose 1 and two
#histogram difference 1 and two
plot <- tmp3
hist_diff12 <- plot %>%
  ggplot(aes(x=age_diff))+
  geom_histogram(binwidth =1, color = "orange", fill = "yellow", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age_diff)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age difference (days)", limits=c(0, 150), 
                     breaks = seq(from= 0, to=150, by= 25))+
  theme_classic()+
  ggtitle("A. Age difference between first and second vaccine dose")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))




###exploring the age difference between first and second dose
tmp3 <- tmp1[Pneu_all_doses >= 3 & (n_dose == 2 | n_dose == 3)]
tmp3[n_dose == 1, age2 := age]
tmp3[n_dose == 2, age3 := age]
tmp3[, age_2 := lapply(.SD, max_NA), .SDcols = c("age2"), by="patid"]
tmp3[, age_3 := lapply(.SD, max_NA), .SDcols = c("age3"), by="patid"]
tmp3[, error := age_2 > age_3]
tmp3[, age_diff := age_3-age_2]
tmp3 <- unique(tmp3[,list(patid, age_2, age_3, age_diff)])
summary(tmp3$age_diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 23.0   119.0   130.0   142.8   148.0  1784.0 

##graph of age difference between dose 2 and 3
#histogram difference one and two
plot <- tmp3
hist_diff23 <- plot %>%
  ggplot(aes(x=age_diff))+
  geom_histogram(binwidth =1, color = "dark green", fill = "green", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age_diff)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age difference (days)", limits=c(0, 150), 
                     breaks = seq(from= 0, to=150, by= 25))+
  theme_classic()+
  ggtitle("B. Age difference between second and third vaccine dose")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))



###calculating how many records will be dropped within every step
###defining the proportion of vaccine events which align
#with those minimum requirements
tmp1[n_dose ==1 & age < 38, timing := FALSE]
tmp1[n_dose ==1 & age > 38, timing := TRUE]
tabyl(tmp1$timing)
# tmp1$timing       n      percent valid_percent
# FALSE    2648 0.0005777601   0.001596494
# TRUE 1655987 0.3613154254   0.998403506
# NA 2924582 0.6381068145            NA
tmp1[, timing := NULL]

tmp1[n_dose ==2 & age < 96, timing := FALSE]
tmp1[n_dose ==2 & age > 96, timing := TRUE]
tabyl(tmp1$timing)
# tmp1$timing       n     percent valid_percent
# FALSE    7974 0.001739826   0.005057244
# TRUE 1568774 0.342286651   0.994942756
# NA 3006469 0.655973523            NA
tmp1[, timing := NULL]

tmp1[n_dose ==3 & age < 347, timing := FALSE]
tmp1[n_dose ==3 & age > 347, timing := TRUE]
tabyl(tmp1$timing)
# tmp1$timing       n     percent valid_percent
# FALSE   21431 0.004675973    0.01619771
# TRUE 1301657 0.284005099    0.98380229
# NA 3260129 0.711318927            NA
tmp1[, timing := NULL]




####final data clean following those thresholds
#first dose given from 38 days
#second dose given from 96 days
#third dose given from 347 days
# age difference between 1&2 ans 2&3 is 26 days

length(unique(tmp1$patid)) #1,660,417
nrow(tmp1) #4,583,217 vaccine records

###---DOSE 1
#define if first dose in time, if not remove
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 38)
# [1] "2648 vaccine records dropped (0.06%)"
# [1] "4580569 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 38)
# [1] "<10 vaccine records dropped (0%)"
# [1] "NA entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 38)
# [1] "<10 vaccine records dropped (0%)"
# [1] "NA entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 38)
# [1] "0 vaccine records dropped (0%)"
# [1] "4579861 entries remaining in df"


###--- DOSE 2 cleaning
#day cleaning first, cut-off day 96
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 96)
# [1] "7074 vaccine records dropped (0.15%)"
# [1] "4573487 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 96)
# [1] "39 vaccine records dropped (0%)"
# [1] "4573448 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 96)
# [1] "0 vaccine records dropped (0%)"
# [1] "4573448 entries remaining in df"

#check the time difference now 
tmp1 <- check_min_diff(df = tmp1, dose_a = 1, dose_b = 2, min_diff = 26)
# [1] "2499 vaccine records dropped (0.05%)"
# [1] "4570949 entries remaining in df"
tmp1 <- check_min_diff(df = tmp1, dose_a = 1, dose_b = 2, min_diff = 26)
# [1] "0 vaccine records dropped (0%)"
# [1] "4570949 entries remaining in df"


###--- DOSE 3 cleaning
#day cleaning first, cut-off day 347
tmp1 <- clean_min_vac_date(df = tmp1, dose = 3, min_age = 347)
# [1] "15951 vaccine records dropped (0.35%)"
# [1] "4554998 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 3, min_age = 347)
# [1] "217 vaccine records dropped (0%)"
# [1] "4554781 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 3, min_age = 347)
# [1] "<10 vaccine records dropped (0%)"
# [1] "4554778 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 3, min_age = 347)
# [1] "0 vaccine records dropped (0%)"
# [1] "4554778 entries remaining in df"

#time difference between 2nd and 3rd dose
#check the time difference now 
tmp1 <- check_min_diff(df = tmp1, dose_a = 2, dose_b = 3, min_diff = 26)
# [1] "1553 vaccine records dropped (0.03%)"
# [1] "4553225 entries remaining in df"
tmp1 <- check_min_diff(df = tmp1, dose_a = 2, dose_b = 3, min_diff = 26)
# [1] "0 vaccine records dropped (0%)"
# [1] "4553225 entries remaining in df"

tmp_id <- unique(tmp1[, list(patid, Pneu_all_doses)])
tabyl(tmp_id$Pneu_all_doses)
# tmp_id$Pneu_all_doses       n      percent valid_percent
# 1   82522 4.963108e-02  4.978306e-02
# 2  267661 1.609790e-01  1.614719e-01
# 3 1300518 7.821686e-01  7.845638e-01
# 4    6666 4.009122e-03  4.021399e-03
# 5     231 1.389300e-04  1.393554e-04
# 6      32 1.924571e-05  1.930465e-05
# 7       2 1.202857e-06  1.206540e-06
# NA    5076 3.052851e-03            NA
tabyl(tmp1$Pneu_all_doses)
# tmp1$Pneu_all_doses       n      percent valid_percent
# 1   82522 1.812386e-02  1.814698e-02
# 2  535322 1.175699e-01  1.177199e-01
# 3 3901554 8.568770e-01  8.579703e-01
# 4   26664 5.856069e-03  5.863541e-03
# 5    1155 2.536664e-04  2.539900e-04
# 6     192 4.216791e-05  4.222172e-05
# 7      14 3.074744e-06  3.078667e-06
# NA    5802 1.274262e-03            NA
tmpc <- tmp1[n_dose <= 3 | is.na(n_dose)]
tabyl(tmpc$n_dose)
tmpc[, vaccine := "Pneu"]
tmpc[, vac_status :=  Pneu_vacc]
tmpc [, n_all_doses := Pneu_all_doses]
tmpc[, Pneu_vacc := NULL]
tmpc[, Pneu_all_doses := NULL]
tmp2[vac_status == "vaccinated", by = .(patid, Pneu_vacc),
     n_all_doses :=  .N]


#saving the final file
tmpc<- tmpc[, list(patid, dob, obsdate, vaccine, vac_status, n_all_doses, n_dose, age)]
write_parquet(tmpc, paste0(data, "Pneu_outcome.parquet"))
