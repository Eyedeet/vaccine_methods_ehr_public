###Script to generate a variable for the MMR vaccine as an outcome

#--- reading in the the code lists for MMR vaccine
csv_files <- list.files(paste0(codelists, "MMR/"))
vacc_codes <- list(rep(data.table(), times = length(csv_files)))

for(i in 1:length(csv_files)){
  csv_files[[i]] <- paste(c(codelists, "MMR/", csv_files[[i]]), collapse="")
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
vacc_codes[[2]] <- vacc_codes[[2]][,prodcodeid := ?..ProdCodeId]
vacc_codes[[2]] <- vacc_codes[[2]][,ProdCodeId := NULL]
vacc_codes[[2]]$prodcodeid <- as.character(vacc_codes[[2]]$prodcodeid)

#attaching the right names to the codelists
MMR_codes <- medcode_list
MMR_codes[[4]] <- vacc_codes[[2]]
names(MMR_codes) <- c("MMR_history", "MMR_terms", "MMR_terms_declinded", 
                       "MMR_products")

MMR_codes$MMR_history <- unique(MMR_codes$MMR_history)
MMR_codes$MMR_terms <- unique(MMR_codes$MMR_terms)
MMR_codes$MMR_terms_declinded <- unique(MMR_codes$MMR_terms_declinded)
MMR_codes$MMR_products <- unique(MMR_codes$MMR_products)

#---extracting the different codes
obs.files<- list.files(path = data, pattern = "\\Observation")
prod.files <- list.files(path = data, pattern = "\\DrugIssue")

setwd(data)
df_terms <- extractor_med(list.files = obs.files, codelist = MMR_codes$MMR_terms)
length(unique(df_terms$patid)) #1,523,541babies with vaccine record
df_terms[given ==1, vac_status := "given"]
df_terms[is.na(given), vac_status := "neutral"]

df_declined <- extractor_med(list.files = obs.files, codelist = MMR_codes$MMR_terms_declinded)
length(unique(df_declined$patid)) # 22,037 babies with declined vaccine record
df_declined[declined ==1, vac_status := "declined"]
df_declined[DNA ==1, vac_status := "declined"]
df_declined[contraindicated ==1, vac_status := "declined"]
df_declined[not_immunised ==1, vac_status := "declined"]

df_products <- extractor_prod(list.files = prod.files, codelist = MMR_codes$MMR_products)
length(unique(df_products$patid)) #19875
df_products <- df_products[, vac_status := "product"]
df_products <- df_products[, obsdate := issuedate]
df_products <- df_products[, issuedate := NULL] #no product codes

#extracting the history and adverse reactions
df_history <- extractor_med(list.files = obs.files, codelist = MMR_codes$MMR_history)
length(unique(df_history$patid))#110 adverse reactions/history



#merging declined and term file
df <- merge(df_terms, df_declined, all = T, by = c("patid", "medcodeid", "dob", "consid",
                                                   "pracid", "obsid", "obsdate", "enterdate",
                                                   "staffid", "parentobsid", "value", "numunitid",
                                                   "obstypeid", "numrangelow", "numrangehigh",
                                                   "probobsid", "Term",  "measles",
                                                   "vac_status", "start_fu", "end_fu"))
#merging with the product terms
colnames(df_products)
df <- merge(df, df_products, all = T, by = c("patid", "dob", "pracid", "enterdate", "obsdate",
                                             "staffid", "measles",
                                             "vac_status", "probobsid", "start_fu", "end_fu"))
#merging with history codes
colnames(df_history)
df_history[, measles := 1]
df <- merge(df, df_history, all = T, by = c("patid", "medcodeid", "dob", "consid",
                                            "pracid", "obsid", "obsdate", "enterdate",
                                            "staffid", "parentobsid", "value", "numunitid",
                                            "obstypeid", "numrangelow", "numrangehigh",
                                            "probobsid", "Term", 
                                            "start_fu", "end_fu", "measles"))

write_parquet(df, paste0(parquet, "filtered/MMR_outcome_uncleaned.parquet"))
df <- read_parquet(paste0(parquet, "filtered/MMR_outcome_uncleaned.parquet"))
#cleaning memory
rm(df_declined)
rm(df_history)
rm(df_products)
rm(df_terms)

###---Vaccine algorithm
#preparation
#only using the pertussis antigen containing vaccines
df <- df[measles == 1]
length(unique(df$patid)) # 1,531,396included in the cohort


#double checking the follow-up period of the study
colnames(df)
df <- df[obsdate <= end_fu &
           obsdate >= start_fu]

length(unique(df$patid)) # 1,521,204 babies remaining

################################################################################
###Step 1: finding out whether a vaccine was given
#flag vaccine as given if there is a product code on a day

tmp <- df[, .(patid, obsdate, enterdate, vac_status, Term, Term.from.EMIS)]
adverse_df <- tmp[is.na(vac_status)] # keep for memory list of adverse effects
length(unique(adverse_df$patid)) # 103
tmp <- tmp[!is.na(vac_status)]
janitor::tabyl(tmp$vac_status)
# tmp$vac_status       n     percent
# declined   23511 0.008297017
# given 1430800 0.504928416
# neutral 1352896 0.477436144
# product   26462 0.009338423

length(unique(tmp$patid)) #1,521,186
length(unique(tmp$obsdate)) #4,627

tmp_wide <- dcast(tmp, patid + obsdate + enterdate+ Term + Term.from.EMIS ~ vac_status, fun.aggregate = length,
                  value.var = "vac_status")
tmp_wide[, .(sum_declined = sum(declined),
             sum_given = sum(given), 
             sum_neutral = sum(neutral),
             sum_prod = sum(product))]

# sum_declined sum_given sum_neutral sum_prod
# 1:        23511   1430800     1352896    26462


#summing up the same codes per day 
tmp_filtered <- tmp_wide[, by=c("patid", "obsdate"), lapply(.SD, sum),
                         .SDcols = c("declined", "given", "neutral", "product")]
#default - all codes count as vaccinated
tmp_filtered[, MMR_vacc := "vaccinated"]
#declined when only declined codes or together with neutral
tmp_filtered[neutral > 0 & declined >0, MMR_vacc:= "declined + neutral"]
tmp_filtered[declined >0 & product == 0 & given == 0 & neutral == 0,
             MMR_vacc:= "declined only"]
#conflict flag when declined and given or declined and product
tmp_filtered[declined >0 & given > 0, MMR_vacc := "declined and given"]
tmp_filtered[declined >0 & product > 0, MMR_vacc := "declined and product"]

#adding for the table
tmp_filtered[declined == 0 & product == 0 & given > 0 & neutral == 0,
             MMR_vacc:= "given only"]
tmp_filtered[declined == 0 & product >0 & given == 0 & neutral == 0,
             MMR_vacc:= "product only"]
tmp_filtered[declined == 0 & product ==0 & given == 0 & neutral > 0,
             MMR_vacc:= "neutal only"]
tmp_filtered[neutral > 0 & given >0, MMR_vacc:= "given + neutral"]
tmp_filtered[neutral > 0 & product >0, MMR_vacc:= "product + neutral"]
tmp_filtered[given > 0 & product >0, MMR_vacc:= "product + given"]

janitor::tabyl(tmp_filtered$MMR_vacc)
# tmp_filtered$MMR_vacc       n      percent
# conflict     337 0.0001237553
# declined   20472 0.0075178610
# vaccinated 2702306 0.9923583837

tmp<- tmp_filtered[MMR_vacc != "conflict"] # drop the conflicted ones 
rm(tmp_filtered)


###Step 2: counting the number of doses
tmp[MMR_vacc == "vaccinated", by = .(patid, MMR_vacc),
    MMR_all_doses :=  .N]
tmp_df <- unique(tmp[, .(patid, MMR_all_doses)])
janitor::tabyl(tmp_df$MMR_all_doses)
# tmp_df$MMR_all_doses       n      percent valid_percent
# 1  347777 2.273097e-01  2.297787e-01
# 2 1144441 7.480153e-01  7.561403e-01
# 3   19770 1.292182e-02  1.306218e-02
# 4    1409 9.209331e-04  9.309363e-04
# 5     109 7.124323e-05  7.201707e-05
# 6      18 1.176494e-05  1.189273e-05
# 7     <10 NA NA
# 8     <10 NA NA
# 11    <10 NA NA
# NA   16440 1.074531e-02            NA
# 


#and just by events
janitor::tabyl(tmp$MMR_all_doses)
# tmp$MMR_all_doses       n      percent valid_percent
# 1  347777 1.277287e-01  1.286964e-01
# 2 2288882 8.406422e-01  8.470107e-01
# 3   59310 2.178290e-02  2.194792e-02
# 4    5636 2.069945e-03  2.085626e-03
# 5     545 2.001632e-04  2.016796e-04
# 6     108 3.966537e-05  3.996587e-05
# 7      21 7.712711e-06  7.771141e-06
# 8      16 5.876351e-06  5.920869e-06
# 11      11 4.039992e-06  4.070597e-06
# NA   20472 7.518791e-03            NA


#checking the distribution over years
hist(year(tmp$obsdate))
tmp <- tmp[order(patid, obsdate)] #increasing towards 2011, then decreasing again

#counting the number of doses 
tmp[MMR_vacc == "vaccinated", n_dose := seq_len(.N), by = "patid"]
tmp[MMR_vacc == "vaccinated", n_dose := rowid(patid)]

#now dob has be added as information again
dob <- unique(df[, list(patid, dob)])
tmp<- merge(tmp, dob, by = "patid")
tmp[, age := obsdate - dob]
tmp[, age := as.integer(age)] # make age an absolute value in days



###---quality checks
length(unique(tmp$patid)) #1,521,107

### Step 1: Excluding all the children who got their first dose of vaccine
###before they were born
q1 <- tmp[n_dose == 1 & age < 16]
tabyl(q1$MMR_all_doses)
length(unique(q1$patid)) #966
dubious_q1 <- q1$patid

#backdating as a cause?
q1vac <- tmp_wide[patid %chin% dubious_q1]
q1vac[, n_dose := seq_len(.N), by = "patid"] #approx for qual check 
q1vac[, n_dose := rowid(patid)]
q1vac[n_dose ==1,late_entry := enterdate >= (obsdate+ years(1))]
tabyl(q1vac$late_entry)
# q1vac$late_entry    n    percent valid_percent
# FALSE  799 0.32178816     0.8271222
# TRUE  167 0.06725735     0.1728778
# NA 1517 0.61095449            NA

#inspecting the records of babies with super early record
q1a <- tmp[patid %chin% dubious_q1]

#checking their date of birth
population <- read_parquet(paste0(data, "study_population.parquet"))
pop_dub <- population[patid %chin% dubious_q1]
hist(pop_dub$yob) ## many recordings for 2006-2010, then drop

population <- read_parquet(paste0(data, "study_population.parquet"))
pop_dub <- population[patid %chin% patids_q1]
write_parquet(pop_dub, paste0(data, "excluded_MMR_clean.parquet"))

rm(population)
rm(pop_dub)
rm(q1vac)
rm(q1)
rm(q1a)


#95 babies had their first vaccine before birth, drop this record
#will be excluded from the final study
dubious_q1 <- c(patids_q1, ex_pat)
tmp1 <- tmp[!patid %chin% dubious_q1]
tabyl(tmp1$MMR_all_doses)
length(unique(tmp1$patid)) #1,517,922

###exploring the data
#at what age was the first MMR dose given
summary(tmp1[n_dose ==1, list(age)])
# Min.   :  16.0  
# 1st Qu.: 387.0  
# Median : 409.0  
# Mean   : 454.5  
# 3rd Qu.: 446.0  
# Max.   :1827.0  


#defining the cut-off age for the first dose
plot <- tmp1[n_dose ==1]
hist_d1 <- plot %>%
  ggplot(aes(x=age))+
  geom_histogram(binwidth =1, color = "darkred", fill = "red", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(0, 1800), 
                     breaks = seq(from= 0, to=1800, by= 100))+
  theme_classic()+
  ggtitle("A. Frequency of first vaccine dose recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))

#at what age was the first MMR dose given
summary(tmp1[n_dose ==2, list(age)])
# age      
# Min.   :  58  
# 1st Qu.:1232  
# Median :1266  
# Mean   :1241  
# 3rd Qu.:1318  
# Max.   :1827  


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


#plotting only the second dose
plot <- tmp1[n_dose ==2]
hist_d2 <- plot %>%
  ggplot(aes(x=age))+
  geom_histogram(binwidth =1, color = "orange", fill = "yellow", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(0,1800), 
                     breaks = seq(from= 0, to=1800, by= 100))+
  theme_classic()+
  ggtitle("B. Frequency of second vaccine dose recorded by age")+
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
  scale_x_continuous(name="age (days)", limits=c(330, 370), 
                     breaks = seq(from= 330, to=370, by= 1))+
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
  scale_x_continuous(name="age (days)", limits=c(50, 75), 
                     breaks = seq(from= 50, to=75, by= 1))+
  theme_classic()+
  ggtitle("B. Second dose : Number of events excluded by age minium")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))

####looking at the difference between the two doses
###making the actual time differences between the  doses
#calculating the individual age difference by vaccine dose
tmp2 <- tmp1[MMR_all_doses >=2 & n_dose <=2]
tmp2[n_dose == 1, age1 := age]
tmp2[n_dose == 2, age2 := age]
tmp2[, age_1 := lapply(.SD, max_NA), .SDcols = c("age1"), by="patid"]
tmp2[, age_2 := lapply(.SD, max_NA), .SDcols = c("age2"), by="patid"]
tmp2[, error := age_1 > age_2]
tmp2[, age_diff := age_2-age_1]
tmp2 <- unique(tmp2[,list(patid, age_1, age_2, age_diff)])
summary(tmp2$age_diff)
nrow(tmp2[age_diff <26])/nrow(tmp2) #0.001790867
nrow(tmp2[age_diff <26]) #2083
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0   805.0   854.0   801.6   900.0  1777.0 

#histogram difference 1 and two
plot <- tmp2
hist_diff12 <- plot %>%
  ggplot(aes(x=age_diff))+
  geom_histogram(binwidth =1, color = "orange", fill = "yellow", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age_diff)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age difference (days)", limits=c(0, 1500), 
                     breaks = seq(from= 0, to=1500, by= 250))+
  theme_classic()+
  ggtitle("Age difference between first and second vaccine dose")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))

theme_classic()



###defining the proportion of vaccine events which align
#with those minimum requirements
tmp1[n_dose ==1 & age < 347, timing := FALSE]
tmp1[n_dose ==1 & age > 347, timing := TRUE]
tabyl(tmp1$timing)
tmp1[, timing := NULL]
# tmp1$timing       n     percent valid_percent
# FALSE   13125 0.004831955   0.008690735
# TRUE 1497104 0.551157239   0.991309265
# NA 1206063 0.444010806            NA


tmp1[n_dose ==2 & age < 530, timing := FALSE]
tmp1[n_dose ==2 & age > 530, timing := TRUE]
tabyl(tmp1$timing)
tmp1[, timing := NULL]
# FALSE   45981 0.01692786    0.03954158
# TRUE 1116871 0.41117487    0.96045842
# NA 1553440 0.57189728            NA



####final data clean following those thresholds
#first dose given from 347 days
#second dose given from 530 days
# age difference between 1&2 is 26 days
length(unique(tmp1$patid)) #1,517,922
nrow(tmp1) # 2,716,292 vaccine records

#
####cleaning the data set
#min age dose one is 347 days
#min age for second dose is 530 days
#time difference is 26 days between two doses
###---DOSE 1
#define if first dose in time, if not remove
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 347)
# [1] "13125 vaccine records dropped (0.48%)"
# [1] "2703167 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 347)
# [1] "347 vaccine records dropped (0.01%)"
# [1] "2702820 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 347)
# [1] "126 vaccine records dropped (0%)"
# [1] "2702694 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 347)
# [1] "<10 vaccine records dropped (0%)"
# [1] "NA entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 347)
# [1] "0 vaccine records dropped (0%)"
# [1] "2702688 entries remaining in df"

#minimum age for the second dose
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 530)
# [1] "39517 vaccine records dropped (1.46%)"
# [1] "2663171 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 530)
# [1] "265 vaccine records dropped (0.01%)"
# [1] "2662906 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 530)
# [1] "<10 vaccine records dropped (0%)"
# [1] "NA entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 530)
# [1] "<10 vaccine records dropped (0%)"
# [1] "NA entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 530)
# [1] "0 vaccine records dropped (0%)"
# [1] "2662897 entries remaining in df"
#check the time difference now 
tmp1 <- check_min_diff(df = tmp1, dose_a = 1, dose_b = 2, min_diff = 26)
# [1] "912 vaccine records dropped (0.03%)"
# [1] "2661985 entries remaining in df"
tmp1 <- check_min_diff(df = tmp1, dose_a = 1, dose_b = 2, min_diff = 26)
# [1] "0 vaccine records dropped (0%)"
# [1] "2661985 entries remaining in df"

tabyl(tmp1$MMR_all_doses)
# tmp1$MMR_all_doses       n      percent valid_percent
# 1  385002 1.446297e-01  1.457492e-01
# 2 2226574 8.364337e-01  8.429085e-01
# 3   27528 1.034116e-02  1.042121e-02
# 4    2360 8.865565e-04  8.934192e-04
# 5      50 1.878298e-05  1.892837e-05
# 6      12 4.507914e-06  4.542810e-06
# 11      11 4.132255e-06  4.164242e-06
# NA   20448 7.681486e-03            NA

#compare when all doses were given
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot <- tmp1
plot$n_dose <- factor(plot$n_dose)
dens_lates <- plot %>%
  ggplot(aes(x=age))+
  geom_density(aes(fill = n_dose), alpha = 0.5) + 
  scale_x_continuous(name="age (days)", limits=c(0, 2000), 
                     breaks = seq(from = 0, to = 2000, by = 250))+
  theme_classic()+
  ggtitle("Density of all MMR vaccine doses recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))+
  scale_fill_manual(name = "vaccine dose", values = cbbPalette[1:6])

colnames(tmp1)
tmp <- tmp1[n_dose >2]
tabyl(tmp, n_dose, MMR_vacc)

####remove all recordings of a fith dose and higher from the data set
###before they were born
tmp2 <- tmp1[n_dose <=2 | is.na(n_dose)]
tabyl(tmp2$n_dose)
# tmp2$n_dose       n     percent valid_percent
# 1 1508068 0.568742735     0.5731628
# 2 1123066 0.423545642     0.4268372
# NA   20448 0.007711623            NA
tmp2[, vaccine := "MMR"]
tmp2[, vac_status :=  MMR_vacc]
tmp2[, n_all_doses := MMR_all_doses]
tmp2[, MMR_vacc := NULL]
tmp2[, MMR_all_doses := NULL]
tmp2[vac_status == "vaccinated", by = .(patid, MMR_vacc),
     n_all_doses :=  .N]



####---save the cleaned data file
colnames(tmp2)
MMR <- tmp2[, list(patid, dob, obsdate, vaccine, vac_status, n_all_doses, n_dose, age)]
write_parquet(MMR, paste0(data, "MMR_outcome.parquet"))
