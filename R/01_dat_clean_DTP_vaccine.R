###Script to generate a variable for the DTP vaccine as an outcome

#--- reading in the the code lists for DTP vaccine
csv_files <- list.files(paste0(codelists, "/DTP/"))
vacc_codes <- list(rep(data.table(), times = length(csv_files)))

for(i in 1:length(csv_files)){
  csv_files[[i]] <- paste(c(codelists, "/DTP/", csv_files[[i]]), collapse="")
  vacc_codes[[i]]<- as.data.table(read.csv(csv_files[[i]]))
}

#prepare the formatting of the code lists
#first using the medcode based ones
medcode_list <- list (vacc_codes[[1]], vacc_codes[[3]], vacc_codes[[4]])
for(i in 1:length(medcode_list)){
  medcode_list[[i]][, medcodeid := MedCodeId]
  medcode_list[[i]][, MedCodeId := NULL]
  medcode_list[[i]]$medcodeid <- as.character(medcode_list[[i]]$medcodeid)
}
#the product code based one
vacc_codes[[2]] <- vacc_codes[[2]][,prodcodeid := ?..ProdCodeId]
vacc_codes[[2]] <- vacc_codes[[2]][,?..ProdCodeId := NULL]
vacc_codes[[2]]$prodcodeid <- as.character(vacc_codes[[2]]$prodcodeid)


#attaching the right names to the codelists
DTP_codes <- medcode_list
DTP_codes[[4]] <- vacc_codes[[2]]
names(DTP_codes) <- c("DTP_history", "DTP_terms", "DTP_terms_declinded", 
                      "DTP_products")
DTP_codes$DTP_terms <- unique(DTP_codes$DTP_terms)
DTP_codes$DTP_history <- unique(DTP_codes$DTP_history)
DTP_codes$DTP_terms_declinded <- unique(DTP_codes$DTP_terms_declinded)
DTP_codes$DTP_products <- unique(DTP_codes$DTP_products)


#---extracting the different codes
obs.files<- list.files(path = data, pattern = "\\Observation")
prod.files <- list.files(path = data, pattern = "\\DrugIssue")

setwd(data)
df_terms <- extractor_med(list.files = obs.files, codelist = DTP_codes$DTP_terms)
length(unique(df_terms$patid)) #1,701,812 babies with vaccine record
df_terms[given ==1, vac_status := "given"]
df_terms[is.na(given), vac_status := "neutral"]

df_declined <- extractor_med(list.files = obs.files, codelist = DTP_codes$DTP_terms_declinded)
length(unique(df_declined$patid)) #9,396 babies with declined vaccine record
df_declined[no_consent ==1, vac_status := "declined"]
df_declined[DNA ==1, vac_status := "declined"]
df_declined[contraindictaed ==1, vac_status := "declined"]
df_declined[not.immunised ==1, vac_status := "declined"]

df_products <- extractor_prod(list.files = prod.files, codelist = DTP_codes$DTP_products)
df_products <- df_products[, vac_status := "product"]
df_products <- df_products[, obsdate := issuedate]
df_products <- df_products[, issuedate := NULL]
length(unique(df_products$patid)) #26330

df_history <- extractor_med(list.files = obs.files, codelist = DTP_codes$DTP_history)
df_history <- df_history[pertussis == 1]
length(unique(df_history$patid)) #218

#merging declined and term file
df <- merge(df_terms, df_declined, all = T, by = c("patid", "medcodeid", "dob", "consid",
                                                   "pracid", "obsid", "obsdate", "enterdate",
                                                   "staffid", "parentobsid", "value", "numunitid",
                                                   "obstypeid", "numrangelow", "numrangehigh",
                                                   "probobsid", "Term", "pertussis", 
                                                   "vac_status", "start_fu", "end_fu"))
#adding the product codes
colnames(df_products)
df <- merge(df, df_products, all = T, by = c("patid", "dob", "pracid", "enterdate", "obsdate",
                                             "staffid", "pertussis",
                                             "vac_status", "probobsid", "start_fu", "end_fu"))
length(unique(df$patid)) # data of 1,702,695 patients contained in data set


#adding the history/adverse reaction code
colnames(df_history)
df <- merge(df, df_history, all = T, by = c("patid", "medcodeid", "dob", "consid",
                                            "pracid", "obsid", "obsdate", "enterdate",
                                            "staffid", "parentobsid", "value", "numunitid",
                                            "obstypeid", "numrangelow", "numrangehigh",
                                            "probobsid", "Term", "pertussis",
                                            "start_fu", "end_fu"))

rm(df_declined)
rm(df_history)
rm(df_products)
rm(df_terms)

write_parquet(df, paste0(data, "DTP_outcome_uncleaned.parquet"))

###---Vaccine algorithm
#preparation
#only using the pertussis antigen containing vaccines
df <- df[pertussis == 1]
length(unique(df$patid)) #1,700,423 babies included in the cohort


#double checking the follow-up period of the study
colnames(df)
df <- df[obsdate <= end_fu &
           obsdate >= start_fu]

length(unique(df$patid)) #1,696,896 babies remaining

################################################################################
###Step 1: finding out whether a vaccine was given
#flag vaccine as given if there is a product code on a day

tmp <- df[, .(patid, obsdate, enterdate, vac_status, Term, Term.from.EMIS)]
adverse_df <- tmp[is.na(vac_status)] # keep for memory list of adverse effects
length(unique(adverse_df$patid)) # 207
tmp <- tmp[!is.na(vac_status)]
janitor::tabyl(tmp$vac_status)
# tmp$vac_status       n     percent
# declined   10408 0.001643882
# given 2812088 0.444152564
# neutral 3451632 0.545164732
# product   57228 0.009038822

length(unique(tmp$patid)) #1,696,893
length(unique(tmp$obsdate)) #4871

tmp_wide <- dcast(tmp, patid + obsdate + enterdate+ Term + Term.from.EMIS ~ vac_status, fun.aggregate = length,
                  value.var = "vac_status")
tmp_wide[, .(sum_declined = sum(declined),
             sum_given = sum(given), 
             sum_neutral = sum(neutral),
             sum_prod = sum(product))]
# sum_declined sum_given sum_neutral sum_prod
# 1:        10408   2812088     3451632    57228

#summing up the same codes per day 
tmp_filtered <- tmp_wide[, by=c("patid", "obsdate"), lapply(.SD, sum),
                         .SDcols = c("declined", "given", "neutral", "product")]
#default - all codes count as vaccinated
tmp_filtered[, DTP_vacc := "vaccinated"]
#declined when only declined codes or together with neutral
tmp_filtered[neutral > 0 & declined >0, DTP_vacc:= "declined"]
tmp_filtered[declined >0 & product == 0 & given == 0 & neutral == 0,
             DTP_vacc:= "declined"]
#conflict flag when declined and given or declined and product
tmp_filtered[declined >0 & given > 0, DTP_vacc := "conflict"]
tmp_filtered[declined >0 & product > 0, DTP_vacc := "conflict"]

janitor::tabyl(tmp_filtered$DTP_vacc)
# tmp_filtered$DTP_vacc       n      percent
# conflict     138 2.258101e-05
# declined    9831 1.608652e-03
# vaccinated 6101360 9.983688e-01

tmp<- tmp_filtered[DTP_vacc != "conflict"] # drop the conflicted ones 
rm(tmp_filtered)
length(unique(tmp$patid)) #1696888

###Step 2: counting the number of doses
tmp[DTP_vacc == "vaccinated", by = .(patid, DTP_vacc),
    DTP_all_doses :=  .N]
tmp_df <- unique(tmp[, .(patid, DTP_all_doses)])
janitor::tabyl(tmp_df$DTP_all_doses)
# tmp_df$DTP_all_doses       n      percent valid_percent
# 1   40616 2.383851e-02  2.394639e-02
# 2   37979 2.229079e-02  2.239167e-02
# 3  506378 2.972054e-01  2.985505e-01
# 4 1091858 6.408377e-01  6.437379e-01
# 5   17898 1.050477e-02  1.055231e-02
# 6    1066 6.256610e-04  6.284925e-04
# 7     285 1.672734e-04  1.680304e-04
# 8      39 2.289004e-05  2.299363e-05
# 9       <10 NA NA
# NA    7676 4.505229e-03            NA

#and just by events
janitor::tabyl(tmp$DTP_all_doses)
# tmp$DTP_all_doses       n      percent valid_percent
# 1   40616 6.646168e-03  6.656876e-03
# 2   75958 1.242933e-02  1.244936e-02
# 3 1519134 2.485823e-01  2.489828e-01
# 4 4367432 7.146613e-01  7.158129e-01
# 5   89490 1.464363e-02  1.466722e-02
# 6    6396 1.046604e-03  1.048291e-03
# 7    1995 3.264503e-04  3.269763e-04
# 8     312 5.105388e-05  5.113614e-05
# 9      27 4.418124e-06  4.425243e-06
# NA    9831 1.608688e-03            NA

#checking the distribution over years
hist(year(tmp$obsdate))
tmp <- tmp[order(patid, obsdate)]



#counting the number of doses 
tmp[DTP_vacc == "vaccinated", n_dose := seq_len(.N), by = "patid"]
tmp[DTP_vacc == "vaccinated", n_dose := rowid(patid)]

#now dob has be added as information again
dob <- unique(df[, list(patid, dob)])
tmp<- merge(tmp, dob, by = "patid")
tmp[, age := obsdate - dob]
tmp[, age := as.integer(age)] # make age an absolute value in days


###---quality checks
length(unique(tmp$patid)) #1,696,888
### Step 1: Excluding all the children who got their first dose of vaccine
###before they were born
q1 <- tmp[n_dose == 1 & age < 16]
tabyl(q1$DTP_all_doses)
length(unique(q1$patid)) #3249

#exploring backdating as a cause?
q1vac <- tmp_wide[patid %chin% dubious_q1]
q1vac[, n_dose := seq_len(.N), by = "patid"] #approx for qual check 
q1vac[, n_dose := rowid(patid)]
q1vac[n_dose ==1,late_entry := enterdate >= (obsdate+ years(1))]
tabyl(q1vac$late_entry)
# q1vac$late_entry     n    percent valid_percent
# FALSE  2520 0.17083588     0.7756233
# TRUE   729 0.04942038     0.2243767
# NA 11502 0.77974375            NA

#checking their date of birth
population <- read_parquet(paste0(data, "study_population.parquet"))
pop_dub <- population[patid %chin% dubious_q1]
hist(pop_dub$yob) #interesting - decline from 2011
write_parquet(pop_dub, paste0(data, "excluded_DTP_cleaning.parquet"))

rm(population)
rm(pop_dub)
rm(dob)
rm(q1)
rm(q1vac)

#2580 babies had their first vaccine before birth, drop this record
#will be excluded from the final study
tmp1 <- tmp[!patid %chin% dubious_q1]
tabyl(tmp1$DTP_all_doses)
length(unique(tmp1$patid)) #1,693,639

###exploring the data
#at what age was the first DTP dose given
summary(tmp1[n_dose ==1, list(age)])
# Min.   :  16.00  
# 1st Qu.:  55.00  
# Median :  64.00  
# Mean   :  80.38  
# 3rd Qu.:  73.00  
# Max.   :1827.00  

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
    ggtitle("Frequency of first pertussis vaccine dose recorded by age")+
    theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 12, hjust = 0.5))

#at what age was the second DTP dose given
summary(tmp1[n_dose ==2, list(age)])
# age        
# Min.   :  18.0  
# 1st Qu.:  87.0  
# Median :  97.0  
# Mean   : 113.5  
# 3rd Qu.: 111.0  
# Max.   :1827.0  

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
  ggtitle("Frequency of second vaccine dose recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 12, hjust = 0.5))


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
  ggtitle("Frequency of third vaccine dose recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 12, hjust = 0.5))



hist_DTP <- plot_grid(hist_d1, hist_d2, hist_d3, align = "h", labels = c("A", "B", "C"), label_size = 12)

ggsave("hist_pertussis_coverage.pdf",
       hist_DTP,
       width = 6,
       height = 12,
       bg = "white",
       path= graphs)

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
  scale_x_continuous(name="age (days)", limits=c(50, 75), 
                     breaks = seq(from= 50, to=75, by= 1))+
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
  scale_x_continuous(name="age (days)", limits=c(85, 110), 
                     breaks = seq(from= 85, to=110, by= 1))+
  theme_classic()+
  ggtitle("C. Third dose : Number of events excluded by age minium")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))


###make a denistiy plot histogram for first three doses
plot <- tmp1[n_dose <=3]
plot$n_dose <- factor(plot$n_dose)
dens_all <- plot %>%
  ggplot(aes(x=age))+
  geom_density(aes(fill = n_dose), alpha = 0.5) + 
  scale_x_continuous(name="age (days)", limits=c(0, 500), 
                     breaks = seq(from = 0, to = 500, by = 25))+
  theme_classic()+
  ggtitle("Density of the first three vaccine doses recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))+
  scale_fill_manual(name = "vaccine dose", values = c("red", "yellow", "green"))

plot$n_dose <- factor(plot$n_dose)
dens_all_b <- plot %>%
  ggplot(aes(x=age))+
  geom_density(aes(fill = n_dose), alpha = 0.5) + 
  scale_x_continuous(name="age (days)", limits=c(0, 1900), 
                     breaks = seq(from = 0, to = 1900, by = 250))+
  theme_classic()+
  ggtitle("Density of the first three vaccine doses recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))+
  scale_fill_manual(name = "vaccine dose", values = c("red", "yellow", "green"))


###calculating how many records will be dropped within every step
###defining the proportion of vaccine events which align
#with those minimum requirements
tmp1[n_dose ==1 & age < 38, timing := FALSE]
tmp1[n_dose ==1 & age > 38, timing := TRUE]
tabyl(tmp1$timing) # 5770 
tmp1[, timing := NULL]

tmp1[n_dose ==2 & age < 66, timing := FALSE]
tmp1[n_dose ==2 & age > 66, timing := TRUE]
tabyl(tmp1$timing) #3522 
tmp1[, timing := NULL]

tmp1[n_dose ==3 & age < 96, timing := FALSE]
tmp1[n_dose ==3 & age > 96, timing := TRUE]
tabyl(tmp1$timing) #3419 
tmp1[, timing := NULL]

tmp1[n_dose ==4 & age < 1077, timing := FALSE]
tmp1[n_dose ==4 & age > 1077, timing := TRUE]
tabyl(tmp1$timing) #30874 
tmp1[, timing := NULL]






###making the actual time differences between the  doses
#calculating the individual age difference by vaccine dose
tmp2 <- tmp1[DTP_all_doses >=2 & n_dose <=2]
tmp2[n_dose == 1, age1 := age]
tmp2[n_dose == 2, age2 := age]


tmp2[, age_1 := lapply(.SD, max_NA), .SDcols = c("age1"), by="patid"]
tmp2[, age_2 := lapply(.SD, max_NA), .SDcols = c("age2"), by="patid"]
tmp2[, error := age_1 > age_2]
tmp2[, age_diff := age_2-age_1]
tmp2 <- unique(tmp2[,list(patid, age_1, age_2, age_diff)])
summary(tmp2$age_diff)
nrow(tmp2[age_diff <26])/nrow(tmp2) #0.01136152
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   28.00   30.00   43.88   39.00 1777.00 
length(which(tmp2$age_diff <28 & tmp2$age_diff!= 0)) #46,850


#histogram difference 1 and two
plot <- tmp2
hist_diff12 <- plot %>%
  ggplot(aes(x=age_diff))+
  geom_histogram(binwidth =1, color = "orange", fill = "yellow", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age_diff)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age difference (days)", limits=c(10, 45), 
                     breaks = seq(from= 10, to=45, by= 1))+
  theme_classic()+
  ggtitle("A. Age difference between first and second vaccine dose")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))




###repeating the same for second and third dose
#calculating the individual age difference by vaccine dose
tmp2 <- tmp1[DTP_all_doses >=3 & n_dose <=3]
tmp2[n_dose == 2, age1 := age]
tmp2[n_dose == 3, age2 := age]
tmp2[, age_1 := lapply(.SD, max_NA), .SDcols = c("age1"), by="patid"]
tmp2[, age_2 := lapply(.SD, max_NA), .SDcols = c("age2"), by="patid"]
tmp2[, error := age_1 > age_2]
tmp2[, age_diff := age_2-age_1]
tmp2 <- unique(tmp2[,list(patid, age_1, age_2, age_diff)])
summary(tmp2$age_diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   28.00   32.00   54.61   42.00 1745.00 
nrow(tmp2[age_diff <26])
nrow(tmp2[age_diff <26])/nrow(tmp2) #0.01005736

#histogram difference two and three
plot <- tmp2
hist_diff23 <- plot %>%
  ggplot(aes(x=age_diff))+
  geom_histogram(binwidth =1, color = "dark green", fill = "green", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age_diff)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age difference (days)", limits=c(10, 45), 
                     breaks = seq(from= 10, to=45, by= 1))+
  theme_classic()+
  ggtitle("B. Age difference between second and third vaccine dose")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))





###exploring the forth dose
#defining the cut-off age for the first dose
summary(tmp1[n_dose ==4, list(age)])
# Min.   :  60  
# 1st Qu.:1239  
# Median :1273  
# Mean   :1288  
# 3rd Qu.:1328  
# Max.   :1827 


plot <- tmp1[n_dose ==4]
hist_d4 <- plot %>%
  ggplot(aes(x=age))+
  geom_histogram(binwidth =1, color = "dark red", fill = "red", alpha = 0.5) + 
  geom_vline(aes(xintercept = median(age)),
             color = "blue", linetype = "dashed", size = 1)+
  scale_x_continuous(name="age (days)", limits=c(0,1800), 
                     breaks = seq(from= 0, to=1800, by= 100))+
  theme_classic()+
  ggtitle("A. Frequency of fourth vaccine dose recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))

#age difference between third and fourth dose
###repeating the same for second and third dose
#calculating the individual age difference by vaccine dose
tmp2 <- tmp1[DTP_all_doses >=4 & n_dose <=4]
tmp2[n_dose == 3, age1 := age]
tmp2[n_dose == 4, age2 := age]
tmp2[, age_1 := lapply(.SD, max_NA), .SDcols = c("age1"), by="patid"]
tmp2[, age_2 := lapply(.SD, max_NA), .SDcols = c("age2"), by="patid"]
tmp2[, error := age_1 > age_2]
tmp2[, age_diff := age_2-age_1]
tmp2 <- unique(tmp2[,list(patid, age_1, age_2, age_diff)])
summary(tmp2$age_diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    1105    1135    1140    1188    1739 
nrow(tmp2[age_diff <363]) #17614
nrow(tmp2[age_diff <363])/nrow(tmp2) #0.01588997





####final data clean following those thresholds
#first dose given from 38 days
#second dose given from 66 days
#third dose given from 96 days
#fourth dose given from 1,077 days
# age difference between 1&2 ans 2&3 is 26 days
#age diff for 3rd and 4th dose is 363 days
length(unique(tmp1$patid)) #1,693,639
nrow(tmp1) #6,097,603vaccine records


###---DOSE 1
#define if first dose in time, if not remove
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 38)
# [1] "5770 vaccine records dropped (0.09%)"
# [1] "6091833 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 38)
# [1] "24 vaccine records dropped (0%)"
# [1] "6091809 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 1, min_age = 38)
# [1] "0 vaccine records dropped (0%)"
# [1] "6091809 entries remaining in df"


###--- DOSE 2 cleaning
#day cleaning first, cut-off day 66
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 66)
# [1] "1635 vaccine records dropped (0.03%)"
# [1] "6090174 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 66)
# [1] "<10 vaccine records dropped (0%)"
# [1] "NA entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 2, min_age = 66)
# [1] "0 vaccine records dropped (0%)"
# [1] "6090168 entries remaining in df"

#check the time difference now 
tmp1 <- check_min_diff(df = tmp1, dose_a = 1, dose_b = 2, min_diff = 26)
# [1] "63788 vaccine records dropped (1.05%)"
# [1] "6026380 entries remaining in df"
tmp1 <- check_min_diff(df = tmp1, dose_a = 1, dose_b = 2, min_diff = 26)
# [1] "0 vaccine records dropped (0%)"
# [1] "6026380 entries remaining in df"


###--- DOSE 3 cleaning
#date cleaning first
tmp1 <- clean_min_vac_date(df = tmp1, dose = 3, min_age = 96)
# [1] "953 vaccine records dropped (0.02%)"
# [1] "6025427 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 3, min_age = 96)
# [1] "<10 vaccine records dropped (0%)"
# [1] "NA entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 3, min_age = 96)
# [1] "0 vaccine records dropped (0%)"
# [1] "6025425 entries remaining in df"

#check the time difference now
tmp1 <- check_min_diff(df = tmp1, dose_a = 2, dose_b = 3, min_diff = 26)
# [1] "54296 vaccine records dropped (0.9%)"
# [1] "5971129 entries remaining in df"
tmp1 <- check_min_diff(df = tmp1, dose_a = 2, dose_b = 3, min_diff = 26)
# [1] "0 vaccine records dropped (0%)"
# [1] "5971129 entries remaining in df"



###---DOSE 4 cleaning
#cleaning date
tmp1 <- clean_min_vac_date(df = tmp1, dose = 4, min_age = 1077)
# [1] "22033 vaccine records dropped (0.37%)"
# [1] "5949096 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 4, min_age = 1077)
# [1] "639 vaccine records dropped (0.01%)"
# [1] "5948457 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 4, min_age = 1077)
# 1] "127 vaccine records dropped (0%)"
# [1] "5948330 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 4, min_age = 1077)
# [1] "12 vaccine records dropped (0%)"
# [1] "5948318 entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 4, min_age = 1077)
# [1] "<10 vaccine records dropped (0%)"
# [1] "NA entries remaining in df"
tmp1 <- clean_min_vac_date(df = tmp1, dose = 4, min_age = 1077)
# [1] "0 vaccine records dropped (0%)"
# [1] "5948317 entries remaining in df"

#min diff for   363 days
tmp1 <- check_min_diff(df = tmp1, dose_a = 3, dose_b = 4, min_diff = 363)
# [1] "5941 vaccine records dropped (0.1%)"
# [1] "5942376 entries remaining in df"
tmp1 <- check_min_diff(df = tmp1, dose_a = 3, dose_b = 4, min_diff = 363)
# [1] "0 vaccine records dropped (0%)"
# [1] "5942376 entries remaining in df"

tabyl(tmp1$DTP_all_doses)
# tmp1$DTP_all_doses       n      percent valid_percent
# 1   40639 6.838847e-03  6.850171e-03
# 2   79168 1.332262e-02  1.334468e-02
# 3 1528578 2.572335e-01  2.576594e-01
# 4 4264772 7.176880e-01  7.188763e-01
# 5   18270 3.074528e-03  3.079619e-03
# 6    1014 1.706388e-04  1.709214e-04
# 7     112 1.884768e-05  1.887889e-05
# NA    9823 1.653042e-03            NA

summary(tmp1[n_dose ==5, age])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1096    1296    1413    1442    1566    1827 

#compare the fourth and the fifth dose
plot <- tmp1[n_dose == 4 | n_dose == 5 | n_dose == 6 |  n_dose == 7]
plot$n_dose <- factor(plot$n_dose)
dens_lates <- plot %>%
  ggplot(aes(x=age))+
  geom_density(aes(fill = n_dose), alpha = 0.5) + 
  scale_x_continuous(name="age (days)", limits=c(1000, 2000), 
                     breaks = seq(from = 1000, to = 2000, by = 100))+
  theme_classic()+
  ggtitle("Density of vaccine doses 4-7 recorded by age")+
  theme(plot.title = element_text(family = "Helvetica", face = "italic", size = 15, hjust = 0.5))+
  scale_fill_manual(name = "vaccine dose", values = c("red", "yellow", "green", "blue"))


####remove all recordings of a fifth dose and higher from the data set
###before they were born
tmp2 <- tmp1[n_dose <=4 | is.na(n_dose)]

tabyl(tmp2$n_dose)
#tmp2$n_dose       n     percent valid_percent
# 1 1659781 0.279502709     0.2799658
# 2 1619142 0.272659210     0.2731110
# 3 1579558 0.265993369     0.2664341
# 4 1070032 0.180190545     0.1804891
# NA    9823 0.001654167            NA

#how many people in latest data set
length(unique(tmp2$patid)) #1,660,666
tmp2[, vaccine := "DTP"]
tmp2[, vac_status :=  DTP_vacc]
tmp2[, n_all_doses := DTP_all_doses]
tmp2[, DTP_vacc := NULL]
tmp2[, DTP_all_doses := NULL]
tmp2[vac_status == "vaccinated", by = .(patid, DTP_vacc),
     n_all_doses :=  .N]

####---save the cleaned data file
colnames(tmp2)
DTP <- tmp2[, list(patid, dob, obsdate, vaccine, vac_status, n_all_doses, n_dose,  age)]
write_parquet(DTP, paste0(data, "DTP_outcome.parquet"))
