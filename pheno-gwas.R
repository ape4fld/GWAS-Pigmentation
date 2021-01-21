################################  11/2019
################################  SCRIPT TO PREPARE PHENOTYPES FOR GWAS
################################

setwd("")

# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("here")
# install.packages("DescTools")
# install.packages("ggpubr")
# install.packages("rlist")
install.packages("brant")

library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(DescTools)
library(ggpubr)
library(rlist)
library(MASS)
library(brant)

here()

pheno <- read.table(here("phenotype-wrangling", "raw-data.csv"), sep = ",", header = TRUE)   # read raw phenotype data
pheno2 <- pheno[,c(1,8,9,11,12)]   # select columns to keep 
colnames(pheno2) <- c("IID","sex", "age", "hair_color", "eye_color")   # change column names

data_ids <- list.files(path = "phenotype-wrangling", pattern = "*.final.ids.fam")   # list plink fam files
IDS = lapply(data_ids, function(x) read.table(here("phenotype-wrangling/", x), header = FALSE) %>%
               .[,2])
names(IDS) <- c("IID","IID","IID","IID","IID")

# PCs for pop structure:
data_pcs <- list.files(path = "phenotype-wrangling", pattern = "*.for-pop-str-control.eigenvec")
PCS = lapply(data_pcs, function(x) read.table(here("phenotype-wrangling", x), header = TRUE) %>%
               .[,c(2,3:10)])

# convert IDs into separate dataframes:
cptp_ids <- as.data.frame(IDS[1])
gsa_4224_ids <- as.data.frame(IDS[2])
gsa_5300_ids <- as.data.frame(IDS[3])
gsa_760_ids <- as.data.frame(IDS[4])
omni_ids <- as.data.frame(IDS[5])

# convert PCs into separate dataframes:
cptp_pcs <- as.data.frame(PCS[1])
gsa_4224_pcs <- as.data.frame(PCS[2])
gsa_5300_pcs <- as.data.frame(PCS[3])
gsa_760_pcs <- as.data.frame(PCS[4])
omni_pcs <- as.data.frame(PCS[5])

# separate phenotypes by chip IDs:
pheno_cptp_1 <- left_join(cptp_ids, pheno2, by = "IID")
pheno_gsa4224_1 <- left_join(gsa_4224_ids, pheno2, by = "IID")
pheno_gsa5300_1 <- left_join(gsa_5300_ids, pheno2, by = "IID")
pheno_gsa760_1 <- left_join(gsa_760_ids, pheno2, by = "IID")
pheno_omni_1 <- left_join(omni_ids, pheno2, by = "IID")

# include PCs on each chip dataframe:
pheno_cptp <- left_join(pheno_cptp_1, cptp_pcs, by = "IID")
pheno_gsa4224 <- left_join(pheno_gsa4224_1, gsa_4224_pcs, by = "IID")
pheno_gsa5300 <- left_join(pheno_gsa5300_1, gsa_5300_pcs, by = "IID")
pheno_gsa760 <- left_join(pheno_gsa760_1, gsa_760_pcs, by = "IID")
pheno_omni <- left_join(pheno_omni_1, omni_pcs, by = "IID")

################################
################################ Binary phenotypes:
################################

####### HAIR:
####### 1. blonde vs. light/dark brown + black
####### 2. red vs. light/dark brown + black
####### 3. light/dark brown vs. black

pheno_cptp$hair_color <- as.character(pheno_cptp$hair_color)
pheno_cptp$eye_color <- as.character(pheno_cptp$eye_color)

# M1 (hair color)
pheno_cptp$hair_m1[pheno_cptp$hair_color == "black"] <- 0
pheno_cptp$hair_m1[pheno_cptp$hair_color == "light_brown"] <- 0
pheno_cptp$hair_m1[pheno_cptp$hair_color == "dark_brown"] <- 0
pheno_cptp$hair_m1[pheno_cptp$hair_color == "blonde"] <- 1
pheno_cptp$hair_m1[pheno_cptp$hair_color == "red"] <- NA

# M2 (hair color)
pheno_cptp$hair_m2[pheno_cptp$hair_color == "black"] <- 0
pheno_cptp$hair_m2[pheno_cptp$hair_color == "light_brown"] <- 0
pheno_cptp$hair_m2[pheno_cptp$hair_color == "dark_brown"] <- 0
pheno_cptp$hair_m2[pheno_cptp$hair_color == "blonde"] <- NA
pheno_cptp$hair_m2[pheno_cptp$hair_color == "red"] <- 1

# M3 (hair color)
pheno_cptp$hair_m3[pheno_cptp$hair_color == "black"] <- 0
pheno_cptp$hair_m3[pheno_cptp$hair_color == "light_brown"] <- 1
pheno_cptp$hair_m3[pheno_cptp$hair_color == "dark_brown"] <- 1
pheno_cptp$hair_m3[pheno_cptp$hair_color == "blonde"] <- NA
pheno_cptp$hair_m3[pheno_cptp$hair_color == "red"] <- NA

####### EYES
####### 1. blue + grey vs. brown
####### 2. blue + grey vs. green
####### 3. blue + grey vs. hazel
####### 4. green vs. brown
####### 5. hazel vs. brown

# M1 (eye color): blue + grey vs. brown
pheno_cptp$eye_m1[pheno_cptp$eye_color == "blue"] <- 1
pheno_cptp$eye_m1[pheno_cptp$eye_color == "grey"] <- 1
pheno_cptp$eye_m1[pheno_cptp$eye_color == "brown"] <- 0
pheno_cptp$eye_m1[pheno_cptp$eye_color == "green"] <- NA
pheno_cptp$eye_m1[pheno_cptp$eye_color == "amber"] <- NA
pheno_cptp$eye_m1[pheno_cptp$eye_color == "hazel"] <- NA

# M2 (eye color): blue + grey vs. green
pheno_cptp$eye_m2[pheno_cptp$eye_color == "blue"] <- 1
pheno_cptp$eye_m2[pheno_cptp$eye_color == "grey"] <- 1
pheno_cptp$eye_m2[pheno_cptp$eye_color == "brown"] <- NA
pheno_cptp$eye_m2[pheno_cptp$eye_color == "green"] <- 0
pheno_cptp$eye_m2[pheno_cptp$eye_color == "amber"] <- NA
pheno_cptp$eye_m2[pheno_cptp$eye_color == "hazel"] <- NA

# M3 (eye color): blue + grey vs. hazel
pheno_cptp$eye_m3[pheno_cptp$eye_color == "blue"] <- 1
pheno_cptp$eye_m3[pheno_cptp$eye_color == "grey"] <- 1
pheno_cptp$eye_m3[pheno_cptp$eye_color == "brown"] <- NA
pheno_cptp$eye_m3[pheno_cptp$eye_color == "green"] <- NA
pheno_cptp$eye_m3[pheno_cptp$eye_color == "amber"] <- NA
pheno_cptp$eye_m3[pheno_cptp$eye_color == "hazel"] <- 0

# M4 (eye color): green vs. brown
pheno_cptp$eye_m4[pheno_cptp$eye_color == "blue"] <- NA
pheno_cptp$eye_m4[pheno_cptp$eye_color == "grey"] <- NA
pheno_cptp$eye_m4[pheno_cptp$eye_color == "brown"] <- 0
pheno_cptp$eye_m4[pheno_cptp$eye_color == "green"] <- 1
pheno_cptp$eye_m4[pheno_cptp$eye_color == "amber"] <- NA
pheno_cptp$eye_m4[pheno_cptp$eye_color == "hazel"] <- NA

# M5 (eye color): hazel vs. brown
pheno_cptp$eye_m5[pheno_cptp$eye_color == "blue"] <- NA
pheno_cptp$eye_m5[pheno_cptp$eye_color == "grey"] <- NA
pheno_cptp$eye_m5[pheno_cptp$eye_color == "brown"] <- 0
pheno_cptp$eye_m5[pheno_cptp$eye_color == "green"] <- NA
pheno_cptp$eye_m5[pheno_cptp$eye_color == "amber"] <- NA
pheno_cptp$eye_m5[pheno_cptp$eye_color == "hazel"] <- 1

# GSA_4224:
pheno_gsa4224$hair_color <- as.character(pheno_gsa4224$hair_color)
# pheno_gsa4224 <- pheno_gsa4224[,c(1:4)]

# M1 (hair color)
pheno_gsa4224$hair_m1[pheno_gsa4224$hair_color == "black"] <- 0
pheno_gsa4224$hair_m1[pheno_gsa4224$hair_color == "light_brown"] <- 0
pheno_gsa4224$hair_m1[pheno_gsa4224$hair_color == "dark_brown"] <- 0
pheno_gsa4224$hair_m1[pheno_gsa4224$hair_color == "blonde"] <- 1
pheno_gsa4224$hair_m1[pheno_gsa4224$hair_color == "red"] <- NA

# M2 (hair color)
pheno_gsa4224$hair_m2[pheno_gsa4224$hair_color == "black"] <- 0
pheno_gsa4224$hair_m2[pheno_gsa4224$hair_color == "light_brown"] <- 0
pheno_gsa4224$hair_m2[pheno_gsa4224$hair_color == "dark_brown"] <- 0
pheno_gsa4224$hair_m2[pheno_gsa4224$hair_color == "blonde"] <- NA
pheno_gsa4224$hair_m2[pheno_gsa4224$hair_color == "red"] <- 1

# M3 (hair color)
pheno_gsa4224$hair_m3[pheno_gsa4224$hair_color == "black"] <- 0
pheno_gsa4224$hair_m3[pheno_gsa4224$hair_color == "light_brown"] <- 1
pheno_gsa4224$hair_m3[pheno_gsa4224$hair_color == "dark_brown"] <- 1
pheno_gsa4224$hair_m3[pheno_gsa4224$hair_color == "blonde"] <- NA
pheno_gsa4224$hair_m3[pheno_gsa4224$hair_color == "red"] <- NA

# GSA5300:
pheno_gsa5300$hair_color <- as.character(pheno_gsa5300$hair_color)
pheno_gsa5300$eye_color <- as.character(pheno_gsa5300$eye_color)

# M1 (hair color)
pheno_gsa5300$hair_m1[pheno_gsa5300$hair_color == "black"] <- 0
pheno_gsa5300$hair_m1[pheno_gsa5300$hair_color == "light_brown"] <- 0
pheno_gsa5300$hair_m1[pheno_gsa5300$hair_color == "dark_brown"] <- 0
pheno_gsa5300$hair_m1[pheno_gsa5300$hair_color == "blonde"] <- 1
pheno_gsa5300$hair_m1[pheno_gsa5300$hair_color == "red"] <- NA
 
# M2 (hair color)
pheno_gsa5300$hair_m2[pheno_gsa5300$hair_color == "black"] <- 0
pheno_gsa5300$hair_m2[pheno_gsa5300$hair_color == "light_brown"] <- 0
pheno_gsa5300$hair_m2[pheno_gsa5300$hair_color == "dark_brown"] <- 0
pheno_gsa5300$hair_m2[pheno_gsa5300$hair_color == "blonde"] <- NA
pheno_gsa5300$hair_m2[pheno_gsa5300$hair_color == "red"] <- 1

# M3 (hair color)
pheno_gsa5300$hair_m3[pheno_gsa5300$hair_color == "black"] <- 0
pheno_gsa5300$hair_m3[pheno_gsa5300$hair_color == "light_brown"] <- 1
pheno_gsa5300$hair_m3[pheno_gsa5300$hair_color == "dark_brown"] <- 1
pheno_gsa5300$hair_m3[pheno_gsa5300$hair_color == "blonde"] <- NA
pheno_gsa5300$hair_m3[pheno_gsa5300$hair_color == "red"] <- NA

########### EYES
####### EYES
####### 1. blue + grey vs. brown
####### 2. blue + grey vs. green
####### 3. blue + grey vs. hazel
####### 4. green vs. brown
####### 5. hazel vs. brown

# M1 (eye color): blue + grey vs. brown
pheno_gsa5300$eye_m1[pheno_gsa5300$eye_color == "blue"] <- 1
pheno_gsa5300$eye_m1[pheno_gsa5300$eye_color == "grey"] <- 1
pheno_gsa5300$eye_m1[pheno_gsa5300$eye_color == "brown"] <- 0
pheno_gsa5300$eye_m1[pheno_gsa5300$eye_color == "green"] <- NA
pheno_gsa5300$eye_m1[pheno_gsa5300$eye_color == "amber"] <- NA
pheno_gsa5300$eye_m1[pheno_gsa5300$eye_color == "hazel"] <- NA

# M2 (eye color): blue + grey vs. green
pheno_gsa5300$eye_m2[pheno_gsa5300$eye_color == "blue"] <- 1
pheno_gsa5300$eye_m2[pheno_gsa5300$eye_color == "grey"] <- 1
pheno_gsa5300$eye_m2[pheno_gsa5300$eye_color == "brown"] <- NA
pheno_gsa5300$eye_m2[pheno_gsa5300$eye_color == "green"] <- 0
pheno_gsa5300$eye_m2[pheno_gsa5300$eye_color == "amber"] <- NA
pheno_gsa5300$eye_m2[pheno_gsa5300$eye_color == "hazel"] <- NA

# M3 (eye color): blue + grey vs. hazel
pheno_gsa5300$eye_m3[pheno_gsa5300$eye_color == "blue"] <- 1
pheno_gsa5300$eye_m3[pheno_gsa5300$eye_color == "grey"] <- 1
pheno_gsa5300$eye_m3[pheno_gsa5300$eye_color == "brown"] <- NA
pheno_gsa5300$eye_m3[pheno_gsa5300$eye_color == "green"] <- NA
pheno_gsa5300$eye_m3[pheno_gsa5300$eye_color == "amber"] <- NA
pheno_gsa5300$eye_m3[pheno_gsa5300$eye_color == "hazel"] <- 0

# M4 (eye color): green vs. brown
pheno_gsa5300$eye_m4[pheno_gsa5300$eye_color == "blue"] <- NA
pheno_gsa5300$eye_m4[pheno_gsa5300$eye_color == "grey"] <- NA
pheno_gsa5300$eye_m4[pheno_gsa5300$eye_color == "brown"] <- 0
pheno_gsa5300$eye_m4[pheno_gsa5300$eye_color == "green"] <- 1
pheno_gsa5300$eye_m4[pheno_gsa5300$eye_color == "amber"] <- NA
pheno_gsa5300$eye_m4[pheno_gsa5300$eye_color == "hazel"] <- NA

# M5 (eye color): hazel vs. brown
pheno_gsa5300$eye_m5[pheno_gsa5300$eye_color == "blue"] <- NA
pheno_gsa5300$eye_m5[pheno_gsa5300$eye_color == "grey"] <- NA
pheno_gsa5300$eye_m5[pheno_gsa5300$eye_color == "brown"] <- 0
pheno_gsa5300$eye_m5[pheno_gsa5300$eye_color == "green"] <- NA
pheno_gsa5300$eye_m5[pheno_gsa5300$eye_color == "amber"] <- NA
pheno_gsa5300$eye_m5[pheno_gsa5300$eye_color == "hazel"] <- 1


# GSA760:
pheno_gsa760$hair_color <- as.character(pheno_gsa760$hair_color)
# pheno_gsa760 <- pheno_gsa760[,c(1:4)]

# M1 (hair color)
pheno_gsa760$hair_m1[pheno_gsa760$hair_color == "black"] <- 0
pheno_gsa760$hair_m1[pheno_gsa760$hair_color == "light_brown"] <- 0
pheno_gsa760$hair_m1[pheno_gsa760$hair_color == "dark_brown"] <- 0
pheno_gsa760$hair_m1[pheno_gsa760$hair_color == "blonde"] <- 1
pheno_gsa760$hair_m1[pheno_gsa760$hair_color == "red"] <- NA

# M2 (hair color)
pheno_gsa760$hair_m2[pheno_gsa760$hair_color == "black"] <- 0
pheno_gsa760$hair_m2[pheno_gsa760$hair_color == "light_brown"] <- 0
pheno_gsa760$hair_m2[pheno_gsa760$hair_color == "dark_brown"] <- 0
pheno_gsa760$hair_m2[pheno_gsa760$hair_color == "blonde"] <- NA
pheno_gsa760$hair_m2[pheno_gsa760$hair_color == "red"] <- 1

# M3 (hair color)
pheno_gsa760$hair_m3[pheno_gsa760$hair_color == "black"] <- 0
pheno_gsa760$hair_m3[pheno_gsa760$hair_color == "light_brown"] <- 1
pheno_gsa760$hair_m3[pheno_gsa760$hair_color == "dark_brown"] <- 1
pheno_gsa760$hair_m3[pheno_gsa760$hair_color == "blonde"] <- NA
pheno_gsa760$hair_m3[pheno_gsa760$hair_color == "red"] <- NA

# OMNI:
pheno_omni$hair_color <- as.character(pheno_omni$hair_color)
# pheno_omni <- pheno_omni[,c(1:4)]

# M1 (hair color)
pheno_omni$hair_m1[pheno_omni$hair_color == "black"] <- 0
pheno_omni$hair_m1[pheno_omni$hair_color == "light_brown"] <- 0
pheno_omni$hair_m1[pheno_omni$hair_color == "dark_brown"] <- 0
pheno_omni$hair_m1[pheno_omni$hair_color == "blonde"] <- 1
pheno_omni$hair_m1[pheno_omni$hair_color == "red"] <- NA

# M2 (hair color)
pheno_omni$hair_m2[pheno_omni$hair_color == "black"] <- 0
pheno_omni$hair_m2[pheno_omni$hair_color == "light_brown"] <- 0
pheno_omni$hair_m2[pheno_omni$hair_color == "dark_brown"] <- 0
pheno_omni$hair_m2[pheno_omni$hair_color == "blonde"] <- NA
pheno_omni$hair_m2[pheno_omni$hair_color == "red"] <- 1

# M3 (hair color)
pheno_omni$hair_m3[pheno_omni$hair_color == "black"] <- 0
pheno_omni$hair_m3[pheno_omni$hair_color == "light_brown"] <- 1
pheno_omni$hair_m3[pheno_omni$hair_color == "dark_brown"] <- 1
pheno_omni$hair_m3[pheno_omni$hair_color == "blonde"] <- NA
pheno_omni$hair_m3[pheno_omni$hair_color == "red"] <- NA


############### Prepare tables for SNPTEST:

# Prepare output tables:
pheno_bin_cptp <- data_frame(ID_1 = pheno_cptp$IID, ID_2 = pheno_cptp$IID, missing = rep(0, nrow(pheno_cptp)),
                               sex = pheno_cptp$sex, age = pheno_cptp$age, PC1 = pheno_cptp$PC1, PC2 = pheno_cptp$PC2, PC3 = pheno_cptp$PC3,
                             PC4 = pheno_cptp$PC4, PC5 = pheno_cptp$PC5, PC6 = pheno_cptp$PC6, PC7 = pheno_cptp$PC7, PC8 = pheno_cptp$PC8,
                             hair_m1 = pheno_cptp$hair_m1, hair_m2 = pheno_cptp$hair_m2,
                             hair_m3 = pheno_cptp$hair_m3, eye_m1 = pheno_cptp$eye_m1, eye_m2 = pheno_cptp$eye_m2, eye_m3 = pheno_cptp$eye_m3, 
                             eye_m4 = pheno_cptp$eye_m4, eye_m5 = pheno_cptp$eye_m5)

pheno_bin_gsa5300 <- data_frame(ID_1 = pheno_gsa5300$IID, ID_2 = pheno_gsa5300$IID, missing = rep(0, nrow(pheno_gsa5300)),
                                sex = pheno_gsa5300$sex, age = pheno_gsa5300$age, PC1 = pheno_gsa5300$PC1, PC2 = pheno_gsa5300$PC2, PC3 = pheno_gsa5300$PC3,
                                PC4 = pheno_gsa5300$PC4, PC5 = pheno_gsa5300$PC5, PC6 = pheno_gsa5300$PC6, PC7 = pheno_gsa5300$PC7, PC8 = pheno_gsa5300$PC8,
                                hair_m1 = pheno_gsa5300$hair_m1, hair_m2 = pheno_gsa5300$hair_m2,
                                hair_m3 = pheno_gsa5300$hair_m3, eye_m1 = pheno_gsa5300$eye_m1, eye_m2 = pheno_gsa5300$eye_m2, eye_m3 = pheno_gsa5300$eye_m3,
                                eye_m4 = pheno_gsa5300$eye_m4, eye_m5 = pheno_gsa5300$eye_m5)
  
pheno_bin_gsa4224 <- data_frame(ID_1 = pheno_gsa4224$IID, ID_2 = pheno_gsa4224$IID, missing = rep(0, nrow(pheno_gsa4224)),
                                sex = pheno_gsa4224$sex, age = pheno_gsa4224$age, PC1 = pheno_gsa4224$PC1, PC2 = pheno_gsa4224$PC2, PC3 = pheno_gsa4224$PC3,
                                PC4 = pheno_gsa4224$PC4, PC5 = pheno_gsa4224$PC5,
                                hair_m1 = pheno_gsa4224$hair_m1, hair_m2 = pheno_gsa4224$hair_m2,
                                hair_m3 = pheno_gsa4224$hair_m3)
  
pheno_bin_gsa760 <- data_frame(ID_1 = pheno_gsa760$IID, ID_2 = pheno_gsa760$IID, missing = rep(0, nrow(pheno_gsa760)),
                              sex = pheno_gsa760$sex, age = pheno_gsa760$age, PC1 = pheno_gsa760$PC1, PC2 = pheno_gsa760$PC2,
                              hair_m1 = pheno_gsa760$hair_m1, hair_m2 = pheno_gsa760$hair_m2,
                              hair_m3 = pheno_gsa760$hair_m3)
  
pheno_bin_omni <- data_frame(ID_1 = pheno_omni$IID, ID_2 = pheno_omni$IID, missing = rep(0, nrow(pheno_omni)),
                             sex = pheno_omni$sex, age = pheno_omni$age, PC1 = pheno_omni$PC1, PC2 = pheno_omni$PC2,
                             hair_m1 = pheno_omni$hair_m1, hair_m2 = pheno_omni$hair_m2,
                             hair_m3 = pheno_omni$hair_m3)


# write tables:

# write.table(pheno_bin_cptp, here("phenotype-wrangling/", "pheno_bin_cptp.sample"), sep = " ", row.names = FALSE, quote = FALSE)
# write.table(pheno_bin_gsa5300, here("phenotype-wrangling/", "pheno_bin_gsa5300.sample"), sep = " ", row.names = FALSE, quote = FALSE)
# write.table(pheno_bin_gsa4224, here("phenotype-wrangling/", "pheno_bin_gsa4224.sample"), sep = " ", row.names = FALSE, quote = FALSE)
# write.table(pheno_bin_gsa760, here("phenotype-wrangling/", "pheno_bin_gsa760.sample"), sep = " ", row.names = FALSE, quote = FALSE)
# write.table(pheno_bin_omni, here("phenotype-wrangling/", "pheno_bin_omni.sample"), sep = " ", row.names = FALSE, quote = FALSE)

# Need to add the second row into each of the sample files after writing files.

###
###
############################################# ORDINARY PHENOTYPES:
###
###

pheno_cptp_or <- pheno_cptp[,1:13]
pheno_gsa5300_or <- pheno_gsa5300[,1:13]
pheno_gsa4224_or <- pheno_gsa4224[,1:13]
pheno_gsa760_or <- pheno_gsa760[,1:13]
pheno_omni_or <- pheno_omni[,1:13]

pheno_cptp_or$hair_color <- as.character(pheno_cptp_or$hair_color)
pheno_cptp_or$eye_color <- as.character(pheno_cptp_or$eye_color)

# Model_ordinal (hair color)
pheno_cptp_or$hair_ord[pheno_cptp_or$hair_color == "black"] <- 3
pheno_cptp_or$hair_ord[pheno_cptp_or$hair_color == "light_brown"] <- 1
pheno_cptp_or$hair_ord[pheno_cptp_or$hair_color == "dark_brown"] <- 2
pheno_cptp_or$hair_ord[pheno_cptp_or$hair_color == "blonde"] <- 0
pheno_cptp_or$hair_ord[pheno_cptp_or$hair_color == "red"] <- NA

# Model_ordinal (eye color)
pheno_cptp_or$eye_ord[pheno_cptp_or$eye_color == "grey"] <- 0
pheno_cptp_or$eye_ord[pheno_cptp_or$eye_color == "blue"] <- 0
pheno_cptp_or$eye_ord[pheno_cptp_or$eye_color == "amber"] <- 2
pheno_cptp_or$eye_ord[pheno_cptp_or$eye_color == "green"] <- 1
pheno_cptp_or$eye_ord[pheno_cptp_or$eye_color == "hazel"] <- 2
pheno_cptp_or$eye_ord[pheno_cptp_or$eye_color == "brown"] <- 3

# GSA5300:
pheno_gsa5300_or$hair_color <- as.character(pheno_gsa5300_or$hair_color)
pheno_gsa5300_or$eye_color <- as.character(pheno_gsa5300_or$eye_color)

# Model_ordinal (hair color)
pheno_gsa5300_or$hair_ord[pheno_gsa5300_or$hair_color == "black"] <- 3
pheno_gsa5300_or$hair_ord[pheno_gsa5300_or$hair_color == "light_brown"] <- 1
pheno_gsa5300_or$hair_ord[pheno_gsa5300_or$hair_color == "dark_brown"] <- 2
pheno_gsa5300_or$hair_ord[pheno_gsa5300_or$hair_color == "blonde"] <- 0
pheno_gsa5300_or$hair_ord[pheno_gsa5300_or$hair_color == "red"] <- NA

# Model_ordinal (eye color)
pheno_gsa5300_or$eye_ord[pheno_gsa5300_or$eye_color == "grey"] <- 0
pheno_gsa5300_or$eye_ord[pheno_gsa5300_or$eye_color == "blue"] <- 0
pheno_gsa5300_or$eye_ord[pheno_gsa5300_or$eye_color == "amber"] <- 2
pheno_gsa5300_or$eye_ord[pheno_gsa5300_or$eye_color == "green"] <- 1
pheno_gsa5300_or$eye_ord[pheno_gsa5300_or$eye_color == "hazel"] <- 2
pheno_gsa5300_or$eye_ord[pheno_gsa5300_or$eye_color == "brown"] <- 3

# GSA4224:

pheno_gsa4224_or$hair_color <- as.character(pheno_gsa4224_or$hair_color)
pheno_gsa4224_or$eye_color <- as.character(pheno_gsa4224_or$eye_color)

# Model_ordinal (hair color)
pheno_gsa4224_or$hair_ord[pheno_gsa4224_or$hair_color == "black"] <- 3
pheno_gsa4224_or$hair_ord[pheno_gsa4224_or$hair_color == "light_brown"] <- 1
pheno_gsa4224_or$hair_ord[pheno_gsa4224_or$hair_color == "dark_brown"] <- 2
pheno_gsa4224_or$hair_ord[pheno_gsa4224_or$hair_color == "blonde"] <- 0
pheno_gsa4224_or$hair_ord[pheno_gsa4224_or$hair_color == "red"] <- NA

# GSA760:

pheno_gsa760_or$hair_color <- as.character(pheno_gsa760_or$hair_color)
pheno_gsa760_or$eye_color <- as.character(pheno_gsa760_or$eye_color)

# Model_ordinal (hair color)
pheno_gsa760_or$hair_ord[pheno_gsa760_or$hair_color == "black"] <- 3
pheno_gsa760_or$hair_ord[pheno_gsa760_or$hair_color == "light_brown"] <- 1
pheno_gsa760_or$hair_ord[pheno_gsa760_or$hair_color == "dark_brown"] <- 2
pheno_gsa760_or$hair_ord[pheno_gsa760_or$hair_color == "blonde"] <- 0
pheno_gsa760_or$hair_ord[pheno_gsa760_or$hair_color == "red"] <- NA

# OMNI:

pheno_omni_or$hair_color <- as.character(pheno_omni_or$hair_color)
pheno_omni_or$eye_color <- as.character(pheno_omni_or$eye_color)

# Model_ordinal (hair color)
pheno_omni_or$hair_ord[pheno_omni_or$hair_color == "black"] <- 3
pheno_omni_or$hair_ord[pheno_omni_or$hair_color == "light_brown"] <- 1
pheno_omni_or$hair_ord[pheno_omni_or$hair_color == "dark_brown"] <- 2
pheno_omni_or$hair_ord[pheno_omni_or$hair_color == "blonde"] <- 0
pheno_omni_or$hair_ord[pheno_omni_or$hair_color == "red"] <- NA

################## test for proportional odds:

pheno_cptp_or$hair_ord <- as.factor(pheno_cptp_or$hair_ord)
pheno_cptp_or$eye_ord <- as.factor(pheno_cptp_or$eye_ord)
cptp_model <- polr(hair_ord ~ sex + age, data = pheno_cptp_or, method = "probit")
brant(cptp_model)

pheno_gsa5300_or$hair_ord <- as.factor(pheno_gsa5300_or$hair_ord)
pheno_gsa5300_or$eye_ord <- as.factor(pheno_gsa5300_or$eye_ord)
gsa5300_model <- polr(hair_ord ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + sex + age, data = pheno_gsa5300_or, method = "probit")
brant(gsa5300_model)

pheno_gsa4224_or$hair_ord <- as.factor(pheno_gsa4224_or$hair_ord)
gsa4224_model <- polr(hair_ord ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + age, data = pheno_gsa4224_or, method = "probit")
brant(gsa4224_model)

# plink:
# The file format is the same as for --pheno (optional header line, FID and IID in first two columns, 
# covariates in remaining columns). By default, the main phenotype is set to missing if any covariate
# is missing; you can disable this with the 'keep-pheno-on-missing-cov' modifier.

# write --pheno file:
pheno_ord_cptp <- data.frame(IID = pheno_cptp_or$IID, hair_ord = pheno_cptp_or$hair_ord, eye_ord = pheno_cptp_or$eye_ord)
pheno_ord_gsa5300 <- data.frame(IID = pheno_gsa5300_or$IID, hair_ord = pheno_gsa5300_or$hair_ord, eye_ord = pheno_gsa5300_or$eye_ord)
pheno_ord_gsa4224 <- data.frame(IID = pheno_gsa4224_or$IID, hair_ord = pheno_gsa4224_or$hair_ord)
pheno_ord_gsa760 <- data.frame(IID = pheno_gsa760_or$IID, hair_ord = pheno_gsa760_or$hair_ord)
pheno_ord_omni <- data.frame(IID = pheno_omni_or$IID, hair_ord = pheno_omni_or$hair_ord)

# write tables with covariates:
cov_ord_cptp <- data.frame(IID = pheno_cptp_or$IID, sex = pheno_cptp_or$sex, age = pheno_cptp_or$age,
                           PC1 = pheno_cptp_or$PC1, PC2 = pheno_cptp_or$PC2, PC3 = pheno_cptp_or$PC3, PC4 = pheno_cptp_or$PC4,
                           PC5 = pheno_cptp_or$PC5, PC6 = pheno_cptp_or$PC6, PC7 = pheno_cptp_or$PC7, PC8 = pheno_cptp_or$PC8)

cov_ord_gsa5300 <- data.frame(IID = pheno_gsa5300_or$IID, sex = pheno_gsa5300_or$sex, age = pheno_gsa5300_or$age,
                              PC1 = pheno_gsa5300_or$PC1, PC2 = pheno_gsa5300_or$PC2, PC3 = pheno_gsa5300_or$PC3, PC4 = pheno_gsa5300_or$PC4,
                              PC5 = pheno_gsa5300_or$PC5, PC6 = pheno_gsa5300_or$PC6, PC7 = pheno_gsa5300_or$PC7, PC8 = pheno_gsa5300_or$PC8)

cov_ord_gsa4224 <- data.frame(IID = pheno_gsa4224_or$IID, sex = pheno_gsa4224_or$sex, age = pheno_gsa4224_or$age,
                              PC1 = pheno_gsa4224_or$PC1, PC2 = pheno_gsa4224_or$PC2, PC3 = pheno_gsa4224_or$PC3, PC4 = pheno_gsa4224_or$PC4,
                              PC5 = pheno_gsa4224_or$PC5)

cov_ord_gsa760 <- data.frame(IID = pheno_gsa760_or$IID, sex = pheno_gsa760_or$sex, age = pheno_gsa760_or$age,
                             PC1 = pheno_gsa760_or$PC1, PC2 = pheno_gsa760_or$PC2)

cov_ord_omni <- data.frame(IID = pheno_omni_or$IID, sex = pheno_omni_or$sex, age = pheno_omni_or$age,
                           PC1 = pheno_omni_or$PC1, PC2 = pheno_omni_or$PC2)

# write tables:

 write.table(pheno_ord_cptp, here("phenotype-wrangling/plink-ord", "pheno_ord_cptp.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 write.table(pheno_ord_gsa5300, here("phenotype-wrangling/plink-ord", "pheno_ord_gsa5300.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 write.table(pheno_ord_gsa4224, here("phenotype-wrangling/plink-ord", "pheno_ord_gsa4224.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 write.table(pheno_ord_gsa760, here("phenotype-wrangling/plink-ord", "pheno_ord_gsa760.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 write.table(pheno_ord_omni, here("phenotype-wrangling/plink-ord", "pheno_ord_omni.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

 write.table(cov_ord_cptp, here("phenotype-wrangling/plink-ord", "cov_ord_cptp.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 write.table(cov_ord_gsa5300, here("phenotype-wrangling/plink-ord", "cov_ord_gsa5300.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 write.table(cov_ord_gsa4224, here("phenotype-wrangling/plink-ord", "cov_ord_gsa4224.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 write.table(cov_ord_gsa760, here("phenotype-wrangling/plink-ord", "cov_ord_gsa760.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 write.table(cov_ord_omni, here("phenotype-wrangling/plink-ord", "cov_ord_omni.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
