library(tidyverse)
library(here)

# read CSV files from GWAS catalogue:
melanoma <- read.table(here("/CPTP-data/skin-cancer-GWAS-catalogue","melanoma_significant_short.csv"),
                       sep=",", header = TRUE, fill = TRUE) %>%
  rename(rsid = ï..Variant)

unique(melanoma$rsid) %>% length(.) -> nsnps
thres_melanoma <- 0.05/nsnps

bsc <- read.table(here("/CPTP-data/skin-cancer-GWAS-catalogue","BSC_significant_short.csv"),
                 sep=",", header = TRUE, fill = TRUE) %>%
  rename(rsid = ï..Variant)

unique(bsc$rsid) %>% length(.) -> nsnps
thres_bsc <- 0.05/nsnps

cscc <- read.table(here("/CPTP-data/skin-cancer-GWAS-catalogue","CSCC_significant_short.csv"),
                   sep=",", header = TRUE, fill = TRUE) %>%
  rename(rsid = ï..Variant)

unique(cscc$rsid) %>% length(.) -> nsnps
thres_cscc <- 0.05/nsnps

suntan <- read.table(here("/CPTP-data/skin-cancer-GWAS-catalogue","sun-tan_significant_short.csv"),
                     sep=",", header = TRUE, fill = TRUE) %>%
  rename(rsid = ï..Variant)

unique(suntan$rsid) %>% length(.) -> nsnps
thres_suntan <- 0.05/nsnps

# read GWAS sumstats and filter by pvalue < 0.01
HC_blonde <- read.table(here("/CPTP-data/Hair-colour-paper","sumstats-HC-blonde.out"), 
                        sep=" ", header=TRUE, fill = TRUE) %>%
  filter(., pvalue_FE <= 0.01)

HC_red <- read.table(here("/CPTP-data/Hair-colour-paper","sumstats-HC-red.out"), 
                     sep=" ", header=TRUE, fill = TRUE) %>%
  filter(., pvalue_FE <= 0.01)

HC_brown <- read.table(here("/CPTP-data/Hair-colour-paper","sumstats-HC-brown.out"), 
                     sep=" ", header=TRUE, fill = TRUE) %>%
  filter(., pvalue_FE <= 0.01)


# merge data:

# melanoma:
filter(HC_blonde, pvalue_FE <= thres_melanoma) %>% 
  inner_join(., melanoma, by = "rsid") -> melanoma_blonde

filter(HC_red, pvalue_FE <= thres_melanoma) %>% 
  inner_join(., melanoma, by = "rsid") -> melanoma_red

filter(HC_brown, pvalue_FE <= thres_melanoma) %>% 
  inner_join(., melanoma, by = "rsid") -> melanoma_brown

rbind(melanoma_blonde, melanoma_brown) %>%
  rbind(., melanoma_red) -> melanoma_out

# bsc:
filter(HC_blonde, pvalue_FE <= thres_bsc) %>% 
  inner_join(., bsc, by = "rsid") -> bsc_blonde

filter(HC_red, pvalue_FE <= thres_bsc) %>% 
  inner_join(., bsc, by = "rsid") -> bsc_red

filter(HC_brown, pvalue_FE <= thres_bsc) %>% 
  inner_join(., bsc, by = "rsid") -> bsc_brown

rbind(bsc_blonde, bsc_brown) %>% rbind(., bsc_red) -> bsc_out

# cscc
filter(HC_blonde, pvalue_FE <= thres_cscc) %>% 
  inner_join(., cscc, by = "rsid") -> cscc_blonde

filter(HC_red, pvalue_FE <= thres_cscc) %>% 
  inner_join(., cscc, by = "rsid") -> cscc_red

filter(HC_brown, pvalue_FE <= thres_cscc) %>% 
  inner_join(., cscc, by = "rsid") -> cscc_brown

rbind(cscc_blonde, cscc_brown) %>% rbind(., cscc_red) -> cscc_out

# suntan
filter(HC_blonde, pvalue_FE <= thres_suntan) %>% 
  inner_join(., suntan, by = "rsid") -> suntan_blonde

filter(HC_red, pvalue_FE <= thres_suntan) %>% 
  inner_join(., suntan, by = "rsid") -> suntan_red

filter(HC_brown, pvalue_FE <= thres_suntan) %>% 
  inner_join(., suntan, by = "rsid") -> suntan_brown

rbind(suntan_blonde, suntan_brown) %>% rbind(., suntan_red) -> suntan_out

write.table(melanoma_out[,c(1:10)], here("/CPTP-data/skin-cancer-GWAS-catalogue/output/","melanoma_follow-up.txt"), row.names = F, sep = "\t", quote = F, na = "NA")
write.table(bsc_out[,c(1:10)], here("/CPTP-data/skin-cancer-GWAS-catalogue/output/","bsc_follow-up.txt"), row.names = F, sep = "\t", quote = F, na = "NA")
write.table(cscc_out[,c(1:10)], here("/CPTP-data/skin-cancer-GWAS-catalogue/output/","cscc_follow-up.txt"), row.names = F, sep = "\t", quote = F, na = "NA")
write.table(suntan_out[,c(1:10)], here("/CPTP-data/skin-cancer-GWAS-catalogue/output/","suntan_follow-up.txt"), row.names = F, sep = "\t", quote = F, na = "NA")

