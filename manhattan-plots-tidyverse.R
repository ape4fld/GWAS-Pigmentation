# 04/2020

# This script runs a function from external source ('ggmanhattan-function.R') to save a Manhattan plot using the Tidyverse package, taken from ________.

# Input files are: 
      # 1) a list (single column, no header) of genome-wide significant SNPs to highlight in red.
      # 2) a table with GWAS summary statistics (rsid, chr, pos, p-value) with header.
      # 3) need to change the directories and paths to where the data is stored.

# It outputs the plot using a chosen name and path in PDF and PNG versions.

setwd("") # set working directory

install.packages("here") # install package if not already in local computer.
library(here)

source(here("GWAS", "ggmanhattan-function.R"))  # where the actual function is stored.
home_dir=here()
mypalette2 <- c("#E69F00", "#56B4E9")  # default colours are light blue and light orange; change as desired.

sig = 5e-8 # significant threshold line to plot.
sugg = 1e-6 # suggestive threshold line to plot.

# Run Function ====
mysnps <- read.table(here("snps-to-hlight.txt"), header=F, colClasses = "character") # read a table with only genome-wide significant SNPs in a column
mysnps <- mysnps[,]

model=M1  # this is a variable used for the prefix name of the plot.

man_filename <- paste("manhattan-",model,"meta", sep="")       # this is the output prefix name of the plot.
gg.manhattan(path = paste(home_dir,"/results", sep = ""),      # this is the path where to store the plot.
             df = "sumstats-table.out",                        # this is the name of the sumstats file.
             threshold = NA, hlight = mysnps, col = mypalette2, ylims = c(0,100),      # change Y-axis limit if needed (-log10(pvalue)).
             title = paste("Manhattan Meta",model,sep=""))     # Title in the plot.
ggsave(file=paste(man_filename,".png", sep=""))
ggsave(file=paste(man_filename,".pdf", sep=""))
