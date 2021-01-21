install.packages("qqman")
install.packages("dplyr")
library(qqman)
library(dplyr)

# FUNCTION TO OBTAIN LAMBDA (GENOMIC CONTROL), QQPLOT AND MANHATTAN PLOT.
# IT USES 2 FILES (with header = TRUE): 
  #file_to_read1 = data frame with p-values only (output from post-gwas script for SNPTEST)
  #file_to_read2 = data frame with 4 columns (SNP, CHR, BP, P) - the clean2 file output from post-gwas script for SNPTEST
# IT NEEDS AS INPUT THE FOLDER WHERE THESE TWO FILES ARE (path_to_folder argument)
# IT NEEDS AS INPUT THE FILE NAME FOR SAVING THE PLOTS (file_name argument)
# The lambda genomic control value will be pasted in the qqplot.
# This function includes the installation of 2 packages: "dplyr" and "qqman".

# example: 
 #post_gwas_plots(path_to_folder = "C:/Users/frida/Dropbox/iris-structures/results-add-new/EUR/", file_to_read1 = "eur-nodules-p-new.out",
 #file_to_read2 = "eur-nodules-clean2-new.txt", file_name = "eur-nodules-fun")

post_gwas_plots <- function(path_to_folder, file_to_read1, file_to_read2, file_name, lambda_pos) {
  
  setwd(path_to_folder)
  
  #read file1:
  p_values <- read.delim(file = file_to_read1, header = TRUE)
  colnames(p_values) <- "frequentist_add_pvalue"

   #genomic control (lambda):
   chisq <- qchisq(1-p_values$frequentist_add_pvalue,1)
   lambda_value <- median(chisq)/qchisq(0.5,1)
   lambda_value <- round(lambda_value, digits = 4)

   qq_filename <- paste("qqplot-", file_name, ".jpeg", sep="")

   #qqplot:
   jpeg(filename = qq_filename, width = 500, height = 500)
   qq(p_values$frequentist_add_pvalue, cex.axis = 1.5)
   text(x = 1, y = lambda_pos, labels = paste("Lambda =", lambda_value, sep=" "))
   dev.off()

   rm(p_values)

  # data wrangling for manhattan plot:
  clean <- read.delim(file = file_to_read2, header = TRUE, sep = " ")
  colnames(clean) <- c("SNP", "CHR", "BP", "P")
  clean$SNP <- as.character(clean$SNP)
  clean2 <- filter(clean, SNP != ".") #remove unnamed snps
  clean3 <- filter(clean2, P <= 0.5)

  rm(clean)
  rm(clean2)

  man_filename <- paste("manhattan-", file_name, ".jpeg", sep="")

  #manhattan plot:
  jpeg(filename = man_filename, width = 1200, height = 1200)
  manhattan(clean3, col=c("black","gray"), chrlabs = c(1:22), cex.axis = 2, ylim=c(0,200))
  dev.off()
}
