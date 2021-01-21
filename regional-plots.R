# 07/2020

# This script will generate regional plots from FINEMAP results, either using the log10BF of the PIP on the Y-axis. 

library(tidyverse)
library(ggplot2)

setwd("")


####################### first extract those SNPs with Bayes Factor >= 2 and save them in a text file:
#######################
#######################

pheno="HC"   # phenotype to read and output as prefix
model=1      # model to read and output as prefix
#sub=4
 
for (sub in 1:10) {
snpfile <- read.csv(paste(pheno,"/",pheno,"-M",model,"-subset",sub,".snp",sep=""), sep = " ", header = TRUE)
filter(snpfile, log10bf >= 2) %>% .[,c(2,3,4,7,11,12)] -> snpfile_bf2
write.table(snpfile_bf2, paste(pheno, "/snp-BF-2/",pheno,"-M",model,"-subset",sub,".snpBF",sep=""), quote = FALSE, row.names = FALSE, sep = " ")
}

######################## plot the PIP results for each region (Y axis; colored by Bayes Factor).
######################## This section takes as input the .snp file from FINEMAP
######################## It also uses the 'snpfile' created in the above section
######################## It outputs a regional plot in PDF and PNG formats.

pheno="EC"
model=1

for (sub in 1:5) {
  snpfile <- read.csv(paste(pheno,"/",pheno,"-M",model,"-subset",sub,".snp",sep=""), sep = " ", header = TRUE)
  chr <- snpfile[1,3]

  ggplot(snpfile, aes(x = position, colour = prob, y = log10bf)) +
  scale_y_continuous() +  geom_point(aes(colour = prob)) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill=FALSE) + labs( x = paste("Chromosome",chr,sep=" "), y = "log10bf") + 
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=14))
  ggsave(file=paste("PLOTS-PIP-BF/",pheno,"/",pheno,"-M",model,"-subset",sub,".png", sep=""))
  ggsave(file=paste("PLOTS-PIP-BF/",pheno,"/",pheno,"-M",model,"-subset",sub,".pdf", sep=""))
}

######################## plot the PIP results for each region (Y axis; colored by LD).
######################## This section takes as input the .snp file from FINEMAP
######################## It also uses the 'snpfile' created in the above section
######################## It also needs the local LD for each top SNP computed with LDStore (r values)
######################## It outputs a regional plot in PDF and PNG formats.

pheno="EC"
model=1

for (sub in 1:5) {
  
  # read SNP file from FINEMAP results for each subset:
  snpfile <- read.csv(paste(pheno,"/",pheno,"-M",model,"-subset",sub,".snp",sep=""), sep = " ", header = TRUE)
  # extract the chromosome #:
  chr <- snpfile[1,3]  
  
  # read the local LD results computed on LDstore:
  ld <- read.csv(paste("results-local-ld/",pheno,"-M",model,"-subset",sub,".ld",sep=""), sep = " ", header = TRUE)
  
  # separate rsids and other columns to compute R^2:
  snps <- colnames(ld[2:ncol(ld)])
  head(ld,n=2)
  name_tmp <- data.frame(ld[,1])
  ld_tmp <- data.frame(ld[,2:ncol(ld)]) 
  ld_tmp <- data.frame(lapply(ld_tmp,"^",2))   # get the R^2

  # bind the tables back again:
  ld <- bind_cols(name_tmp, ld_tmp)
  
  # rename columns:
  colnames(ld) <- c("rsid", snps)  # change colnames for LD

  # join snp table and LD table:
  snp_plot <- left_join(snpfile, ld, by = "rsid")  # join LD data with finemap data
  
  # get the label for the legend:
  rsq=expression(paste(r^2))

  # plot and save png:
  for (x in 1:length(snps)) {
   index <- snps[x]
  
   top <- filter(snp_plot, rsid == !!index)   # filter to retain only index SNP
   ggplot(snp_plot, aes_string(x = "position", colour = index, y = "log10bf")) +
    scale_y_continuous() +  
     geom_point(colour = "black", size = 1.5) +
     geom_point(aes_string(colour = index), size = 1, shape = 19) +
    #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
    scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
    theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
    guides(fill=FALSE) + labs( x = paste("Chromosome",chr,sep=" "), y = "log10 (Bayes Factor)") +
    theme(axis.text.x = element_text(face="bold", size=12), axis.text.y = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
    theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) + #xlim(88079191,89512533) # define precise zlim
    geom_point(data=top, colour = "black", size = 5.1, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4.6) +
     geom_hline(yintercept=2, colour = "black", linetype="dashed")
    ggsave(file=paste("PLOTS-PIP-LD/",pheno,"/",pheno,"-M",model,"-subset",sub,"-",snps[x],"-logBF.png", sep=""))
    ggsave(file=paste("PLOTS-PIP-LD/",pheno,"/",pheno,"-M",model,"-subset",sub,"-",snps[x],"-logBF.pdf", sep=""))
  }
}
