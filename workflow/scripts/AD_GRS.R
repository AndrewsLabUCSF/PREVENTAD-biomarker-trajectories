#===========================================
#   AD Genetic Risk Scores
#   Data set constuction
#
#===========================================

#Che, R., & Motsinger-Reif, A. A. (2012). A new explained-variance based genetic risk score for predictive modeling of disease risk. Statistical Applications in Genetics and Molecular Biology, 11(4), Article 15.

#Che, R., & Motsinger-Reif, A. A. (2013). Evaluation of genetic risk score models in the presence of interaction and linkage disequilibrium. Frontiers in Genetics, 4, 138. doi:10.3389/fgene.2013.00138

#-------------------------------------------------------------
#  Required Packages
#-------------------------------------------------------------

library(foreign)
library(car)
library(reshape)
library(memisc)
library(plyr)
library(dplyr)

#-------------------------------------------------------------
#  Functions
#-------------------------------------------------------------
#compares curents risk scores to those previously calculated
#check(variable1, variabl2)
check <- function(x, y){
  thetable <- table(x, y, useNA = 'always')
  list(table = thetable)
}

#-------------------------------------------------------------
#  Importing data & Subset Data
#-------------------------------------------------------------

#Import
snps.raw <- as.data.frame(as.data.set(spss.system.file("~/Dropbox/Research/PhD/Analysis/0 Dataset Construction/1) Raw Data/PATH_OpenArray_Genotypes.sav")))
newload.raw <- read.csv("~/Dropbox/Research/Data/Original (PATH)/Genetics/PATH_OpenArray_load16.csv", header = TRUE)
apoe.raw <- read.spss("~/Dropbox/Research/PhD/Analysis/0 Dataset Construction/1) Raw Data/PATH_genes_total_120511.sav", use.value.labels = FALSE, to.data.frame = TRUE)
GRSweights <- read.csv("~/Dropbox/Research/PhD/Analysis/0 Dataset Construction/1) Raw Data/GRS weights.csv", header = TRUE)

#Susbet
snps <- subset(snps.raw, cohort == 60, select = c(pathid, rs3764650, rs744373,  rs9296559, rs34813869, rs11136000,  rs3818361,  rs11767557, rs4938933,  rs670139,  rs610932,  rs3851179))
newload <- subset(newload.raw, select = c(pathid, rs9271100, rs28834970, rs11218343, rs10498633, rs8093731, rs35349669, rs304132, rs2718058, rs1476679, rs7933019, rs17125944, rs927174))

snps <- merge(snps, newload, by = 'pathid', all = TRUE)

apoe <- subset(apoe.raw, AGEGRP == 60, select = c(pathid, apoecode))

#-------------------------------------------------------------
#  Recode SNP data
#-------------------------------------------------------------
#0 = homozygous protective allele 
#1 = heterozygours risk allele 
#2 = homozygours risk allele 

snps$ABCA7 <- with(snps, car:::recode(rs3764650, "'T/T' = 0; 'G/T' = 1; 'G/G' = 2; else = NA", as.factor.result = FALSE))
snps$BIN1 <- with(snps, car:::recode(rs744373, "'A/A' = 0; 'A/G' = 1; 'G/G' = 2; else = NA", as.factor.result = FALSE))
snps$CD2AP <- with(snps, car:::recode(rs9296559, "'T/T' = 0; 'C/T' = 1; 'C/C' = 2; else = NA", as.factor.result = FALSE))
snps$CD33 <- with(snps, car:::recode(rs34813869, "'G/G' = 0; 'A/G' = 1; 'A/A' = 2; else = NA", as.factor.result = FALSE))
snps$CLU <- with(snps, car:::recode(rs11136000, "'T/T' = 0; 'C/T' = 1; 'C/C' = 2; else = NA", as.factor.result = FALSE))
snps$CR1 <- with(snps, car:::recode(rs3818361, "'G/G' = 0; 'A/G' = 1; 'A/A' = 2; else = NA", as.factor.result = FALSE))
snps$EPHA1 <- with(snps, car:::recode(rs11767557, "'C/C' = 0; 'C/T' = 1; 'T/T' = 2; else = NA", as.factor.result = FALSE))
snps$MS4A4A <- with(snps, car:::recode(rs4938933, "'C/C' = 0; 'C/T' = 1; 'T/T' = 2; else = NA", as.factor.result = FALSE))
snps$MS4A4E <- with(snps, car:::recode(rs670139, "'G/G' = 0; 'G/T' = 1; 'T/T' = 2; else = NA", as.factor.result = FALSE))
snps$MS4A6A <- with(snps, car:::recode(rs610932, "'T/T' = 0;'G/T' = 1; 'G/G' = 2; else = NA", as.factor.result = FALSE))
snps$PICALM <- with(snps, car:::recode(rs3851179, "'T/T' = 0; 'C/T' = 1; 'C/C' = 2; else = NA", as.factor.result = FALSE))

snps$HLA <- with(snps, car:::recode(rs9271100, "'C/C' = 0; 'C/T' = 1; 'T/T' = 2; else = NA", as.factor.result = FALSE))
snps$PTK2B <- with(snps, car:::recode(rs28834970, "'T/T' = 0; 'C/T' = 1; 'C/C' = 2; else = NA", as.factor.result = FALSE))
snps$SORL1 <- with(snps, car:::recode(rs11218343, "'C/C' = 0; 'C/T' = 1; 'T/T' = 2; else = NA", as.factor.result = FALSE))
snps$SLC24A4_RIN3 <- with(snps, car:::recode(rs10498633, "'T/T' = 0; 'G/T' = 1; 'G/G' = 2; else = NA", as.factor.result = FALSE))
snps$DSG2 <- with(snps, car:::recode(rs8093731, "'T/T' = 0; 'C/T' = 1; 'C/C' = 2; else = NA", as.factor.result = FALSE))
snps$INPP5D <- with(snps, car:::recode(rs35349669, "'C/C' = 0; 'C/T' = 1; 'T/T' = 2; else = NA", as.factor.result = FALSE))
snps$MEF2C <- with(snps, car:::recode(rs304132, "'A/A' = 0; 'A/G' = 1; 'G/G' = 2; else = NA", as.factor.result = FALSE))
snps$NME8 <- with(snps, car:::recode(rs2718058, "'G/G' = 0; 'A/G' = 1; 'A/A' = 2; else = NA", as.factor.result = FALSE))
snps$ZCWPW1 <- with(snps, car:::recode(rs1476679, "'C/C' = 0; 'C/T' = 1; 'T/T' = 2; else = NA", as.factor.result = FALSE))
snps$CELF1 <- with(snps, car:::recode(rs7933019, "'G/G' = 0; 'C/G' = 1; 'C/C' = 2; else = NA", as.factor.result = FALSE))
snps$FERMT2 <- with(snps, car:::recode(rs17125944, "'T/T' = 0;'C/T' = 1; 'C/C' = 2; else = NA", as.factor.result = FALSE))
snps$CASS4 <- with(snps, car:::recode(rs927174, "'C/C' = 0; 'A/C' = 1; 'A/A' = 2; else = NA", as.factor.result = FALSE))

#0 = non carrier 
#1 = risk allele carrier
snps$ABCA7_carrier <- with(snps, car:::recode(ABCA7, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$BIN1_carrier <- with(snps, car:::recode(BIN1, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$CD2AP_carrier <- with(snps, car:::recode(CD2AP, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$CD33_carrier <- with(snps, car:::recode(CD33, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$CLU_carrier <- with(snps, car:::recode(CLU, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$CR1_carrier <- with(snps, car:::recode(CR1, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$EPHA1_carrier <- with(snps, car:::recode(EPHA1, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$MS4A4A_carrier <- with(snps, car:::recode(MS4A4A, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$MS4A4E_carrier <- with(snps, car:::recode(MS4A4E, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$MS4A6A_carrier <- with(snps, car:::recode(MS4A6A, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$PICALM_carrier <- with(snps, car:::recode(PICALM, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))

snps$HLA_carrier <- with(snps, car:::recode(HLA, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$PTK2B_carrier <- with(snps, car:::recode(PTK2B, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$SORL1_carrier <- with(snps, car:::recode(SORL1, "0 = 0; 1 = 0; 2 = 1", as.factor.result = FALSE))  #recessive model
snps$SLC24A4_RIN3_carrier <- with(snps, car:::recode(SLC24A4_RIN3, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$DSG2_carrier <- with(snps, car:::recode(DSG2, "0 = 0; 1 = 0; 2 = 1", as.factor.result = FALSE))    #recessive model
snps$INPP5D_carrier <- with(snps, car:::recode(INPP5D, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$MEF2C_carrier <- with(snps, car:::recode(MEF2C, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$NME8_carrier <- with(snps, car:::recode(NME8, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$ZCWPW1_carrier <- with(snps, car:::recode(ZCWPW1, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$CELF1_carrier <- with(snps, car:::recode(CELF1, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$FERMT2_carrier <- with(snps, car:::recode(FERMT2, "0 = 0; 1 = 1; 2 = 1", as.factor.result = FALSE))
snps$CASS4_carrier <- with(snps, car:::recode(CASS4, "0 = 0; 1 = 0; 2 = 1", as.factor.result = FALSE))  #recessive model

#-------------------------------------------------------------
#  Recode APOE SNP data
#-------------------------------------------------------------
#apoecode: e2/e2 = 1; e3/e3 = 2; e4/e4 = 3; e2/e3 = 4; e2/e4 = 5; e3/e4 = 6

#coding e2 carriers, e4 carriers and e3/e3 (reference)
#apoegeno: e3/e3 = 0; e2+ = 1; e4+ = 2
apoe$apoegeno <- with(apoe, recode(apoecode, "1 = 1; 2 = 0; 3 = 2; 4 = 1; 5 = NA; 6 = 2", as.factor.result = FALSE))

#codeing APOE e4 genotypes
#2 = e3/e3; 1 = e3/e4; 2 = e4/e4
apoe$apoee4_geno <- with(apoe, recode(apoecode, "2 = 0; 6 = 1; 3 = 2; else = 0", as.factor.result = FALSE))
apoe $apoe4 <- with(apoe, car:::recode(apoecode, '2 = 0; 3 = 1; 6 = 1; else = NA'))

#codeing APOE e2 genotypes
#2 = e3/e3; 1 = e3/e2; 2 = e2/e2
apoe$apoee2_geno <- with(apoe, recode(apoecode, "2 = 0; 4 = 1; 1 = 2; else = 0", as.factor.result = FALSE))
apoe $apoe2 <- with(apoe, car:::recode(apoecode, '2 = 0; 1 = 1; 4 = 1; else = NA'))

#merge APOE and OpenArray genotype data
geno <- merge(apoe, snps[, c(1, 25:70)], by = "pathid", all = TRUE)
apply(geno[,8:30], 2, table, useNA = 'always')
apply(geno[,31:53], 2, table, useNA = 'always')

#-------------------------------------------------------------
#  Genetic Risk Scores 
#-------------------------------------------------------------

#Calculating genotype weights
GRSweights$Wor <- log(GRSweights$OR) #weight of OR GRS 
GRSweights$Wev <- log(GRSweights$OR)*sqrt((2*GRSweights$MAF)*(1-GRSweights$MAF)) #weights for EV GRS

#caculating the weight for the genotypes at each SNP for each individual
#weights based on OR 
geno$worapoee4 <- geno$apoee4_geno*GRSweights[GRSweights$Gene == "APOE e4","Wor"]
geno$worapoee2 <- geno$apoee2_geno*GRSweights[GRSweights$Gene == "APOE e2","Wor"]
geno$WorABCA7 <- geno$ABCA7*GRSweights[GRSweights$Gene == "ABCA7","Wor"]
geno$worBIN1 <- geno$BIN1*GRSweights[GRSweights$Gene == "BIN1","Wor"]
geno$worCD2AP <- geno$CD2AP*GRSweights[GRSweights$Gene == "CD2AP","Wor"]
geno$worCD33 <- geno$CD33*GRSweights[GRSweights$Gene == "CD33","Wor"]
geno$worCLU <- geno$CLU*GRSweights[GRSweights$Gene == "CLU","Wor"]
geno$worCR1 <- geno$CR1*GRSweights[GRSweights$Gene == "CR1","Wor"]
geno$worEPHA1 <- geno$EPHA1*GRSweights[GRSweights$Gene == "EPHA1","Wor"]
geno$worMS4A4A <- geno$MS4A4A*GRSweights[GRSweights$Gene == "MS4A4A","Wor"]
geno$worMS4A4E <- geno$MS4A4E*GRSweights[GRSweights$Gene == "MS4A4E","Wor"]
geno$worMS4A6A <- geno$MS4A6A*GRSweights[GRSweights$Gene == "MS4A6A","Wor"]
geno$worPICALM <- geno$PICALM*GRSweights[GRSweights$Gene == "PICALM","Wor"]

geno$worHLA <- geno$HLA*GRSweights[GRSweights$Gene == "HLA-DRB5","Wor"]
geno$worPTK2B <- geno$PTK2B*GRSweights[GRSweights$Gene == "PTK2B","Wor"]
geno$worSORL1 <- geno$SORL1*GRSweights[GRSweights$Gene == "SORL1","Wor"]
geno$worSLC24A4_RIN3 <- geno$SLC24A4_RIN3*GRSweights[GRSweights$Gene == "SLC24A4-RIN3","Wor"]
geno$worDSG2 <- geno$DSG2*GRSweights[GRSweights$Gene == "DSG2","Wor"]
geno$worINPP5D <- geno$INPP5D*GRSweights[GRSweights$Gene == "INPP5D","Wor"]
geno$worMEF2C <- geno$MEF2C*GRSweights[GRSweights$Gene == "MEF2C","Wor"]
geno$worNME8 <- geno$NME8*GRSweights[GRSweights$Gene == "NME8","Wor"]
geno$worZCWPW1 <- geno$ZCWPW1*GRSweights[GRSweights$Gene == "ZCWPW1","Wor"]
geno$worCELF1 <- geno$CELF1*GRSweights[GRSweights$Gene == "CELF1","Wor"]
geno$worFERMT2 <- geno$FERMT2*GRSweights[GRSweights$Gene == "FERMT2","Wor"]
geno$worCASS4 <- geno$CASS4*GRSweights[GRSweights$Gene == "CASS4","Wor"]


#weights based on OR and MAF (explained variance)
geno$Wevapoee4 <- geno$apoee4_geno*GRSweights[GRSweights$Gene == "APOE e4","Wev"]
geno$Wevapoee2 <- geno$apoee2_geno*GRSweights[GRSweights$Gene == "APOE e2","Wev"]
geno$WevABCA7 <- geno$ABCA7*GRSweights[GRSweights$Gene == "ABCA7","Wev"]
geno$WevBIN1 <- geno$BIN1*GRSweights[GRSweights$Gene == "BIN1","Wev"]
geno$WevCD2AP <- geno$CD2AP*GRSweights[GRSweights$Gene == "CD2AP","Wev"]
geno$WevCD33 <- geno$CD33*GRSweights[GRSweights$Gene == "CD33","Wev"]
geno$WevCLU <- geno$CLU*GRSweights[GRSweights$Gene == "CLU","Wev"]
geno$WevCR1 <- geno$CR1*GRSweights[GRSweights$Gene == "CR1","Wev"]
geno$WevEPHA1 <- geno$EPHA1*GRSweights[GRSweights$Gene == "EPHA1","Wev"]
geno$WevMS4A4A <- geno$MS4A4A*GRSweights[GRSweights$Gene == "MS4A4A","Wev"]
geno$WevMS4A4E <- geno$MS4A4E*GRSweights[GRSweights$Gene == "MS4A4E","Wev"]
geno$WevMS4A6A <- geno$MS4A6A*GRSweights[GRSweights$Gene == "MS4A6A","Wev"]
geno$WevPICALM <- geno$PICALM*GRSweights[GRSweights$Gene == "PICALM","Wev"]

geno$WevHLA <- geno$HLA*GRSweights[GRSweights$Gene == "HLA-DRB5","Wev"]
geno$WevPTK2B <- geno$PTK2B*GRSweights[GRSweights$Gene == "PTK2B","Wev"]
geno$WevSORL1 <- geno$SORL1*GRSweights[GRSweights$Gene == "SORL1","Wev"]
geno$WevSLC24A4_RIN3 <- geno$SLC24A4_RIN3*GRSweights[GRSweights$Gene == "SLC24A4-RIN3","Wev"]
geno$WevDSG2 <- geno$DSG2*GRSweights[GRSweights$Gene == "DSG2","Wev"]
geno$WevINPP5D <- geno$INPP5D*GRSweights[GRSweights$Gene == "INPP5D","Wev"]
geno$WevMEF2C <- geno$MEF2C*GRSweights[GRSweights$Gene == "MEF2C","Wev"]
geno$WevNME8 <- geno$NME8*GRSweights[GRSweights$Gene == "NME8","Wev"]
geno$WevZCWPW1 <- geno$ZCWPW1*GRSweights[GRSweights$Gene == "ZCWPW1","Wev"]
geno$WevCELF1 <- geno$CELF1*GRSweights[GRSweights$Gene == "CELF1","Wev"]
geno$WevFERMT2 <- geno$FERMT2*GRSweights[GRSweights$Gene == "FERMT2","Wev"]
geno$WevCASS4 <- geno$CASS4*GRSweights[GRSweights$Gene == "CASS4","Wev"]

#-------------------------------------------------------------
#  Genetic risk scores
#-------------------------------------------------------------

#summing the weights for each individuels SNPs (including APOE)
geno$SCgrs_apoe <- rowSums(geno[,c('apoee4_geno', 'ABCA7', 'BIN1', 'CD2AP', 'CD33', 'CLU', 'CR1', 'EPHA1', 'MS4A4A', 'MS4A4E', 'MS4A6A', 'PICALM')], na.rm = FALSE) #effect of apoe e2 is not included
geno$ORgrs_apoe <- rowSums(geno[54:66], na.rm = FALSE)
geno$EVgrs_apoe <- rowSums(geno[79:91], na.rm = FALSE)

#summing the weights for each individuels SNPs (APOE excluded)
geno$SCgrs <- rowSums(geno[8:18], na.rm = FALSE)
geno$ORgrs <- rowSums(geno[56:66], na.rm = FALSE)
geno$EVgrs <- rowSums(geno[81:91], na.rm = FALSE)

#Overall: All IGAP SNPs plus APOE.
#including new load variants
#summing the weights for each individuels SNPs (including APOE)
geno$SCgrs_all_apoe <- rowSums(geno[,c(4, 8:30)], na.rm = FALSE) #effect of apoe e2 is not included
geno$ORgrs_all_apoe <- rowSums(geno[54:78], na.rm = FALSE)
geno$EVgrs_all_apoe <- rowSums(geno[79:103], na.rm = FALSE)

#summing the weights for each individuels SNPs (APOE excluded)
geno$SCgrs_all <- rowSums(geno[,c(8:30)], na.rm = FALSE)
geno$ORgrs_all <- rowSums(geno[56:78], na.rm = FALSE)
geno$EVgrs_all <- rowSums(geno[81:103], na.rm = FALSE)

##See Darst 2016 Pathway-Specific Polygenic Risk Scores as Predictors of Amyloid-beta Deposition and Cognitive Function in a Sample at Increased Risk for Alzheimer's Disease. 
#Immune Response Pathway: INPP5D, CLU, CR1, HLA-DRB1, MEF2C, PTK2B 
geno$IRPgrs <- rowSums(geno[,c('worINPP5D', 'worCLU', 'worCR1', 'worHLA', 'worPTK2B')])
#Aβ Clearance Pathway: APOE, CLU, CR1, PICALM.
geno$ACPgrs <- rowSums(geno[,c('worapoee4', 'worCLU', 'worCR1', 'worPICALM')])
#Cholesterol Pathway: APOE, ABCA7, CLU.
geno$CPgrs <- rowSums(geno[,c('worapoee4', 'WorABCA7', 'worCLU')])


#-------------------------------------------------------------
#  Export Data
#-------------------------------------------------------------
#export
write.csv(geno, "~/Dropbox/Research/PhD/Analysis/0 Dataset Construction/2) Derived Data/AD_GRS.csv")
grs <- read.csv("~/Dropbox/Research/PhD/Analysis/0 Dataset Construction/2) Derived Data/AD_GRS.csv", header = TRUE, row.names = 1)



#-------------------------------------------------------------
#  Plots
#-------------------------------------------------------------


library(ggplot2)
library(gridExtra)
p1 <- ggplot(data = geno, aes(x = factor(0), y = SCgrs_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Blue") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Simple Count GRS", y = "Genetic Risk Scores") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

p2 <- ggplot(data = geno, aes(x = factor(0), y = ORgrs_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Red") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Odds Ratio GRS", y = "") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

p3 <- ggplot(data = geno, aes(x = factor(0), y = EVgrs_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Green") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Explained Variance GRS", y = "") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

grid.arrange(p1, p2, p3, ncol = 3)

p4 <- ggplot(data = geno, aes(x = factor(0), y = SCgrs_all_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Blue") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Simple Count GRS", y = "Genetic Risk Scores") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

p5 <- ggplot(data = geno, aes(x = factor(0), y = ORgrs_all_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Red") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Odds Ratio GRS", y = "") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

p6 <- ggplot(data = geno, aes(x = factor(0), y = EVgrs_all_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Green") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Explained Variance GRS", y = "") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

grid.arrange(p4, p5, p6, ncol = 3)


#-------------------------------------------------------------
#  normalized

df <- dplyr:::select(geno, pathid, SCgrs_apoe, ORgrs_apoe, EVgrs_apoe, SCgrs_all_apoe, ORgrs_all_apoe, EVgrs_all_apoe)
zscores <- function(x){
  out <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  out
}

df$zSCgrs_apoe <- zscores(df$SCgrs_apoe)
df$zORgrs_apoe <- zscores(df$ORgrs_apoe)
df$zEVgrs_apoe <- zscores(df$EVgrs_apoe)
df$zSCgrs_all_apoe <- zscores(df$SCgrs_all_apoe)
df$zORgrs_all_apoe <- zscores(df$ORgrs_all_apoe)
df$zEVgrs_all_apoe <- zscores(df$EVgrs_all_apoe)

p1 <- ggplot(data = df, aes(x = factor(0), y = zSCgrs_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Blue") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Simple Count GRS", y = "Genetic Risk Scores") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

p2 <- ggplot(data = df, aes(x = factor(0), y = zORgrs_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Red") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Odds Ratio GRS", y = "") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

p3 <- ggplot(data = df, aes(x = factor(0), y = zEVgrs_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Green") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Explained Variance GRS", y = "") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

grid.arrange(p1, p2, p3, ncol = 3)

p4 <- ggplot(data = df, aes(x = factor(0), y = zSCgrs_all_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Blue") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Simple Count GRS", y = "Genetic Risk Scores") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

p5 <- ggplot(data = df, aes(x = factor(0), y = zORgrs_all_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Red") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Odds Ratio GRS", y = "") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

p6 <- ggplot(data = df, aes(x = factor(0), y = EVgrs_all_apoe)) + 
  geom_violin(scale = "width", alpha = .5, fill = "Green") +  
  geom_boxplot(width = .15) + 
  theme_bw() + labs( x = "Explained Variance GRS", y = "") + 
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=10), axis.text.x=element_blank())

grid.arrange(p4, p5, p6, ncol = 3)





























