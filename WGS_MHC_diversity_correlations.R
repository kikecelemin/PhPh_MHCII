
setwd("C:/Users/kikec/Desktop/KIKE/PhD_Potsdam_Universität/Harbour_Porpoise/Whole_genome_resequencing/WGS_Data_Analysis/3rd_paper/05_Genome-wide_MHC_relationship/")

library(tidyverse)
library(stats)

###Load data
WGS_MHC <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")

###WGS_Het_DQQ_genetic_distance
pdf("WGS_Het_DQB_genetic_distance_correlation.pdf")
ggplot(data= WGS_MHC, aes(x=Heterozygosity,y=DQB_genetic_distance))+ geom_point(color="#024C81")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
dev.off()

lm_WGS_DQB_Het=lm(Heterozygosity~DQB_genetic_distance, data = WGS_MHC)
summary(lm_WGS_DQB_Het)

###WGS_Het_DQQ_functional_distance
pdf("WGS_Het_DQB_functional_distance_correlation.pdf")
ggplot(data= WGS_MHC, aes(x=Heterozygosity,y=DQB_functional_distance))+ geom_point(color="#024C81")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
dev.off()

lm_WGS_DQB_Het=lm(Heterozygosity~DQB_functional_distance, data = WGS_MHC)
summary(lm_WGS_DQB_Het)

###WGS_Het_DRB2_genetic_distance
pdf("WGS_Het_DRB2_genetic_distance_correlation.pdf")
ggplot(data= WGS_MHC, aes(x=Heterozygosity,y=DRB2_genetic_distance))+ geom_point(color="#024C81")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
dev.off()

lm_WGS_DRB2_Het=lm(Heterozygosity~DRB2_genetic_distance, data = WGS_MHC)
summary(lm_WGS_DRB2_Het)

###WGS_Het_DRB2_functional_distance
pdf("WGS_Het_DRB2_functional_distance_correlation.pdf")
ggplot(data= WGS_MHC, aes(x=Heterozygosity,y=DRB2_functional_distance))+ geom_point(color="#024C81")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
dev.off()

lm_WGS_DQB_Het=lm(Heterozygosity~DRB2_functional_distance, data = WGS_MHC)
summary(lm_WGS_DQB_Het)

###Correlation functional and genetic distance
pdf("Functional_genetic_correlation_DRB2.pdf")
ggplot(data= WGS_MHC, aes(x=DRB2_genetic_distance,y=DRB2_functional_distance))+ geom_point(color="#024C81")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
dev.off()

lm_WGS_DQB_Het=lm(DRB2_genetic_distance~DRB2_functional_distance, data = WGS_MHC)
summary(lm_WGS_DQB_Het)

pdf("Functional_genetic_correlation_DQB.pdf")
ggplot(data= WGS_MHC, aes(x=DQB_genetic_distance,y=DQB_functional_distance))+ geom_point(color="#024C81")+ theme_classic()+geom_smooth(method=lm,color="#a70000")
dev.off()

lm_WGS_DQB_Het=lm(DQB_genetic_distance~DQB_functional_distance, data = WGS_MHC)
summary(lm_WGS_DQB_Het)


###Logistic regression MHC zygosity and genome-wide heterozygosity
##DRB1
log_data <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")
glm_DRB1<-glm(DRB1_Zygosity ~ Heterozygosity, data = log_data, family = "binomial")
summary(glm_DRB1)

glm_DRB1<-glm(DRB1_Supertype_Zygosity ~ Heterozygosity, data = log_data, family = "binomial")
summary(glm_DRB1)

log_data_plot <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")
p5<-ggplot(log_data_plot,aes(x=DRB1_Zygosity,y=Heterozygosity,fill=DRB1_Zygosity))
p5.box<-p5+geom_boxplot(notch=FALSE,outlier.colour="black", outlier.shape=16,outlier.size=2,)
pbox5<- p5.box+scale_fill_manual(values = c("Homo"="#A70000","Hetero"="#024C81"))+theme_classic()
pbox5

p7<-ggplot(log_data_plot,aes(x=DRB1_Supertype_Zygosity,y=Heterozygosity,fill=DRB1_Supertype_Zygosity))
p7.box<-p7+geom_boxplot(notch=FALSE,outlier.colour="black", outlier.shape=16,outlier.size=2,)
pbox7<- p7.box+scale_fill_manual(values = c("Homo"="#A70000","Hetero"="#024C81"))+theme_classic()
pbox7

pdf("DRB1_zygosity_genome_heterozygosity.pdf")
pbox5
dev.off()

pdf("DRB1_supertype_zygosity_genome_heterozygosity.pdf")
pbox7
dev.off()

##DQB
glm_DQB<-glm(DQB_Zygosity ~ Heterozygosity, data = log_data, family = "binomial")
summary(glm_DQB)

glm_DQB<-glm(DQB_Supertype_Zygosity ~ Heterozygosity, data = log_data, family = "binomial")
summary(glm_DQB)

p6<-ggplot(log_data_plot,aes(x=DQB_Zygosity,y=Heterozygosity,fill=DQB_Zygosity))
p6.box<-p6+geom_boxplot(notch=FALSE,outlier.colour="black", outlier.shape=16,outlier.size=2,)
pbox6<- p6.box+scale_fill_manual(values = c("Homo"="#A70000","Hetero"="#024C81"))+theme_classic()
pbox6

p8<-ggplot(log_data_plot,aes(x=DQB_Supertype_Zygosity,y=Heterozygosity,fill=DQB_Supertype_Zygosity))
p8.box<-p8+geom_boxplot(notch=FALSE,outlier.colour="black", outlier.shape=16,outlier.size=2,)
pbox8<- p8.box+scale_fill_manual(values = c("Homo"="#A70000","Hetero"="#024C81"))+theme_classic()
pbox8

pdf("DQB_zygosity_genome_heterozygosity.pdf")
pbox6
dev.off()

pdf("DQB_supertype_zygosity_genome_heterozygosity.pdf")
pbox8
dev.off()

