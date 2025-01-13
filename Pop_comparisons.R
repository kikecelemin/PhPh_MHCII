setwd("C:/Users/kikec/Desktop/KIKE/PhD_Potsdam_Universität/Harbour_Porpoise/Whole_genome_resequencing/WGS_Data_Analysis/3rd_paper/06_Pop_comparisons/")

library(ggplot2)

Divergence <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")

###DRB1 Functional divergence
p1<-ggplot(Divergence,aes(x=Region,y=DRB1_functional_distance,fill=Region))
p1.box<-p1+geom_boxplot(notch=FALSE,outlier.colour="black", outlier.shape=16,outlier.size=2,)
pbox1<- p1.box+scale_fill_manual(values = c("1_BLS"="black", "2_P.p.p"="white","3_CA"="#FFC000", "4_GRE"="#3B7D23","5_ICE"="#A6A6A6", "6_NOS / SKA"="#A70000","7_BES"="#024C81", "8_PBS"="#82AEC0","9_IBE"="#DAB600"))+theme_classic()
pbox1

pdf("Intra_individual_DRB1_functional_divergence_per_pop.pdf", 7, 5)
pbox1
dev.off()

###DRB1 Genetic divergence
p2<-ggplot(Divergence,aes(x=Region,y=DRB1_genetic_distance,fill=Region))
p2.box<-p2+geom_boxplot(notch=FALSE,outlier.colour="black", outlier.shape=16,outlier.size=2,)
pbox2<- p2.box+scale_fill_manual(values = c("1_BLS"="black", "2_P.p.p"="white","3_CA"="#FFC000", "4_GRE"="#3B7D23","5_ICE"="#A6A6A6", "6_NOS / SKA"="#A70000","7_BES"="#024C81", "8_PBS"="#82AEC0","9_IBE"="#DAB600"))+theme_classic()
pbox2

pdf("Intra_individual_DRB1_genetic_divergence_per_pop.pdf", 7, 5)
pbox2
dev.off()

###DQB Functional divergence
p3<-ggplot(Divergence,aes(x=Region,y=DQB_functional_distance,fill=Region))
p3.box<-p3+geom_boxplot(notch=FALSE,outlier.colour="black", outlier.shape=16,outlier.size=2,)
pbox3<- p3.box+scale_fill_manual(values = c("1_BLS"="black", "2_P.p.p"="white","3_CA"="#FFC000", "4_GRE"="#3B7D23","5_ICE"="#A6A6A6", "6_NOS / SKA"="#A70000","7_BES"="#024C81", "8_PBS"="#82AEC0","9_IBE"="#DAB600"))+theme_classic()
pbox3

pdf("Intra_individual_DQB_functional_divergence_per_pop.pdf", 7, 5)
pbox3
dev.off()

###DQB Genetic divergence
p4<-ggplot(Divergence,aes(x=Region,y=DQB_genetic_distance,fill=Region))
p4.box<-p4+geom_boxplot(notch=FALSE,outlier.colour="black", outlier.shape=16,outlier.size=2,)
pbox4<- p4.box+scale_fill_manual(values = c("1_BLS"="black", "2_P.p.p"="white","3_CA"="#FFC000", "4_GRE"="#3B7D23","5_ICE"="#A6A6A6", "6_NOS / SKA"="#A70000","7_BES"="#024C81", "8_PBS"="#82AEC0","9_IBE"="#DAB600"))+theme_classic()
pbox4

pdf("Intra_individual_DQB_genetic_divergence_per_pop.pdf", 7, 5)
pbox4
dev.off()

##Subspecies
Subspecies_comparison <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")
shapiro.test(Subspecies_comparison$DQB_functional_distance)

wilcox.test(DQB_genetic_distance ~ Region, data = Subspecies_comparison)
wilcox.test(DQB_functional_distance ~ Region, data = Subspecies_comparison)
wilcox.test(DRB1_genetic_distance ~ Region, data = Subspecies_comparison)
wilcox.test(DRB1_functional_distance ~ Region, data = Subspecies_comparison)

##P.p.p regions
Region_comparison <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")

kruskal.test(DQB_genetic_distance ~ Region, data = Region_comparison)
kruskal.test(DQB_functional_distance ~ Region, data = Region_comparison)
kruskal.test(DRB1_genetic_distance ~ Region, data = Region_comparison)
kruskal.test(DRB1_functional_distance ~ Region, data = Region_comparison)

##CAN
CAN <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")
shapiro.test(CAN$DQB_functional_distance)

wilcox.test(DQB_genetic_distance ~ Region, data = CAN)
wilcox.test(DQB_functional_distance ~ Region, data = CAN)
wilcox.test(DRB1_genetic_distance ~ Region, data = CAN)
wilcox.test(DRB1_functional_distance ~ Region, data = CAN)

##GRE
GRE <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")
shapiro.test(GRE$DRB1_functional_distance)

wilcox.test(DQB_genetic_distance ~ Region, data = GRE)
wilcox.test(DQB_functional_distance ~ Region, data = GRE)
wilcox.test(DRB1_genetic_distance ~ Region, data = GRE)
wilcox.test(DRB1_functional_distance ~ Region, data = GRE)

##ICE
ICE <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")

wilcox.test(DQB_genetic_distance ~ Region, data = ICE)
wilcox.test(DQB_functional_distance ~ Region, data = ICE)
wilcox.test(DRB1_genetic_distance ~ Region, data = ICE)
wilcox.test(DRB1_functional_distance ~ Region, data = ICE)

##NOS
NOS <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")

wilcox.test(DQB_genetic_distance ~ Region, data = NOS)
wilcox.test(DQB_functional_distance ~ Region, data = NOS)
wilcox.test(DRB1_genetic_distance ~ Region, data = NOS)
wilcox.test(DRB1_functional_distance ~ Region, data = NOS)

##BES
BES <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")

wilcox.test(DQB_genetic_distance ~ Region, data = BES)
wilcox.test(DQB_functional_distance ~ Region, data = BES)
wilcox.test(DRB1_genetic_distance ~ Region, data = BES)
wilcox.test(DRB1_functional_distance ~ Region, data = BES)

##PBS
PBS <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")

wilcox.test(DQB_genetic_distance ~ Region, data = PBS)
wilcox.test(DQB_functional_distance ~ Region, data = PBS)
wilcox.test(DRB1_genetic_distance ~ Region, data = PBS)
wilcox.test(DRB1_functional_distance ~ Region, data = PBS)

##True_PBS
True_PBS <- read.table(file= file.choose(), header=TRUE, sep = ",", dec = ".")

wilcox.test(DQB_genetic_distance ~ Region, data = True_PBS)
wilcox.test(DQB_functional_distance ~ Region, data = True_PBS)
wilcox.test(DRB1_genetic_distance ~ Region, data = True_PBS)
wilcox.test(DRB1_functional_distance ~ Region, data = True_PBS)

