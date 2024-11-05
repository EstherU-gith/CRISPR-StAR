###------------------------------------------------ Load packages --------------------------------------------------------

library(tidyverse)
library(data.table)
library(ggrepel)
library(xlsx)
library(dplyr)
#library(seqinr) # write fasta files
library(stringr)
library(grid)
library(gridExtra)
library("gtools")
library(ggExtra)
#library(ggstatsplot)
library("ggpubr")
library(patchwork) #plot_spacer
# library(ggblend)
#conflict_prefer("select", "dplyr")
library(scales)
# install.packages('pROC')
library(pROC)

# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
### ------------------------------------------------ Load control lists --------------------------------------------------------
Essentials<-read.table('/Volumes/groups/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/coreessentials.tsv', header=T )
Controls<- read.table('/Volumes/groups/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/nonessentials.tsv', header=T )
Vitro.A.B.E.depleting.min3.RRAsc10 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Vitro.A.B.E.depleting.min3.RRAsc10.txt', sep = "\t", header = T)
#Vitro.A.B.E.depleting.min3.RRAsc10 = 566 genes
Vitro.A.B.E.FC2_depleting.min3.RRAsc10 <- read.table('/Volumes/Esther/NGS_results/Scripts_gene_lists/Gene_lists_ABE_screen/Vitro.A.B.E.FC2_depleting.min3.RRAsc10.txt', sep = "\t", header = T)
#Vitro.A.B.E.FC2_depleting.min3.RRAsc10 = 566 genes
Mouse.ortholog_Controls_472 <- read.table("/Volumes/groups/elling/Esther/NGS_results/Scripts_gene_lists/Mouse.ortholog_Controls_472.txt", sep = "\t", header = T)
Mouse.ortholog_Controls_550 <- read.table("/Volumes/groups/elling/Esther/NGS_results/Scripts_gene_lists/Mouse.ortholog_Controls_550.txt", sep = "\t", header = T)
Mouse.ortholog_Controls_650 <- read.table("/Volumes/groups/elling/Esther/NGS_results/Scripts_gene_lists/Mouse.ortholog_Controls_650.txt", sep = "\t", header = T)
Mouse.ortholog_CoreEssentials_934 <- read.table("/Volumes/groups/elling/Esther/NGS_results/Scripts_gene_lists/Mouse.ortholog_CoreEssentials_934.txt", sep = "\t", header = T)

#

# VITRO # ---------- A.B.E_FC2  guide - load data ------------------ 
FC1.reads.i.a.A <-read.table("/Volumes/groups/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/FC1.reads.i.a.txt",header = T, sep = "\t")
FC1.reads.i.a.vitro.A <- FC1.reads.i.a.A %>%
  filter(vitro_vivo == "in_vitro") %>% group_by(guide) %>% mutate(Gene = str_split_fixed(guide, '_', 2)[,1])
FC1.reads.i.a.vitro.A <- FC1.reads.i.a.vitro.A %>% group_by(guide, Gene, Replicate) %>% summarise(active=sum(active),inactive=sum(inactive))
#write.table(FC1.reads.i.a.vitro.A, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/FC1.reads.i.a.vitro.A.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
FC1.reads.i.a.vitro.A <- read.table("/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/FC1.reads.i.a.vitro.A.txt", header = T, sep = "\t")

FC2.unfiltered.reads.vitro.B.E<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/FC2.unfiltered.reads.vitro.B.E.txt", header=T, sep = '\t')
FC2.only.guides.vitro.B.E <- FC2.unfiltered.reads.vitro.B.E %>% group_by(guide, Gene, Replicate, inactive_active) %>% summarise(reads=sum(reads))
FC2_reads.i.a.vitro.B.E <- reshape2::dcast(FC2.only.guides.vitro.B.E %>% dplyr::select(guide, Gene, Replicate, inactive_active, reads),
                                           guide + Gene + Replicate ~ inactive_active, fill = 0 , value.var = 'reads')
#write.table(FC2_reads.i.a.vitro.B.E, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/FC2_reads.i.a.vitro.B.E.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
FC2_reads.i.a.vitro.B.E <- read.table("/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/FC2_reads.i.a.vitro.B.E.txt", header = T, sep = "\t")

##### Vitro # ---------- A.B.E - Merge ------------------
A.B.E.FC2_guide.reads.vitro.Repl <- merge(FC1.reads.i.a.vitro.A, FC2_reads.i.a.vitro.B.E, by = c('guide', 'Gene','Replicate'), all = T, suffixes = c(".A",".B.E"), fill = 0)

##### Vitro # ---------- A.B.E - Pool ------------------
A.B.E.FC2_guide.reads.vitro.Repl$active <- rowSums(A.B.E.FC2_guide.reads.vitro.Repl[, c("active.A", "active.B.E")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.Repl$inactive <- rowSums(A.B.E.FC2_guide.reads.vitro.Repl[, c("inactive.A", "inactive.B.E")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.Repl$sumReads <- rowSums(A.B.E.FC2_guide.reads.vitro.Repl[, c("active", "inactive")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.Repl_wide <- dcast(setDT(A.B.E.FC2_guide.reads.vitro.Repl), guide + Gene ~ Replicate, value.var = c('active','inactive'), fill = 0) %>% data.frame()
#write.table(A.B.E.FC2_guide.reads.vitro.Repl_wide, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.Repl_wide.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E.FC2_guide.reads.vitro.Repl_wide<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.Repl_wide.txt", header=T, sep = '\t')

##### Vitro # ---------- A.B.E - Add pseudo count ------------------
#A.B.E library 106.601 guides, 21.348 genes
A.B.E.FC2_guide.reads.vitro.Repl_wide_pseudo <- A.B.E.FC2_guide.reads.vitro.Repl_wide %>%
  mutate_at(vars(active_Blue, active_Green, active_Red, inactive_Blue, inactive_Green, inactive_Red), ~ . + 0.5)  # 106,638 guides
#write.table(A.B.E.FC2_guide.reads.vitro.Repl_wide_pseudo, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.Repl_wide_pseudo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# VITRO # ---------- RATIO A screen and B.E screen --------------------------------------------------------

vitro.read_counts.A <- FC1.reads.i.a.vitro.A %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

vitro.read_counts.A <- vitro.read_counts.A %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
vitro.read_counts.A <- vitro.read_counts.A %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
vitro.read_counts.A <- vitro.read_counts.A %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
vitro.read_counts.A <- vitro.read_counts.A %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
vitro.read_counts.A <- vitro.read_counts.A %>%
  mutate(Data = "Vitro_counts_A")


vitro.read_counts.B.E <- FC2_reads.i.a.vitro.B.E %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

vitro.read_counts.B.E <- vitro.read_counts.B.E %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
vitro.read_counts.B.E <- vitro.read_counts.B.E %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
vitro.read_counts.B.E <- vitro.read_counts.B.E %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
vitro.read_counts.B.E <- vitro.read_counts.B.E %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
vitro.read_counts.B.E <- vitro.read_counts.B.E %>%
  mutate(Data = "Vitro_counts_B.E")


read_counts.A.B.E <- A.B.E_FC2.pooled_most_filters.ia0 %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
vitro.read_counts.B.E <- vitro.read_counts.B.E %>%
  mutate(Data = "Vitro_counts_B.E")


read_counts.A.B.E.sumfilter <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3 %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
vitro.read_counts.B.E <- vitro.read_counts.B.E %>%
  mutate(Data = "Vitro_counts_B.E")

# Ratio.read.counts <- rbind(read_counts.A, read_counts.B.E, read_counts.A.B.E,read_counts.A.B.E.sumfilter)
Ratio.read.counts.vitro <- rbind(vitro.read_counts.A, vitro.read_counts.B.E)

# VITRO # ---------- Ratio A screen and B.E screen separate - normalisation -------------------------------------------------------
head(A.B.E.FC2_guide.reads.vitro.Repl)

Vitro.NonEss_A.B.E_guide_UMI <- A.B.E.FC2_guide.reads.vitro.Repl %>%
  filter((Gene) %in% Mouse.ortholog_Controls_550$gene)
head(Vitro.NonEss_A.B.E_guide_UMI)
read_counts.Ctrl.A.B.E.vitro <- Vitro.NonEss_A.B.E_guide_UMI %>%
  summarise(
    Reads_Active.A = sum(active.A, na.rm = TRUE),
    Reads_Inactive.A = sum(inactive.A, na.rm = TRUE),
    Reads_Active.B.E = sum(active.B.E, na.rm = TRUE),
    Reads_Inactive.B.E = sum(inactive.B.E, na.rm = TRUE))

read_counts.Ctrl.A.B.E.vitro <- read_counts.Ctrl.A.B.E.vitro %>%
  mutate(Ratio_Active_to_Inactive.A = round(Reads_Active.A / Reads_Inactive.A, 1),
         Ratio_Active_to_Inactive.B.E = round(Reads_Active.B.E / Reads_Inactive.B.E, 1))
read_counts.Ctrl.A.B.E.vitro <- read_counts.Ctrl.A.B.E.vitro %>%
  mutate(SumReads.A = Reads_Active.A + Reads_Inactive.A,
         SumReads.B.E = Reads_Active.B.E + Reads_Inactive.B.E)
read_counts.Ctrl.A.B.E.vitro <- read_counts.Ctrl.A.B.E.vitro %>%
  mutate(Perc.active.A = round((Reads_Active.A / SumReads.A)*100, 2),
         Perc.inactive.A = round((Reads_Inactive.A / SumReads.A)*100, 2),
         Perc.active.B.E = round((Reads_Active.B.E / SumReads.B.E)*100, 2),
         Perc.inactive.B.E = round((Reads_Inactive.B.E / SumReads.B.E)*100, 2))
#read_counts.Ctrl.A.B.E.vitro <- read_counts.Ctrl.A.B.E.vitro %>% mutate(Perc.total = )

head(read_counts.Ctrl.A.B.E.vitro)
read_counts.Ctrl.A.B.E.vitro.long <- reshape2::melt(read_counts.Ctrl.A.B.E.vitro, id.vars = c("Sample", "Cells", "StAR.vector", "Virus", "StAR.vector.New"), measure.vars = c("Perc.total", "Perc.inactive"))                                  

Vitro.NonEss_A.B.E_guide_UMI <- Vitro.NonEss_A.B.E_guide_UMI %>%
  mutate(sumReads.A = active.A + inactive.A,
         sumReads.B.E = active.B.E + inactive.B.E)
Vitro.NonEss_A.B.E_guide_UMI.Ratio <- Vitro.NonEss_A.B.E_guide_UMI %>% mutate(Perc.active.A = round((active.A / sumReads.A)*100, 2),
                                                                              Perc.inactive.A = round((inactive.A / sumReads.A)*100, 2),
                                                                              Perc.active.B.E = round((active.B.E / sumReads.B.E)*100, 2),
                                                                              Perc.inactive.B.E = round((inactive.B.E / sumReads.B.E)*100, 2))
read_counts.Ctrl.A.B.E.vitro.median <- Vitro.NonEss_A.B.E_guide_UMI.Ratio %>%
  summarise(
    Perc.active.A = median(Perc.active.A, na.rm = TRUE),
    Perc.inactive.A = median(Perc.inactive.A, na.rm = TRUE),
    Perc.active.B.E = median(Perc.active.B.E, na.rm = TRUE),
    Perc.inactive.B.E = median(Perc.inactive.B.E, na.rm = TRUE))

head(read_counts.Ctrl.A.B.E.vitro.median)

Normalisation.factor.for.activeB.E.vitro <- read_counts.Ctrl.A.B.E.vitro.median$Perc.active.A/read_counts.Ctrl.A.B.E.vitro.median$Perc.active.B.E
Normalisation.factor.for.inactiveB.E.vitro <- read_counts.Ctrl.A.B.E.vitro.median$Perc.inactive.A/read_counts.Ctrl.A.B.E.vitro.median$Perc.inactive.B.E


A.B.E.FC2_guide.reads.vitro.Repl.norm <- A.B.E.FC2_guide.reads.vitro.Repl %>% 
  mutate(active.B.E.norm = round((active.B.E *Normalisation.factor.for.activeB.E.vitro),2),
         inactive.B.E.norm = round((inactive.B.E *Normalisation.factor.for.inactiveB.E.vitro),2))
head(A.B.E.FC2_guide.reads.vitro.Repl.norm)


##### Vitro # ---------- A.B.E - Pool & add pseudo count after normalization ------------------
A.B.E.FC2_guide.reads.vitro.Repl.norm$active <- rowSums(A.B.E.FC2_guide.reads.vitro.Repl.norm[, c("active.A", "active.B.E.norm")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.Repl.norm$inactive <- rowSums(A.B.E.FC2_guide.reads.vitro.Repl.norm[, c("inactive.A", "inactive.B.E.norm")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.Repl.norm$sumReads <- rowSums(A.B.E.FC2_guide.reads.vitro.Repl.norm[, c("active", "inactive")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.Repl_wide.norm <- dcast(setDT(A.B.E.FC2_guide.reads.vitro.Repl.norm), guide + Gene ~ Replicate, value.var = c('active','inactive'), fill = 0) %>% data.frame()
#write.table(A.B.E.FC2_guide.reads.vitro.Repl_wide.norm, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.Repl_wide.norm.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#A.B.E.FC2_guide.reads.vitro.Repl_wide.norm <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.Repl_wide.norm.txt", header=T, sep = '\t')

#PSEUDO
#A.B.E library 106.601 guides, 21.348 genes
A.B.E.FC2_guide.reads.vitro.Repl_wide.norm_pseudo <- A.B.E.FC2_guide.reads.vitro.Repl_wide.norm %>%
  mutate_at(vars(active_Blue, active_Green, active_Red, inactive_Blue, inactive_Green, inactive_Red), ~ . + 0.5)  # 106,638 guides
write.table(A.B.E.FC2_guide.reads.vitro.Repl_wide.norm_pseudo, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.Repl_wide.norm_pseudo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



##### Vitro # ---------- A.B.E - Pre-Mageck ------------------
# For guide only Paired analysis
#write.table(A.B.E.FC2_guide.reads.vitro.Repl_wide_pseudo, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E.FC2_guide.reads.vitro.Repl_wide_pseudo.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Pre_Mageck_A.B.E.FC2_guide.reads.vitro.Repl_wide_pseudo <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E.FC2_guide.reads.vitro.Repl_wide_pseudo.tsv', header = T,  sep = "\t")

# sbatch mageck test -k A.B.E.FC2_guide.reads.vitro.Repl_wide_pseudo.tsv -t active_Blue,active_Green,active_Red -c inactive_Blue,inactive_Green,inactive_Red -n Mageck.A.B.E_FC2.invitro.pseudo.paired.standard --remove-zero any --paired 

# Vitro # ---------- A.B.E - Pre-Mageck - after normalization ------------------
# For guide only Paired analysis
#write.table(A.B.E.FC2_guide.reads.vitro.Repl_wide.norm_pseudo, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E.FC2_guide.reads.vitro.Repl_wide.norm_pseudo.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Pre_Mageck_A.B.E.FC2_guide.reads.vitro.Repl_wide.norm_pseudo <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E.FC2_guide.reads.vitro.Repl_wide.norm_pseudo.tsv', header = T,  sep = "\t")

# sbatch mageck test -k Pre_Mageck_A.B.E.FC2_guide.reads.vitro.Repl_wide.norm_pseudo.tsv -t active_Blue,active_Green,active_Red -c inactive_Blue,inactive_Green,inactive_Red -n Mageck.A.B.E_FC2.invitro.norm.pseudo.paired.standard --remove-zero any --paired 

# VITRO # ---------- MAGECK ------------------
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck.A.B.E_FC2.invitro.pseudo.paired.standard.gene_summary.txt', header = T)
#mageck_A.B.E.FC2_invitro_guides_pseudo.paired.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck.A.B.E_FC2.invitro.pseudo.paired.standard.sgrna_summary.txt', header = T)
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old$minRRA.fdr <- -log10(apply(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old$minRRA.score <- -log10(apply(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old[,c("neg.score","pos.score")], 1, min))
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old$minpvalue <- -log10(apply(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old$gene <- str_split_fixed(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old$id, '_', 2)[,1] 

# VITRO # ---------- MAGECK - after normalization ------------------------------------ 

#from here the dataframes are labeled the same as without normalization to keep it easier to rerun certain parts of the script

mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck.A.B.E_FC2.invitro.norm.pseudo.paired.standard.gene_summary.txt', header = T)
#mageck_A.B.E.FC2_invitro_guides_pseudo.paired.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck.A.B.E_FC2.invitro.norm.pseudo.paired.standard.sgrna_summary.txt', header = T)
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard$minRRA.fdr <- -log10(apply(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard$minRRA.score <- -log10(apply(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard[,c("neg.score","pos.score")], 1, min))
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard$minpvalue <- -log10(apply(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard$gene <- str_split_fixed(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard$id, '_', 2)[,1] 



# VITRO # ---------- Number of guides or UMIS per gene ------------------
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>%
  ggplot(aes(x = num)) +
  geom_histogram(binwidth = 1, color = "black", fill= "grey") +
  theme_bw()+
  xlab('Number of guide_UMIs per gene') + ylab('Count')

nrGene_invitro <- nrow(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard) #A.B.E. pools 21348 genes, this screen 21361 genes

# VITRO # ---------- Defining in vitro essentials ------------------
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard
Vitro.A.B.E.FC2_depleting.min3.RRAsc10 <- mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(neg.lfc <= -3 & minRRA.score >= 10)
#write.table(Vitro.A.B.E.FC2_depleting.min3.RRAsc10, '/Volumes/Esther/NGS_results/Scripts_gene_lists/Gene_lists_ABE_screen/Vitro.A.B.E.FC2_depleting.min3.RRAsc10.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Vitro.A.B.E.FC2_depleting.min2.5.RRAsc10 <- mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(neg.lfc <= -2.5 & minRRA.score >= 10)
#write.table(Vitro.A.B.E.FC2_depleting.min2.5.RRAsc10, '/Volumes/Esther/NGS_results/Scripts_gene_lists/Gene_lists_ABE_screen/Vitro.A.B.E.FC2_depleting.min2.5.RRAsc10.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Vitro.A.B.E.FC2_depleting.min2.RRAsc10 <- mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(neg.lfc <= -2 & minRRA.score >= 10)
#write.table(Vitro.A.B.E.FC2_depleting.min2.RRAsc10, '/Volumes/Esther/NGS_results/Scripts_gene_lists/Gene_lists_ABE_screen/Vitro.A.B.E.FC2_depleting.min2.RRAsc10.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Vitro_old vs Vitro.norm # ---------- Correlation ---------- 
mageck_A.B.E.FC2_invitro.old_norm_genes.standard <- merge(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard_old [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
                                                         mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
                                                         by = 'id', suffixes = c(".invitro_old",".invitro"))
head(mageck_A.B.E.FC2_invitro.old_norm_genes.standard)
colors4 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'Depleting validation library' = 'orange', 
             'Enriching validation library' = 'limegreen', 'mito_inner_membr_GO_0005743' = 'Blue', 'Vivo.validation.enriching.vivo1x' = 'magenta3',
             'MicroLibrary' = 'cyan')
mageck_A.B.E.FC2_invitro.old_norm_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro_old, y = neg.lfc.invitro)) +
  geom_hline(yintercept=0, color="gray20", linewidth = 0.5) +
  geom_vline(xintercept=0, color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other')) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invitro.old_norm_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invitro.old_norm_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) +
  xlab(bquote('Log'[2]~'fold change in vitro before recalculation')) + ylab(bquote('Log'[2]~'fold change in vitro'))+ labs(color = ' ') + 
  theme_pubr(base_size = 14) + scale_color_manual(values = colors4)  +
  xlim(-8.2,5.2) + 
  ylim(-8.2,5.2) +
  theme(plot.title = element_text(hjust = 0.5))

#
# VITRO # ---------- StAR Volcano plots - minRRA.score ------------------
colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(toupper(gene) %in% Controls$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'StAR analysis \nIn vitro gene LFC + RRA score \nno filter on sumReads') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,32) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')

mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter((gene) %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'StAR analysis \nIn vitro gene LFC + RRA score \nno filter on sumReads') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,33) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')

#
#!!! VITRO # ---------- StAR Volcano plots - minRRA.score ------------------ RRA.score max 15 ------------------ With density plots ---------------
VPlot1 <- 
  mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'),  size = 1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter((gene) %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ 
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,33) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")  #+ #guides(colour = 'none')

Plot.Dens.V1x <- mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>%
  ggplot(aes(x = neg.lfc)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  xlab(NULL)  +
  xlim(-7.5,7.5) +  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  theme(legend.position = "none", axis.text.x = element_blank())

Plot.Dens.V1y <- mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>%
  ggplot(aes(x = minRRA.score)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  #geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), alpha=0.5, linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  xlab(NULL) +
  xlim(0,33) + theme_pubr(base_size = 12) + scale_color_manual(values = colors) + 
  theme(legend.position = "none", axis.text.y = element_blank()) + coord_flip()

Plot.Dens.V1x + plot_spacer() + VPlot1 + theme(aspect.ratio = 1)+ Plot.Dens.V1y +
  plot_layout(ncol = 2, nrow = 2, widths = c(2, 0.5), heights = c(0.5, 2), guides = 'collect') & theme(legend.position = 'top')

#
#!!! VITRO # dAUC # ---------- StAR ---------- ---------- no Others! different Essential groups ---------- 
#A.B.E.FC2_guide.reads.vitro.noRepl <- A.B.E.FC2_guide.reads.vitro.Repl %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
A.B.E.FC2_guide.reads.vitro.noRepl <- A.B.E.FC2_guide.reads.vitro.Repl.norm %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo <- A.B.E.FC2_guide.reads.vitro.noRepl %>% mutate_at(vars(active, inactive), ~ . + 0.5)
A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$StAR.LFC <- log2(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$active / 
                                                             A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$inactive)

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo<- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo %>% 
  mutate(group= ifelse(Gene %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(Gene %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo %>% filter(group == "Essential" | group == "Non-Essential")
nrow( A.B.E.FC2_guide.reads.vitro.noRepl.pseudo %>% filter(group == "Non-Essential"))

A.B.E.FC2_guide.reads.vitro.noRepl.noOthers.cumulative <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers %>% 
  dplyr::select(guide, Gene, StAR.LFC, group) %>% na.omit()

A.B.E.vitro.guide.dAUC.rank.LFCStAR <- A.B.E.FC2_guide.reads.vitro.noRepl.noOthers.cumulative %>% arrange(StAR.LFC)
A.B.E.vitro.guide.dAUC.rank.LFCStAR$all_rank <- 1:nrow(A.B.E.vitro.guide.dAUC.rank.LFCStAR)
A.B.E.vitro.guide.dAUC.rank.LFCStAR$all_rank_perc_sgRNA <- (A.B.E.vitro.guide.dAUC.rank.LFCStAR$all_rank*1)/nrow(A.B.E.vitro.guide.dAUC.rank.LFCStAR) 

A.B.E.vitro.guide.dAUC.rank.LFCStAR <- A.B.E.vitro.guide.dAUC.rank.LFCStAR %>% group_by(group) %>% mutate(rank_group = 1:n(), rank_group_prct = 1:n()/n(), rank_sum = cumsum(rank_group_prct)) 
#head(A.B.E.guide.dAUC.rank.LFCStAR)

# Group the data by "group"
grouped_data.vitro_StAR <- A.B.E.vitro.guide.dAUC.rank.LFCStAR %>%
  group_by(group) %>% arrange(all_rank_perc_sgRNA)
# Function to calculate AUC using the trapezoidal rule
calculate_auc.vitro_StAR <- function(x, y) {
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2) }
# Calculate AUC for each group
auc_results.vitro_StAR <- grouped_data.vitro_StAR %>%
  summarize(AUC = calculate_auc.vitro_StAR(all_rank_perc_sgRNA, rank_group_prct))
print(auc_results.vitro_StAR)


Essential_auc.vitro_StAR <- auc_results.vitro_StAR$AUC[auc_results.vitro_StAR$group == "Essential"]
NonEssential_auc.vitro_StAR <- auc_results.vitro_StAR$AUC[auc_results.vitro_StAR$group == "Non-Essential"]
dAUC.vitro_Ess_nonEss_StAR <- Essential_auc.vitro_StAR - NonEssential_auc.vitro_StAR
# Print the dAUC value
cat("dAUC vitro =", dAUC.vitro_Ess_nonEss_StAR, "\n")

ggplot(A.B.E.vitro.guide.dAUC.rank.LFCStAR, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group)) + 
  geom_line(size = 0.75) + 
  geom_hline(yintercept = max(A.B.E.vitro.guide.dAUC.rank.LFCStAR$rank_group_prct), linetype = "solid", color = "black") +
  geom_vline(xintercept = max(A.B.E.vitro.guide.dAUC.rank.LFCStAR$all_rank_perc_sgRNA), linetype = "solid", color = "black") +
  theme_pubr() + 
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "black", "Other" = "gray")) +
  labs(title = "StAR AUC in Vitro \nEssentials: Vitro.depleting.min3.RRAsc10 \nControls: Mouse.ortholog_Controls_550", y = "Cumulative fraction", x = "Percentage rank of sgRNAs") + 
  annotate("text", x = 0.85, y = 0.1, label = paste("dAUC =", round(dAUC.vitro_Ess_nonEss_StAR, 2)), color = "black")

ggplot(A.B.E.vitro.guide.dAUC.rank.LFCStAR, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group)) + 
  geom_line(size = 0.75) + 
  geom_hline(yintercept = max(A.B.E.vitro.guide.dAUC.rank.LFCStAR$rank_group_prct), linetype = "solid", color = "black") +
  geom_vline(xintercept = max(A.B.E.vitro.guide.dAUC.rank.LFCStAR$all_rank_perc_sgRNA), linetype = "solid", color = "black") +
  theme_pubr() + 
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "black", "Other" = "gray")) +
  labs(title = "StAR AUC in Vitro \nEssentials: Vitro.depleting.min3.RRAsc10 \nControls: Mouse.ortholog_Controls_550", y = "Cumulative fraction", x = "Percentage rank of sgRNAs") + 
  annotate("text", x = 0.85, y = 0.1, label = paste("dAUC =", round(dAUC.vitro_Ess_nonEss_StAR, 2)), color = "black")


#



#VITRO # dAUC # TRIALS ---------- StAR ---------- ---------- 
A.B.E.FC2_guide.reads.vitro.noRepl <- A.B.E.FC2_guide.reads.vitro.Repl %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
A.B.E.FC2_guide.reads.vitro.noRepl.pseudo <- A.B.E.FC2_guide.reads.vitro.noRepl %>% mutate_at(vars(active, inactive), ~ . + 0.5)
A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$StAR.LFC <- log2(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$active / 
                                                             A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$inactive)

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo<- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo %>% 
  mutate(group= ifelse(Gene %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(Gene %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo %>% filter(group == "Essential" | group == "Non-Essential")
nrow( A.B.E.FC2_guide.reads.vitro.noRepl.pseudo %>% filter(group == "Non-Essential"))

essential_indices <- which(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss$group == "Essential")

# Randomly select 250 indices to remove
indices_to_remove <- sample(essential_indices, 250)

# Remove the selected rows
A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss[-indices_to_remove, ]
nrow( A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss %>% filter(group == "Non-Essential"))
nrow( A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss %>% filter(group == "Essential"))


A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss.cumulative <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss %>% 
  dplyr::select(guide, Gene, StAR.LFC, group) %>% na.omit()

A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.noOthers.lessEss.cumulative %>% arrange(StAR.LFC)
A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$all_rank <- 1:nrow(A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss)
A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$all_rank_perc_sgRNA <- (A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$all_rank*1)/nrow(A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss) 

A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss <- A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss %>% group_by(group) %>% mutate(rank_group = 1:n(), rank_group_prct = 1:n()/n(), rank_sum = cumsum(rank_group_prct)) 
#head(A.B.E.guide.dAUC.rank.LFCStAR)

# Group the data by "group"
grouped_data.vitro_StAR.lessEss <- A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss %>%
  group_by(group) %>% arrange(all_rank_perc_sgRNA)
# Function to calculate AUC using the trapezoidal rule
calculate_auc.vitro_StAR <- function(x, y) {
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2) }
# Calculate AUC for each group
auc_results.vitro_StAR.lessEss <- grouped_data.vitro_StAR.lessEss %>%
  summarize(AUC = calculate_auc.vitro_StAR(all_rank_perc_sgRNA, rank_group_prct))
print(auc_results.vitro_StAR.lessEss)

Essential_auc.vitro_StAR.lessEss <- auc_results.vitro_StAR.lessEss$AUC[auc_results.vitro_StAR.lessEss$group == "Essential"]
NonEssential_auc.vitro_StAR.lessEss <- auc_results.vitro_StAR.lessEss$AUC[auc_results.vitro_StAR.lessEss$group == "Non-Essential"]
dAUC.vitro_Ess_nonEss_StAR.lessEss <- Essential_auc.vitro_StAR.lessEss - NonEssential_auc.vitro_StAR.lessEss
# Print the dAUC value
cat("dAUC vitro =", dAUC.vitro_Ess_nonEss_StAR.lessEss, "\n")

ggplot(A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group)) + 
  geom_line(size = 0.75) + 
  geom_hline(yintercept = max(A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$rank_group_prct), linetype = "solid", color = "black") +
  geom_vline(xintercept = max(A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$all_rank_perc_sgRNA), linetype = "solid", color = "black") +
  theme_pubr() + 
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "black", "Other" = "gray")) +
  labs(title = "StAR AUC in Vitro \nEssentials: Vitro.depleting.min3.RRAsc10 \nControls: Mouse.ortholog_Controls_550", y = "Cumulative fraction", x = "Percentage rank of sgRNAs") + 
  annotate("text", x = 0.85, y = 0.1, label = paste("dAUC =", round(dAUC.vitro_Ess_nonEss_StAR, 2)), color = "black")


library(zoo)

x <- 1:10
y <- 3*x+25
id <- order(x)

AUC <- sum(diff(x[id])*rollmean(y[id],2))


library(pracma)

A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess <- A.B.E.vitro.guide.dAUC.rank.LFCStAR %>% filter(group == "Essential")
A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss <- A.B.E.vitro.guide.dAUC.rank.LFCStAR %>% filter(group == "Non-Essential")

AUC.vitrolessEss <- trapz(A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$rank_group_prct)
AUC.vitroEss <- trapz(A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$rank_group_prct)
AUC.vitroNonEss <- trapz(A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$rank_group_prct)
dAUC.vitro <- AUC.vitroEss - AUC.vitroNonEss


#library(MESS)
auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$rank_group_prct, type = 'spline')
auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$rank_group_prct, type = 'spline')
#dAUC.vitro <- 
(auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$rank_group_prct, to = max(0.5),type = 'spline')
  - auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$rank_group_prct,from = max(0.5), type = 'spline'))

(auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$rank_group_prct, type = 'spline')
  - auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$rank_group_prct, type = 'spline'))

(auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.lessEss$rank_group_prct, type = 'spline')
  - auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.NonEss$rank_group_prct, type = 'spline'))

AUC(A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$rank_group_prct,  method = 'spline')



# Create two data frames for the rows you want to add
row1 <- data.frame(guide = "x", 
                   Gene = "x", 
                   StAR.LFC = -11, 
                   group = "Non-Essential", 
                   all_rank = 0, 
                   all_rank_perc_sgRNA = 0, 
                   rank_group = 0, 
                   rank_group_prct = 0, 
                   rank_sum = 0)

row2 <- data.frame(guide = "x", 
                   Gene = "x", 
                   StAR.LFC = 2, 
                   group = "Essential", 
                   all_rank = 5395, 
                   all_rank_perc_sgRNA = 1, 
                   rank_group = 2827, 
                   rank_group_prct = 1, 
                   rank_sum = 1414.5)

# Add the rows to your dataframe
A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted <- rbind(A.B.E.vitro.guide.dAUC.rank.LFCStAR, row1, row2)

#A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted <- A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted %>% group_by(group) %>% mutate(rank_group = 1:n(), rank_group_prct = 1:n()/n(), rank_sum = cumsum(rank_group_prct)) 
#head(A.B.E.guide.dAUC.rank.LFCStAR)

# Group the data by "group"
grouped_data.vitro_StAR.adjusted <- A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted %>%
  group_by(group) %>% arrange(all_rank_perc_sgRNA)


(auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adj.Ess$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.adj.Ess$rank_group_prct, type = 'spline')
  - auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adj.NonEss$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.adj.NonEss$rank_group_prct, type = 'spline'))

auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adj.Ess$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.adj.Ess$rank_group_prct)

ggplot(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group)) + 
  geom_line(size = 0.75) + 
  geom_hline(yintercept = max(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$rank_group_prct), linetype = "solid", color = "black") +
  geom_vline(xintercept = max(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$all_rank_perc_sgRNA), linetype = "solid", color = "black") +
  theme_pubr() + 
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "black", "Other" = "gray")) +
  labs(title = "StAR AUC in Vitro \nEssentials: Vitro.depleting.min3.RRAsc10 \nControls: Mouse.ortholog_Controls_550", y = "Cumulative fraction", x = "Percentage rank of sgRNAs") #+ 
# annotate("text", x = 0.85, y = 0.1, label = paste("dAUC =", round(dAUC.vitro_Ess_nonEss_StAR, 2)), color = "black")



# Assuming 'rank_group_prct' as the response and 'all_rank_perc_sgRNA' as the predictor
auc_value <- auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$rank_group_prct, A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$all_rank_perc_sgRNA)

# Calculate AUC for the "Essential" group
auc_essential <- auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$rank_group_prct[A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$group == "Essential"], 
                     A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$all_rank_perc_sgRNA[A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$group == "Essential"])

# Calculate AUC for the "Non-Essential" group
auc_non_essential <- auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$rank_group_prct[A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$group == "Non-Essential"], 
                         A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$all_rank_perc_sgRNA[A.B.E.vitro.guide.dAUC.rank.LFCStAR.adjusted$group == "Non-Essential"])

# Calculate the difference in AUC
auc_difference <- auc_essential - auc_non_essential

# Print the calculated AUC values and their difference
print(auc_essential)
print(auc_non_essential)
print(auc_difference)



#
# VITRO # rank # ---------- StAR ---------- ---------- no Others! different Essential groups ---------- 
# A.B.E.FC2_guide.reads.vitro.noRepl <- A.B.E.FC2_guide.reads.vitro.Repl %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
# A.B.E.FC2_guide.reads.vitro.noRepl.pseudo <- A.B.E.FC2_guide.reads.vitro.noRepl %>% mutate_at(vars(active, inactive), ~ . + 0.5)
# A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$StAR.LFC <- log2(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$active / 
#                                                              A.B.E.FC2_guide.reads.vitro.noRepl.pseudo$inactive)
# 
# A.B.E.FC2_guide.reads.vitro.noRepl.pseudo<- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo %>% 
#   mutate(group= ifelse(Gene %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(Gene %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo %>% arrange(StAR.LFC)
head(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank)

gene_counts <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank %>% group_by(Gene) %>% summarize(count = n())
A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank %>%
  left_join(gene_counts, by = "Gene") %>% filter(count <= 5) %>% select(-count)

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5$Rank <- 1:nrow(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5)
A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 %>% select(-active, -inactive)
head(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5)

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 %>%
  mutate(Unique_Genes_Count = cumsum(!duplicated(Gene)))

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 %>% mutate(average_guides = Rank / Unique_Genes_Count)

# A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 %>%
#   arrange(Rank) %>%
#   group_by(Gene) %>%
#   mutate(sgRNA.gene = row_number()) %>%
#   ungroup()


ggplot(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5, aes(x = Rank, y = average_guides)) +
  geom_point() + 
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene") +
  #xlim(0,5000) +
  scale_x_continuous(labels = label_scientific()) +
  labs(title = 'In vitro - StAR') +
  theme_pubr()



head(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5)
A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 %>% filter(group == "Essential") %>%
  distinct(Gene) %>%
  nrow()
#

A.B.E.FC2_guide.reads.vitro.noRepl.random <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5 %>%  mutate(random_order = runif(nrow(.))) %>%  # create a column of random numbers
  arrange(random_order) %>%  # sort the dataframe by the random column
  select(-random_order)  # remove the random column

A.B.E.FC2_guide.reads.vitro.noRepl.random$Rank <- 1:nrow(A.B.E.FC2_guide.reads.vitro.noRepl.random)


A.B.E.FC2_guide.reads.vitro.noRepl.random <- A.B.E.FC2_guide.reads.vitro.noRepl.random %>%
  mutate(Unique_Genes_Count = cumsum(!duplicated(Gene)))

A.B.E.FC2_guide.reads.vitro.noRepl.random <- A.B.E.FC2_guide.reads.vitro.noRepl.random %>% mutate(average_guides = Rank / Unique_Genes_Count)


A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5_select <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5[c("Rank", "average_guides")]
A.B.E.FC2_guide.reads.vitro.noRepl.random_select <- A.B.E.FC2_guide.reads.vitro.noRepl.random[c("Rank", "average_guides")]
combined.vitro <- merge(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5_select, A.B.E.FC2_guide.reads.vitro.noRepl.random_select, by = "Rank", suffixes = c('.vitro', '.vitro.rdm'))

colnames(combined.vitro) 

combined_long_vitro <- reshape2::melt(combined.vitro, id.vars = "Rank", variable.name = "Method", value.name = "Avg_guides")

ggplot(combined_long_vitro, aes(x = Rank, y = Avg_guides, color = Method)) +
  geom_point(data = combined_long_vitro[combined_long_vitro$Method == "average_guides.vitro.rdm",], aes(x = Rank, y = Avg_guides)) +
  geom_point(data = combined_long_vitro[combined_long_vitro$Method == "average_guides.vitro",], aes(x = Rank, y = Avg_guides)) +
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene",
       title = 'Reproducibility of sgRNAs') +
  xlim(0,3e+03)+
  scale_x_continuous(labels = scales::label_scientific()) +
  scale_color_manual(values =c('average_guides.vitro' = "#000000", 'average_guides.vitro.rdm' = "#999999")) +
  theme_pubr()

#
#

# VITRO # ---------- minRRA.score, minRRA.fdr, min.pvalue plot in VITRO ------------------
colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
colors1 <- c('Non-Essential' = 'black', 'Other' = 'grey', 'In_Vitro_Depleting' = 'red')
colors2 <- c('Non-Essential' = 'black', 'Essential' = 'red', 'Other' = 'grey', 'In_Vitro_Depleting' = 'red', 'In_Vivo_Enriching' = 'limegreen', 
             'In_Vivo_Depleting' = 'gold')

plot_VT1 <- mageck_A.B.E_invitro_genes_pseudo.standard %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_hline(aes(yintercept=0), color="gray90", linewidth = 2) +
  geom_vline(aes(xintercept=0), color="gray90", linewidth = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E_invitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E_invitro_genes_pseudo.standard %>% filter(id %in% Vitro.A.B.E.depleting.min2.RRAsc10$id), alpha=0.5, size=2.5) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'In vitro gene LFC + RRA score') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  guides(colour = 'none') 
plot_VT2 <- mageck_A.B.E_invitro_genes_pseudo.standard %>%
  ggplot(aes(x = neg.lfc, y = minRRA.fdr)) +
  geom_hline(aes(yintercept=0), color="gray90", linewidth = 2) +
  geom_vline(aes(xintercept=0), color="gray90", linewidth = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E_invitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E_invitro_genes_pseudo.standard %>% filter(id %in% Vitro.A.B.E.depleting.min2.RRAsc10$id), alpha=0.5, size=2.5) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA fdr'))+ labs(color = ' ', title = 'In vitro gene LFC + RRA FDR') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  guides(colour = 'none') 
plot_VT3 <- mageck_A.B.E_invitro_genes_pseudo.standard %>%
  ggplot(aes(x = neg.lfc, y = minpvalue)) +
  geom_hline(aes(yintercept=0), color="gray90", linewidth = 2) +
  geom_vline(aes(xintercept=0), color="gray90", linewidth = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E_invitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'Essential'), data=mageck_invitro_genes_pseudo.standard %>% filter(id %in% Vitro.A.B.E.depleting.min2.RRAsc10$id), alpha=0.5, size=2.5) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('pValue'))+ labs(color = ' ', title = 'In vitro gene LFC + pValue') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  guides(colour = 'none') 

grid.arrange(plot_VT1, plot_VT2, plot_VT3, ncol = 3, top = textGrob('In vitro LFC + significance scores \nVitro.A.B.E.depleting.min2.RRAsc10 = Essential group (Red)', gp = gpar(fontsize = 14)))
#  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots



# VITRO # ---------- Conventional ---------- 
#A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0 <- read.table("/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.txt", sep = "\t", header = T)

A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.pseudo <- A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0 %>% mutate_at(vars(pooled_Calculated_reads, active, inactive), ~ . + 0.5)  

names(A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.pseudo)

# VITRO # ---------- Conventional - after normalization ---------- 
A.B.E.FC2_guide.reads.vitro.Repl.norm

A.B.E.FC2_guide.reads.vitro.norm <- A.B.E.FC2_guide.reads.vitro.Repl.norm %>% group_by(guide, Gene) %>% summarise(active=sum(active),inactive=sum(inactive))
#write.table(A.B.E.FC2_guide.reads.vitro.norm, "/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.norm.txt", sep = "\t", row.names = F, col.names = T, quote = F)
A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.norm <- combined_data_pooled %>% left_join(A.B.E.FC2_guide.reads.vitro.norm, by = "guide") %>% replace(is.na(.), 0) #104.735 rows

A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.norm$Gene <- str_split_fixed(A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.norm$guide, '_', 2)[, 1]

A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.norm.pseudo <- A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.norm %>% mutate_at(vars(pooled_Calculated_reads, active, inactive), ~ . + 0.5)  

names(A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.norm.pseudo)

# VITRO # ---------- CONVENTIONAL - Pre-Mageck ------------------------------------ 
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt

Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt <- A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.pseudo %>% select(guide, Gene, active, pooled_Calculated_reads)
#write.table(Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt.tsv', header = T)

#sbatch mageck test -k Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt.tsv -t active -c pooled_Calculated_reads -n Mageck_A.B.E_FC2.CONV.vitro.guide.level.standard --remove-zero any 

# VITRO # ---------- CONVENTIONAL - Pre-Mageck - after normalization ------------------------------------ 

Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt.norm <- A.B.E_FC2.pooled_guides.vitro.plasmidLibr.w0.ia.0.norm.pseudo %>% select(guide, Gene, active, pooled_Calculated_reads)
#write.table(Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt.norm, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt.norm.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt.tsv', header = T)

#sbatch mageck test -k Pre_Mageck_A.B.E_F2.guide.vitro.plasmid.nofilt.norm.tsv -t active -c pooled_Calculated_reads -n Mageck_A.B.E_FC2.CONV.vitro.norm.guide.level.standard --remove-zero any 

# VITRO # ---------- CONVENTIONAL - Mageck ------------------------------------ 
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.CONV.vitro.guide.level.standard.gene_summary.txt', header = T)
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minRRA.fdr <- -log10(apply(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minRRA.score <- -log10(apply(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt[,c("neg.score","pos.score")], 1, min))
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minpvalue <- -log10(apply(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$gene <- str_split_fixed(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$id, '_', 2)[,1] # 21382 genes
length(unique(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$id)) # 21242 genes

mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$adj.minRRA.score <- ifelse(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minRRA.score > 15, 15, 
                                                                             mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minRRA.score)
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt <-  mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>%
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>%
  ggplot(aes(x = num)) +
  geom_histogram(binwidth = 1, color = "black", fill= "grey") +
  theme_bw()+
  xlim(0, 30) +
  xlab('Number of guide_UMIs per gene') + ylab('Count')

#write.table(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# VITRO # ---------- CONVENTIONAL - Mageck - after normalization------------------------------------ 

#from here the dataframes are labeled the same as without normalization to keep it easier to rerun certain parts of the script

mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.CONV.vitro.norm.guide.level.standard.gene_summary.txt', header = T)
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minRRA.fdr <- -log10(apply(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minRRA.score <- -log10(apply(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt[,c("neg.score","pos.score")], 1, min))
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minpvalue <- -log10(apply(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$gene <- str_split_fixed(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$id, '_', 2)[,1] # 21382 genes
length(unique(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$id)) # 21242 genes

mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$adj.minRRA.score <- ifelse(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minRRA.score > 15, 15, 
                                                                             mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt$minRRA.score)
mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt <-  mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>%
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

#write.table(mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt.norm.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#!!! VITRO # ---------- CONVENTIONAL Volcano plot - minRRA.score ------------------ With density plots ---------------

VPlot1 <- 
  mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'),  size = 1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>% filter((gene) %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  #geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ 
  #labs(color = ' ', title = 'Conventional analysis \nIn vitro gene LFC + RRA score \nno filters') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-10,5) + ylim(0,20) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")  #+ #guides(colour = 'none')
#theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank())  #+ #guides(colour = 'none')

Plot.Dens.V1x <- mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>%
  ggplot(aes(x = neg.lfc)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL)  +
  xlim(-10,5) +  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  theme(legend.position = "none", axis.text.x = element_blank())

Plot.Dens.V1y <- mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>%
  ggplot(aes(x = adj.minRRA.score)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  #geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), alpha=0.5, linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vitro_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL) +
  xlim(0,20) + theme_pubr(base_size = 12) + scale_color_manual(values = colors) + 
  theme(legend.position = "none", axis.text.y = element_blank()) + coord_flip()

# Plot.Dens.V1x + plot_spacer() + VPlot1 + theme(aspect.ratio = 1)+ Plot.Dens.V1y +
#   plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4), guides = 'collect') & theme(legend.position = 'top')
Plot.Dens.V1x + plot_spacer() + VPlot1 + theme(aspect.ratio = 1)+ Plot.Dens.V1y +
  plot_layout(ncol = 2, nrow = 2, widths = c(2, 0.5), heights = c(0.5, 2), guides = 'collect') & theme(legend.position = 'top')
# VPlot1 + theme(aspect.ratio = 1)+ Plot.Dens.V1y + Plot.Dens.V1x + plot_spacer() + 
#   plot_layout(ncol = 2, nrow = 2, widths = c(3, 1), heights = c(3, 1))

#


# VIVO # ---------- A screen and B.E. screen - load data --------------------------------------------------------
FC1_FC2_poly_hop.filtered.reads.vivo_0.001.A<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/FC1_FC2_poly_hop.filtered.reads.vivo_0.001.A.txt", header=T, sep = '\t')
FC2_poly_hop.filtered.reads.vivo_0.001.B.E<- read.table("/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/FC2_poly_hop.filtered.reads.vivo_0.001.B.E.txt", header=T, sep = '\t')

# This dataframes are separately the A and the B.E. screen ONLY VIVO reads, poly UMI and Hopping filtered, not yet filtered for sumReads or inactive/active reads

# Vivo # ---------- A.B.E - Merge ------------------
A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001 <- merge(FC1_FC2_poly_hop.filtered.reads.vivo_0.001.A, FC2_poly_hop.filtered.reads.vivo_0.001.B.E, by = c('guide', 'UMI', 'guide_UMI', 'Gene', 'Sample', 'Replicate'), all = T, suffixes = c(".A",".B.E"), fill = 0)
names(A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001) # "guide"        "UMI"          "guide_UMI"    "Gene"         "Sample"       "Replicate"    "active.A"     "inactive.A"   "sumReads.A"   "active.B.E"   "inactive.B.E" "sumReads.B.E"

#write.table(A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.txt", header=T, sep = '\t')

# Vivo # ---------- A.B.E -  Pool ------------------
A.B.E_FC2.pooled_most_filters.ia0 <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001
A.B.E_FC2.pooled_most_filters.ia0$active <- rowSums(A.B.E_FC2.pooled_most_filters.ia0[, c("active.A", "active.B.E")], na.rm = TRUE)
A.B.E_FC2.pooled_most_filters.ia0$inactive <- rowSums(A.B.E_FC2.pooled_most_filters.ia0[, c("inactive.A", "inactive.B.E")], na.rm = TRUE)
A.B.E_FC2.pooled_most_filters.ia0$sumReads <- rowSums(A.B.E_FC2.pooled_most_filters.ia0[, c("sumReads.A", "sumReads.B.E")], na.rm = TRUE)
A.B.E_FC2.pooled_most_filters.ia0 <- A.B.E_FC2.pooled_most_filters.ia0 %>% select(guide, UMI, guide_UMI, Gene, Sample, Replicate, active, inactive, sumReads)
#write.table(A.B.E_FC2.pooled_most_filters.ia0, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2.pooled_most_filters.ia0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E_FC2.pooled_most_filters.ia0<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2.pooled_most_filters.ia0.txt", header=T, sep = '\t')
A.B.E_FC2.pooled_most_filters.ia0 <- A.B.E_FC2.pooled_most_filters.ia0 %>% mutate(guide_UMI_Repl = paste(guide_UMI, Replicate, sep = '_')) # 1,940,058 guide_UMI_Repl
A.B.E_FC2.pooled_most_filters.ia0.guide_UMI_Repl <- A.B.E_FC2.pooled_most_filters.ia0 %>% select("guide", "UMI", "guide_UMI", "guide_UMI_Repl", "Gene", "Sample", "Replicate", "active", "inactive","sumReads") 
#write.table(A.B.E_FC2.pooled_most_filters.ia0.guide_UMI_Repl, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_most_filters.ia0.guide_UMI_Repl.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

names(A.B.E_FC2.pooled_most_filters.ia0)
head(A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001)
# Vivo # ---------- A.B.E - filter raw data - STAR ------------------ SR20.UMIs3 ------------------ 
A.B.E_FC2.pooled.guide_UMI_Repl.SR20 <- A.B.E_FC2.pooled_most_filters.ia0 %>% 
  filter(sumReads >= 20) # 185,515 guide_UMI_Repl
nrow(A.B.E_FC2.pooled.guide_UMI_Repl.SR20) # 185.515 guide_UMI_Repl
length(unique(A.B.E_FC2.pooled.guide_UMI_Repl.SR20$Gene)) # 21,309 genes
#write.table(A.B.E_FC2.pooled.guide_UMI_Repl.SR20, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled.guide_UMI_Repl.SR20.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E_FC2.pooled.guide_UMI_Repl.SR20<- read.table("/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled.guide_UMI_Repl.SR20.txt", header=T, sep = '\t')


A.B.E_FC2.pooled.guide_UMI_Repl.SR20$to_delete_UMIs3 <- TRUE # Create a new column in the original dataframe to mark rows to be deleted
A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20 %>% group_by(Gene) %>% filter(n() >= 3) %>% ungroup()
A.B.E_FC2.pooled.guide_UMI_Repl.SR20$to_delete_UMIs3[A.B.E_FC2.pooled.guide_UMI_Repl.SR20$Gene %in% A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3$Gene] <- FALSE # Mark the rows to be kept as FALSE in the original dataframe
deleted_rows.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20[
  A.B.E_FC2.pooled.guide_UMI_Repl.SR20$to_delete_UMIs3, ] # Identify rows that were deleted  # 1406 guide_UMI_Repl
A.B.E_FC2.pooled.guide_UMI_Repl.SR20 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20[, !names(A.B.E_FC2.pooled.guide_UMI_Repl.SR20) %in% "to_delete_UMIs3"] # Remove the temporary column from the original dataframe
A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3[, !names(A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3) 
                                                                                         %in% "to_delete_UMIs3"] # Remove the temporary column from the original dataframe  # 184,109 guide_UMI_Repl
#write.table(A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3 <- read.table("/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3.txt", header=T, sep = '\t')

length(unique(A.B.E_FC2.pooled.guide_UMI_Repl.SR20$Gene)) # 21309 genes
length(unique(A.B.E_FC2.pooled.guide_UMI_Repl.SR20$guide)) # 81217 genes
length(unique(A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3$Gene)) # 20490 genes
length(unique(A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3$guide)) # 79938 genes

# VIVO # ---------- Ratio A screen and B.E screen --------------------------------------------------------

read_counts.A <- FC1_FC2_poly_hop.filtered.reads.vivo_0.001.A %>%
  # group_by(Sample) %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.A <- read_counts.A %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.A <- read_counts.A %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.A <- read_counts.A %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.A <- read_counts.A %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.A <- read_counts.A %>%
  mutate(Data = "Vivo_counts_A")


read_counts.B.E <- FC2_poly_hop.filtered.reads.vivo_0.001.B.E %>%
  # group_by(Sample) %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.B.E <- read_counts.B.E %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.B.E <- read_counts.B.E %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.B.E <- read_counts.B.E %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.B.E <- read_counts.B.E %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.B.E <- read_counts.B.E %>%
  mutate(Data = "Vivo_counts_B.E")


read_counts.A.B.E <- A.B.E_FC2.pooled_most_filters.ia0 %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Data = "Vivo_counts_A.B.E")


read_counts.A.B.E.sumfilter <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3 %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Data = "Vivo_counts_A.B.E.sumfilter")

Ratio.read.counts <- rbind(read_counts.A, read_counts.B.E, read_counts.A.B.E,read_counts.A.B.E.sumfilter)
Ratio.read.counts.vt.vv <- rbind(read_counts.A, read_counts.B.E, read_counts.A.B.E,read_counts.A.B.E.sumfilter, vitro.read_counts.A, vitro.read_counts.B.E)


# VIVO # ---------- Ratio A screen and B.E screen separate - normalisation -------------------------------------------------------

A.B.E_FC2_filtered.forRatio <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001 %>% 
  mutate(sumReads.A.B.E = rowSums(across(c(sumReads.A, sumReads.B.E)), na.rm = TRUE))
A.B.E_FC2_filtered.forRatio <- A.B.E_FC2_filtered.forRatio %>% 
  filter(sumReads.A.B.E >= 20) 

A.B.E_FC2_filtered.forRatio$to_delete_UMIs3 <- TRUE # Create a new column in the original dataframe to mark rows to be deleted
A.B.E_FC2_filtered.forRatio.SR20.UMIs3 <- A.B.E_FC2_filtered.forRatio %>% group_by(Gene) %>% filter(n() >= 3) %>% ungroup()
A.B.E_FC2_filtered.forRatio$to_delete_UMIs3[A.B.E_FC2_filtered.forRatio$Gene %in% A.B.E_FC2_filtered.forRatio$Gene] <- FALSE # Mark the rows to be kept as FALSE in the original dataframe
Ratio.deleted_rows.SR20.UMIs3 <- A.B.E_FC2_filtered.forRatio[
  A.B.E_FC2_filtered.forRatio$to_delete_UMIs3, ] # Identify rows that were deleted  # 1406 guide_UMI_Repl
A.B.E_FC2_filtered.forRatio <- A.B.E_FC2_filtered.forRatio[, !names(A.B.E_FC2_filtered.forRatio) %in% "to_delete_UMIs3"] # Remove the temporary column from the original dataframe
A.B.E_FC2_filtered.forRatio.SR20.UMIs3 <- A.B.E_FC2_filtered.forRatio.SR20.UMIs3[, !names(A.B.E_FC2_filtered.forRatio.SR20.UMIs3) 
                                                                                 %in% "to_delete_UMIs3"] # Remove the temporary column from the original dataframe  # 184,109 guide_UMI_Repl

nrow(A.B.E_FC2_filtered.forRatio)
nrow(A.B.E_FC2_filtered.forRatio.SR20.UMIs3)
head(A.B.E_FC2_filtered.forRatio)
head(A.B.E_FC2_filtered.forRatio.SR20.UMIs3)

Vivo.NonEss_A.B.E_guide_UMI <- A.B.E_FC2_filtered.forRatio.SR20.UMIs3 %>%
  filter((Gene) %in% Mouse.ortholog_Controls_550$gene)
head(Vivo.NonEss_A.B.E_guide_UMI)
read_counts.Ctrl.A.B.E <- Vivo.NonEss_A.B.E_guide_UMI %>%
  summarise(
    Reads_Active.A = sum(active.A, na.rm = TRUE),
    Reads_Inactive.A = sum(inactive.A, na.rm = TRUE),
    Reads_Active.B.E = sum(active.B.E, na.rm = TRUE),
    Reads_Inactive.B.E = sum(inactive.B.E, na.rm = TRUE))

read_counts.Ctrl.A.B.E <- read_counts.Ctrl.A.B.E %>%
  mutate(Ratio_Active_to_Inactive.A = round(Reads_Active.A / Reads_Inactive.A, 1),
         Ratio_Active_to_Inactive.B.E = round(Reads_Active.B.E / Reads_Inactive.B.E, 1))
read_counts.Ctrl.A.B.E <- read_counts.Ctrl.A.B.E %>%
  mutate(SumReads.A = Reads_Active.A + Reads_Inactive.A,
         SumReads.B.E = Reads_Active.B.E + Reads_Inactive.B.E)
read_counts.Ctrl.A.B.E <- read_counts.Ctrl.A.B.E %>%
  mutate(Perc.active.A = round((Reads_Active.A / SumReads.A)*100, 2),
         Perc.inactive.A = round((Reads_Inactive.A / SumReads.A)*100, 2),
         Perc.active.B.E = round((Reads_Active.B.E / SumReads.B.E)*100, 2),
         Perc.inactive.B.E = round((Reads_Inactive.B.E / SumReads.B.E)*100, 2))

head(read_counts.Ctrl.A.B.E)
Vivo.NonEss_A.B.E_guide_UMI.Ratio <- Vivo.NonEss_A.B.E_guide_UMI %>% mutate(Perc.active.A = round((active.A / sumReads.A)*100, 2),
                                                                            Perc.inactive.A = round((inactive.A / sumReads.A)*100, 2),
                                                                            Perc.active.B.E = round((active.B.E / sumReads.B.E)*100, 2),
                                                                            Perc.inactive.B.E = round((inactive.B.E / sumReads.B.E)*100, 2))
read_counts.Ctrl.A.B.E.median <- Vivo.NonEss_A.B.E_guide_UMI.Ratio %>%
  summarise(
    Perc.active.A = median(Perc.active.A, na.rm = TRUE),
    Perc.inactive.A = median(Perc.inactive.A, na.rm = TRUE),
    Perc.active.B.E = median(Perc.active.B.E, na.rm = TRUE),
    Perc.inactive.B.E = median(Perc.inactive.B.E, na.rm = TRUE))

head(read_counts.Ctrl.A.B.E.median)

Normalisation.factor.for.activeB.E <- read_counts.Ctrl.A.B.E.median$Perc.active.A/read_counts.Ctrl.A.B.E.median$Perc.active.B.E
Normalisation.factor.for.inactiveB.E <- read_counts.Ctrl.A.B.E.median$Perc.inactive.A/read_counts.Ctrl.A.B.E.median$Perc.inactive.B.E


A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001 %>% 
  mutate(active.B.E.norm = round((active.B.E *Normalisation.factor.for.activeB.E),2),
         inactive.B.E.norm = round((inactive.B.E *Normalisation.factor.for.inactiveB.E),2))
A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm %>%
  mutate(sumReads.B.E.norm = active.B.E.norm + inactive.B.E.norm)

head(A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm)

# Vivo # ---------- A.B.E -  Pool & filter (SR20.UMIs3) after normalisation------------------
A.B.E_FC2.pooled_most_filters.ia0.norm <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm
A.B.E_FC2.pooled_most_filters.ia0.norm$active <- rowSums(A.B.E_FC2.pooled_most_filters.ia0.norm[, c("active.A", "active.B.E.norm")], na.rm = TRUE)
A.B.E_FC2.pooled_most_filters.ia0.norm$inactive <- rowSums(A.B.E_FC2.pooled_most_filters.ia0.norm[, c("inactive.A", "inactive.B.E.norm")], na.rm = TRUE)
A.B.E_FC2.pooled_most_filters.ia0.norm$sumReads <- rowSums(A.B.E_FC2.pooled_most_filters.ia0.norm[, c("sumReads.A", "sumReads.B.E.norm")], na.rm = TRUE)
A.B.E_FC2.pooled_most_filters.ia0.norm <- A.B.E_FC2.pooled_most_filters.ia0.norm %>% select(guide, UMI, guide_UMI, Gene, Sample, Replicate, active, inactive, sumReads)
#write.table(A.B.E_FC2.pooled_most_filters.ia0.norm, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2.pooled_most_filters.ia0.norm.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E_FC2.pooled_most_filters.ia0.norm<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2.pooled_most_filters.ia0.norm.txt", header=T, sep = '\t')
A.B.E_FC2.pooled_most_filters.ia0.norm <- A.B.E_FC2.pooled_most_filters.ia0.norm %>% mutate(guide_UMI_Repl = paste(guide_UMI, Replicate, sep = '_')) # 1,940,058 guide_UMI_Repl
A.B.E_FC2.pooled_most_filters.ia0.norm.guide_UMI_Repl <- A.B.E_FC2.pooled_most_filters.ia0.norm %>% select("guide", "UMI", "guide_UMI", "guide_UMI_Repl", "Gene", "Sample", "Replicate", "active", "inactive","sumReads") 
#write.table(A.B.E_FC2.pooled_most_filters.ia0.norm.guide_UMI_Repl, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_most_filters.ia0.norm.guide_UMI_Repl.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
head(A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm)
head(A.B.E_FC2.pooled_most_filters.ia0.norm)

A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20 <- A.B.E_FC2.pooled_most_filters.ia0.norm %>% 
  filter(sumReads >= 20) # 185,515 guide_UMI_Repl

head(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20)

A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20$to_delete_UMIs3 <- TRUE # Create a new column in the original dataframe to mark rows to be deleted
A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20 %>% group_by(Gene) %>% filter(n() >= 3) %>% ungroup()
A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20$to_delete_UMIs3[A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20$Gene %in% A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20$Gene] <- FALSE # Mark the rows to be kept as FALSE in the original dataframe
Ratio.deleted_rows.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20[
  A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20$to_delete_UMIs3, ] # Identify rows that were deleted  # 1406 guide_UMI_Repl
A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20 <- A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20[, !names(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20) %in% "to_delete_UMIs3"] # Remove the temporary column from the original dataframe
A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20.UMIs3[, !names(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20.UMIs3) 
                                                                                                   %in% "to_delete_UMIs3"] # Remove the temporary column from the original dataframe  # 184,109 guide_UMI_Repl
head(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20.UMIs3)
nrow(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20.UMIs3)
nrow(A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3)
# Vivo # ---------- StAR - Pre-Mageck ------------------------------------ 
Pre_Mageck_A.B.E_F2.guide_UMI_Repl.SR20 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20 %>% select(guide_UMI_Repl, Gene, active, inactive)
Pre_Mageck_A.B.E_F2.guide_UMI_Repl.SR20 <- Pre_Mageck_A.B.E_F2.guide_UMI_Repl.SR20 %>% mutate_at(vars(active, inactive), ~ . + 0.5)
#write.table(Pre_Mageck_A.B.E_F2.guide_UMI_Repl.SR20, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide_UMI_Repl.SR20.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Pre_Mageck_A.B.E_F2.guide_UMI_Repl.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide_UMI_Repl.SR20.tsv', header = T)

#sbatch mageck test -k Pre_Mageck_A.B.E_F2.guide_UMI_Repl.SR20.tsv -t active -c inactive -n Mageck_A.B.E_FC2.SR20.vivo.UMI.level.standard --remove-zero any 

# Vivo # ---------- StAR - Pre-Mageck - after normalization------------------------------------ 
Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20 <- A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20 %>% select(guide_UMI_Repl, Gene, active, inactive)
Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20 <- Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20 %>% mutate_at(vars(active, inactive), ~ . + 0.5)
#write.table(Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20.tsv', header = T)

#sbatch mageck test -k Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20.tsv -t active -c inactive -n Mageck_A.B.E_FC2.norm.SR20.vivo.UMI.level.standard --remove-zero any 

head(Pre_Mageck_A.B.E_F2.guide_UMI_Repl.norm.SR20)

# Vivo # ---------- StAR - Mageck ------------------------------------ 
mageck_A.B.E.FC2_vivo_genes.standard.SR20_old <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.SR20.vivo.UMI.level.standard.gene_summary.txt', header = T)
#mageck_A.B.E.FC2_vivo_guides.standard.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.SR20.vivo.UMI.level.standard.sgrna_summary.txt', header = T)
mageck_A.B.E.FC2_vivo_genes.standard.SR20_old$minRRA.fdr <- -log10(apply(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.B.E.FC2_vivo_genes.standard.SR20_old$minRRA.score <- -log10(apply(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old[,c("neg.score","pos.score")], 1, min))
mageck_A.B.E.FC2_vivo_genes.standard.SR20_old$minpvalue <- -log10(apply(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.B.E.FC2_vivo_genes.standard.SR20_old$gene <- str_split_fixed(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old$id, '_', 2)[,1] # 20,896 genes
mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3 <- mageck_A.B.E.FC2_vivo_genes.standard.SR20_old %>% filter(num >=3) #19,907 genes
nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3)
length(unique(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3$id)) # 20490 genes
length(unique(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old$id)) # 21309 genes

mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3$adj.minRRA.score <- ifelse(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3$minRRA.score > 15, 15, 
                                                                           mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3$minRRA.score)
mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3 <-  mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3 %>%
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

#write.table(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables//mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3.txt', header = T)

# Vivo # ---------- StAR - Mageck - after normalization ------------------------------------ 

#from here the dataframes are labeled the same as without normalization to keep it easier to rerun certain parts of the script

mageck_A.B.E.FC2_vivo_genes.standard.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.norm.SR20.vivo.UMI.level.standard.gene_summary.txt', header = T)
mageck_A.B.E.FC2_vivo_guides.standard.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.norm.SR20.vivo.UMI.level.standard.sgrna_summary.txt', header = T)
mageck_A.B.E.FC2_vivo_genes.standard.SR20$minRRA.fdr <- -log10(apply(mageck_A.B.E.FC2_vivo_genes.standard.SR20[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.B.E.FC2_vivo_genes.standard.SR20$minRRA.score <- -log10(apply(mageck_A.B.E.FC2_vivo_genes.standard.SR20[,c("neg.score","pos.score")], 1, min))
mageck_A.B.E.FC2_vivo_genes.standard.SR20$minpvalue <- -log10(apply(mageck_A.B.E.FC2_vivo_genes.standard.SR20[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.B.E.FC2_vivo_genes.standard.SR20$gene <- str_split_fixed(mageck_A.B.E.FC2_vivo_genes.standard.SR20$id, '_', 2)[,1] # 20,896 genes
mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 <- mageck_A.B.E.FC2_vivo_genes.standard.SR20 %>% filter(num >=3) #19,907 genes
nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3)
length(unique(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3$id)) # 20490 genes
length(unique(mageck_A.B.E.FC2_vivo_genes.standard.SR20$id)) # 21309 genes

mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3$adj.minRRA.score <- ifelse(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3$minRRA.score > 15, 15, 
                                                                           mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3$minRRA.score)
mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 <-  mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

#write.table(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables//mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3.txt', header = T)


# VIVO # ---------- Number of guides or UMIS per gene ------------------
mageck_A.B.E.FC2_vivo_genes.standard.SR20$adj.num.vivo <- ifelse(mageck_A.B.E.FC2_vivo_genes.standard.SR20$num > 50, 50, 
                                                                 mageck_A.B.E.FC2_vivo_genes.standard.SR20$num) 
mageck_A.B.E.FC2_vivo_genes.standard.SR20 %>%
  ggplot(aes(x = num)) +
  geom_histogram(binwidth = 1, color = "black", fill= "grey") +
  theme_pubr()+
  xlim(0.5, 40) +
  xlab('Number of guide_UMI_Repl per gene') + ylab('Count')

mean(mageck_A.B.E.FC2_vivo_genes.standard.SR20$num)

A.B.E_FC2.pooled.guide_UMI_Repl.SR20 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20 %>% mutate(UMI_Repl = paste(UMI, Replicate, sep = '_'))
ABE.screen_unique_umi_counts <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20 %>%
  group_by(guide) %>%
  summarise(Unique_UMI_Count = n_distinct(UMI_Repl))
ABE.screen_unique_umi_counts$adj.num.vivo <- ifelse(ABE.screen_unique_umi_counts$Unique_UMI_Count > 15, 15, 
                                                    ABE.screen_unique_umi_counts$Unique_UMI_Count) 

ABE.screen_unique_umi_counts %>%
  ggplot(aes(x = adj.num.vivo)) +
  geom_histogram(binwidth = 1, color = "black", fill= "grey") +
  theme_pubr()+
  #xlim(-0.5, 50.5) +
  xlab('Number of UMI_Repl per guide') + ylab('Count')

ABE.screen_unique_umi_counts %>%
  ggplot(aes(x = adj.num.vivo)) +
  geom_histogram(binwidth = 1, color = "black", fill= "grey") +
  theme_pubr()+
  xlim(0.5, 16) +
  xlab('Number of UMI_Repl per guide') + ylab('Count')

mean(ABE.screen_unique_umi_counts$Unique_UMI_Count)


# Vivo # ---------- Number of guides or UMIS per replicate / take rate ---------- 

replicate_counts <- A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20 %>%
  group_by(Replicate) %>%
  summarise(count = n())
#write.table(replicate_counts, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/replicate_counts.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

Vivo_replicate_mice_numbers <- read.table("/Volumes/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Vivo_replicate_mice_numbers.txt", header = TRUE)

Vivo_replicate_counts <- replicate_counts %>%
  left_join(Vivo_replicate_mice_numbers, by = "Replicate")  %>% mutate(TakeRate = (count/Mice))

Total.mice <- sum(Vivo_replicate_counts$Mice)
Total.UMIs <- sum(Vivo_replicate_counts$count)
Total.average.TR <- Total.UMIs/Total.mice
Total.median.TR <- median(Vivo_replicate_counts$TakeRate)

head(Vivo_replicate_counts)

# ggplot(Vivo_replicate_counts, aes(x = Replicate.nr, y = TakeRate, group = Batch)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(title = "Ratio active vs inactive reads",
#        x = "Batch",
#        y = "Engrafted UMIs") +
#   theme_pubr()


ggplot(Vivo_replicate_counts, aes(x = factor(Replicate.nr), y = TakeRate, fill = factor(Batch))) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_grid(~ Batch, scales = "free_x", space = "free_x") +  # Separate plots for each Batch
  labs(title = "TakeRate by Replicate Number and Batch",
       x = "Tumor sample",
       y = "Engrafted UMIs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Grey50", "grey50"), 
                    name = "Batch",
                    labels = c("Batch 1", "Batch 2")) + # Customize legend labels
  theme_pubr()

# Vivo # ---------- Number of unique UMIS  ---------- 

head(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20)
A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm
num.uniqueUMIs <- length(unique(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20$UMI))
uniqueUMIs <- unique(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20$UMI)
umi_counts <- A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20 %>%
  count(UMI, name = "Frequency")

head(A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm)
head(uniqueUMIs)

umi_counts %>%
  ggplot(aes(x = Frequency)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey") +
  scale_y_log10() +  # Apply logarithmic scale to y-axis
  #scale_y_continuous( breaks = scales::pretty_breaks(n = 3), expand = expansion(mult = c(0, 0.05)) ) +
  theme_pubr() +
  xlim(0, 30) +
  xlab('Frequency of each UMI') + ylab('Log10(Number of UMIs)') 

#xlim(-0.5, max(umi_counts$Frequency) + 0.5)

sum(umi_counts$Frequency == 10)



# Vivo # ---------- Number of sgRNAs from whole library  ---------- 

head(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20)
unique.sgRNAs.vivo <- length(unique(A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20$guide))

mageck_A.B.E.FC2_vivo_genes.standard.SR20
head(combined_data_pooled)
unique.sgRNAs.plasmid <- length(unique(combined_data_pooled$guide_name))


#
#Vivo.Vitro # ---------- RATIO PLOT ---------- 
A.B.E_FC2.average.ratio <- read.table("/Volumes/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2.average.ratio.txt", header = TRUE)
A.B.E_FC2.average.ratio.total.inactive <- A.B.E_FC2.average.ratio[A.B.E_FC2.average.ratio$Variable %in% c("Perc.inactive", "Perc.total"), ]
A.B.E_FC2.average.ratio.total.inactive<- A.B.E_FC2.average.ratio.total.inactive[order(A.B.E_FC2.average.ratio.total.inactive$Variable, decreasing = TRUE), ]

head(A.B.E_FC2.average.ratio)
ggplot(A.B.E_FC2.average.ratio.total.inactive, aes(x = vivo.vitro, y = value, fill = Variable, group = Batch)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Ratio active vs inactive reads",
       x = "Batch",
       y = "Percentage",
       fill = "Variable") +
  scale_fill_manual(name = NULL, values = c("Perc.total" = "#2F3093", "Perc.inactive" = "#DB1B5C"),
                    labels = c( "Perc.total" = "Percentage active reads", "Perc.inactive" = "Percentage inactive reads")) +
  theme_pubr()

A.B.E_FC2.median.ratio <- read.table("/Volumes/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2.median.ratio.txt", header = TRUE)
A.B.E_FC2.median.ratio.total.inactive <- A.B.E_FC2.median.ratio[A.B.E_FC2.median.ratio$Variable %in% c("Perc.inactive", "Perc.total"), ]
A.B.E_FC2.median.ratio.total.inactive<- A.B.E_FC2.median.ratio.total.inactive[order(A.B.E_FC2.median.ratio.total.inactive$Variable, decreasing = TRUE), ]

head(A.B.E_FC2.median.ratio)
ggplot(A.B.E_FC2.median.ratio.total.inactive, aes(x = vivo.vitro, y = value, fill = Variable, group = Batch)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Ratio active vs inactive reads",
       x = "Batch",
       y = "Percentage",
       fill = "Variable") +
  scale_fill_manual(name = NULL, values = c("Perc.total" = "#2F3093", "Perc.inactive" = "#DB1B5C"),
                    labels = c( "Perc.total" = "Percentage active reads", "Perc.inactive" = "Percentage inactive reads")) +
  theme_pubr()

# Control-Essential lists # ---------- Nr. genes that match ------------------

nrow(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id)) #566 genes that match
nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id)) #526 genes that match
nrow(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard %>% filter(toupper(gene) %in% Controls$gene)) #486 genes that match
nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(toupper(gene) %in% Controls$gene)) #469 genes that match
nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter((gene) %in% Mouse.ortholog_CoreEssentials_934$gene)) #844 genes that match


# VIVO # ---------- StAR Volcano plots - minRRA.score ------------------
colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(toupper(gene) %in% Controls$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo gene LFC + RRA score \ncut off SR20.UMIs3') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,38) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')

nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id)) #529 genes
nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene)) # 498 genes
#filter(toupper(gene) %in% Controls$gene)

#!!! VIVO # ---------- StAR Volcano plots - minRRA.score ------------------ RRA.score max 15 ------------------

colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = adj.minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo gene LFC + RRA score \ncut off SR20.UMIs3') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,15) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')

mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = adj.minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  geom_text_repel(aes(label=ifelse((neg.lfc < -3 & adj.minRRA.score > 6 & group == 'Other') | (neg.lfc > 2 & adj.minRRA.score > 6 & group == 'Other'), as.character(id),'')),box.padding = 2, max.overlaps = Inf,
                  segment.color = "blue", segment.size = 0.1, size = 3) +
  
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo gene LFC + RRA score \ncut off SR20.UMIs3') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,15) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')
#geom_point(aes(color = 'Blue'), data = mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id == "Mapk1"), size = 2) +

toLabel <- mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter((neg.lfc < -3 & adj.minRRA.score > 5 & group == 'Other') | (neg.lfc > 2 & adj.minRRA.score > 5 & group == 'Other')) #121 genes
toLabel2 <- mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter((neg.lfc < -3 & adj.minRRA.score > 6 & group == 'Other') | (neg.lfc > 2 & adj.minRRA.score > 6 & group == 'Other')) #82 genes

#!!! VIVO # ---------- StAR Volcano plots - minRRA.score ------------------ RRA.score max 15 ------------------ With density plots ---------------

VPlot1 <- 
  mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = adj.minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'),  size = 1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter((gene) %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ 
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,15) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")  #+ #guides(colour = 'none')

Plot.Dens.V1x <- mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  xlab(NULL)  +
  xlim(-7.5,7.5) +  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  theme(legend.position = "none", axis.text.x = element_blank())

Plot.Dens.V1y <- mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = adj.minRRA.score)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  #geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), alpha=0.5, linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  xlab(NULL) +
  xlim(0,15) + theme_pubr(base_size = 12) + scale_color_manual(values = colors) + 
  theme(legend.position = "none", axis.text.y = element_blank()) + coord_flip()

Plot.Dens.V1x + plot_spacer() + VPlot1 + theme(aspect.ratio = 1)+ Plot.Dens.V1y +
  plot_layout(ncol = 2, nrow = 2, widths = c(2, 0.5), heights = c(0.5, 2), guides = 'collect') & theme(legend.position = 'top')

#
# VIVO # ---------- StAR Volcano plots - minRRA.score ------------------ RRA.score max 15 & highlighting groups ------------------

colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = adj.minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo gene LFC + RRA score \ncut off SR20.UMIs3') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,15) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')

A.B.E.FC2_vivo_spec.depleting.vivo1xRRA3
colors5 <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black', 'A.B.E.FC2_vivo_spec.depleting.vivo1x' = 'Gold')
mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = adj.minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  geom_point(aes(color = 'A.B.E.FC2_vivo_spec.depleting.vivo1x'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% A.B.E.FC2_vivo_spec.depleting.vivo1x$id),  size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo gene LFC + RRA score \ncut off SR20.UMIs3') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors5) +
  xlim(-7.5,7.5) + ylim(0,15) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')

Vitro.A.B.E.FC2_depleting.min3.RRAsc10
# VIVO # ---------- StAR Volcano plots - minRRA.score, minRRA.fdr, min.pvalue plot in VIVO ------------------
colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
plot_VV1 <- mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(toupper(gene) %in% Controls$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo gene LFC + RRA score \ncut off SR20.UMIs3') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,38) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')
plot_VV2 <- mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = minRRA.fdr)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(toupper(gene) %in% Controls$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA fdr'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo gene LFC + RRA fdr \ncut off SR20.UMIs3') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,38) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')
plot_VV3 <- mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc, y = minpvalue)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(toupper(gene) %in% Controls$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('pValue'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo gene LFC + pValue \ncut off SR20.UMIs3') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-7.5,7.5) + ylim(0,38) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')
grid.arrange(plot_VV1, plot_VV2, plot_VV3, ncol = 3, top = textGrob('In vivo LFC + significance scores \nVitro.A.B.E.FC2_depleting.min3.RRAsc10 = Essential group (Red)', gp = gpar(fontsize = 14)))

# Vivo_old vs Vivo.norm # ---------- Correlation ---------- 
mageck_A.B.E.FC2_invivo.old_norm_genes.standard <- merge(mageck_A.B.E.FC2_vivo_genes.standard.SR20_old.UMIs3 [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
                                                     mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
                                                     by = 'id', suffixes = c(".invivo_old",".invivo"))
head(mageck_A.B.E.FC2_invivo.old_norm_genes.standard)
colors4 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'Depleting validation library' = 'orange', 
             'Enriching validation library' = 'limegreen', 'mito_inner_membr_GO_0005743' = 'Blue', 'Vivo.validation.enriching.vivo1x' = 'magenta3',
             'MicroLibrary' = 'cyan')
mageck_A.B.E.FC2_invivo.old_norm_genes.standard %>%
  ggplot(aes(x = neg.lfc.invivo_old, y = neg.lfc.invivo)) +
  geom_hline(yintercept=0, color="gray20", linewidth = 0.5) +
  geom_vline(xintercept=0, color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other')) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivo.old_norm_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivo.old_norm_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) +
  xlab(bquote('Log'[2]~'fold change in vivo before recalculation')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ') + 
  theme_pubr(base_size = 14) + scale_color_manual(values = colors4)  +
  xlim(-8.2,5.2) + 
  ylim(-8.2,5.2) +
  theme(plot.title = element_text(hjust = 0.5))


#
# Vitro.Vivo # ---------- Correlation ---------- 
mageck_A.B.E.FC2_invivovitro_genes.standard <- merge(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
                                                     mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
                                                     by = 'id', suffixes = c(".invitro",".invivo"))

mageck_A.B.E.FC2_invivovitro_genes.standard.noUMIrule <- merge(mageck_A.B.E.FC2_invitro_genes_pseudo.paired.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
                                                               mageck_A.B.E.FC2_vivo_genes.standard.SR20 [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
                                                               by = 'id', suffixes = c(".invitro",".invivo"))

#write.table(mageck_A.B.E.FC2_invivovitro_genes.standard, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_invivovitro_genes.standard.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(mageck_A.B.E.FC2_invivovitro_genes.standard.noUMIrule, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_invivovitro_genes.standard.noUMIrule.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#mageck_A.B.E.FC2_invivovitro_genes.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_invivovitro_genes.standard.txt', header = T)
names(mageck_A.B.E.FC2_invivovitro_genes.standard)

mageck_A.B.E.FC2_invivovitro_genes.standard <- mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

mageck_A.B.E.FC2_invivovitro_genes.standard$dot_size_RRAvivo <- ifelse(mageck_A.B.E.FC2_invivovitro_genes.standard$minRRA.score.invivo > 15, 15, mageck_A.B.E.FC2_invivovitro_genes.standard$minRRA.score.invivo)

#write.table(mageck_A.B.E.FC2_invivovitro_genes.standard, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_invivovitro_genes.standard.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#mageck_A.B.E.FC2_invivovitro_genes.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_invivovitro_genes.standard.txt', header = T)

common.vivo.vitro.CoreEssential.screen.ABE <- mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(neg.lfc.invivo < -2 & neg.lfc.invitro < -2)
#write.table(common.vivo.vitro.CoreEssential.screen.ABE,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/common.vivo.vitro.CoreEssential.screen.ABE.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#common.vivo.vitro.CoreEssential.screen.ABE <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/common.vivo.vitro.CoreEssential.screen.ABE.txt', header = T)


colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  theme(plot.title = element_text(hjust = 0.5))

nrow(mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id)) #526 genes that match
nrow(mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(toupper(id) %in% Controls$gene)) #469 genes that match
nrow(mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) #472 genes that match
nrow(mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) #844 genes that match

colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene), size=1) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors2) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  theme(plot.title = element_text(hjust = 0.5))

# Vitro.Vivo # ---------- Correlation ---------- Bubble-size RRA score--------
library(scales)  # Load the scales library for breaks_width
desired_breaks <- 5  # Change this to the desired number of breaks

colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo, size = dot_size_RRAvivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other')) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ', 
                                                                                                     title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene', size = 'RRA score') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors2) + scale_size_continuous(range = c(0.2, 4), breaks = seq(0, 15, by = 3))  +
  xlim(-7.5,5) + ylim(-7.5,5) +
  theme(plot.title = element_text(hjust = 0.5))

#!!!! Vitro.Vivo # ---------- Correlation ---------- Bubble-size RRA score + GO-terms --------
mito_inner_membr_GO_0005743 <- read.table('/Volumes/groups/elling/Esther/NGS_results/Scripts_gene_lists/GOterm_Gene_lists/mito_inner_membr_GO_0005743.txt', header = T, sep = '\t', quote = "")
library(scales)  # Load the scales library for breaks_width
desired_breaks <- 5  # Change this to the desired number of breaks

colors4 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'Depleting validation library' = 'orange', 
             'Enriching validation library' = 'limegreen', 'mito_inner_membr_GO_0005743' = 'Blue', 'Vivo.validation.enriching.vivo1x' = 'magenta3',
             'MicroLibrary' = 'cyan')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo, size = dot_size_RRAvivo)) +
  geom_hline(yintercept=0, color="gray20", linewidth = 0.5) +
  geom_vline(xintercept=0, color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other')) +
  # geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) +
  geom_point(aes(color = 'mito_inner_membr_GO_0005743'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% mito_inner_membr_GO_0005743$gene & minRRA.score.invivo > 2)) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) +
  geom_hline(yintercept=-1, linetype = 'dashed', color="black", linewidth = 0.5) +
  geom_vline(xintercept=-1, linetype = 'dashed', color="black", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dotdash', color="black", linewidth = 0.5) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ') + 
  #labs(title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene', size = 'RRA score') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors4) + scale_size_continuous(range = c(0.2, 4), breaks = seq(0, 15, by = 3))  +
  guides(color = guide_legend(title = " ", nrow = 2), size = guide_legend(title = "RRAscore in vivo", title.position = "top")) +
  xlim(-8.2,5.2) + 
  ylim(-8.2,5.2) +
  theme(plot.title = element_text(hjust = 0.5))

#write.table(mageck_A.B.E.FC2_invivovitro_genes.standard, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_invivovitro_genes.standard.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(mageck_A.B.E.FC2_invivovitro_genes.standard, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_invivovitro_genes.standard.norm.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Vitro.Vivo # ---------- Correlation ---------- Vivo depleting & mitochondrial & oxphos highlighted --------
Oxphos <- read.table('/Volumes/groups/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Gene lists/Oxphos pool A (KEGG).txt', header=T, sep = '\t')
GO_0005753_mito.proton.trans.ATP <- read.table('/Volumes/groups/elling/Esther/NGS_results/Scripts_gene_lists/GOterm_Gene_lists/GO_0005753_mito.proton.trans.ATP.txt', header=T, sep = '\t')
GO_TCA <- read.table('/Volumes/groups/elling/Esther/NGS_results/Scripts_gene_lists/GOterm_Gene_lists/tricarboxylic acid cycle TCA GO_0006099.txt', header=F, sep = '\t')
GO_TCA <- GO_TCA %>% distinct(V2) %>% rename(gene=V2)
nrow(Oxphos)
colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black')
colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'GO_0005753_mito.proton.trans.ATP' = 'Blue')
colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'Oxphos' = 'Blue')
colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'TCA' = 'Blue')

mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$id), alpha=0.5, size=1) +
  #geom_point(aes(color = 'GO_0005753_mito.proton.trans.ATP'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% GO_0005753_mito.proton.trans.ATP$Gene), alpha=0.5, size=1.8) +
  #geom_point(aes(color = 'Oxphos'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Oxphos$gene), alpha=0.5, size=1.8) +
  geom_point(aes(color = 'TCA'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% GO_TCA$gene), alpha=0.5, size=1.8) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors2) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  theme(plot.title = element_text(hjust = 0.5))

CoreEssentials.ABE.screen <- mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)
#CoreEssentials.ABE.screen.498 <- CoreEssentials.ABE.screen[sample(nrow(CoreEssentials.ABE.screen), 498), "id", drop = FALSE]
nrow(CoreEssentials.ABE.screen) #844 genes
nrow(mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene))#498 genes
nrow(mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) #844 genes
##

# Vitro.Vivo # ---------- Correlation ---------- All genes that show movement highlighted --------

mageck_A.B.E.FC2_invivovitro_genes.standard
A.B.E.FC2_vitro.vivo.moving.genes.1.5 <- mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(neg.lfc.invitro < -1.5 | neg.lfc.invivo < -1.5) #2387 genes
A.B.E.FC2_vitro.vivo.moving.genes.1 <- mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(neg.lfc.invitro < -1 | neg.lfc.invivo < -1) #3545 genes
nrow(A.B.E.FC2_vitro.vivo.moving.genes.1)

#write.table(A.B.E.FC2_vitro.vivo.moving.genes.1.5$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vitro.vivo.moving.genes.1.5_genenames.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(A.B.E.FC2_vitro.vivo.moving.genes.1$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vitro.vivo.moving.genes.1_genenames.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

colors3 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'All geness with movement' = 'Aquamarine')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$id), size=1) +
  geom_point(aes(color = 'All geness with movement'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% A.B.E.FC2_vitro.vivo.moving.genes.1$id), size=1) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors3) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  theme(plot.title = element_text(hjust = 0.5))

# Vitro.Vivo # ---------- Correlation ---------- In VIVO specific depleting highlighted --------
A.B.E.FC2_vivo_spec.depleting1 <- mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(neg.lfc.invivo < -1) # <-2 = 1490 genes, <-1.5 = 2094 genes, <-1 = 3227 genes
A.B.E.FC2_vivo_spec.depleting1$vivo1.3x <- ifelse(A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1.3 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro), TRUE, "x") #2313 genes (<-2), 468 genes (<-1.5)
A.B.E.FC2_vivo_spec.depleting.vivo1.3 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1.3x == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting.vivo1.3$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting.vivo1.3_genenames.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo1.5x <- ifelse(A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1.5 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro), TRUE, "x") #2121 genes (<-2), 414 genes (<-1.5)
A.B.E.FC2_vivo_spec.depleting.vivo1.5 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1.5x == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting.vivo1.5$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting.vivo1.5_genenames.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo1x <- ifelse(A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo < (1 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro), TRUE, "x") #2696 genes (<-2), 414 genes (<-1.5)
A.B.E.FC2_vivo_spec.depleting.vivo1x <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1x == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting.vivo1x$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting.vivo1x_genenames.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting2 <- mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(neg.lfc.invivo < -2) # <-2 = 1490 genes, <-1.5 = 2094 genes, <-1 = 3227 genes
nrow(A.B.E.FC2_vivo_spec.depleting2 %>% filter(vivo1.3xRRA2 == TRUE))

A.B.E.FC2_vivo_spec.depleting1$vivo1.5xRRA5 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1.5 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                      & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 5, TRUE, "x") 
A.B.E.FC2_vivo_spec.depleting.vivo1.5xRRA5 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1.5xRRA5 == TRUE)

A.B.E.FC2_vivo_spec.depleting1$vivo1xRRA5 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                    & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 5, TRUE, "x") 
A.B.E.FC2_vivo_spec.depleting1$vivo1xRRA3 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                    & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 3, TRUE, "x")
A.B.E.FC2_vivo_spec.depleting.vivo1xRRA3 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1xRRA3 == TRUE)

A.B.E.FC2_vivo_spec.depleting2$vivo1.3xRRA2 <- ifelse((A.B.E.FC2_vivo_spec.depleting2$neg.lfc.invivo <= (1.3 * A.B.E.FC2_vivo_spec.depleting2$neg.lfc.invitro)) 
                                                      & A.B.E.FC2_vivo_spec.depleting2$minRRA.score.invivo > 2, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting2.vivo1.3xRRA2 <- A.B.E.FC2_vivo_spec.depleting2 %>% filter(vivo1.3xRRA2 == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting2.vivo1.3xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting2.vivo1.3xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo1.3xRRA2 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1.3 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                      & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 2, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting1.vivo1.3xRRA2 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1.3xRRA2 == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting1.vivo1.3xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting1.vivo1.3xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo1.5xRRA2 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1.5 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                      & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 2, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting.vivo1.5xRRA2 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1.5xRRA2 == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting.vivo1.5xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting.vivo1.5xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo1.75xRRA2 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1.75 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                       & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 2, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting1.vivo1.75xRRA2 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1.75xRRA2 == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting1.vivo1.75xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting1.vivo1.75xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo1.5xRRA3 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1.5 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                      & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 3, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting1$vivo2xRRA2 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (2 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                    & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 2, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting1.vivo2xRRA2 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo2xRRA2 == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting1.vivo2xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting1.vivo2xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo2xRRA2 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (2 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro)) 
                                                    & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 2, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting1.vivo2xRRA2 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo2xRRA2 == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting1.vivo2xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting1.vivo2xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo1x.interc1.RRA2 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro) - 1) 
                                                             & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 2, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting1.vivo1x.interc1.RRA2 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1x.interc1.RRA2 == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting1.vivo1x.interc1.RRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting1.vivo1x.interc1.RRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_vivo_spec.depleting1$vivo1x.interc0.75.RRA2 <- ifelse((A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invivo <= (1 * A.B.E.FC2_vivo_spec.depleting1$neg.lfc.invitro) - 0.75) 
                                                                & A.B.E.FC2_vivo_spec.depleting1$minRRA.score.invivo > 2, TRUE, "x") #663 genes
A.B.E.FC2_vivo_spec.depleting1.vivo1x.interc0.75.RRA2 <- A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1x.interc0.75.RRA2 == TRUE)
#write.table(A.B.E.FC2_vivo_spec.depleting1.vivo1x.interc0.75.RRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_vivo_spec.depleting1.vivo1x.interc0.75.RRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

colors4 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'LFC In vivo <1x in vitro' = 'Gold', 
             'Validating A & B' = 'Blue', 'extra' = 'magenta', 'mito_inner_membr_GO_0005743' = 'Blue')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo, size = dot_size_RRAvivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other')) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) +
  geom_point(aes(color = 'LFC In vivo <1x in vitro'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% A.B.E.FC2_vivo_spec.depleting.vivo1x$id)) +
  geom_point(aes(color = 'extra'), data=A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo2xRRA2 == TRUE)) +
  geom_point(aes(color = 'mito_inner_membr_GO_0005743'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% mito_inner_membr_GO_0005743$gene & minRRA.score.invivo > 2)) +
  # geom_text_repel(data = mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% GOterm.interest), aes(label = id), box.padding = 1, max.overlaps = Inf,
  #                 segment.color = "black", segment.size = 0.2,  size = 4) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ') +
  #labs(title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors4) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  #annotate("text", x = 3.5, y = -6, label = paste("LFC.vivo < 1.5x LFC.vitro & \nRRAscore >3"), color = "black") +
  guides(color = guide_legend(title = " ", nrow = 2), size = guide_legend(title = "RRAscore in vivo", title.position = "top")) +
  #annotate("text", x = 3.5, y = -6, label = paste("LFC.vivo < 1x LFC.vitro"), color = "black") +
  theme(plot.title = element_text(hjust = 0.5))

GOterm.interest <- c('Terf1', 'Ago1', 'Ago2', 'Ago3', 'Ago4', 'Drosha', 'Tcf3', 'Sox9')

colors4 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'In vivo depleting' = 'Gold')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene), size=1) +
  geom_point(aes(color = 'In vivo depleting'), data=A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1.75xRRA2 == TRUE),size=1) +
  #geom_point(aes(color = 'In vivo depleting'), data=A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1xRRA3 == TRUE), size=1) +
  # geom_text_repel(data = A.B.E.FC2_vivo_spec.depleting %>% filter(vivo1.5xRRA5 == TRUE), aes(label = id), box.padding = 1, max.overlaps = Inf,
  #                 segment.color = "gray20", segment.size = 0.1,  size = 3) +  
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors4) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  annotate("text", x = 3.5, y = -6, label = paste("LFC.vivo < 1x LFC.vitro & \nRRAscore >3"), color = "black") +
  theme(plot.title = element_text(hjust = 0.5))

Grey.validating.AB <- mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(neg.lfc.A <= -2 & neg.lfc.B <= -2 & group == 'Other') 

#
# Vitro.Vivo # ---------- Correlation ---------- In VITRO specific depleting highlighted --------
A.B.E.FC2_VITRO_spec.depleting1 <- mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(neg.lfc.invitro < -1) # <-2 = 1070 genes, <-1.5 = 1388 genes, <-1 = 1754 genes
A.B.E.FC2_VITRO_spec.depleting1$vitro1x <- ifelse((A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invitro <= (1 * A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invivo)), TRUE, "x") 
A.B.E.FC2_VITRO_spec.depleting1.vitro1x <- A.B.E.FC2_VITRO_spec.depleting1 %>% filter(vitro1x == TRUE)


A.B.E.FC2_VITRO_spec.depleting1$vitro1.5xRRA2 <- ifelse((A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invitro <= (1.5 * A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invivo)) 
                                                        & A.B.E.FC2_VITRO_spec.depleting1$minRRA.score.invitro > 2, TRUE, "x") 
A.B.E.FC2_VITRO_spec.depleting1.vitro1.5xRRA2 <- A.B.E.FC2_VITRO_spec.depleting1 %>% filter(vitro1.5xRRA2 == TRUE)
#write.table(A.B.E.FC2_VITRO_spec.depleting1.vitro1.5xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_VITRO_spec.depleting1.vitro1.5xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_VITRO_spec.depleting1$vitro1.75xRRA2 <- ifelse((A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invitro <= (1.75 * A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invivo)) 
                                                         & A.B.E.FC2_VITRO_spec.depleting1$minRRA.score.invitro > 2, TRUE, "x") 
A.B.E.FC2_VITRO_spec.depleting1.vitro1.75xRRA2 <- A.B.E.FC2_VITRO_spec.depleting1 %>% filter(vitro1.75xRRA2 == TRUE)
#write.table(A.B.E.FC2_VITRO_spec.depleting1.vitro1.75xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_VITRO_spec.depleting1.vitro1.75xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_VITRO_spec.depleting1$vitro1.25xRRA2 <- ifelse((A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invitro <= (1.25 * A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invivo)) 
                                                         & A.B.E.FC2_VITRO_spec.depleting1$minRRA.score.invitro > 2, TRUE, "x") 
A.B.E.FC2_VITRO_spec.depleting1.vitro1.25xRRA2 <- A.B.E.FC2_VITRO_spec.depleting1 %>% filter(vitro1.25xRRA2 == TRUE)
#write.table(A.B.E.FC2_VITRO_spec.depleting1.vitro1.25xRRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_VITRO_spec.depleting1.vitro1.25xRRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_VITRO_spec.depleting1$vitro1x.interc1.25.RRA2 <- ifelse((A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invitro <= (1 * A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invivo) - 1.25) 
                                                                  & A.B.E.FC2_VITRO_spec.depleting1$minRRA.score.invitro > 2, TRUE, "x") 
A.B.E.FC2_VITRO_spec.depleting1.vitro1x.interc1.25.RRA2 <- A.B.E.FC2_VITRO_spec.depleting1 %>% filter(vitro1x.interc1.25.RRA2 == TRUE)
#write.table(A.B.E.FC2_VITRO_spec.depleting1.vitro1x.interc1.25.RRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_VITRO_spec.depleting1.vitro1x.interc1.25.RRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

A.B.E.FC2_VITRO_spec.depleting1$vitro1.75x.interc0.5.RRA2 <- ifelse((A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invitro <= (1.75 * A.B.E.FC2_VITRO_spec.depleting1$neg.lfc.invivo) - 0.5) 
                                                                    & A.B.E.FC2_VITRO_spec.depleting1$minRRA.score.invitro > 2, TRUE, "x") 
A.B.E.FC2_VITRO_spec.depleting1.vitro1.75x.interc0.5.RRA2 <- A.B.E.FC2_VITRO_spec.depleting1 %>% filter(vitro1.75x.interc0.5.RRA2 == TRUE)
#write.table(A.B.E.FC2_VITRO_spec.depleting1.vitro1.75x.interc0.5.RRA2$id,'/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_04_GOterm_lists/A.B.E.FC2_VITRO_spec.depleting1.vitro1.75x.interc0.5.RRA2_genenames.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

colors8 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'LFC In vitro <1x in vivo' = 'plum1', 
             'Validating A & B' = 'Blue', 'extra' = 'magenta', 'mito_inner_membr_GO_0005743' = 'Blue', 'In vitro depleting' = 'aquamarine' )
#colors4 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'In vivo depleting' = 'Gold', 'In vitro depleting' = 'Magenta')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo, size = dot_size_RRAvivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other')) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) +
  geom_point(aes(color = 'LFC In vitro <1x in vivo'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% A.B.E.FC2_VITRO_spec.depleting1.vitro1x$id)) +
  #geom_point(aes(color = 'extra'), data=A.B.E.FC2_vivo_spec.depleting1 %>% filter(vivo1x.interc0.75.RRA2 == TRUE)) +
  geom_point(aes(color = 'In vitro depleting'), data=A.B.E.FC2_VITRO_spec.depleting1 %>% filter(vitro1.75x.interc0.5.RRA2 == TRUE)) +
  geom_point(aes(color = 'mito_inner_membr_GO_0005743'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% mito_inner_membr_GO_0005743$gene & minRRA.score.invivo > 2)) +
  # geom_text_repel(data = mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% GOterm.interest), aes(label = id), box.padding = 1, max.overlaps = Inf,
  #                 segment.color = "black", segment.size = 0.2,  size = 4) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ') + 
  #labs(title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors8) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  guides(color = guide_legend(title = " ", nrow = 2), size = guide_legend(title = "RRAscore in vivo", title.position = "top")) +
  #annotate("text", x = 3.5, y = -6, label = paste("LFC.vivo < 1x LFC.vitro"), color = "black") +
  theme(plot.title = element_text(hjust = 0.5))

#
regression_data <- vitro_D14_vivo_P1_merged %>%
  filter(neg.lfc.x < 0)  # Consider only x < 0 for regression calculation

slope <- cov(mageck_A.B.E.FC2_invivovitro_genes.standard$neg.lfc.invitro, mageck_A.B.E.FC2_invivovitro_genes.standard$neg.lfc.vivo) / var(mageck_A.B.E.FC2_invivovitro_genes.standard$neg.lfc.invitro)
intercept <- mean(regression_data$neg.lfc.y) - slope * mean(regression_data$neg.lfc.x)



# Your existing ggplot code
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  
  geom_point(aes(color = 'Other')) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene)) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) +
  geom_point(aes(color = 'LFC In vitro <1x in vivo'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% A.B.E.FC2_VITRO_spec.depleting1.vitro1x$id)) +
  geom_point(aes(color = 'In vitro depleting'), data=A.B.E.FC2_VITRO_spec.depleting1 %>% filter(vitro1x.interc1.RRA2 == TRUE)) +
  geom_point(aes(color = 'mito_inner_membr_GO_0005743'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% mito_inner_membr_GO_0005743$gene & minRRA.score.invivo > 2)) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ') + 
  theme_pubr(base_size = 14) + scale_color_manual(values = colors8) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  guides(color = guide_legend(title = " ", nrow = 2), size = guide_legend(title = "RRAscore in vivo", title.position = "top")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  
  # Adding regression line and equation
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) + 
  geom_text(x = -6, y = 2.5, label = paste("y = ", round(coef(lm(neg.lfc.invivo ~ neg.lfc.invitro, data = mageck_A.B.E.FC2_invivovitro_genes.standard))[1], 2), " + ",
                                           round(coef(lm(neg.lfc.invivo ~ neg.lfc.invitro, data = mageck_A.B.E.FC2_invivovitro_genes.standard))[2], 2), "x"))


#


# Vitro.Vivo # ---------- Correlation ---------- highlighting genes ---------- 
colors1 <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black', "Blue" = 'Blue')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  # geom_text_repel(aes(label=ifelse(neg.lfc.invitro > 0 & neg.lfc.invivo > 2, as.character(id),'')),box.padding = 2, max.overlaps = Inf,
  #                 segment.color = "blue", segment.size = 0.1, size = 3) +
  geom_text_repel(aes(label=ifelse(id == "Mapk1", as.character(id),'')), max.overlaps = Inf, size = 3) +
  geom_point(aes(color = 'Blue'), data = mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id == "Mapk1"), size = 2) +
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors1) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  theme(plot.title = element_text(hjust = 0.5))

# Vitro.Vivo # ---------- Correlation ---------- highlighting genes ---------- 
Vivo.validation.depleting.vivo1x.noEss <- read.table('/Volumes/Esther/NGS_results/Screen_results/20220314_Validation_Screen_Yumm450R/Analysis_YummR_Validation/5_GOterm_lists/Vivo.validation.depleting.vivo1x.noEss.txt', sep = "\t", header = T)

colors1 <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black', "Blue" = 'Blue')
colors6 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black', 'In vivo depleting' = 'Gold', 'Vivo.validation.depleting.vivo1x.noEss' = 'Blue')
mageck_A.B.E.FC2_invivovitro_genes.standard %>%
  ggplot(aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter((id) %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene), size=1) +
  geom_point(aes(color = 'Vivo.validation.depleting.vivo1x.noEss'), data=mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Vivo.validation.depleting.vivo1x.noEss$id), size=1) +
  geom_text_repel(data = mageck_A.B.E.FC2_invivovitro_genes.standard %>% filter(id %in% Vivo.validation.depleting.vivo1x.noEss$id), aes(label = id), box.padding = 1, max.overlaps = Inf,
                  segment.color = "gray20", segment.size = 0.1,  size = 3) +  
  xlab(bquote('Log'[2]~'fold change in vitro')) + ylab(bquote('Log'[2]~'fold change in vivo'))+ labs(color = ' ', title = 'StAR analysis \nIn vivo vs In vitro gene LFC \nIn vivo 3 or more UMIs per gene') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors6) +
  xlim(-7.5,5) + ylim(-7.5,5) +
  #annotate("text", x = 3.5, y = -6, label = paste("LFC.vivo < 1x LFC.vitro"), color = "black") +
  theme(plot.title = element_text(hjust = 0.5))


#!!! dAUC # ---------- StAR ---------- SR20.UMIs3 ---------- no Others! different Essential groups ---------- 
A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20
#A.B.E_FC2.pooled.guideLevel.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20 %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
A.B.E_FC2.pooled.guideLevel.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.norm.SR20 %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3 %>% mutate_at(vars(active, inactive), ~ . + 0.5)
A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo$StAR.LFC <- log2(A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo$active / 
                                                                 A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo$inactive)

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo<- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo %>% 
  mutate(group= ifelse(Gene %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(Gene %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.noOthers <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo %>% filter(group == "Essential" | group == "Non-Essential")

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.cumulative <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.noOthers %>% 
  dplyr::select(guide, Gene, StAR.LFC, group) %>% na.omit()

A.B.E.guide.dAUC.rank.LFCStAR <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.cumulative %>% arrange(StAR.LFC)
A.B.E.guide.dAUC.rank.LFCStAR$all_rank <- 1:nrow(A.B.E.guide.dAUC.rank.LFCStAR)
A.B.E.guide.dAUC.rank.LFCStAR$all_rank_perc_sgRNA <- (A.B.E.guide.dAUC.rank.LFCStAR$all_rank*1)/nrow(A.B.E.guide.dAUC.rank.LFCStAR) 

A.B.E.guide.dAUC.rank.LFCStAR <- A.B.E.guide.dAUC.rank.LFCStAR %>% group_by(group) %>% mutate(rank_group = 1:n(), rank_group_prct = 1:n()/n(), rank_sum = cumsum(rank_group_prct)) 
#head(A.B.E.guide.dAUC.rank.LFCStAR)

# Group the data by "group"
grouped_data_StAR <- A.B.E.guide.dAUC.rank.LFCStAR %>%
  group_by(group) %>% arrange(all_rank_perc_sgRNA)
# Function to calculate AUC using the trapezoidal rule
calculate_auc_StAR <- function(x, y) {
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2) }
# Calculate AUC for each group
auc_results_StAR <- grouped_data_StAR %>%
  summarize(AUC = calculate_auc_StAR(all_rank_perc_sgRNA, rank_group_prct))
print(auc_results_StAR)

Essential_auc_StAR <- auc_results_StAR$AUC[auc_results_StAR$group == "Essential"]
NonEssential_auc_StAR <- auc_results_StAR$AUC[auc_results_StAR$group == "Non-Essential"]
dAUC_Ess_nonEss_StAR <- Essential_auc_StAR - NonEssential_auc_StAR
# Print the dAUC value
cat("dAUC SR20.UMIs3 =", dAUC_Ess_nonEss_StAR, "\n")

ggplot(A.B.E.guide.dAUC.rank.LFCStAR, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group)) + 
  geom_line(size = 0.75) + 
  geom_hline(yintercept = max(A.B.E.guide.dAUC.rank.LFCStAR$rank_group_prct), linetype = "solid", color = "black") +
  geom_vline(xintercept = max(A.B.E.guide.dAUC.rank.LFCStAR$all_rank_perc_sgRNA), linetype = "solid", color = "black") +
  theme_pubr() + 
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "black", "Other" = "gray")) +
  labs(title = "StAR AUC - SR20.UMIs3 \nEssentials: Vitro.depleting.min3.RRAsc10 \nControls: Mouse.ortholog_Controls_550", y = "Cumulative fraction", x = "Percentage rank of sgRNAs") + 
  annotate("text", x = 0.8, y = 0.1, label = paste("dAUC =", round(dAUC_Ess_nonEss_StAR, 2)), color = "black")

# VIVO # rank # ---------- StAR ---------- ---------- no Others! different Essential groups ---------- 
# A.B.E_FC2.pooled.guideLevel.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20 %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
# A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3 %>% mutate_at(vars(active, inactive), ~ . + 0.5)
# A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo$StAR.LFC <- log2(A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo$active / 
#                                                                  A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo$inactive)
# 
# A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo<- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo %>% 
#   mutate(group= ifelse(Gene %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(Gene %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo %>% arrange(StAR.LFC)
head(A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank)

gene_counts <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank %>% group_by(Gene) %>% summarize(count = n())
A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank %>%
  left_join(gene_counts, by = "Gene") %>% filter(count <= 5) %>% select(-count)

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5$Rank <- 1:nrow(A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5)
A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 %>% select(-active, -inactive)
head(A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5)


A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 %>%
  mutate(Unique_Genes_Count = cumsum(!duplicated(Gene)))

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 %>% mutate(average_guides = Rank / Unique_Genes_Count)

# A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 %>%
#   arrange(Rank) %>%
#   group_by(Gene) %>%
#   mutate(sgRNA.gene = row_number()) %>%
#   ungroup()


ggplot(A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5, aes(x = Rank, y = average_guides)) +
  geom_point() + 
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene") +
  #xlim(0,5000) +
  scale_x_continuous(labels = label_scientific()) +
  labs(title = 'StAR') +
  theme_pubr()



A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 %>% filter(group == "Essential") %>%
  distinct(Gene) %>%
  nrow()

#

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5 %>%  mutate(random_order = runif(nrow(.))) %>%  # create a column of random numbers
  arrange(random_order) %>%  # sort the dataframe by the random column
  select(-random_order)  # remove the random column

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random$Rank <- 1:nrow(A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random)


A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random %>%
  mutate(Unique_Genes_Count = cumsum(!duplicated(Gene)))

A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random %>% mutate(average_guides = Rank / Unique_Genes_Count)


A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5_select <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5[c("Rank", "average_guides")]
A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random_select <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random[c("Rank", "average_guides")]

combined.vivo.StAR <- merge(A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5_select, A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.random_select, by = "Rank", suffixes = c('.vvStAR', '.vvStAR.rdm'))

colnames(combined.vivo.StAR) 

combined_long_vivo.StAR <- reshape2::melt(combined.vivo.StAR, id.vars = "Rank", variable.name = "Method", value.name = "Avg_guides")

ggplot(combined_long_vivo.StAR, aes(x = Rank, y = Avg_guides, color = Method)) +
  geom_point(data = combined_long_vivo.StAR[combined_long_vivo.StAR$Method == "average_guides.vvStAR.rdm",], aes(x = Rank, y = Avg_guides)) +
  geom_point(data = combined_long_vivo.StAR[combined_long_vivo.StAR$Method == "average_guides.vvStAR",], aes(x = Rank, y = Avg_guides)) +
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene",
       title = 'Reproducibility of sgRNAs') +
  xlim(0,3e+03)+
  scale_x_continuous(labels = scales::label_scientific()) +
  scale_color_manual(values =c('average_guides.vvStAR' = "#000000", 'average_guides.vvStAR.rdm' = "#999999")) +
  theme_pubr()

#


#!!! Vivo AB overlap# ---------- Correlation of A & B OVERLAP genes ---------- SR20.UMIs3 only B.E_FC2 ---------- MAGECK ---------- 

#after normalization
head(A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm)

Reads.i.a.vivo.A.norm <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm %>% select(guide, Gene, guide_UMI, Replicate, active.A, inactive.A)
Reads.i.a.vivo.A.norm <- na.omit(Reads.i.a.vivo.A.norm)
colnames(Reads.i.a.vivo.A.norm) <- c("guide", "Gene", "guide_UMI","Replicate", "active", "inactive")
Reads.i.a.vivo.A.norm <- Reads.i.a.vivo.A.norm %>% mutate(sumReads = active + inactive)

A.vivo.guide_UMI_Repl <- Reads.i.a.vivo.A.norm %>% mutate(guide_UMI_Repl = paste(guide_UMI, Replicate, sep = '_')) # 1,042,721 guide_UMI_Repl
A.vivo.guide_UMI_Repl.SR20 <- A.vivo.guide_UMI_Repl %>% filter(sumReads >= 20) # 82,079 guide_UMI_Repl
#A.vivo.guide_UMI_Repl.SR30 <- A.vivo.guide_UMI_Repl %>% filter(sumReads >= 30) # 78714 guide_UMI_Repl

A.vivo.guide_UMI_Repl.SR20.pseudo <- A.vivo.guide_UMI_Repl.SR20 %>% mutate_at(vars(active, inactive), ~ . + 0.5)
Pre_Mageck_A.vivo.guide_UMI_Repl.SR20.pseudo <- A.vivo.guide_UMI_Repl.SR20.pseudo %>% select(guide_UMI_Repl, Gene, active, inactive)
#write.table(Pre_Mageck_A.vivo.guide_UMI_Repl.SR20.pseudo, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.vivo.guide_UMI_Repl.SR20.pseudo.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#Pre_Mageck_A.vivo.guide_UMI_Repl.SR20.pseudo <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/Pre_Mageck_A.vivo.guide_UMI_Repl.SR20.pseudo.tsv', header = T)
#sbatch mageck test -k Pre_Mageck_A.vivo.guide_UMI_Repl.SR20.pseudo.tsv -t active -c inactive -n Mageck.A.vivo.SR20.UMI.level.standard --remove-zero any 

# A.vivo.guide_UMI_Repl.SR30.pseudo <- A.vivo.guide_UMI_Repl.SR30 %>% mutate_at(vars(active, inactive), ~ . + 0.5)
# Pre_Mageck_A.vivo.guide_UMI_Repl.SR30.pseudo <- A.vivo.guide_UMI_Repl.SR30.pseudo %>% select(guide_UMI_Repl, Gene, active, inactive)
# #write.table(Pre_Mageck_A.vivo.guide_UMI_Repl.SR30.pseudo, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/Pre_Mageck_A.vivo.guide_UMI_Repl.SR30.pseudo.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# #sbatch mageck test -k Pre_Mageck_A.vivo.guide_UMI_Repl.SR30.pseudo.tsv -t active -c inactive -n Mageck.A.vivo.SR30.UMI.level.standard --remove-zero any


Reads.i.a.vivo.B.E.norm <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm %>% select(guide, Gene, guide_UMI, Replicate, active.B.E.norm, inactive.B.E.norm)
Reads.i.a.vivo.B.E.norm <- na.omit(Reads.i.a.vivo.B.E.norm)
colnames(Reads.i.a.vivo.B.E.norm) <- c("guide", "Gene", "guide_UMI","Replicate", "active", "inactive")
Reads.i.a.vivo.B.E.norm <- Reads.i.a.vivo.B.E.norm %>% mutate(sumReads = active + inactive)

B.E.FC2.vivo.guide_UMI_Repl <- Reads.i.a.vivo.B.E.norm %>% mutate(guide_UMI_Repl = paste(guide_UMI, Replicate, sep = '_')) # 1,042,721 guide_UMI_Repl

B.E.FC2.vivo.guide_UMI_Repl.SR20 <- B.E.FC2.vivo.guide_UMI_Repl %>% filter(sumReads >= 20) # 92,081 guide_UMI_Repl
B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo <- B.E.FC2.vivo.guide_UMI_Repl.SR20 %>% mutate_at(vars(active, inactive), ~ . + 0.5)
Pre_Mageck_B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo <- B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo %>% select(guide_UMI_Repl, Gene, active, inactive)
#write.table(Pre_Mageck_B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#Pre_Mageck_B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/Pre_Mageck_B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo.tsv', header = T)
#sbatch mageck test -k Pre_Mageck_B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo.tsv -t active -c inactive -n Mageck.B.E.FC2.vivo.SR20.UMI.level.standard --remove-zero any 

head(Pre_Mageck_A.vivo.guide_UMI_Repl.SR20.pseudo)
#A.B.E_FC2_guide_UMI_Repl <- rbind(Pre_Mageck_A.vivo.guide_UMI_Repl.SR20.pseudo,Pre_Mageck_B.E.FC2.vivo.guide_UMI_Repl.SR20.pseudo)
#write.table(A.B.E_FC2_guide_UMI_Repl, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2_guide_UMI_Repl.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

mageck_A_invivo_genes_pseudo.standard.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/Mageck.A.vivo.SR20.UMI.level.standard.gene_summary.txt', header = T)
#mageck_A_invivo_guides_pseudo.standard.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/Mageck.A.vivo.SR20.UMI.level.standard.sgrna_summary.txt', header = T)
mageck_A_invivo_genes_pseudo.standard.SR20$minRRA.fdr <- -log10(apply(mageck_A_invivo_genes_pseudo.standard.SR20[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A_invivo_genes_pseudo.standard.SR20$minRRA.score <- -log10(apply(mageck_A_invivo_genes_pseudo.standard.SR20[,c("neg.score","pos.score")], 1, min))
mageck_A_invivo_genes_pseudo.standard.SR20$minpvalue <- -log10(apply(mageck_A_invivo_genes_pseudo.standard.SR20[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A_invivo_genes_pseudo.standard.SR20$gene <- str_split_fixed(mageck_A_invivo_genes_pseudo.standard.SR20$id, '_', 2)[,1] # 21,267 genes
mageck_A_invivo_genes_pseudo.standard.SR20.UMIs3 <- mageck_A_invivo_genes_pseudo.standard.SR20 %>% filter(num >=3) #19,907 genes
#write.table(mageck_A_invivo_genes_pseudo.standard.SR20.UMIs3, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/mageck_A_invivo_genes_pseudo.standard.SR20.UMIs3.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
nrow(mageck_A_invivo_genes_pseudo.standard.SR20.UMIs3)

mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/Mageck.B.E.FC2.vivo.SR20.UMI.level.standard.gene_summary.txt', header = T)
#mageck_B.E_FC2.invivo_guides_pseudo.standard.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/Mageck.B.E.FC2.vivo.SR20.UMI.level.standard.sgrna_summary.txt', header = T)
mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20$minRRA.fdr <- -log10(apply(mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20[,c("neg.fdr","pos.fdr")], 1, min))
mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20$minRRA.score <- -log10(apply(mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20[,c("neg.score","pos.score")], 1, min))
mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20$minpvalue <- -log10(apply(mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20[,c("neg.p.value","pos.p.value")], 1, min))
mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20$gene <- str_split_fixed(mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20$id, '_', 2)[,1] # 21,267 genes
mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3 <- mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20 %>% filter(num >=3) #19,907 genes
#write.table(mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_05_Mageck_A.B.E/mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
nrow(mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3)


#mageck_A.B_overlap.genes.vivo.SR20.UMIs3 <- merge(mageck_A_invivo_genes_pseudo.standard.SR20 [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
# mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20 [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
# by = 'id', suffixes = c(".A",".B"))
mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 <- merge(mageck_A_invivo_genes_pseudo.standard.SR20.UMIs3 [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
                                                      mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3 [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
                                                      by = 'id', suffixes = c(".A",".B"))
#write.table(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3.txt', sep = "\t", header = T)

mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 <- mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% 
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

names(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3)


mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-9,7) + ylim(-9,7) +
  theme(plot.title = element_text(hjust = 0.5))

# Vivo AB overlap # ---------- Correlation of A & B OVERLAP genes ---------- regression line ---------- 

linear_model <- lm(neg.lfc.A ~ neg.lfc.B, data = mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3)
coefficients <- coef(linear_model)
slope <- coefficients[2] 
print(slope)

lm_output <- summary(lm(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3$neg.lfc.A ~ mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3$neg.lfc.B))
slope = lm_output$coefficients[2,1]
cor(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3$neg.lfc.A, mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3$neg.lfc.B)

colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')

mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE) + 
  stat_cor(method = 'pearson', p.accuracy = 0.001, r.accuracy = 0.01)+
  
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-9,7) + ylim(-9,7) +
  theme(plot.title = element_text(hjust = 0.5))

mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
  stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-9,8.5) + ylim(-9,8.5) +
  theme(plot.title = element_text(hjust = 0.5))
#

# Vivo AB overlap # ---------- Correlation of A & B OVERLAP genes ---------- different plots ---------- 

colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red', 'Non-Essential' = 'black')
mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), alpha=0.5, size=1) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'StAR analysis - In vivo pool A vs pool B \nSR20.UMIs3 & \nCore-Essentials: Mouse.ortholog_CoreEssentials_934 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors2) +
  xlim(-9,7) + ylim(-9,7) +
  theme(plot.title = element_text(hjust = 0.5))

nrow(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene)) #99 genes
nrow(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_650$gene)) #113 genes
nrow(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene)) #308 genes,   208 genes > depmap.Core-Ess.min1.counts400
nrow(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$gene)) #187 genes 

colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-9,7) + ylim(-9,7) +
  theme(plot.title = element_text(hjust = 0.5))

Grey.validating.AB <- mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(neg.lfc.A <= -2 & neg.lfc.B <= -2 & group == 'Other') 

colors5 <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black', 'Validating A & B' = 'Blue')
mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id),  size=1) +
  geom_point(aes(color = 'Validating A & B'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Grey.validating.AB$id),  size=1) +
  
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors5) +
  xlim(-9,7) + ylim(-9,7) +
  theme(plot.title = element_text(hjust = 0.5))

A.B.E.FC2_vivo_spec.depleting.vivo1xRRA3
colors5 <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black', 'A.B.E.FC2_vivo_spec.depleting.vivo1x' = 'Gold')
mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id),  size=1) +
  geom_point(aes(color = 'A.B.E.FC2_vivo_spec.depleting.vivo1x'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% A.B.E.FC2_vivo_spec.depleting.vivo1x$id),  size=1) +
  
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors5) +
  xlim(-9,7) + ylim(-9,7) +
  theme(plot.title = element_text(hjust = 0.5))



# Vivo AB overlap # ---------- Correlation of A & B OVERLAP genes ---  With density plots ---------------
mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl <- mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id)
mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3

# Plot1 <- mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
#   ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
#   geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
#   geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
#   #geom_point(aes(color = 'Other'), size = 1) +
#   geom_point(aes(color = 'Other'), data = mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl$id), size = 1) +
#   geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
#   geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
#   geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
#   #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
#   stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
#   xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B')) +
#   #labs(color = ' ', title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
#   theme_pubr(base_size = 14) + scale_color_manual(values = colors) + 
#   theme(panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = 'none') +
#   xlim(-9,8.5) + ylim(-9,8.5) +   
#   coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
#   theme(plot.title = element_text(hjust = 0.5))

Plot1 <- mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
  stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B')) +
  #labs(color = ' ', title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) + 
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = 'none') +
  xlim(-9,8.5) + ylim(-9,8.5) +   
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  theme(plot.title = element_text(hjust = 0.5))

Plot.Dens.A <- mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.A)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL)  +
  xlim(-9,8.5) +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  theme(legend.position = "none", axis.text.x = element_blank())

Plot.Dens.B <- mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.B)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), alpha=0.5, linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL) +
  xlim(-9,8.5) +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) + 
  theme(legend.position = "none", axis.text.y = element_blank()) + coord_flip()

Plot.Dens.A + plot_spacer() + Plot1 + theme(aspect.ratio = 1)+ Plot.Dens.B + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

#

#!!! Vitro AB overlap# ---------- Correlation of A & B OVERLAP genes ---------- MAGECK ---------- 

# FC1.reads.i.a.vitro.A.wide <- dcast(setDT(FC1.reads.i.a.vitro.A), guide + Gene ~ Replicate, value.var = c('active','inactive'), fill = 0) %>% data.frame()
# FC2_reads.i.a.vitro.B.E.wide <- dcast(setDT(FC2_reads.i.a.vitro.B.E), guide + Gene ~ Replicate, value.var = c('active','inactive'), fill = 0) %>% data.frame()

#after normalization
FC1.reads.i.a.vitro.A.wide <- A.B.E.FC2_guide.reads.vitro.Repl.norm %>% select(guide, Gene, Replicate, active.A, inactive.A)
FC1.reads.i.a.vitro.A.wide <- na.omit(FC1.reads.i.a.vitro.A.wide)
colnames(FC1.reads.i.a.vitro.A.wide) <- c("guide", "Gene", "Replicate", "active", "inactive")
FC1.reads.i.a.vitro.A.wide <- dcast(setDT(FC1.reads.i.a.vitro.A.wide), guide + Gene ~ Replicate, value.var = c('active','inactive'), fill = 0, ) %>% data.frame()

FC2_reads.i.a.vitro.B.E.wide <- A.B.E.FC2_guide.reads.vitro.Repl.norm %>% select(guide, Gene, Replicate, active.B.E.norm, inactive.B.E.norm)
FC2_reads.i.a.vitro.B.E.wide <- na.omit(FC2_reads.i.a.vitro.B.E.wide)
colnames(FC2_reads.i.a.vitro.B.E.wide) <- c("guide", "Gene", "Replicate", "active", "inactive")
FC2_reads.i.a.vitro.B.E.wide <- dcast(setDT(FC2_reads.i.a.vitro.B.E.wide), guide + Gene ~ Replicate, value.var = c('active','inactive'), fill = 0, ) %>% data.frame()

Pre_Mageck_FC1.reads.i.a.vitro.A.wide <- FC1.reads.i.a.vitro.A.wide %>% mutate_at(vars(active_Blue, active_Green, active_Red, inactive_Blue, inactive_Green, inactive_Red), ~ . + 0.5)  # 61993 guides
Pre_Mageck_FC2_reads.i.a.vitro.B.E.wide <- FC2_reads.i.a.vitro.B.E.wide %>% mutate_at(vars(active_Blue, active_Green, active_Red, inactive_Blue, inactive_Green, inactive_Red), ~ . + 0.5)  # 68718 guides

#write.table(Pre_Mageck_FC1.reads.i.a.vitro.A.wide, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_FC1.reads.i.a.vitro.A.wide.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(Pre_Mageck_FC2_reads.i.a.vitro.B.E.wide, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_FC2_reads.i.a.vitro.B.E.wide.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# sbatch mageck test -k Pre_Mageck_FC1.reads.i.a.vitro.A.wide.tsv -t active_Blue,active_Green,active_Red -c inactive_Blue,inactive_Green,inactive_Red -n Mageck.A.invitro.pseudo.paired.standard --remove-zero any --paired 
# sbatch mageck test -k Pre_Mageck_FC2_reads.i.a.vitro.B.E.wide.tsv -t active_Blue,active_Green,active_Red -c inactive_Blue,inactive_Green,inactive_Red -n Mageck.B.E_FC2.invitro.pseudo.paired.standard --remove-zero any --paired 

mageck_A_invitro_genes_pseudo.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck.A.invitro.pseudo.paired.standard.gene_summary.txt', header = T)
#mageck_A_invivo_guides_pseudo.standard.SR20 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck.A.invitro.pseudo.paired.standard.sgrna_summary.txt', header = T)
mageck_A_invitro_genes_pseudo.standard$minRRA.fdr <- -log10(apply(mageck_A_invitro_genes_pseudo.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A_invitro_genes_pseudo.standard$minRRA.score <- -log10(apply(mageck_A_invitro_genes_pseudo.standard[,c("neg.score","pos.score")], 1, min))
mageck_A_invitro_genes_pseudo.standard$minpvalue <- -log10(apply(mageck_A_invitro_genes_pseudo.standard[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A_invitro_genes_pseudo.standard$gene <- str_split_fixed(mageck_A_invitro_genes_pseudo.standard$id, '_', 2)[,1] # 12395 genes
nrow(mageck_A_invitro_genes_pseudo.standard)

mageck_B.E_FC2.invitro_genes_pseudo.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck.B.E_FC2.invitro.pseudo.paired.standard.gene_summary.txt', header = T)
#mageck_B.E_FC2.invitro_guides_pseudo.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck.B.E_FC2.invitro.pseudo.paired.standard.sgrna_summary.txt', header = T)
mageck_B.E_FC2.invitro_genes_pseudo.standard$minRRA.fdr <- -log10(apply(mageck_B.E_FC2.invitro_genes_pseudo.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_B.E_FC2.invitro_genes_pseudo.standard$minRRA.score <- -log10(apply(mageck_B.E_FC2.invitro_genes_pseudo.standard[,c("neg.score","pos.score")], 1, min))
mageck_B.E_FC2.invitro_genes_pseudo.standard$minpvalue <- -log10(apply(mageck_B.E_FC2.invitro_genes_pseudo.standard[,c("neg.p.value","pos.p.value")], 1, min))
mageck_B.E_FC2.invitro_genes_pseudo.standard$gene <- str_split_fixed(mageck_B.E_FC2.invitro_genes_pseudo.standard$id, '_', 2)[,1] # 12395 genes
nrow(mageck_B.E_FC2.invitro_genes_pseudo.standard)


mageck_A.B.FC2_overlap.genes.vitro <- merge(mageck_A_invitro_genes_pseudo.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
                                            mageck_B.E_FC2.invitro_genes_pseudo.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
                                            by = 'id', suffixes = c(".A",".B"))
# #write.table(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3.txt', sep = "\t", header = T)

mageck_A.B.FC2_overlap.genes.vitro <- mageck_A.B.FC2_overlap.genes.vitro %>% 
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

names(mageck_A.B.FC2_overlap.genes.vitro)

mageck_A.B.FC2_overlap.genes.vitro.other.nondepl <- mageck_A.B.FC2_overlap.genes.vitro %>% filter(group == "Other" & neg.lfc.A < 1 & neg.lfc.A > -1 & neg.lfc.B < 1 & neg.lfc.B > -1)
mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl <- mageck_A.B.FC2_overlap.genes.vitro %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id)


colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')

mageck_A.B.FC2_overlap.genes.vitro %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  xlab(bquote('Log'[2]~'fold change in vitro pool A')) + ylab(bquote('Log'[2]~'fold change in vitro pool B'))+ labs(color = ' ', 
                                                                                                                    title = 'StAR analysis - In vitro pool A vs pool B \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-9,7) + ylim(-9,7) +
  theme(plot.title = element_text(hjust = 0.5))


colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black')
mageck_A.B.FC2_overlap.genes.vitro %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene),  size=1) +
  xlab(bquote('Log'[2]~'fold change in vitro pool A')) + ylab(bquote('Log'[2]~'fold change in vitro pool B'))+ labs(color = ' ', 
                                                                                                                    title = 'StAR analysis - In vitro pool A vs pool B \nSR20.UMIs3 & \nCore-Essentials: Mouse.ortholog_CoreEssentials_934 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors2) +
  xlim(-9,7) + ylim(-9,7) +
  theme(plot.title = element_text(hjust = 0.5))

mageck_A.B.FC2_overlap.genes.vitro %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
  stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
  xlab(bquote('Log'[2]~'fold change in vitro pool A')) + ylab(bquote('Log'[2]~'fold change in vitro pool B'))+ labs(color = ' ', 
                                                                                                                    title = 'StAR analysis - In vitro pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-9,8.5) + ylim(-9,8.5) +
  theme(plot.title = element_text(hjust = 0.5))
#
# Vitro AB overlap # ---------- Correlation of A & B OVERLAP genes ---  With density plots ---------------
mageck_A.B.FC2_overlap.genes.vitro.other.nondepl <- mageck_A.B.FC2_overlap.genes.vitro %>% filter(group == "Other" & neg.lfc.A < 1 & neg.lfc.A > -1 & neg.lfc.B < 1 & neg.lfc.B > -1)
mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl <- mageck_A.B.FC2_overlap.genes.vitro %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id)
mageck_A.B.FC2_overlap.genes.vitro

colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
colors10 <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black', 'Blue' = 'Blue')

Plot1 <- mageck_A.B.FC2_overlap.genes.vitro %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  # geom_point(aes(color = 'Other'), data = mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl$id), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') +
  #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
  stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
  xlab(bquote('Log'[2]~'fold change in vitro pool A')) + ylab(bquote('Log'[2]~'fold change in vitro pool B')) +
  #labs(color = ' ', title = 'StAR analysis - In vitro pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors10) +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = 'none') +
  xlim(-9,8.5) + ylim(-9,8.5) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  theme(plot.title = element_text(hjust = 0.5))

Plot1 <- mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
  stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
  xlab(bquote('Log'[2]~'fold change in vitro pool A')) + ylab(bquote('Log'[2]~'fold change in vitro pool B')) +
  #labs(color = ' ', title = 'StAR analysis - In vitro pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors10) + 
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = 'none') +
  xlim(-9,8.5) + ylim(-9,8.5) +   
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  theme(plot.title = element_text(hjust = 0.5))

Plot.Dens.A <- mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.A)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'),  linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id),linewidth=1) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL)  +
  xlim(-9,8.5) +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  theme(legend.position = "none", axis.text.x = element_blank())

Plot.Dens.B <- mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.B)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL) +
  xlim(-9,8.5) +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) + 
  theme(legend.position = "none", axis.text.y = element_blank()) + coord_flip()

Plot.Dens.A + plot_spacer() + Plot1 + theme(aspect.ratio = 1)+ Plot.Dens.B + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
#
##### Vitro # ---------- A.B.E_FC2 - guide level - Merge, Pool ------------------
A.B.E.FC2_guide.reads.vitro.pooled <- A.B.E.FC2_guide.reads.vitro.Repl_wide
A.B.E.FC2_guide.reads.vitro.pooled$active <- rowSums(A.B.E.FC2_guide.reads.vitro.pooled[, c("active_Blue", "active_Green", "active_Red")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.pooled$inactive <- rowSums(A.B.E.FC2_guide.reads.vitro.pooled[, c("inactive_Blue", "inactive_Green", "inactive_Red")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.pooled$sumReads <- rowSums(A.B.E.FC2_guide.reads.vitro.pooled[, c("active", "inactive")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.pooled <- A.B.E.FC2_guide.reads.vitro.pooled %>% select(guide, Gene, active, inactive, sumReads)
#write.table(A.B.E.FC2_guide.reads.vitro.pooled, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.pooled.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E.FC2_guide.reads.vitro.pooled<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.pooled.txt", header=T, sep = '\t')

##### Vitro # ---------- A.B.E_FC2 - guide level - Merge, Pool - after normalization ------------------
A.B.E.FC2_guide.reads.vitro.norm.pooled <- A.B.E.FC2_guide.reads.vitro.Repl_wide.norm
A.B.E.FC2_guide.reads.vitro.norm.pooled$active <- rowSums(A.B.E.FC2_guide.reads.vitro.norm.pooled[, c("active_Blue", "active_Green", "active_Red")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.norm.pooled$inactive <- rowSums(A.B.E.FC2_guide.reads.vitro.norm.pooled[, c("inactive_Blue", "inactive_Green", "inactive_Red")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.norm.pooled$sumReads <- rowSums(A.B.E.FC2_guide.reads.vitro.norm.pooled[, c("active", "inactive")], na.rm = TRUE)
A.B.E.FC2_guide.reads.vitro.norm.pooled <- A.B.E.FC2_guide.reads.vitro.norm.pooled %>% select(guide, Gene, active, inactive, sumReads)
write.table(A.B.E.FC2_guide.reads.vitro.norm.pooled, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.norm.pooled.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E.FC2_guide.reads.vitro.norm.pooled<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.norm.pooled.txt", header=T, sep = '\t')


# Vivo # ---------- A.B.E_F2 -  Pool guide level------------------
A.B.E_FC2.pooled_guides.vivo.ia0 <- A.B.E_FC2.pooled_most_filters.ia0 %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
#write.table(A.B.E_FC2.pooled_guides.vivo.ia0, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.ia0.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E_FC2.pooled_guides.vivo.ia0<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.ia0.txt", header=T, sep = '\t')

A.B.E_FC2.pooled_guides.vivo.ia0.norm <- A.B.E_FC2.pooled_most_filters.ia0.norm %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))

A.B.E_FC2.pooled_guides.vivo.SR20 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20 %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
#write.table(A.B.E_FC2.pooled_guides.vivo.SR20, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.SR20.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E_FC2.pooled_guides.vivo.SR20<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.SR20.txt", header=T, sep = '\t')

A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3 <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3 %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
#write.table(A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3, "/Volumes/groups/elling/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3<- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3.txt", header=T, sep = '\t')


##### A.B.E_FC2 pooled with PlasmidLibr data # ----- SR20 And/OR UMIs and nonfiltered------------
#combined_data_pooled <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Calculated.Plasmid.reads_combined_data_pooled.txt", sep = "\t", header = T)
#combined_data_pooled <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Calculated.Plasmid.reads_combined_data_pooled.norm.txt", sep = "\t", header = T)

A.B.E_FC2.pooled_guides.vivo.SR20.plasmidLibr <- merge(combined_data_pooled, A.B.E_FC2.pooled_guides.vivo.SR20, by = "guide")
names(A.B.E_FC2.pooled_guides.vivo.SR20.plasmidLibr)
#write.table(A.B.E_FC2.pooled_guides.vivo.SR20.plasmidLibr, "/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.SR20.plasmidLibr.txt", sep = "\t", row.names = F, col.names = T, quote = F)
A.B.E_FC2.pooled_guides.vivo.SR20.plasmidLibr <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.SR20.plasmidLibr.txt", sep = "\t", header = T)

A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3.plasmidLibr <- merge(combined_data_pooled, A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3, by = "guide")
names(A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3.plasmidLibr)
#write.table(A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3.plasmidLibr, "/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3.plasmidLibr.txt", sep = "\t", row.names = F, col.names = T, quote = F)
A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3.plasmidLibr <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.SR20.UMIs3.plasmidLibr.txt", sep = "\t", header = T)

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt <- bind_rows(left_join(combined_data_pooled, A.B.E_FC2.pooled_guides.vivo.ia0, by = "guide"),
                                                             A.B.E_FC2.pooled_guides.vivo.ia0 %>% anti_join(combined_data_pooled, by = "guide")) %>% replace(is.na(.), 0) #106,598 rows
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt$Gene <- str_split_fixed(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt$guide, '_', 2)[, 1]
names(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt)
#write.table(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt, "/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.txt", sep = "\t", row.names = F, col.names = T, quote = F)
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.txt", sep = "\t", header = T)

A.B.E.FC2_guide.reads.vitro <- A.B.E.FC2_guide.reads.vitro.Repl %>% group_by(guide, Gene) %>% summarise(active=sum(active),inactive=sum(inactive))
#write.table(A.B.E.FC2_guide.reads.vitro, "/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E.FC2_guide.reads.vitro.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Vivo # ---------- A.B.E_F2 -  Pool guide level + pool with plasmidlibr data - after normalization------------------

A.B.E_FC2.pooled_guides.vivo.norm.ia0 <- A.B.E_FC2.pooled_most_filters.ia0.norm %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))

#combined_data_pooled <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Calculated.Plasmid.reads_combined_data_pooled.norm.txt", sep = "\t", header = T)
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm <- bind_rows(left_join(combined_data_pooled, A.B.E_FC2.pooled_guides.vivo.norm.ia0, by = "guide"),
                                                                  A.B.E_FC2.pooled_guides.vivo.norm.ia0 %>% anti_join(combined_data_pooled, by = "guide")) %>% replace(is.na(.), 0) #106,598 rows
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm$Gene <- str_split_fixed(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm$guide, '_', 2)[, 1]
names(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm)
#write.table(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm, "/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm.txt", sep = "\t", row.names = F, col.names = T, quote = F)
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm.txt", sep = "\t", header = T)


# Vivo # ---------- CONVENTIONAL - Pre-Mageck ------------------------------------ 
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt

Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt %>% select(guide, Gene, active, pooled_Calculated_reads)
Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt <- Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt %>% mutate_at(vars(active, pooled_Calculated_reads), ~ . + 0.5)
#write.table(Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.tsv', header = T)

#sbatch mageck test -k Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.tsv -t active -c pooled_Calculated_reads -n Mageck_A.B.E_FC2.CONV.vivo.guide.level.standard --remove-zero any 

# Vivo # ---------- CONVENTIONAL - Pre-Mageck - after normalization ------------------------------------ 

Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.norm <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm %>% select(guide, Gene, active, pooled_Calculated_reads)
Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.norm <- Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.norm %>% mutate_at(vars(active, pooled_Calculated_reads), ~ . + 0.5)
#write.table(Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.norm, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.norm.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.norm <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.norm.tsv', header = T)

#sbatch mageck test -k Pre_Mageck_A.B.E_F2.guide.vivo.plasmid.nofilt.norm.tsv -t active -c pooled_Calculated_reads -n Mageck_A.B.E_FC2.CONV.vivo.norm.guide.level.standard --remove-zero any 

# Vivo # ---------- CONVENTIONAL - Mageck ------------------------------------ 
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.CONV.vivo.guide.level.standard.gene_summary.txt', header = T)
#mageck_A.B.E.FC2_vivo_Conv_guides.standard.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.CONV.vivo.guide.level.standard.sgrna_summary.txt', header = T)
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minRRA.fdr <- -log10(apply(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minRRA.score <- -log10(apply(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt[,c("neg.score","pos.score")], 1, min))
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minpvalue <- -log10(apply(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$gene <- str_split_fixed(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$id, '_', 2)[,1] # 21382 genes
length(unique(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$id)) # 21382 genes

mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$adj.minRRA.score <- ifelse(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minRRA.score > 15, 15, 
                                                                            mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minRRA.score)
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt <-  mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>%
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

#write.table(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3.txt', header = T)

# Vivo # ---------- CONVENTIONAL - Mageck - after normalization ------------------------------------ 

#from here the dataframes are labeled the same as without normalization to keep it easier to rerun certain parts of the script

mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.CONV.vivo.norm.guide.level.standard.gene_summary.txt', header = T)
#mageck_A.B.E.FC2_vivo_Conv_guides.standard.nofilt <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.B.E_FC2.CONV.vivo.norm.guide.level.standard.sgrna_summary.txt', header = T)
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minRRA.fdr <- -log10(apply(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minRRA.score <- -log10(apply(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt[,c("neg.score","pos.score")], 1, min))
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minpvalue <- -log10(apply(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$gene <- str_split_fixed(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$id, '_', 2)[,1] # 21382 genes
length(unique(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$id)) # 21382 genes

mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$adj.minRRA.score <- ifelse(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minRRA.score > 15, 15, 
                                                                            mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt$minRRA.score)
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt <-  mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>%
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

#write.table(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3.txt', header = T)
#write.table(mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt.norm.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#!!! VIVO # ---------- CONVENTIONAL Volcano plot - minRRA.score ------------------
colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')
mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>%
  ggplot(aes(x = neg.lfc, y = adj.minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'),  size = 1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter((gene) %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  #geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ labs(color = ' ', title = 'Conventional analysis \nIn vitro gene LFC + RRA score \nno filters') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-8,14) + ylim(0,15) +
  theme(plot.title = element_text(hjust = 0.5)) # + guides(colour = 'none')

nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id)) #529 genes
nrow(mageck_A.B.E.FC2_vivo_genes.standard.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene)) # 498 genes
#filter(toupper(gene) %in% Controls$gene)
#!!! VIVO # ---------- CONVENTIONAL Volcano plot - minRRA.score ------------------ With density plots ---------------

VPlot1 <- 
  mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>%
  ggplot(aes(x = neg.lfc, y = adj.minRRA.score)) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'),  size = 1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter((gene) %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  #geom_point(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), alpha=0.5, size=1) +
  xlab(bquote('Log'[2]~'fold change')) + ylab(bquote('RRA score'))+ 
  #labs(color = ' ', title = 'Conventional analysis \nIn vitro gene LFC + RRA score \nno filters') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-8,14) + ylim(0,15) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")  #+ #guides(colour = 'none')
#theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank())  #+ #guides(colour = 'none')

Plot.Dens.V1x <- mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>%
  ggplot(aes(x = neg.lfc)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5, linetype = 'dashed') +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL)  +
  xlim(-8,14) +  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  theme(legend.position = "none", axis.text.x = element_blank())

Plot.Dens.V1y <- mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>%
  ggplot(aes(x = adj.minRRA.score)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  #geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), alpha=0.5, linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.E.FC2_vivo_Conv_genes.standard.nofilt %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL) +
  xlim(0,15) + theme_pubr(base_size = 12) + scale_color_manual(values = colors) + 
  theme(legend.position = "none", axis.text.y = element_blank()) + coord_flip()

# Plot.Dens.V1x + plot_spacer() + VPlot1 + theme(aspect.ratio = 1)+ Plot.Dens.V1y +
#   plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4), guides = 'collect') & theme(legend.position = 'top')
Plot.Dens.V1x + plot_spacer() + VPlot1 + theme(aspect.ratio = 1)+ Plot.Dens.V1y +
  plot_layout(ncol = 2, nrow = 2, widths = c(2, 0.5), heights = c(0.5, 2), guides = 'collect') & theme(legend.position = 'top')
# VPlot1 + theme(aspect.ratio = 1)+ Plot.Dens.V1y + Plot.Dens.V1x + plot_spacer() + 
#   plot_layout(ncol = 2, nrow = 2, widths = c(3, 1), heights = c(3, 1))

#


#!!! dAUC # ---------- CONVENTIONAL ---------- no Others! different Essential groups ---------- 

# A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt %>%
#   mutate(group= ifelse(Gene %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(Gene %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.norm %>%
  mutate(group= ifelse(Gene %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(Gene %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups %>% mutate_at(vars(pooled_Calculated_reads, active, inactive), ~ . + 0.5)
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups$Conv.LFC <- log2(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups$active/A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups$pooled_Calculated_reads)

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.noOthers <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups %>% filter(group == "Essential" | group == "Non-Essential")

names(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.noOthers)
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.cumulative <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.noOthers %>% 
  select(guide, Gene, pooled_Calculated_reads, active, Conv.LFC, group) %>% na.omit()


A.B.E.guide.dAUC.rank.LFCConv <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.cumulative %>% arrange(Conv.LFC)
A.B.E.guide.dAUC.rank.LFCConv$all_rank <- 1:nrow(A.B.E.guide.dAUC.rank.LFCConv)
A.B.E.guide.dAUC.rank.LFCConv$all_rank_perc_sgRNA <- (A.B.E.guide.dAUC.rank.LFCConv$all_rank*1)/nrow(A.B.E.guide.dAUC.rank.LFCConv) 

A.B.E.guide.dAUC.rank.LFCConv <- A.B.E.guide.dAUC.rank.LFCConv %>% group_by(group) %>% mutate(rank_group = 1:n(), 
                                                                                              rank_group_prct = 1:n()/n(), 
                                                                                              rank_sum = cumsum(rank_group_prct)) 
# Group the data by "group"
grouped_data_Conv <- A.B.E.guide.dAUC.rank.LFCConv %>%
  group_by(group) %>% arrange(all_rank_perc_sgRNA)
# Function to calculate AUC using the trapezoidal rule
calculate_auc_Conv <- function(x, y) {
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2) }
# Calculate AUC for each group
auc_results_Conv <- grouped_data_Conv %>%
  summarize(AUC = calculate_auc_Conv(all_rank_perc_sgRNA, rank_group_prct))

Essential_auc_Conv <- auc_results_Conv$AUC[auc_results_Conv$group == "Essential"]
NonEssential_auc_Conv <- auc_results_Conv$AUC[auc_results_Conv$group == "Non-Essential"]
dAUC_Ess_nonEss_Conv <- Essential_auc_Conv - NonEssential_auc_Conv
# Print the dAUC value
print(dAUC_Ess_nonEss_Conv)

ggplot(A.B.E.guide.dAUC.rank.LFCConv, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group)) + 
  geom_line(size = 0.75) + 
  geom_hline(yintercept = max(A.B.E.guide.dAUC.rank.LFCConv$rank_group_prct), linetype = "solid", color = "black") +
  geom_vline(xintercept = max(A.B.E.guide.dAUC.rank.LFCConv$all_rank_perc_sgRNA), linetype = "solid", color = "black") +
  theme_pubr() + 
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "blue", "Other" = "gray")) +
  labs(title = "Conventional AUC - nofilters \nEssentials: Vitro.depleting.min3.RRAsc10 \nControls: Mouse.ortholog_Controls_550", y = "Cumulative fraction", x = "Percentage rank of sgRNAs") + 
  annotate("text", x = 0.8, y = 0.1, label = paste("dAUC =", round(dAUC_Ess_nonEss_Conv, 2)), color = "black")

# VIVO # rank # ---------- CONVENTIONAL ---------- ----------  ---------- 
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt %>%
  mutate(group= ifelse(Gene %in% Mouse.ortholog_Controls_550$gene, 'Non-Essential', ifelse(Gene %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups %>% mutate_at(vars(pooled_Calculated_reads, active, inactive), ~ . + 0.5)
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups$Conv.LFC <- log2(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups$active/A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups$pooled_Calculated_reads)

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups %>% arrange(Conv.LFC)
head(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank)

gene_counts <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank %>% group_by(Gene) %>% summarize(count = n())
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank %>%
  left_join(gene_counts, by = "Gene") %>% filter(count <= 5) %>% select(-count)

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5$Rank <- 1:nrow(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5)
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 %>% select(-active, -inactive, -guide_name, -pooled_Calculated_reads)
head(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5)


A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 %>%
  mutate(Unique_Genes_Count = cumsum(!duplicated(Gene)))

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 %>% mutate(average_guides = Rank / Unique_Genes_Count)

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 %>%
  arrange(Rank) %>%
  group_by(Gene) %>%
  mutate(sgRNA.gene = row_number()) %>%
  ungroup()


ggplot(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5, aes(x = Rank, y = average_guides)) +
  geom_point() + 
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene") +
  #xlim(0,5000) +
  theme_pubr()


ggplot(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5, aes(x = Rank, y = average_guides)) +
  geom_point() + 
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene") +
  #xlim(0,5000) +
  scale_x_continuous(labels = label_scientific()) +
  labs(title = 'Conventional') +
  theme_pubr()


A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 %>% filter(group == "Essential") %>%
  distinct(Gene) %>%
  nrow()
#

#

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5 %>%  mutate(random_order = runif(nrow(.))) %>%  # create a column of random numbers
  arrange(random_order) %>%  # sort the dataframe by the random column
  select(-random_order)  # remove the random column

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random$Rank <- 1:nrow(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random)


A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random %>%
  mutate(Unique_Genes_Count = cumsum(!duplicated(Gene)))

A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random %>% mutate(average_guides = Rank / Unique_Genes_Count)


A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5_select <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5[c("Rank", "average_guides")]
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random_select <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random[c("Rank", "average_guides")]

combined.vivo.Conv <- merge(A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5_select, A.B.E_FC2.pooled_guides.vivo.plasmidLibr.random_select, by = "Rank", suffixes = c('.vvConv', '.vvConv.rdm'))

colnames(combined.vivo.Conv) 

combined_long_vivo.Conv <- reshape2::melt(combined.vivo.Conv, id.vars = "Rank", variable.name = "Method", value.name = "Avg_guides")


ggplot(combined_long_vivo.Conv, aes(x = Rank, y = Avg_guides, color = Method)) +
  geom_point(data = combined_long_vivo.Conv[combined_long_vivo.Conv$Method == "average_guides.vvConv.rdm",], aes(x = Rank, y = Avg_guides)) +
  geom_point(data = combined_long_vivo.Conv[combined_long_vivo.Conv$Method == "average_guides.vvConv",], aes(x = Rank, y = Avg_guides)) +
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene",
       title = 'Reproducibility of sgRNAs - Conventional') +
  xlim(0,3e+03)+
  scale_x_continuous(labels = scales::label_scientific()) +
  scale_color_manual(values =c('average_guides.vvConv' = "#000000", 'average_guides.vvConv.rdm' = "#999999")) +
  theme_pubr()



#!!! Vivo AB overlap # ---------- CONVENTIONAL ---------- Correlation of A & B OVERLAP genes ---------- MAGECK ---------- 


FC1.FC2.reads.i.a.vivo.A.guide <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm %>% select(guide, Gene, active.A, inactive.A)
FC1.FC2.reads.i.a.vivo.A.guide <- na.omit(FC1.FC2.reads.i.a.vivo.A.guide)
colnames(FC1.FC2.reads.i.a.vivo.A.guide) <- c("guide", "Gene", "active", "inactive")

FC2_reads.i.a.vivo.B.E.guide <- A.B.E_FC2_poly_hop.filtered.reads.ia0.vivo_0.001.norm %>% select(guide, Gene, active.B.E.norm, inactive.B.E.norm)
FC2_reads.i.a.vivo.B.E.guide <- na.omit(FC2_reads.i.a.vivo.B.E.guide)
colnames(FC2_reads.i.a.vivo.B.E.guide) <- c("guide", "Gene", "active", "inactive")


FC1.FC2.reads.i.a.vivo.A.guide <- FC1.FC2.reads.i.a.vivo.A.guide %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
FC2_reads.i.a.vivo.B.E.guide <- FC2_reads.i.a.vivo.B.E.guide %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))


#FC1.FC2.reads.i.a.vivo.A.guide <- FC1_FC2_poly_hop.filtered.reads.vivo_0.001.A %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))
#FC2_reads.i.a.vivo.B.E.guide <- FC2_poly_hop.filtered.reads.vivo_0.001.B.E %>% group_by(guide, Gene) %>% summarise(active = sum(active), inactive = sum(inactive))

# write.table(FC1.FC2.reads.i.a.vivo.A.guide, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/FC1.FC2.reads.i.a.vivo.A.guide.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(FC2_reads.i.a.vivo.B.E.guide, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_02_Data_lists/FC2_reads.i.a.vivo.B.E.guide.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

combined_data_pooled_A <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Calculated.Plasmid.reads_combined_data_pooled_A.txt", sep = "\t", header = T)

A.pooled_guides.vivo.plasmidLibr.nofilt <- bind_rows(left_join(combined_data_pooled_A, FC1.FC2.reads.i.a.vivo.A.guide, by = "guide"),
                                                     FC1.FC2.reads.i.a.vivo.A.guide %>% anti_join(combined_data_pooled_A, by = "guide")) %>% replace(is.na(.), 0) #105234 rows
A.pooled_guides.vivo.plasmidLibr.nofilt$sumReads <- rowSums(A.pooled_guides.vivo.plasmidLibr.nofilt[, c("active", "inactive")], na.rm = TRUE)

A.pooled_guides.vivo.plasmidLibr.nofilt <- A.pooled_guides.vivo.plasmidLibr.nofilt %>% group_by(Gene) %>% filter(!(pooled_Calculated_reads == 0| sum(sumReads == 0) >= 5))
A.pooled_guides.vivo.plasmidLibr.nofilt$Gene <- str_split_fixed(A.pooled_guides.vivo.plasmidLibr.nofilt$guide, '_', 2)[, 1]
Pre_Mageck_A.guide.vivo.plasmid.nofilt <- A.pooled_guides.vivo.plasmidLibr.nofilt %>% select(guide, Gene, active, pooled_Calculated_reads)
Pre_Mageck_A.guide.vivo.plasmid.nofilt <- Pre_Mageck_A.guide.vivo.plasmid.nofilt %>% mutate_at(vars(active, pooled_Calculated_reads), ~ . + 0.5)
#write.table(Pre_Mageck_A.guide.vivo.plasmid.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.guide.vivo.plasmid.nofilt.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#sbatch mageck test -k Pre_Mageck_A.guide.vivo.plasmid.nofilt.tsv -t active -c pooled_Calculated_reads -n Mageck_A.CONV.vivo.guide.level.standard --remove-zero any 

combined_data_pooled_B <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Calculated.Plasmid.reads_combined_data_pooled_B.txt", sep = "\t", header = T)

B.E.pooled_guides.vivo.plasmidLibr.nofilt <- bind_rows(left_join(combined_data_pooled_B, FC2_reads.i.a.vivo.B.E.guide, by = "guide"),
                                                       FC2_reads.i.a.vivo.B.E.guide %>% anti_join(combined_data_pooled_B, by = "guide")) %>% replace(is.na(.), 0) #106466 rows
B.E.pooled_guides.vivo.plasmidLibr.nofilt$sumReads <- rowSums(B.E.pooled_guides.vivo.plasmidLibr.nofilt[, c("active", "inactive")], na.rm = TRUE)
B.E.pooled_guides.vivo.plasmidLibr.nofilt <- B.E.pooled_guides.vivo.plasmidLibr.nofilt %>% group_by(Gene) %>% filter(!(pooled_Calculated_reads == 0| sum(sumReads == 0) >= 5))
B.E.pooled_guides.vivo.plasmidLibr.nofilt$Gene <- str_split_fixed(B.E.pooled_guides.vivo.plasmidLibr.nofilt$guide, '_', 2)[, 1]
Pre_Mageck_B.E.guide.vivo.plasmid.nofilt <- B.E.pooled_guides.vivo.plasmidLibr.nofilt %>% select(guide, Gene, active, pooled_Calculated_reads)
Pre_Mageck_B.E.guide.vivo.plasmid.nofilt <- Pre_Mageck_B.E.guide.vivo.plasmid.nofilt %>% mutate_at(vars(active, pooled_Calculated_reads), ~ . + 0.5)
#write.table(Pre_Mageck_B.E.guide.vivo.plasmid.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_B.E.guide.vivo.plasmid.nofilt.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#sbatch mageck test -k Pre_Mageck_B.E.guide.vivo.plasmid.nofilt.tsv -t active -c pooled_Calculated_reads -n Mageck_B.E.CONV.vivo.guide.level.standard --remove-zero any 

mageck_A.CONV.vivo_genes.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.CONV.vivo.guide.level.standard.gene_summary.txt', header = T)
mageck_A.CONV.vivo_genes.standard$minRRA.fdr <- -log10(apply(mageck_A.CONV.vivo_genes.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.CONV.vivo_genes.standard$minRRA.score <- -log10(apply(mageck_A.CONV.vivo_genes.standard[,c("neg.score","pos.score")], 1, min))
mageck_A.CONV.vivo_genes.standard$minpvalue <- -log10(apply(mageck_A.CONV.vivo_genes.standard[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.CONV.vivo_genes.standard$gene <- str_split_fixed(mageck_A.CONV.vivo_genes.standard$id, '_', 2)[,1] # 12349 genes
nrow(mageck_A.CONV.vivo_genes.standard)

mageck_B.E.CONV.vivo_genes.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_B.E.CONV.vivo.guide.level.standard.gene_summary.txt', header = T)
mageck_B.E.CONV.vivo_genes.standard$minRRA.fdr <- -log10(apply(mageck_B.E.CONV.vivo_genes.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_B.E.CONV.vivo_genes.standard$minRRA.score <- -log10(apply(mageck_B.E.CONV.vivo_genes.standard[,c("neg.score","pos.score")], 1, min))
mageck_B.E.CONV.vivo_genes.standard$minpvalue <- -log10(apply(mageck_B.E.CONV.vivo_genes.standard[,c("neg.p.value","pos.p.value")], 1, min))
mageck_B.E.CONV.vivo_genes.standard$gene <- str_split_fixed(mageck_B.E.CONV.vivo_genes.standard$id, '_', 2)[,1] # 13650 genes
nrow(mageck_B.E.CONV.vivo_genes.standard)


mageck_A.B.CONV_overlap.genes.vivo <- merge(mageck_A.CONV.vivo_genes.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
                                            mageck_B.E.CONV.vivo_genes.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
                                            by = 'id', suffixes = c(".A",".B"))
# #write.table(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_B.E_FC2.invivo_genes_pseudo.standard.SR20.UMIs3.txt', sep = "\t", header = T)

mageck_A.B.CONV_overlap.genes.vivo <- mageck_A.B.CONV_overlap.genes.vivo %>% 
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

names(mageck_A.B.CONV_overlap.genes.vivo)

colors <- c('Other' = 'gray', 'Essential' = 'red', 'Non-Essential' = 'black')

mageck_A.B.CONV_overlap.genes.vivo %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Mouse.ortholog_Controls_550$gene),  size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'Conventional analysis - In vivo pool A vs pool B \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
  xlim(-9,8) + ylim(-9,8) +
  theme(plot.title = element_text(hjust = 0.5))


colors2 <- c('Other' = 'gray', 'Core-Essential' = 'red3', 'Non-Essential' = 'black')
mageck_A.B.CONV_overlap.genes.vivo %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), alpha=0.5, size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Core-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Mouse.ortholog_CoreEssentials_934$gene),  size=1) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'Conventional analysis - In vivo pool A vs pool B \nCore-Essentials: Mouse.ortholog_CoreEssentials_934 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors2) +
  xlim(-9,8.5) + ylim(-9,8.5) +
  theme(plot.title = element_text(hjust = 0.5))

mageck_A.B.CONV_overlap.genes.vivo %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
  stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B'))+ labs(color = ' ', 
                                                                                                                  title = 'Conventional analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  xlim(-9,8.5) + ylim(-9,8.5) +
  theme(plot.title = element_text(hjust = 0.5))

# Vivo AB overlap # ---------- CONVENTIONAL ---------- Correlation of A & B OVERLAP genes ---  With density plots ---------------
mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl <- mageck_A.B.CONV_overlap.genes.vivo %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id)
mageck_A.B.CONV_overlap.genes.vivo

# Plot1 <- mageck_A.B.CONV_overlap.genes.vivo %>%
#   ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
#   geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
#   geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
#   geom_point(aes(color = 'Other'), size = 1) +
#   # geom_point(aes(color = 'Other'), data = mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl$id), size = 1) +
#   geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
#   geom_point(aes(color = 'Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
#   geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') +
#   #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
#   stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
#   xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B')) +
#   #labs(color = ' ', title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
#   theme_pubr(base_size = 14) + scale_color_manual(values = colors) +
#   theme(panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = 'none') +
#   xlim(-9,8.5) + ylim(-9,8.5) +
#   coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
#   theme(plot.title = element_text(hjust = 0.5))

Plot1 <- mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Other'), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'spearman', p.accuracy = 0.001, r.accuracy = 0.01)+
  stat_cor(method = 'pearson',aes(label = paste(..rr.label.., sep = "~`,`~")), label.x = 3) +
  xlab(bquote('Log'[2]~'fold change in vivo pool A')) + ylab(bquote('Log'[2]~'fold change in vivo pool B')) +
  #labs(color = ' ', title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 14) + scale_color_manual(values = colors) + 
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = 'none') +
  xlim(-9,8.5) + ylim(-9,8.5) +   
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  theme(plot.title = element_text(hjust = 0.5))

Plot.Dens.A <- mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.A)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL)  +
  xlim(-9,8.5) +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) +
  theme(legend.position = "none", axis.text.x = element_blank())

Plot.Dens.B <- mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>%
  ggplot(aes(x = neg.lfc.B)) +
  #geom_hline(aes(yintercept=0), color="gray90", linewidth = 1) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_density(aes(color = 'Other'), alpha=0.5, linewidth = 1) +
  geom_density(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>% filter(id %in% Mouse.ortholog_Controls_550$gene), linewidth=1) +
  geom_density(aes(color = 'Essential'), data=mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), linewidth=1) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab(NULL) +
  xlim(-9,8.5) +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors) + 
  theme(legend.position = "none", axis.text.y = element_blank()) + coord_flip()

Plot.Dens.A + plot_spacer() + Plot1 + theme(aspect.ratio = 1)+ Plot.Dens.B + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

#

#!!! Vitro AB overlap # ---------- CONVENTIONAL ---------- Correlation of A & B OVERLAP genes ---------- MAGECK ---------- 

#after normalization
FC1.reads.i.a.vitro.A.norm <- A.B.E.FC2_guide.reads.vitro.Repl.norm %>% select(guide, Gene, Replicate, active.A, inactive.A)
FC1.reads.i.a.vitro.A.norm <- na.omit(FC1.reads.i.a.vitro.A.norm)
colnames(FC1.reads.i.a.vitro.A.norm) <- c("guide", "Gene", "Replicate", "active", "inactive")

FC2_reads.i.a.vitro.B.E.norm <- A.B.E.FC2_guide.reads.vitro.Repl.norm %>% select(guide, Gene, Replicate, active.B.E.norm, inactive.B.E.norm)
FC2_reads.i.a.vitro.B.E.norm <- na.omit(FC2_reads.i.a.vitro.B.E.norm)
colnames(FC2_reads.i.a.vitro.B.E.norm) <- c("guide", "Gene", "Replicate", "active", "inactive")

FC1.reads.i.a.vitro.A.noRepl <- FC1.reads.i.a.vitro.A.norm %>% group_by(guide, Gene) %>% summarise(active=sum(active),inactive=sum(inactive))

FC2_reads.i.a.vitro.B.E.noRepl <- FC2_reads.i.a.vitro.B.E.norm %>% group_by(guide, Gene) %>% summarise(active=sum(active),inactive=sum(inactive))

combined_data_pooled_A <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Calculated.Plasmid.reads_combined_data_pooled_A.txt", sep = "\t", header = T)

A.pooled_guides.vitro.plasmidLibr.nofilt <- bind_rows(left_join(combined_data_pooled_A, FC1.reads.i.a.vitro.A.noRepl, by = "guide"),
                                                      FC1.reads.i.a.vitro.A.noRepl %>% anti_join(combined_data_pooled_A, by = "guide")) %>% replace(is.na(.), 0) #105234 rows
A.pooled_guides.vitro.plasmidLibr.nofilt$sumReads <- rowSums(A.pooled_guides.vitro.plasmidLibr.nofilt[, c("active", "inactive")], na.rm = TRUE)

A.pooled_guides.vitro.plasmidLibr.nofilt <- A.pooled_guides.vitro.plasmidLibr.nofilt %>% group_by(Gene) %>% filter(!(pooled_Calculated_reads == 0| sum(sumReads == 0) >= 5))
A.pooled_guides.vitro.plasmidLibr.nofilt$Gene <- str_split_fixed(A.pooled_guides.vitro.plasmidLibr.nofilt$guide, '_', 2)[, 1]
A.pooled_guides.vitro.plasmidLibr.nofilt <- A.pooled_guides.vitro.plasmidLibr.nofilt %>% select(guide, Gene, active, pooled_Calculated_reads)
Pre_Mageck_A.guide.vitro.plasmid.nofilt <- A.pooled_guides.vitro.plasmidLibr.nofilt %>% mutate_at(vars(active, pooled_Calculated_reads), ~ . + 0.5)
#write.table(Pre_Mageck_A.guide.vitro.plasmid.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_A.guide.vitro.plasmid.nofilt.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#sbatch mageck test -k Pre_Mageck_A.guide.vitro.plasmid.nofilt.tsv -t active -c pooled_Calculated_reads -n Mageck_A.CONV.vitro.guide.level.standard --remove-zero any 

combined_data_pooled_B <- read.table("/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_B.E.screen.2ndTime/_02_Data_lists/Calculated.Plasmid.reads_combined_data_pooled_B.txt", sep = "\t", header = T)

B.E.pooled_guides.vitro.plasmidLibr.nofilt <- bind_rows(left_join(combined_data_pooled_B, FC2_reads.i.a.vitro.B.E.noRepl, by = "guide"),
                                                        FC2_reads.i.a.vitro.B.E.noRepl %>% anti_join(combined_data_pooled_B, by = "guide")) %>% replace(is.na(.), 0) #106466 rows
B.E.pooled_guides.vitro.plasmidLibr.nofilt$sumReads <- rowSums(B.E.pooled_guides.vitro.plasmidLibr.nofilt[, c("active", "inactive")], na.rm = TRUE)
B.E.pooled_guides.vitro.plasmidLibr.nofilt <- B.E.pooled_guides.vitro.plasmidLibr.nofilt %>% group_by(Gene) %>% filter(!(pooled_Calculated_reads == 0| sum(sumReads == 0) >= 5))
B.E.pooled_guides.vitro.plasmidLibr.nofilt$Gene <- str_split_fixed(B.E.pooled_guides.vitro.plasmidLibr.nofilt$guide, '_', 2)[, 1]
Pre_Mageck_B.E.guide.vitro.plasmid.nofilt <- B.E.pooled_guides.vitro.plasmidLibr.nofilt %>% select(guide, Gene, active, pooled_Calculated_reads)
Pre_Mageck_B.E.guide.vitro.plasmid.nofilt <- Pre_Mageck_B.E.guide.vitro.plasmid.nofilt %>% mutate_at(vars(active, pooled_Calculated_reads), ~ . + 0.5)
#write.table(Pre_Mageck_B.E.guide.vitro.plasmid.nofilt, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Pre_Mageck_B.E.guide.vitro.plasmid.nofilt.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#sbatch mageck test -k Pre_Mageck_B.E.guide.vitro.plasmid.nofilt.tsv -t active -c pooled_Calculated_reads -n Mageck_B.E.CONV.vitro.guide.level.standard --remove-zero any 

mageck_A.CONV.vitro_genes.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_A.CONV.vitro.guide.level.standard.gene_summary.txt', header = T)
mageck_A.CONV.vitro_genes.standard$minRRA.fdr <- -log10(apply(mageck_A.CONV.vitro_genes.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_A.CONV.vitro_genes.standard$minRRA.score <- -log10(apply(mageck_A.CONV.vitro_genes.standard[,c("neg.score","pos.score")], 1, min))
mageck_A.CONV.vitro_genes.standard$minpvalue <- -log10(apply(mageck_A.CONV.vitro_genes.standard[,c("neg.p.value","pos.p.value")], 1, min))
mageck_A.CONV.vitro_genes.standard$gene <- str_split_fixed(mageck_A.CONV.vitro_genes.standard$id, '_', 2)[,1] # 12351 genes
nrow(mageck_A.CONV.vitro_genes.standard)

mageck_B.E.CONV.vitro_genes.standard <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_03_Mageck/Mageck_B.E.CONV.vitro.guide.level.standard.gene_summary.txt', header = T)
mageck_B.E.CONV.vitro_genes.standard$minRRA.fdr <- -log10(apply(mageck_B.E.CONV.vitro_genes.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_B.E.CONV.vitro_genes.standard$minRRA.score <- -log10(apply(mageck_B.E.CONV.vitro_genes.standard[,c("neg.score","pos.score")], 1, min))
mageck_B.E.CONV.vitro_genes.standard$minpvalue <- -log10(apply(mageck_B.E.CONV.vitro_genes.standard[,c("neg.p.value","pos.p.value")], 1, min))
mageck_B.E.CONV.vitro_genes.standard$gene <- str_split_fixed(mageck_B.E.CONV.vitro_genes.standard$id, '_', 2)[,1] # 13661 genes
nrow(mageck_B.E.CONV.vitro_genes.standard)


mageck_A.B.CONV_overlap.genes.vitro <- merge(mageck_A.CONV.vitro_genes.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")], 
                                             mageck_B.E.CONV.vitro_genes.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num", "gene")],
                                             by = 'id', suffixes = c(".A",".B"))
# write.table(mageck_A.B.CONV_overlap.genes.vitro, '/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.CONV_overlap.genes.vitro.txt', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# mageck_A.B.CONV_overlap.genes.vitro <- read.table('/Volumes/Esther/NGS_results/Screen_results/20230810_Yumm450R_B.E.screen.2ndTime/Analysis_Yumm450R_A.B.E/_05_Mageck_R_tables/mageck_A.B.CONV_overlap.genes.vitro.txt', sep = "\t", header = T)

mageck_A.B.CONV_overlap.genes.vitro <- mageck_A.B.CONV_overlap.genes.vitro %>% 
  mutate(group= ifelse(id %in% Mouse.ortholog_Controls_550$gene, 'None-Essential', ifelse(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id, 'Essential', 'Other')))

names(mageck_A.B.CONV_overlap.genes.vitro)

mageck_A.B.CONV_overlap.genes.vitro.noOther.nondepl <- mageck_A.B.CONV_overlap.genes.vitro %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id)

#!!! ALL AB overlap plots # ---------- Correlation of A & B OVERLAP genes ---------- different plots ---------- 

colors10 <- c('Other' = 'gold', 'Essential' = 'red', 'Non-Essential' = 'black', 'Neutral Others' = 'gray40')

#previously neutral others lightblue1, alpha = 0.5 and others in gray, alpha 0.7
#colors10 <- c('Other' = 'gray40', 'Essential' = 'red', 'Non-Essential' = 'black', 'Neutral Others' = 'gray40')

# Vitro StAR 
mageck_A.B.FC2_overlap.genes.vitro %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Neutral Others'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id), size=1, alpha = 0.1) +
  geom_point(aes(color = 'Other'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vitro %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  #geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  # stat_cor(method = 'pearson') +
  xlab(bquote('Log'[2]~'fold change in vivo batch 1')) + ylab(bquote('Log'[2]~'fold change in vivo batch 2'))+ 
  #labs(color = ' ', title = 'StAR analysis - In vitro pool A vs pool B (FC2) & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors10) +
  xlim(-9,8.5) + ylim(-9,8.5) +
  theme(plot.title = element_text(hjust = 0.5))

cor.test(mageck_A.B.FC2_overlap.genes.vitro$neg.lfc.A, mageck_A.B.FC2_overlap.genes.vitro$neg.lfc.B)
cor.test(mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl$neg.lfc.A, mageck_A.B.FC2_overlap.genes.vitro.noOther.nondepl$neg.lfc.B)

# Vivo Conventional 
mageck_A.B.CONV_overlap.genes.vivo %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Neutral Others'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id), size=1, alpha = 0.1) +
  geom_point(aes(color = 'Other'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.CONV_overlap.genes.vivo %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  #geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'pearson') +
  xlab(bquote('Log'[2]~'fold change in vivo batch 1')) + ylab(bquote('Log'[2]~'fold change in vivo batch 2'))+ 
  #labs(color = ' ', title = 'Conventional analysis - In vivo pool A vs pool B (FC2) & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors10) +
  xlim(-9,8.5) + ylim(-9,8.5) +
  theme(plot.title = element_text(hjust = 0.5))

cor.test(mageck_A.B.CONV_overlap.genes.vivo$neg.lfc.A, mageck_A.B.CONV_overlap.genes.vivo$neg.lfc.B)
cor.test(mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl$neg.lfc.A, mageck_A.B.CONV_overlap.genes.vivo.noOther.nondepl$neg.lfc.B)

# Vivo StAR 
mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Neutral Others'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id), size=1, alpha = 0.1) +
  geom_point(aes(color = 'Other'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3 %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  #geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'pearson') +
  xlab(bquote('Log'[2]~'fold change in vivo batch 1')) + ylab(bquote('Log'[2]~'fold change in vivo batch 2'))+ 
  #labs(color = ' ', title = 'StAR analysis - In vivo pool A vs pool B (FC2) \nSR20.UMIs3 & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors10) +
  xlim(-9,8.5) + ylim(-9,8.5) +
  theme(plot.title = element_text(hjust = 0.5))

cor.test(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3$neg.lfc.A, mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3$neg.lfc.B)
cor.test(mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl$neg.lfc.A, mageck_A.B.FC2_overlap.genes.vivo.SR20.UMIs3.noOther.nondepl$neg.lfc.B)

# Vitro Conventional 
mageck_A.B.CONV_overlap.genes.vitro %>%
  ggplot(aes(x = neg.lfc.A, y = neg.lfc.B)) +
  geom_hline(aes(yintercept=0), color="gray20", linewidth = 0.5) +
  geom_vline(aes(xintercept=0), color="gray20", linewidth = 0.5) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'Neutral Others'), data=mageck_A.B.CONV_overlap.genes.vitro %>% filter(id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id), size=1, alpha = 0.1) +
  geom_point(aes(color = 'Other'), data=mageck_A.B.CONV_overlap.genes.vitro %>% filter(!id %in% mageck_A.B.FC2_overlap.genes.vitro.other.nondepl$id), size = 1) +
  geom_point(aes(color = 'Non-Essential'), data=mageck_A.B.CONV_overlap.genes.vitro %>% filter(id %in% Mouse.ortholog_Controls_550$gene), size=1) +
  geom_point(aes(color = 'Essential'), data=mageck_A.B.CONV_overlap.genes.vitro %>% filter(id %in% Vitro.A.B.E.FC2_depleting.min3.RRAsc10$id), size=1) +
  #geom_smooth(method = lm, formula = y ~ x, se = FALSE, size = 0.75, color = 'black') + 
  #stat_cor(method = 'pearson') +
  xlab(bquote('Log'[2]~'fold change in vitro batch 1')) + ylab(bquote('Log'[2]~'fold change in vitro batch 2'))+ 
  #labs(color = ' ', title = 'Conventional analysis - In vitro pool A vs pool B (FC2) & \nEssentials: Vitro.depl.min3.RRAsc10 \nNon-Essential: Mouse.ortholog_Controls_550') +
  theme_pubr(base_size = 12) + scale_color_manual(values = colors10) +
  xlim(-9,8.5) + ylim(-9,8.5) +
  theme(plot.title = element_text(hjust = 0.5))

cor.test(mageck_A.B.CONV_overlap.genes.vitro$neg.lfc.A, mageck_A.B.CONV_overlap.genes.vitro$neg.lfc.B)
cor.test(mageck_A.B.CONV_overlap.genes.vitro.noOther.nondepl$neg.lfc.A, mageck_A.B.CONV_overlap.genes.vitro.noOther.nondepl$neg.lfc.B)


#
#!!! rank - ALL lines of sgRNA reproducibility plots # ---------- ---------- 
A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5
A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5

A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5_select <- A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5[c("Rank", "average_guides")]
A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5_select <- A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5[c("Rank", "average_guides")]
A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5_select <- A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5[c("Rank", "average_guides")]

# Merge the selected columns from all three data frames into one
combined <- merge(merge(A.B.E.FC2_guide.reads.vitro.noRepl.pseudo.rank.5_select, A.B.E_FC2.pooled_guides.vivo.plasmidLibr.nofilt.groups.rank.5_select, by = "Rank"), A.B.E_FC2.pooled.guideLevel.SR20.UMIs3.pseudo.rank.5_select, by = "Rank")
colnames(combined) <- c('Rank', 'Vitro.avg.guides' , 'Conv.avg.guides', 'StAR.avg.guides' )

combined_long <- reshape2::melt(combined, id.vars = "Rank", variable.name = "Method", value.name = "Avg_guides")

# All lines together
ggplot(combined_long, aes(x = Rank, y = Avg_guides, color = Method)) +
  geom_point(data = combined_long[combined_long$Method == "Vitro.avg.guides",], aes(x = Rank, y = Avg_guides)) +
  geom_point(data = combined_long[combined_long$Method == "Conv.avg.guides",], aes(x = Rank, y = Avg_guides)) +
  geom_point(data = combined_long[combined_long$Method == "StAR.avg.guides",], aes(x = Rank, y = Avg_guides)) +
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene",
       title = 'Reproducibility of sgRNAs') +
  xlim(0,3e+03)+
  scale_x_continuous(labels = scales::label_scientific()) +
  scale_color_manual(values =c('Vitro.avg.guides' = "#999999", 'Conv.avg.guides' = "#E69F00", 'StAR.avg.guides' = "#56B4E9")) +
  theme_pubr()

## ## or all plots separate

# Vitro StAR 
ggplot(combined_long_vitro, aes(x = Rank, y = Avg_guides, color = Method)) +
  geom_point(data = combined_long_vitro[combined_long_vitro$Method == "average_guides.vitro.rdm",], aes(x = Rank, y = Avg_guides)) +
  geom_point(data = combined_long_vitro[combined_long_vitro$Method == "average_guides.vitro",], aes(x = Rank, y = Avg_guides)) +
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene",
       title = 'Reproducibility of sgRNAs - Vitro StAR') +
  scale_x_continuous(labels = scales::label_scientific()) +
  xlim(0,3e+04) + ylim(1,2.95) +
  scale_color_manual(values =c('average_guides.vitro' = "#000000", 'average_guides.vitro.rdm' = "#999999")) +
  theme_pubr()

# Vivo Conventional
ggplot(combined_long_vivo.Conv, aes(x = Rank, y = Avg_guides, color = Method)) +
  geom_point(data = combined_long_vivo.Conv[combined_long_vivo.Conv$Method == "average_guides.vvConv.rdm",], aes(x = Rank, y = Avg_guides)) +
  geom_point(data = combined_long_vivo.Conv[combined_long_vivo.Conv$Method == "average_guides.vvConv",], aes(x = Rank, y = Avg_guides)) +
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene",
       title = 'Reproducibility of sgRNAs - Vivo Conventional') +
  scale_x_continuous(labels = scales::label_scientific()) +
  xlim(0,3e+04) + ylim(1,1.86) +
  scale_color_manual(values =c('average_guides.vvConv' = "#000000", 'average_guides.vvConv.rdm' = "#999999")) +
  theme_pubr()

# Vivo StAR
ggplot(combined_long_vivo.StAR, aes(x = Rank, y = Avg_guides, color = Method)) +
  geom_point(data = combined_long_vivo.StAR[combined_long_vivo.StAR$Method == "average_guides.vvStAR.rdm",], aes(x = Rank, y = Avg_guides)) +
  geom_point(data = combined_long_vivo.StAR[combined_long_vivo.StAR$Method == "average_guides.vvStAR",], aes(x = Rank, y = Avg_guides)) +
  labs(x = "sgRNAs sorted by rank (x)",
       y = "#sgRNAs/gene",
       title = 'Reproducibility of sgRNAs - Vivo StAR') +
  scale_x_continuous(labels = scales::label_scientific()) +
  xlim(0,3e+04) + ylim(1,1.86) +
  scale_color_manual(values =c('average_guides.vvStAR' = "#000000", 'average_guides.vvStAR.rdm' = "#999999")) +
  theme_pubr()

#!!! ROC # ---------- Vitro.vivo.Conv/StAR ---------- ---------- 
pseudo_roc = grouped_data.vitro_StAR %>% roc(group , StAR.LFC)
star_roc = grouped_data_StAR %>% roc(group , StAR.LFC)
conv_roc = grouped_data_Conv %>% roc(group , Conv.LFC)
pseudo_roc_random =  grouped_data.vitro_StAR %>% ungroup() %>%
  mutate(group = sample(group) ,
         StAR.LFC= sample(StAR.LFC)) %>%
  roc(group , StAR.LFC)

# roc-auc values
auc(pseudo_roc)
auc(star_roc)
auc(pseudo_roc_random)
auc(conv_roc)

roc_curves = rbind(coords(pseudo_roc, ret = c('tpr', 'fpr','precision','recall')) %>% mutate(sample = 'Vitro_Theoretical max') ,
                   coords(star_roc, ret = c( 'tpr', 'fpr','precision','recall')) %>% mutate(sample = 'StAR') ,
                   coords(pseudo_roc_random, ret = c( 'tpr', 'fpr','precision','recall')) %>% mutate(sample = 'Random'),
                   coords(conv_roc, ret = c( 'tpr', 'fpr','precision','recall')) %>% mutate(sample = 'Conventional'))

ROC_plot <- 
  roc_curves %>%
  ggplot(aes(x = fpr , y = tpr, color = sample)) +
  geom_line(size = 0.75) +
  theme_pubr()+
  scale_color_manual(values = c("StAR" = "#DA1B5D", "Vitro_Theoretical max" = "black", 'Random' = 'gray', 'Conventional' = 'Blue')) +
  labs(x = 'False positive rate', y = 'True positive rate') +
  theme(legend.position = 'bottom')

ROC_precision_plot <- 
  roc_curves %>%
  ggplot(aes(x = recall, y = precision,  color = sample)) +
  geom_line(size = 0.75) +
  theme_pubr()+
  scale_color_manual(values = c("StAR" = "#DA1B5D", "Vitro_Theoretical max" = "black", 'Random' = 'gray', 'Conventional' = 'Blue')) +
  labs(x = 'Recall', y = 'Precision') +
  lims(y = c(0,1)) +
  theme(legend.position = 'bottom')


ROC_plot + ROC_precision_plot  +
  plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = 'top')

#

#!!! AUC # ---------- Vitro.vivo.Conv/StAR ---------- ---------- 




ggplot() +
  geom_line(data = A.B.E.vitro.guide.dAUC.rank.LFCStAR, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group), linetype = 'dashed', size = 0.75) +
  geom_line(data = A.B.E.guide.dAUC.rank.LFCStAR, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group), linetype = "solid", size = 0.75) +
  geom_line(data = A.B.E.guide.dAUC.rank.LFCConv, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group), linetype = "dotted", size = 0.75) +
  geom_hline(yintercept = max(A.B.E.guide.dAUC.rank.LFCConv$rank_group_prct), linetype = "solid", color = "black") +
  geom_vline(xintercept = max(A.B.E.guide.dAUC.rank.LFCConv$all_rank_perc_sgRNA), linetype = "solid", color = "black") +
  theme_pubr()+
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "black", "Other" = "gray")) +
  labs(y = "Cumulative fraction", x = "Percentage rank of sgRNAs") +
  theme(legend.position = 'bottom') +
  labs(title = "dAUC Vitro, Vivo-Conventional, Vivo-StAR") 


ggplot() +
  geom_line(data = A.B.E.vitro.guide.dAUC.rank.LFCStAR, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group, linetype = 'dashed'), size = 0.75) +
  geom_line(data = A.B.E.guide.dAUC.rank.LFCStAR, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group, linetype = "solid"), linetype = "solid", size = 0.75) +
  geom_line(data = A.B.E.guide.dAUC.rank.LFCConv, aes(x = all_rank_perc_sgRNA, y  = rank_group_prct, color = group, linetype = "dotted"), size = 0.75) +
  geom_hline(yintercept = max(A.B.E.guide.dAUC.rank.LFCConv$rank_group_prct), linetype = "solid", color = "black") +
  geom_vline(xintercept = max(A.B.E.guide.dAUC.rank.LFCConv$all_rank_perc_sgRNA), linetype = "solid", color = "black") +
  theme_pubr()+
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "black", "Other" = "gray")) +
  labs(y = "Cumulative fraction", x = "Percentage rank of sgRNAs") +
  theme(legend.position = 'bottom') +
  labs(title = "dAUC Vitro, Vivo-Conventional, Vivo-StAR") 






ggplot() +
  geom_line(data = A.B.E.vitro.guide.dAUC.rank.LFCStAR, aes(x = ifelse(rank_group_prct == 0, 0, ifelse(rank_group_prct == 100, 1, all_rank_perc_sgRNA)), 
                                                            y = ifelse(all_rank_perc_sgRNA == 0, 0, ifelse(all_rank_perc_sgRNA == 100, 1, rank_group_prct)), 
                                                            color = group), 
            linetype = 'dashed', size = 0.75) +
  # geom_line(data = A.B.E.guide.dAUC.rank.LFCStAR, aes(x = ifelse(rank_group_prct == 0, 0, ifelse(rank_group_prct == 100, 1, all_rank_perc_sgRNA)), 
  #                                                     y = ifelse(all_rank_perc_sgRNA == 0, 0, ifelse(all_rank_perc_sgRNA == 100, 1, rank_group_prct)), 
  #                                                     color = group), 
  #           linetype = "solid", size = 0.75) +
  # geom_line(data = A.B.E.guide.dAUC.rank.LFCConv, aes(x = ifelse(rank_group_prct == 0, 0, ifelse(rank_group_prct == 100, 1, all_rank_perc_sgRNA)), 
  #                                                     y = ifelse(all_rank_perc_sgRNA == 0, 0, ifelse(all_rank_perc_sgRNA == 100, 1, rank_group_prct)), 
  #                                                     color = group), 
  #           linetype = "dotted", size = 0.75) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  theme_pubr() +
  scale_color_manual(values = c("Essential" = "red", "Non-Essential" = "black", "Other" = "gray")) +
  labs(y = "Cumulative fraction", x = "Percentage rank of sgRNAs") +
  theme(legend.position = 'bottom') +
  labs(title = "dAUC Vitro, Vivo-Conventional, Vivo-StAR")


auc(A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$all_rank_perc_sgRNA,A.B.E.vitro.guide.dAUC.rank.LFCStAR.Ess$rank_group_prct, type = 'spline')

head(A.B.E.vitro.guide.dAUC.rank.LFCStAR)



# VIVO # ---------- Ratio A screen and B.E screen --------------------------------------------------------

read_counts.A <- FC1_FC2_poly_hop.filtered.reads.vivo_0.001.A %>%
  # group_by(Sample) %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.A <- read_counts.A %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.A <- read_counts.A %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.A <- read_counts.A %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.A <- read_counts.A %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.A <- read_counts.A %>%
  mutate(Data = "Vivo_counts_A")


Vivo.NonEss_A_guide_UMI <- FC1_FC2_poly_hop.filtered.reads.vivo_0.001.A %>%
  filter((Gene) %in% Mouse.ortholog_Controls_550$gene)
read_counts.Ctrl.A <- Vivo.NonEss_A_guide_UMI %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))
read_counts.Ctrl.A <- read_counts.Ctrl.A %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.Ctrl.A <- read_counts.Ctrl.A %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.Ctrl.A <- read_counts.Ctrl.A %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.Ctrl.A <- read_counts.Ctrl.A %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.Ctrl.A <- read_counts.Ctrl.A %>%
  mutate(Data = "Vivo_counts_ctrls_A")

Vivo.NonEss_A_guide_UMI.Ratio <- Vivo.NonEss_A_guide_UMI %>% mutate(Perc.active = round((active / sumReads)*100, 2),
                                                                    Perc.inactive = round((inactive / sumReads)*100, 2))
read_counts.Ctrl.A.median <- Vivo.NonEss_A_guide_UMI.Ratio %>%
  summarise(
    Perc.active = median(Perc.active),
    Perc.inactive = median(Perc.inactive), 
    Reads_Active_0 = sum(active == 0),
    Reads_Inactive_0 = sum(inactive == 0))
read_counts.Ctrl.A.median <-read_counts.Ctrl.A.median %>% mutate(Total_UMIs = nrow(Vivo.NonEss_A_guide_UMI.Ratio))
read_counts.Ctrl.A.median <- read_counts.Ctrl.A.median %>% mutate(Data = "Vivo_counts_ctrls_A.median")
head(read_counts.Ctrl.A.median)

read_counts.B.E <- FC2_poly_hop.filtered.reads.vivo_0.001.B.E %>%
  # group_by(Sample) %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.B.E <- read_counts.B.E %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.B.E <- read_counts.B.E %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.B.E <- read_counts.B.E %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.B.E <- read_counts.B.E %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.B.E <- read_counts.B.E %>%
  mutate(Data = "Vivo_counts_B.E")


Vivo.NonEss_B.E_guide_UMI <- FC2_poly_hop.filtered.reads.vivo_0.001.B.E %>%
  filter((Gene) %in% Mouse.ortholog_Controls_550$gene)
read_counts.Ctrl.B.E <- Vivo.NonEss_B.E_guide_UMI %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))
read_counts.Ctrl.B.E <- read_counts.Ctrl.B.E %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.Ctrl.B.E <- read_counts.Ctrl.B.E %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.Ctrl.B.E <- read_counts.Ctrl.B.E %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.Ctrl.B.E <- read_counts.Ctrl.B.E %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.Ctrl.B.E <- read_counts.Ctrl.B.E %>%
  mutate(Data = "Vivo_counts_ctrls_B.E")

Vivo.NonEss_B.E_guide_UMI.Ratio <- Vivo.NonEss_B.E_guide_UMI %>% mutate(Perc.active = round((active / sumReads)*100, 2),
                                                                        Perc.inactive = round((inactive / sumReads)*100, 2))
read_counts.Ctrl.B.E.median <- Vivo.NonEss_B.E_guide_UMI.Ratio %>%
  summarise(
    Perc.active = median(Perc.active),
    Perc.inactive = median(Perc.inactive),
    Reads_Active_0 = sum(active == 0),
    Reads_Inactive_0 = sum(inactive == 0)) 
read_counts.Ctrl.B.E.median <-read_counts.Ctrl.B.E.median %>% mutate(Total_UMIs = nrow(Vivo.NonEss_B.E_guide_UMI.Ratio))
read_counts.Ctrl.B.E.median <- read_counts.Ctrl.B.E.median %>% mutate(Data = "Vivo_counts_ctrls_B.E.median")

head(read_counts.Ctrl.B.E.median)



read_counts.A.B.E <- A.B.E_FC2.pooled_most_filters.ia0 %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.A.B.E <- read_counts.A.B.E %>%
  mutate(Data = "Vivo_counts_A.B.E")


read_counts.A.B.E.sumfilter <- A.B.E_FC2.pooled.guide_UMI_Repl.SR20.UMIs3 %>%
  summarise(
    Reads_Active = sum(active),
    Reads_Inactive = sum(inactive))

read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Ratio_Active_to_Inactive = round(Reads_Active / Reads_Inactive, 1))
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(SumUMIs = Reads_Active + Reads_Inactive)
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Perc.active = round((Reads_Active / SumUMIs)*100, 2),
         Perc.inactive = round((Reads_Inactive / SumUMIs)*100, 2))
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Perc.total = Perc.active + Perc.inactive)
read_counts.A.B.E.sumfilter <- read_counts.A.B.E.sumfilter %>%
  mutate(Data = "Vivo_counts_A.B.E.sumfilter")

Ratio.read.counts <- rbind(read_counts.A, read_counts.B.E, read_counts.A.B.E,read_counts.A.B.E.sumfilter)
Ratio.read.counts.vt.vv <- rbind(read_counts.A, read_counts.B.E, read_counts.A.B.E,read_counts.A.B.E.sumfilter, vitro.read_counts.A, vitro.read_counts.B.E)

Ratio.read.counts <- rbind(read_counts.A, read_counts.B.E, read_counts.A.B.E,read_counts.A.B.E.sumfilter)



