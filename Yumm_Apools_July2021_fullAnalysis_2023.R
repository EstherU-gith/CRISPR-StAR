###------------------------------------------------ load packages --------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggrepel)
library(xlsx)
library(dplyr)
#library(biomaRt)

###------------------------------------------------ Load control lists --------------------------------------------------------
#setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/')
Essentials<-read.table('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/coreessentials.tsv', header=T )
Controls<- read.table('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/nonessentials.tsv', header=T )

###------------------------------------------------ load data FC1 --------------------------------------------------------
#setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/mapping/MDs_noShadows/')
unfiltered.reads.FC1 <- list.files('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/mapping/MDs_noShadows/') %>%
  map_df(~fread(.x, header = TRUE) %>% 
           mutate(vitro_vivo = ifelse(str_detect(index, 'vitro'), 'in_vitro','in_vivo'),
                  Replicate = str_split_fixed(index, pattern='_', n=4)[,3],
                  inactive_active = ifelse(str_detect(index, 'inactive'), 'inactive','active'),
                  guide_UMI = paste(guide, UMI, sep = '_'),
                  Gene = str_split_fixed(guide, '_', 2)[,1]))
#write.table(unfiltered.reads.FC1, "/Volumes/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/unfiltered.reads.FC1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
unfiltered.reads.FC1 <- read.table("/Volumes/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/unfiltered.reads.FC1.txt", header = T)

###------------------------------------------------ load data FC2 --------------------------------------------------------
#setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/ngs_raw_flowcell2/mapping/MDs_noShadows_shadows_discarded/')
unfiltered.reads.FC2 <- list.files('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/ngs_raw_flowcell2/mapping/MDs_noShadows_shadows_discarded/') %>%
  map_df(~fread(.x, header = TRUE) %>% 
           mutate(vitro_vivo = ifelse(str_detect(index, 'vitro'), 'in_vitro','in_vivo'),
                  Replicate = str_split_fixed(index, pattern='_', n=4)[,3],
                  inactive_active = ifelse(str_detect(index, 'inactive'), 'inactive','active'),
                  guide_UMI = paste(guide, UMI, sep = '_'),
                  Gene = str_split_fixed(guide, '_', 2)[,1]))
#write.table(unfiltered.reads.FC2, "/Volumes/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/unfiltered.reads.FC2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
unfiltered.reads.FC2 <- read.table("/Volumes/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/unfiltered.reads.FC2.txt", header = T)


#IN VITRO#--------------------------------- IN VITRO----------------------- COLLAPSE TO GUIDE LEVEL ------------ + pseudo counts -------------------------------------
#no UMIs / guide level dataset
only.guides.invitro <- unfiltered.reads.FC1 %>%
  group_by(Replicate, guide, inactive_active, vitro_vivo = if_else(str_detect(index, "vitro"), "in_vitro", "in_vivo")) %>%
  summarise(reads = sum(reads)) %>%
  filter(vitro_vivo == "in_vitro") %>% mutate(Gene = str_split_fixed(guide, '_', 2)[,1])

a.i.guides.only.invitro <- reshape2::dcast(only.guides.invitro %>% dplyr::select(guide, Replicate, inactive_active, vitro_vivo, reads, Gene),
                                           guide + Replicate + vitro_vivo + Gene ~ inactive_active, fill = 0 , value.var = 'reads')
a.i.guides.only.invitro_wide <- dcast(setDT(a.i.guides.only.invitro), guide + Gene ~ Replicate, value.var = c('active','inactive'), fill = 0) %>% data.frame()

a.i.guides.only.invitro.noRepl <- a.i.guides.only.invitro %>%
  group_by(guide, Gene, active, inactive) %>%
  summarise(active = sum(active), inactive = sum(inactive))

a.i.guides.only.invitro.noRepl <- a.i.guides.only.invitro %>%
  group_by(guide, Gene) %>%
  summarise(active = sum(active), inactive = sum(inactive))
#write.table(a.i.guides.only.invitro.noRepl,  "/Volumes/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/a.i.guides.only.invitro.noRepl.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# a.i.guides.only.invitro_wide$sumReads_Blue <- a.i.guides.only.invitro_wide$active_Blue + a.i.guides.only.invitro_wide$inactive_Blue
# a.i.guides.only.invitro_wide$sumReads_Green <- a.i.guides.only.invitro_wide$active_Green + a.i.guides.only.invitro_wide$inactive_Green
# a.i.guides.only.invitro_wide$sumReads_Red <- a.i.guides.only.invitro_wide$active_Red + a.i.guides.only.invitro_wide$inactive_Red
# a.i.guides.only.invitro_wide.filtered <- a.i.guides.only.invitro_wide %>% filter(sumReads_Blue >10, sumReads_Green >10, sumReads_Red >10)
# 
#             # ---------------------------------- ADDITION OF PSEUDO COUNT --------------------------
# 
# a.i.guides.only.invitro.filtered.pseudo_wide <- a.i.guides.only.invitro_wide.filtered
# a.i.guides.only.invitro.filtered.pseudo_wide$reads.a_Blue <- a.i.guides.only.invitro.filtered.pseudo_wide$reads.a_Blue + 0.5 
# a.i.guides.only.invitro.filtered.pseudo_wide$reads.a_Green <- a.i.guides.only.invitro.filtered.pseudo_wide$reads.a_Green + 0.5 
# a.i.guides.only.invitro.filtered.pseudo_wide$reads.a_Red <- a.i.guides.only.invitro.filtered.pseudo_wide$reads.a_Red + 0.5 
# a.i.guides.only.invitro.filtered.pseudo_wide$reads.i_Blue <- a.i.guides.only.invitro.filtered.pseudo_wide$reads.i_Blue + 0.5 
# a.i.guides.only.invitro.filtered.pseudo_wide$reads.i_Green <- a.i.guides.only.invitro.filtered.pseudo_wide$reads.i_Green + 0.5 
# a.i.guides.only.invitro.filtered.pseudo_wide$reads.i_Red <- a.i.guides.only.invitro.filtered.pseudo_wide$reads.i_Red + 0.5 

            # ---------------------------------- ADDITION OF PSEUDO COUNT ------------ Replicates together--------------
# a.i.guides.only.invitro.pseudo.noRepl <- a.i.guides.only.invitro.noRepl
# a.i.guides.only.invitro.pseudo.noRepl$reads.a <- a.i.guides.only.invitro.pseudo.noRepl$reads.a + 0.5 
# a.i.guides.only.invitro.pseudo.noRepl$reads.i <- a.i.guides.only.invitro.pseudo.noRepl$reads.i + 0.5

#IN VIVO#--------------------------------- IN VIVO----------------------- Combine Flowcells ------------ + pseudo counts -------------------------------------

FC1.reads.i.a <- reshape2::dcast(unfiltered.reads.FC1 %>% select(guide, UMI, Replicate, inactive_active, vitro_vivo, reads),
                       guide + UMI + Replicate + vitro_vivo ~ inactive_active, fill = 0 , value.var = 'reads')
FC2.reads.i.a <- reshape2::dcast(unfiltered.reads.FC2 %>% select(guide, UMI, Replicate, inactive_active, vitro_vivo, reads),
                       guide + UMI + Replicate + vitro_vivo ~ inactive_active, fill = 0 , value.var = 'reads')
FC1_FC2_reads.i.a <- merge(FC1.reads.i.a, FC2.reads.i.a, by = c('guide', 'UMI', 'Replicate', 'vitro_vivo'), all = T, suffixes = c(".FC1",".FC2"), fill = 0)
FC1_FC2_reads.i.a$guide_UMI <- paste(FC1_FC2_reads.i.a$guide, FC1_FC2_reads.i.a$UMI, sep = '_')
FC1_FC2_reads.i.a$Gene<- str_split_fixed(FC1_FC2_reads.i.a$guide, '_', 2)[,1] #FC1_FC2_reads$guide %>% str_sub(1, -3) 
FC1_FC2_reads.i.a <- FC1_FC2_reads.i.a[,c('guide', 'UMI', 'guide_UMI', 'Gene', 'active.FC1', 'inactive.FC1', 'active.FC2', 'inactive.FC2', 'Replicate', 'vitro_vivo') ]
setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/')
write.table(FC1_FC2_reads.i.a, "/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/FC1_FC2_reads.i.a.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(FC1.reads.i.a, "/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/FC1.reads.i.a.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(FC2.reads.i.a, "/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/FC2.reads.i.a.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


FC1_FC2_reads.invivo <- FC1_FC2_reads.i.a %>% filter(vitro_vivo == 'in_vivo')
FC1_FC2_reads.invivo$sumreadsFC1 <- FC1_FC2_reads.invivo$active.FC1 + FC1_FC2_reads.invivo$inactive.FC1
FC1_FC2_reads.invivo$sumreadsFC2 <- FC1_FC2_reads.invivo$active.FC2 + FC1_FC2_reads.invivo$inactive.FC2
FC1_FC2_reads.invivo$sumReads <- FC1_FC2_reads.invivo$sumreadsFC1 + FC1_FC2_reads.invivo$sumreadsFC2

FC1_FC2_reads.invivo$active <- FC1_FC2_reads.invivo$active.FC1 + FC1_FC2_reads.invivo$active.FC2 
FC1_FC2_reads.invivo$inactive <- FC1_FC2_reads.invivo$inactive.FC1 + FC1_FC2_reads.invivo$inactive.FC2 
FC1_FC2_reads.invivo.pooled <- FC1_FC2_reads.invivo %>% select(guide, UMI, guide_UMI, Gene, Replicate, active, inactive)
FC1_FC2_reads.invivo.pooled.nopolyUMIs <- FC1_FC2_reads.invivo.pooled %>% filter(!grepl('GGGGG|CCCCCCC|AAAAAAA|TTTTTTT', UMI))
FC1_FC2_reads.invivo.pooled.nopolyUMIs.noRepl <- FC1_FC2_reads.invivo.pooled.nopolyUMIs %>%
  group_by(guide, UMI, guide_UMI, Gene) %>%
  summarise(active = sum(active), inactive = sum(inactive))

FC1_FC2_reads.invivo.pooled.nopolyUMIs.noRepl.guideLevel <- FC1_FC2_reads.invivo.pooled.nopolyUMIs %>%
  group_by(guide, Gene) %>%
  summarise(active = sum(active), inactive = sum(inactive))

FC1_FC2_reads.invivo.pooled.nopolyUMIs.noRepl.guideLevel
#merged.Libr.vitro.Screen <- read.table('/Volumes/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/merged.Libr.vitro.Screen.txt', header = T)

Vivo_Libr_vitro_Reads <- merge(FC1_FC2_reads.invivo.pooled.nopolyUMIs.noRepl.guideLevel, merged.Libr.vitro.Screen, by = 'guide')
Vivo_Libr_vitro_Reads = Vivo_Libr_vitro_Reads %>% select("guide", "Gene.x", "active", "inactive.x", "pooled_reads", "inactive.y", "ratio", "Calculated.Reads")
colnames(Vivo_Libr_vitro_Reads) <- c('guide', 'Gene', 'active.vivo', 'inactive.vivo', 'reads.libr', 'inactive.vitro', 'ratio.vitro', 'Calculated.Libr.Reads')

#write.table(Vivo_Libr_vitro_Reads, '/Volumes/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/Vivo_Libr_vitro_Reads.txt', sep = "\t", row.names = F, col.names = T, quote = F)


Vivo_Libr_vitro_Reads$Conv.LFC <- Vivo_Libr_vitro_Reads$active.vivo/ Vivo_Libr_vitro_Reads$inactive.vitro
Vivo_Libr_vitro_Reads$StAR.LFC <- Vivo_Libr_vitro_Reads$active.vivo/ Vivo_Libr_vitro_Reads$inactive.vivo

Vivo_Libr_vitro_Reads.filtered <- Vivo_Libr_vitro_Reads %>% filter(active.vivo >5, inactive.vivo >5)
Vivo_Libr_vitro_Reads.filtered.Ctrls <- Vivo_Libr_vitro_Reads.filtered %>% filter(toupper(Gene) %in% Controls$gene)
Vivo_Libr_vitro_Reads.filtered.Ctrls <- arrange(Vivo_Libr_vitro_Reads.filtered.Ctrls,  by =  Conv.LFC)
Vivo_Libr_vitro_Reads.filtered.Ess <- Vivo_Libr_vitro_Reads.filtered %>% filter((Gene) %in% Essentials$gene)
Vivo_Libr_vitro_Reads.filtered.Ess <- arrange(Vivo_Libr_vitro_Reads.filtered.Ess,  by =  Conv.LFC)

Vivo_Libr_vitro_all <- merge(a.i.guides.only.invitro.noRepl, Vivo_Libr_vitro_Reads, by = 'guide')
colnames(Vivo_Libr_vitro_all) <- c("guide", "Gene", "active.vitro", "inactive.vitro", "Gene.y", "active.vivo", "inactive.vivo", "reads.libr",  
                                   "inactive.vitro", "ratio.vitro", "Calculated.Libr.Reads", "Conv.LFC", "StAR.LFC")
Vivo_Libr_vitro_all <- Vivo_Libr_vitro_all %>% select("guide", "Gene", "active.vitro", "inactive.vitro", "active.vivo", "inactive.vivo", "reads.libr",  
                               "Conv.LFC", "StAR.LFC")

Vivo_Libr_vitro_all$Vitro.LFC <- Vivo_Libr_vitro_all$active.vitro/Vivo_Libr_vitro_all$inactive.vitro
Vivo_Libr_vitro_all.unfiltered <- Vivo_Libr_vitro_all

Vivo_Libr_vitro_all <- Vivo_Libr_vitro_all %>% filter(active.vivo >5, inactive.vivo >5)

Vivo_Libr_vitro_all$vitro.depl0.2 <- ifelse((Vivo_Libr_vitro_all$Vitro.LFC < 0.2), TRUE, '-')
Vivo_Libr_vitro_all$vitro.depl0.1 <- ifelse((Vivo_Libr_vitro_all$Vitro.LFC < 0.1), TRUE, '-')
Vivo_Libr_vitro_all.unfiltered$vitro.depl0.2 <- ifelse((Vivo_Libr_vitro_all.unfiltered$Vitro.LFC < 0.2), TRUE, '-')
Vivo_Libr_vitro_all.unfiltered$vitro.depl0.1 <- ifelse((Vivo_Libr_vitro_all.unfiltered$Vitro.LFC < 0.1), TRUE, '-')
Vitro.depleting0.2 <- Vivo_Libr_vitro_all.unfiltered %>% filter(vitro.depl0.2 == TRUE)
Vitro.depleting0.1 <- Vivo_Libr_vitro_all %>% filter(vitro.depl0.1 == TRUE)
Vitro.depleting0.2 <- arrange(Vitro.depleting0.2,  by =  Conv.LFC)
write.table(Vitro.depleting0.2, '/Volumes/Esther/NGS_results/Scripts_gene_lists/Yumm450R.Vitro.depleting0.2.txt',quote = F, row.names = F, col.names = T, sep = '\t')


Vivo_Libr_vitro_all.unfiltered$Controls <- ifelse((toupper(Vivo_Libr_vitro_all.unfiltered$Gene) %in% Controls$gene), TRUE, '-')

Vivo_Libr_vitro_all.reads.0.5 <- Vivo_Libr_vitro_all.unfiltered %>% filter(active.vivo <6)
Vivo_Libr_vitro_all.reads.0 <- Vivo_Libr_vitro_all.unfiltered %>% filter(active.vivo <1)
Vivo_Libr_vitro_all.reads.ai.0.5 <- Vivo_Libr_vitro_all.unfiltered %>% filter(active.vivo <6 ,inactive.vivo <6)

count




             # ---------------------------------- FILTERING sumreads and poly UMIs ------------------
FC1_FC2_reads.invivo.30 <- FC1_FC2_reads.invivo %>% filter(sumReads >= 30)
FC1_FC2_reads.invivo.30.nopolyUMIs <- FC1_FC2_reads.invivo.30 %>% filter(!grepl('GGGGG|CCCCCCC|AAAAAAA|TTTTTTT', UMI))

            # ---------------------------------- FILTERING UMI hopping -----------------------------
FC1_FC2_reads.invivo.30.nopolyUMIs <- FC1_FC2_reads.invivo.30.nopolyUMIs %>% arrange(desc(sumReads))
temp <- FC1_FC2_reads.invivo.30.nopolyUMIs
FC1_FC2_filtered_UMI <- data.frame()

while(temp$sumReads[1] > 100000){
  jumped_UMI <- temp$UMI[1]
  all_jumped_UMI <- temp %>% filter(UMI == jumped_UMI) %>% mutate(ratio = sumReads/max(sumReads))
  
  FC1_FC2_filtered_UMI <- rbind(FC1_FC2_filtered_UMI, all_jumped_UMI %>% filter(ratio > 0.001) %>% select(-ratio))
  temp <- temp %>% filter(UMI != jumped_UMI)
}
FC1_FC2_filtered_UMI <- rbind(FC1_FC2_filtered_UMI, temp)

           # ---------------------------------- ADDITION OF PSEUDO COUNT ----------------------------
FC1_FC2_filtered_UMI.pseudo <- FC1_FC2_filtered_UMI
FC1_FC2_filtered_UMI.pseudo$reads.a <- FC1_FC2_filtered_UMI$active.FC1 + FC1_FC2_filtered_UMI$active.FC2 + 0.5
FC1_FC2_filtered_UMI.pseudo$reads.i <- FC1_FC2_filtered_UMI$inactive.FC1 + FC1_FC2_filtered_UMI$inactive.FC2 + 0.5

#IN VIVO#--------------------------------- IN VIVO----------------------- on Guide level -------------------------------------
  # FC1_FC2_reads.invivo.guidelevel <- FC1_FC2_reads.invivo.30.nopolyUMIs
  # FC1_FC2_reads.invivo.guidelevel$reads.a <- FC1_FC2_reads.invivo.30.nopolyUMIs$active.FC1 + FC1_FC2_reads.invivo.30.nopolyUMIs$active.FC2
  # FC1_FC2_reads.invivo.guidelevel$reads.i <- FC1_FC2_reads.invivo.30.nopolyUMIs$inactive.FC1 + FC1_FC2_reads.invivo.30.nopolyUMIs$inactive.FC2
  # #no UMIs / guide level dataset
  # FC1_FC2_reads.invivo.guides<- FC1_FC2_reads.invivo.guidelevel %>% 
  #   group_by(Replicate, guide, Gene, `in_vitro/in_vivo`) %>% 
  #   summarise(reads.i=sum(reads.i), reads.a=sum(reads.a))
  # FC1_FC2_reads.invivo.guides.pseudo <- FC1_FC2_reads.invivo.guides
  # FC1_FC2_reads.invivo.guides.pseudo$reads.a <- FC1_FC2_reads.invivo.guides.pseudo$reads.a + 0.5
  # FC1_FC2_reads.invivo.guides.pseudo$reads.i <- FC1_FC2_reads.invivo.guides.pseudo$reads.i + 0.5
  # 
  # FC1_FC2_reads.invivo.guides.pseudo_sep <- FC1_FC2_reads.invivo.guides.pseudo %>% 
  #   group_by(guide, Gene, reads.a, reads.i) %>% 
  #   summarise(reads.a = sum(reads.a), reads.i = sum(reads.i))
  # 
  # setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Mageck_bothFC/')
  # write.table(FC1_FC2_reads.invivo.guides.pseudo_sep, 'FC1_FC2.guide.invivo.pseudo.sep.tsv', row.names = F, quote = F, sep ='\t')

#BOTH#--------------------------------- WRITE TABLES  of pseudo counts----------------------------------------------------------
setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/')
write.table(a.i.guides.only.invitro.pseudo.noRepl, "a.i.guides.only.invitro.pseudo.noRepl.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(FC1_FC2_filtered_UMI.pseudo, "FC1_FC2_filtered_UMI.pseudo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#IN VITRO#------------------------------------ PRE-MAGECK ----------------------------------------------------------------
setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Mageck_bothFC/')
write.table(a.i.guides.only.invitro.filtered.pseudo_wide, 'a.i.guides.only.invitro.filtered.pseudo_wide.tsv', row.names = F, quote = F, sep ='\t')

#IN VIVO#------------------------------------ PRE-MAGECK ----------------------------------------------------------------
FC1_FC2.guideUMI.invivo.pseudo.sep <- FC1_FC2_filtered_UMI.pseudo %>% 
  dplyr::select(guide_UMI, Gene, reads.a, reads.i)

setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Mageck_bothFC/')
write.table(FC1_FC2.guideUMI.invivo.pseudo.sep, 'FC1_FC2.guideUMI.invivo.pseudo.sep.tsv', row.names = F, quote = F, sep ='\t')

#IN VITRO#------------------------------------ MAGECK ----------------------------------------------------------------
setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Mageck_bothFC/')
mageck_invitro_genes_pseudo.standard <- read.table('mageck.invitro.filteredsum10.pseudo.wide.standard.gene_summary.txt', header = T)
mageck_invitro_guides_pseudo.standard <- read.table('mageck.invitro.filteredsum10.pseudo.wide.standard.sgrna_summary.txt', header = T)
mageck_invitro_genes_pseudo.standard$minRRA.fdr <- -log10(apply(mageck_invitro_genes_pseudo.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_invitro_genes_pseudo.standard$minRRA.score <- -log10(apply(mageck_invitro_genes_pseudo.standard[,c("neg.score","pos.score")], 1, min))
mageck_invitro_genes_pseudo.standard$minpvalue <- -log10(apply(mageck_invitro_genes_pseudo.standard[,c("neg.p.value","pos.p.value")], 1, min))

#IN VIVO#------------------------------------ MAGECK ----------------------------------------------------------------
setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Mageck_bothFC/')
mageck_invivo_genes_pseudo.standard <- read.table('mageck.invivo.bothFC.filtered.pseudo.standard.gene_summary.txt', header = T)
mageck_invivo_guides_pseudo.standard <- read.table('mageck.invivo.bothFC.filtered.pseudo.standard.sgrna_summary.txt', header = T)
mageck_invivo_genes_pseudo.standard$minRRA.fdr <- -log10(apply(mageck_invivo_genes_pseudo.standard[,c("neg.fdr","pos.fdr")], 1, min))
mageck_invivo_genes_pseudo.standard$minRRA.score <- -log10(apply(mageck_invivo_genes_pseudo.standard[,c("neg.score","pos.score")], 1, min))
mageck_invivo_genes_pseudo.standard$minpvalue <- -log10(apply(mageck_invivo_genes_pseudo.standard[,c("neg.p.value","pos.p.value")], 1, min))

#IN VITRO#------------------------------------ Number of guides or UMIS per gene ----------------------------------------------------------------
mageck_invitro_genes_pseudo.standard %>%
  ggplot(aes(x = num)) +
  geom_histogram(binwidth = 1, color = "black", fill= "grey") +
  theme_bw()+
  xlab('Number of guide_UMIs per gene') + ylab('Count')

nrGene_invitro <- nrow(mageck_invitro_genes_pseudo.standard)

#IN VIVO#------------------------------------ Number of guides or UMIS per gene ----------------------------------------------------------------
mageck_invivo_genes_pseudo.standard %>%
  ggplot(aes(x = num)) +
  geom_histogram(binwidth = 1, color = "black", fill= "grey") +
  theme_bw()+
  xlim(0,25) +
  xlab('Number of UMIs per gene') + ylab('Count') 
 

#IN VIVO#------------------------------------ Number of UNIQUE UMIS  / TAKE RATE----------------------------------------------------------------
nr.replicates <- count(unique(FC1_FC2_filtered_UMI.pseudo[, "Replicate"]))
unique.UMIs <- nrow(mageck_invivo_guides_pseudo.standard)
Take_rate_perReplicate <- unique.UMIs/nr.replicates
Take_rate_perMouse <- unique.UMIs/55 #nr.mice in complete screen

nrGene_invivo <- nrow(mageck_invivo_genes_pseudo.standard)

#IN VITRO#------------------------------------ minRRA.score plot in VITRO ----------------------------------------------------------------
colors <- c('other' = 'gray', 'essential' = 'red', 'control' = 'black')

mageck_invitro_genes_pseudo.standard %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'essential'), data=mageck_invitro_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'control'), data=mageck_invitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('LFC') + ylab('RRA score')+ labs(color = ' ', title = 'In vitro gene LFC + RRA score') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

#IN VITRO#------------------------------------ minRRA.fdr plot in VITRO ----------------------------------------------------------------
mageck_invitro_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc, y = minRRA.fdr)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'control'), data=mageck_invitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'essential'), data=mageck_invitro_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5, size=2.5) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('LFC') + ylab('RRA fdr')+ labs(color = ' ', title = 'In vitro gene LFC + RRA fdr') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

#IN VITRO#------------------------------------ p Value plot in VITRO ----------------------------------------------------------------
mageck_invitro_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc, y = minpvalue)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'control'), data=mageck_invitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'essential'), data=mageck_invitro_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5, size=2.5) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('LFC') + ylab('p Value')+ labs(color = ' ', title = 'In vitro gene LFC + p value') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

#IN VIVO#------------------------------------ minRRA.score plot in VIVO ----------------------------------------------------------------
colors <- c('other' = 'gray', 'essential' = 'red', 'control' = 'black')

mageck_invivo_genes_pseudo.standard %>%
  ggplot(aes(x = neg.lfc, y = minRRA.score)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'essential'), data=mageck_invivo_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'control'), data=mageck_invivo_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  # coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('LFC') + ylab('RRA score')+ labs(color = ' ', title = 'In vivo gene LFC + RRA score') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

#IN VIVO#------------------------------------ minRRA.fdr plot in VIVO ----------------------------------------------------------------
mageck_invivo_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc, y = minRRA.fdr)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'control'), data=mageck_invivo_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'essential'), data=mageck_invivo_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5, size=2.5) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('LFC') + ylab('RRA fdr')+ labs(color = ' ', title = 'In vivo gene LFC + RRA fdr') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

#IN VIVO#------------------------------------ p Value plot in VIVO ----------------------------------------------------------------
mageck_invivo_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc, y = minpvalue)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'control'), data=mageck_invivo_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'essential'), data=mageck_invivo_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5, size=2.5) +
  #coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('LFC') + ylab('p Value')+ labs(color = ' ', title = 'In vivo gene LFC + p value') +
  #xlim(-c,c) + 
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

## ---------- correlation plots  ---------------------------------------------------------------------------
mageck_invivovitro_genes_pseudo.standard <- merge(mageck_invitro_genes_pseudo.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num")], 
                                                  mageck_invivo_genes_pseudo.standard [,c("id","neg.lfc","minRRA.score", "minRRA.fdr", "minpvalue", "num")],
                                                  by = 'id', suffixes = c(".invitro",".invivo"))
# write.table(mageck_invivovitro_genes_pseudo.standard,"/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Mageck_bothFC/mageck_invivovitro_genes_pseudo.standard.txt",
#             sep = '\t', row.names = F, quote = F)
mageck_invivovitro_genes_pseudo.standard <- read.table('mageck_invivovitro_genes_pseudo.standard,"/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Mageck_bothFC/mageck_invivovitro_genes_pseudo.standard.txt', header = T)

colors <- c('other' = 'gray', 'essential' = 'red', 'control' = 'black')

# this is to calculate the maximum value that we need to plot, I then use this value to set the limits of the x & y axes so that the plot looks nicer
g<-max(abs(mageck_invivovitro_genes_pseudo.standard$neg.lfc.invivo))
h<-max(abs(mageck_invivovitro_genes_pseudo.standard$neg.lfc.invitro))
i<-round(max(g,h),1) + 0.5

mageck_invivovitro_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc.invitro, y = neg.lfc.invivo)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'essential'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'control'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5, size=2.5) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('In vitro LFC') + ylab('In vivo LFC')+ labs(color = ' ',
                                                   title = 'In vivo versus in vitro LFC correlation') +
  #xlim(-i,5) + ylim(-i,5) +
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

## ---------- correlation plots  -------------------------------- + BUBBLE SIZE ------------------------------------
mageck_invivovitro_genes_pseudo.standard$bubble_size <- ifelse(
  mageck_invivovitro_genes_pseudo.standard$num.invivo > 15, 15, 
  mageck_invivovitro_genes_pseudo.standard$num.invivo)
#write.table(mageck_invivovitro_genes_pseudo.standard, "/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/mageck_invivovitro_genes_pseudo.standard.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#more complete version is more below 
mageck_invivovitro_genes_pseudo.standard <- read.table("/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Mageck_bothFC/mageck_invivovitro_genes_pseudo.standard.txt", header = T)

colors <- c('other' = 'gray', 'essential' = 'red', 'control' = 'black')

mageck_invivovitro_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc.invitro, y = neg.lfc.invivo, size = bubble_size)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5) +
  geom_point(aes(color = 'essential'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5) +
  geom_point(aes(color = 'control'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('In vitro LFC') + ylab('In vivo LFC')+ labs(color = ' ',
                                                   title = 'In vivo versus in vitro LFC correlation', size = 'Number of UMIs per \ngene in vivo') +
  xlim(-10,2.5) + ylim(-10,5.5) +
  theme_bw() + scale_color_manual(values = colors) + scale_size_continuous(range = c(0.5,7)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

## ---------- MA plots  --------------------------------------------------------------------
mageck_invitro_guides_pseudo.standard$sumI.A <- mageck_invitro_guides_pseudo.standard$treatment_count + mageck_invitro_guides_pseudo.standard$control_count
ncol(mageck_invitro_guides_pseudo.standard)

mageck_invitro_guides_pseudo.standard %>%
  ggplot( aes(x = log10(control_count), y = LFC)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5, size = 2.5) +
  geom_point(aes(color = 'essential'), data=mageck_invitro_guides_pseudo.standard %>% filter(Gene %in% Essentials$gene), alpha=0.5, size=2.5) +
  geom_point(aes(color = 'control'), data=mageck_invitro_guides_pseudo.standard %>% filter(toupper(Gene) %in% Controls$gene), alpha=0.5, size=2.5) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('log10(sumI.A)') + ylab('In vitro LFC per guide')+ labs(color = ' ',
                                                               title = 'LFC in vitro vs sum reads') +
  #xlim(0,1000000) + ylim(-5,10) +
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

## ---------- correlation plots  -------------------------------- HIGHLIGHTING DIFFERENT GENES ------------------------------------

mageck_invivovitro_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc.invitro, y = neg.lfc.invivo,)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5) +
  geom_point(aes(color = 'essential'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5) +
  geom_point(aes(color = 'control'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5) +
  #geom_point(aes(color = ''), data=mageck_invivovitro_genes_pseudo.standard %>% filter(id == 'Aldoa'), alpha=0.5, color = 'blue') +
  #geom_point(aes(color = 'Nduf'), data=mageck_invitro_guides_pseudo.standard %>% filter(Gene %in% Oxphos$gene), alpha=0.5) +
  #geom_text_repel(aes(label=ifelse(id %in% ATPase$gene, as.character(id),'')),box.padding = 2, max.overlaps = Inf, 
  #                segment.color = "grey50", segment.size = 0.1, size = 3) +
  geom_text_repel(aes(label=ifelse(neg.lfc.invivo < -5 & neg.lfc.invitro > -2, as.character(id),'')),box.padding = 2, max.overlaps = Inf, 
                  segment.color = "grey50", segment.size = 0.1, size = 3) +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('In vitro LFC') + ylab('In vivo LFC')+ labs(color = ' ', title = 'In vivo versus in vitro LFC correlation \npseudo.standard', 
                                                   size = 'Guide_UMIs  \nper gene in vivo') +
  xlim(-i,5.5) + ylim(-i,5.5) +
  theme_bw() + scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

## ---------- line plots for guides per gene  --------------------------------------------------------------------

mageck_invitro_guides_pseudo.standard$guide_number <- str_split_fixed(mageck_invitro_guides_pseudo.standard$sgrna, '_',3)[,2] 
mageck_invivo_guides_pseudo.standard$guide_number <- str_split_fixed(mageck_invivo_guides_pseudo.standard$sgrna, '_',3)[,2]

mageck_invitro_guides_pseudo.standard %>% filter(Gene %in% 
    c('Aldoa', 'Sod1',  'Lamtor4', 'Ndufb5', 'Atp5a1', 'Atp5k', 'Pigr', 'Pkn2', 'Birc6', 'Defa2','Rpl9-ps1','Tpx2','Prkar1a', 'Rad17', 'Csnk2b')) %>% 
  ggplot(aes(x = LFC, y = Gene, color = guide_number)) + 
  geom_point(size = 2) +
  geom_vline(xintercept= 0, linetype = 'dashed') +
  #scale_color_manual(values = c('darkgreen', 'gray60', 'red', 'cornflowerblue', 'gold2')) +
  theme_bw() +
  theme(legend.position = 'none')

mageck_invitro_guides_pseudo.standard %>% filter(Gene %in% invitro.depletion$id) %>% 
  ggplot(aes(x = LFC, y = Gene, color = guide_number)) + 
  geom_point(size = 1) +
  geom_vline(xintercept= 0, linetype = 'dashed') +
  theme_bw() 

## ----------  curved lines hit retrieval ------------ with lm line--------------------------------------------------------
# lm_output <- summary(lm(mageck_invivovitro_genes_pseudo.standard$neg.lfc.invivo ~ mageck_invivovitro_genes_pseudo.standard$neg.lfc.invitro))
# 
# sd.vivoContr <- sd(mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene) %>% pull(neg.lfc.invivo) ) * 1.5 #sd.vivoContr
# sd.vitroContr <- sd(mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene) %>% pull(neg.lfc.invitro) ) * 2
# 
# intercept = sd.vivoContr
# slope= lm_output$coefficients[2,1]
# 
# data_p=mageck_invivovitro_genes_pseudo.standard %>% filter(id == 'Aldoa')
# data_p$model_value <-  slope * data_p$neg.lfc.invitro + intercept
# data_p$above <- data_p$model_value < data_p$neg.lfc.invivo
# 
# mageck_invivovitro_genes_pseudo.standard <-  mageck_invivovitro_genes_pseudo.standard %>% 
#   mutate(model_high = slope * neg.lfc.invitro + 1.160818,
#          model_low = slope * neg.lfc.invitro - 1.160818, 
#          type = 'other')
# mageck_invivovitro_genes_pseudo.standard <- mageck_invivovitro_genes_pseudo.standard %>% 
#   mutate(type = case_when((model_high < neg.lfc.invivo & neg.lfc.invitro < -0.4160973) | (neg.lfc.invitro < -1 & neg.lfc.invivo > -0.5 )  ~ 'in vitro specific',
#                           model_low > neg.lfc.invivo & neg.lfc.invitro < -0.4160973 ~ 'in vivo lowish', 
#                           model_high < neg.lfc.invivo & between(neg.lfc.invitro, -1,2) ~ 'in vivo specific up',
#                           model_low > neg.lfc.invivo & between(neg.lfc.invitro, -1, 2) ~ 'in vivo specific down', 
#                           TRUE ~ 'other'))

mageck_invivovitro_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc.invitro, y = neg.lfc.invivo,)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5) +
  geom_point(aes(color = type), alpha=0.5) +
  geom_abline(slope = lm_output$coefficients[2,1], intercept = 1.160818, linetype = 'dashed') +
  geom_abline(slope = lm_output$coefficients[2,1], intercept = -1.160818, linetype = 'dashed') +
  geom_vline(xintercept = 0.4160973, linetype = 'dashed') +
  geom_vline(xintercept = -0.4160973, linetype = 'dashed') +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('In vitro LFC') + ylab('In vivo LFC')+ labs(color = ' ', title = 'In vivo versus in vitro LFC correlation \npseudo.standard', 
                                                   size = 'Guide_UMIs  \nper gene in vivo') +
  xlim(-i,5.5) + ylim(-i,5.5) +
  theme_bw() #+ scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

## ----------  hit retrieval ------------ with diagonal line through 0,0 --------------------------------------------------------
slope = 1
sd.vivoContr <- sd(mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene) %>% pull(neg.lfc.invivo) ) * 1.5 #sd.vivoContr
sd.vitroContr <- sd(mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene) %>% pull(neg.lfc.invitro) ) * 2

mageck_invivovitro_genes_pseudo.standard <-  mageck_invivovitro_genes_pseudo.standard %>% 
  mutate(model_high = slope * neg.lfc.invitro + sd.vivoContr,
         model_low = slope * neg.lfc.invitro - sd.vivoContr, 
         type = 'other')
mageck_invivovitro_genes_pseudo.standard <- mageck_invivovitro_genes_pseudo.standard %>% 
  mutate(type = case_when(neg.lfc.invivo > 2.5  ~ 'in vivo specific up', 
                          (model_high < neg.lfc.invivo & neg.lfc.invitro < -sd.vitroContr) 
                          | (neg.lfc.invitro < -1 & neg.lfc.invivo > -0.5 )  ~ 'in vitro specific',
                          model_low > neg.lfc.invivo & neg.lfc.invitro < -sd.vitroContr ~ 'in vivo lowish', 
                          model_low > neg.lfc.invivo & between(neg.lfc.invitro, -1, 2) ~ 'in vivo specific down', 
                          TRUE ~ 'other'))

mageck_invivovitro_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc.invitro, y = neg.lfc.invivo,)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.5) +
  geom_point(aes(color = type), alpha=0.5) +
  geom_point(aes(color = 'essential'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.5) +
  geom_point(aes(color = 'control'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.5) +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = sd.vivoContr, linetype = 'dashed') +
  geom_abline(slope = 1, intercept = -sd.vivoContr, linetype = 'dashed') +
  geom_vline(xintercept = sd.vitroContr, linetype = 'dashed') +
  geom_vline(xintercept = -sd.vitroContr, linetype = 'dashed') +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('In vitro LFC') + ylab('In vivo LFC')+ labs(color = ' ', title = 'In vivo versus in vitro LFC correlation \npseudo.standard mageck', 
                                                   size = 'Guide_UMIs  \nper gene in vivo') +
  xlim(-i,5.5) + ylim(-i,5.5) +
  theme_bw() + scale_color_brewer(palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom') +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

## --------- MITO list ----------------
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("biomaRt")
library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mouse ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('external_gene_name','ensembl_gene_id', 'go_id', 'name_1006', 'description'),
                   filters = 'go', values = 'GO:0005739', mart = ensembl)
listFilters(mart=ensembl)
mito.genes <- gene.data[gene.data$go_id == 'GO:0005739',]
#write.table(mito.genes, "/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/GO analysis/GO term list. matched Apools/mito.genes.GO.0005739.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

detach("package:biomaRt", unload = TRUE)

mito.genes <- read.table("/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/GO analysis/GO term list. matched Apools/mito.genes.GO.0005739.txt", header = TRUE)

## --------- Validation list Naomi ----------------
setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis')
validation.list.Naomi <- read.table("Validation guide list Naomi.txt", header = TRUE)

validation.list.Naomi$gene <- str_split_fixed(validation.list.Naomi$guide, '_', 2)[,1]
validation.list.Naomi.vivo <- validation.list.Naomi[validation.list.Naomi$group == 'vivo',]
validation.list.Naomi.both <- validation.list.Naomi[validation.list.Naomi$group == 'both',]

## ----------  tables hit retrieval --------------------------------------------------------------------
table(mageck_invivovitro_genes_pseudo.standard$type)
invivo.depletion <-  mageck_invivovitro_genes_pseudo.standard %>% filter(type == 'in vivo specific down')
invivo.enriched <-  mageck_invivovitro_genes_pseudo.standard %>% filter(type == 'in vivo specific up')
invivo.low <-  mageck_invivovitro_genes_pseudo.standard %>% filter(type == 'in vivo lowish')  
invitro.depletion <- mageck_invivovitro_genes_pseudo.standard %>% filter(type == 'in vitro specific')  

mageck_invivovitro_genes_pseudo.standard$Gene.group <- "Other"
mageck_invivovitro_genes_pseudo.standard[mageck_invivovitro_genes_pseudo.standard$id %in% Essentials$gene, "Gene.group"] <- "Essentials"
mageck_invivovitro_genes_pseudo.standard[toupper(mageck_invivovitro_genes_pseudo.standard$id) %in% Controls$gene, "Gene.group"] <- "Controls"
mageck_invivovitro_genes_pseudo.standard[mageck_invivovitro_genes_pseudo.standard$id %in% mito.genes$external_gene_name, "Gene.group"] <- "Mito"

#write.table(mageck_invivovitro_genes_pseudo.standard, "/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/mageck_invivovitro_genes_pseudo.standard.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Complete list with gene groups
mageck_invivovitro_genes_pseudo.standard <- read.table( "/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists/Yumm_Apools_July2021_Analysis/02_Data_lists/mageck_invivovitro_genes_pseudo.standard.txt", header = TRUE)

## ----------  pre-GO term --------------------------------------------------------------------
setwd("/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Scripts and gene lists")
pool.A1 <- read.table('A1.mm.FINAL.VBC_sgRNA.guide_list (per pool).txt', header = T, sep = "\t", quote = "") # I have saved the excel sheets as tab-delimited text files
pool.A2 <- read.table('A2.mm.FINAL.VBC_sgRNA.guide_list (per pool).txt', header = T, sep = "\t", quote = "")
pool.A3 <- read.table('A3.mm.FINAL.VBC_sgRNA.guide_list (per pool).txt', header = T, sep = "\t", quote = "")
pool.A4 <- read.table('A4.mm.FINAL.VBC_sgRNA.guide_list (per pool).txt', header = T, sep = "\t", quote = "")
pool.A1$cloning.FW_oligo
all.pools <- rbind(pool.A1[,c('ensembl_gene_id', 'gene', 'sgRNA', 'Auto.pick.top.sgRNAs', 'cloning.FW_oligo')], pool.A2[,c('ensembl_gene_id', 'gene', 'sgRNA', 'Auto.pick.top.sgRNAs', 'cloning.FW_oligo')], 
                   pool.A3[,c('ensembl_gene_id', 'gene', 'sgRNA', 'Auto.pick.top.sgRNAs', 'cloning.FW_oligo')], pool.A4[,c('ensembl_gene_id', 'gene', 'sgRNA', 'Auto.pick.top.sgRNAs', 'cloning.FW_oligo')])  # this has all gene names & all ensembl IDs for pool A
# remove all the duplicated rows (because 5 gRNAs per gene in the original file so 5 rows per gene)
all.pools <- distinct(all.pools)
mageck_invivovitro_genes_pseudo.standard <- merge(mageck_invivovitro_genes_pseudo.standard, all.pools, by.x = 'id', by.y = 'gene', all.x = TRUE)
names(mageck_invivovitro_genes_pseudo.standard)[names(mageck_invivovitro_genes_pseudo.standard) == 'ensembl_gene_id'] <- 'ensembl_id'
#mageck_invivovitro_genes_pseudo.standard <- subset(mageck_invivovitro_genes_pseudo.standard, select = -c(12,16,18))

all.pools <- read.table("/Volumes/elling/Esther/NGS_results/Scripts_gene_lists/Mouse_library/all.Apools.sgRNAs.txt", header = TRUE)

## ----------  LOAD GO term lists --------------------------------------------------------------------

######## Go to document Gene lists generation (from GO term lists).R or Script_Gene_list_libraries.R

## ----------  ADD column for GO term lists --------------------------------------------------------------------

vivovitro.LFC.genegroups <- mageck_invivovitro_genes_pseudo.standard %>% 
  dplyr::select(id, neg.lfc.invitro, minRRA.score.invitro, num.invitro, neg.lfc.invivo, minRRA.fdr.invivo, num.invivo, type, ensembl_id)
vivovitro.LFC.genegroups$Controls <- ifelse((toupper(vivovitro.LFC.genegroups$id) %in% Controls$gene), TRUE, '-')
vivovitro.LFC.genegroups$Oxphos <- ifelse(vivovitro.LFC.genegroups$id %in% Oxphos$gene, TRUE, '-')
vivovitro.LFC.genegroups$Ars2 <- ifelse(vivovitro.LFC.genegroups$id %in% Ars2$gene, TRUE, '-')
vivovitro.LFC.genegroups$mito <- ifelse(vivovitro.LFC.genegroups$id %in% mito.genes$external_gene_name, TRUE, '-')
vivovitro.LFC.genegroups$Nduf <- ifelse(vivovitro.LFC.genegroups$id %in% Nduf$gene, TRUE, '-')
vivovitro.LFC.genegroups$Essentials <- ifelse(vivovitro.LFC.genegroups$id %in% Essentials$gene, TRUE, '-')

## ---------- hit retrieval and lists ------------------------------------------------------------------
names(vivovitro.LFC.genegroups) [9] <- 'ensembl_id'

Allgenes.Apools.inscreen <- vivovitro.LFC.genegroups %>% dplyr::select(ensembl_id)
invivo.depletion<- vivovitro.LFC.genegroups %>% filter(id %in% invivo.depletion$id)
invivo.enriched<- vivovitro.LFC.genegroups %>% filter(id %in% invivo.enriched$id)
invivo.low<- vivovitro.LFC.genegroups %>% filter(id %in% invivo.low$id)
invitro.depletion<- vivovitro.LFC.genegroups %>% filter(id %in% invitro.depletion$id)
invivo.alldepleting.filtered <- rbind(invivo.low,invivo.depletion) %>% filter((neg.lfc.invivo < -2.5) & (neg.lfc.invitro > -2.5))

setwd("/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/GO analysis/")
# write.table(invivo.depletion, "hits.invivo.depletion.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(invivo.enriched, "hits.invivo.enriched.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(invivo.low, "hits.invivo.low.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(invitro.depletion, "hits.invitro.depletion.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(vivovitro.LFC.genegroups, "vivovitro.LFC.genegroups.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(Allgenes.Apools.inscreen, "Allgenes.Apools.inscreen.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(invivo.alldepleting.filtered, "invivo.alldepleting.filtered.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
## -------------- GO term analysis EXAMPLE IN VIVO All DEPLETING ------- combined - Biological Process and Cellular Component Go term analysis ------
invivo.alldepleting.filtered.id <- invivo.alldepleting.filtered$id
invivo.alldepleting.filtered.id <- dplyr::data_frame(invivo.alldepleting.filtered.id) #or invivo.alldepleting.filtered.id <- tibble(invivo.alldepleting.filtered.id)
names(invivo.alldepleting.filtered.id) <- 'id'
setdiff(invivo.alldepleting.filtered.id$id, all.pools$gene)

invivo.alldepleting.filtered.id.ensembl <- all.pools %>% filter(gene %in% invivo.alldepleting.filtered.id$id)
setwd("/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/GO analysis/")
#write.table(invivo.alldepleting.filtered.id.ensembl$ensembl_gene_id, "invivo.alldepleting.filtered.ensembl.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## ----------  HIT retrieval ------------ with diagonal line through 0,0 --------------------------------------------------------

slope = 1
sd.vivoContr <- sd(mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene) %>% pull(neg.lfc.invivo) ) * 1.5 #sd.vivoContr
sd.vitroContr <- sd(mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene) %>% pull(neg.lfc.invitro) ) * 2

mageck_invivovitro_genes_pseudo.standard <-  mageck_invivovitro_genes_pseudo.standard %>% 
  mutate(model_high = slope * neg.lfc.invitro + sd.vivoContr,
         model_low = slope * neg.lfc.invitro - sd.vivoContr, 
         type = 'other')
mageck_invivovitro_genes_pseudo.standard <- mageck_invivovitro_genes_pseudo.standard %>% 
  mutate(type = case_when(neg.lfc.invivo > 2.5  ~ 'in vivo specific up', 
                          (model_high < neg.lfc.invivo & neg.lfc.invitro < -sd.vitroContr) 
                          | (neg.lfc.invitro < -1 & neg.lfc.invivo > -0.5 )  ~ 'in vitro specific',
                          model_low > neg.lfc.invivo & neg.lfc.invitro < -sd.vitroContr ~ 'in vivo lowish', 
                          model_low > neg.lfc.invivo & between(neg.lfc.invitro, -1, 2) ~ 'in vivo specific down', 
                          TRUE ~ 'other'))

invivo.alldepleting.filtered.2 <- rbind(invivo.low,invivo.depletion) %>% filter((neg.lfc.invivo < -2) & (neg.lfc.invitro > -1.5))
invivo.alldepleting.filtered.noEssentials.2 <- invivo.alldepleting.filtered.2 %>% filter(!Essentials == TRUE)

colors3 <- c('other' = 'gray', 'essential' = 'red', 'control' = 'black', 'invivo.depleting' = 'limegreen')

mageck_invivovitro_genes_pseudo.standard %>%
  ggplot( aes(x = neg.lfc.invitro, y = neg.lfc.invivo,)) +
  geom_hline(aes(yintercept=0), color="gray90", size = 2) +
  geom_vline(aes(xintercept=0), color="gray90", size = 2) +  # add thick gray lines on the x & y axis
  geom_point(aes(color = 'other'), alpha=0.2) +
  # geom_point(aes(color = ''), data=mageck_invivovitro_genes_pseudo.standard %>% 
  #              filter(model_high < neg.lfc.invivo, neg.lfc.invitro < -1), alpha=0.5, color = 'blue') +
  # geom_point(aes(color = ''), data=mageck_invivovitro_genes_pseudo.standard %>% 
  # filter(model_low > neg.lfc.invivo, neg.lfc.invitro < -1), alpha=0.5, color = 'red') +
  #geom_point(aes(color = type), alpha=0.5) +
  geom_point(aes(color = 'essential'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% Essentials$gene), alpha=0.2) +
  geom_point(aes(color = 'control'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene), alpha=0.2) +
  geom_point(aes(color = 'invivo.depleting'), data=mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% invivo.alldepleting.filtered.2$id), alpha=0.2) +
  #geom_point(aes(color = ''), data=mageck_invivovitro_genes_pseudo.standard %>% filter(Gene %in% Nduf$gene), alpha=0.5) +
  #geom_point(aes(color = ''), data=mageck_invivovitro_genes_pseudo.standard %>% filter(id == 'Hars2'), alpha=0.5, color = 'blue') +
  # geom_text_repel(aes(label=ifelse((neg.lfc.invivo > 2.4 & neg.lfc.invivo < 3), as.character(id),'')),box.padding = 2, max.overlaps = Inf,
  #                 segment.color = "grey50", segment.size = 0.1, size = 3) +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = sd.vivoContr, linetype = 'dashed') +
  geom_abline(slope = 1, intercept = -sd.vivoContr, linetype = 'dashed') +
  geom_vline(xintercept = sd.vitroContr, linetype = 'dashed') +
  geom_vline(xintercept = -sd.vitroContr, linetype = 'dashed') +
  #geom_smooth(stat = 'lm') +
  coord_fixed(ratio = 1) +  # this makes sure the x & y axes have the same scale
  xlab('In vitro LFC') + ylab('In vivo LFC')+ labs(color = ' ', title = 'In vivo versus in vitro LFC correlation \npseudo.standard mageck', 
                                                   size = 'Guide_UMIs  \nper gene in vivo') +
  xlim(-10,5.1) + ylim(-10,5.1) +
  theme_bw() + scale_color_manual(values = colors3) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom', text = element_text(size=13) ) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.25)))  # this sets the opacity of the legend dots

## ----------  DESIGN VALIDATION LIST --------------------------------------------------------------------
#invivo.alldepleting.filtered.2 -- (invivo.low,invivo.depletion) %>% filter((neg.lfc.invivo < -2) & (neg.lfc.invitro > -1.5)) -- 207 genes
#Nduf  -- 42 genes
#Merge_aerobic_resp_transp.screen  -- 88 genes
#validation.list.Naomi.screen  -- 21 genes
#Controls  -- 98 genes after filtering -- 495 genes after filtering
#invivo.enriched.screen -- 97 genes
invivo.alldepleting.filtered.2 <- rbind(invivo.low,invivo.depletion) %>% filter((neg.lfc.invivo < -2) & (neg.lfc.invitro > -1.5))
invivo.alldepleting.filtered.noEssentials.2 <- invivo.alldepleting.filtered.2 %>% filter(!Essentials == TRUE)
Merge_aerobic_resp_transp <- read.table("/Volumes/elling/Esther/NGS_results/Scripts_gene_lists/GOterm_Gene_lists/Merge_aerobic_resp_transp.txt", header =TRUE, sep = "\t")

invivo.alldepleting.filtered.2.screen <- mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% invivo.alldepleting.filtered.2$id)
Nduf.screen <- mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% Nduf$gene)
Merge_aerobic_resp_transp.screen <- mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% Merge_aerobic_resp_transp$Gene)
validation.list.Naomi.screen <- mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% validation.list.Naomi$gene)
Controls.screen <- mageck_invivovitro_genes_pseudo.standard %>% filter(toupper(id) %in% Controls$gene)  %>% 
  filter((between(neg.lfc.invivo, -1,1)) & (between(neg.lfc.invitro, -0.25,0.25))) %>% filter(num.invivo > 8) #98 genes
invivo.enriched.screen <- mageck_invivovitro_genes_pseudo.standard %>% filter(id %in% invivo.enriched$id)

Validation.list.Esther <- rbind(invivo.alldepleting.filtered.2.screen, Nduf.screen, Merge_aerobic_resp_transp.screen, validation.list.Naomi.screen,Controls.screen, invivo.enriched.screen)
Validation.list.Esther.noDupl <- Validation.list.Esther[!duplicated(Validation.list.Esther),]
#write.table(Validation.list.Esther.noDupl, "Validation.list.Esther.noDupl.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
Validations.Esther.woEnriched <- rbind(invivo.alldepleting.filtered.2.screen, Nduf.screen, Merge_aerobic_resp_transp.screen, validation.list.Naomi.screen)
Validations.Esther.woEnriched.noDupl <- Validations.Esther.woEnriched[!duplicated(Validations.Esther.woEnriched),]

#n_distinct(Validations.Esther.Controls.guides$id)

Validation.list.Esther.noDupl.genesOnly <- Validation.list.Esther.noDupl %>% dplyr::select(id,ensembl_id)
Validations.Esther.woEnriched.genesOnly <- Validations.Esther.woEnriched.noDupl %>% dplyr::select(id,ensembl_id)
Validations.Esther.vivoEnriched <- invivo.enriched.screen[!duplicated(invivo.enriched.screen),]
Validations.Esther.vivoEnriched.genesOnly <- Validations.Esther.vivoEnriched %>% dplyr::select(id,ensembl_id)
Validations.Esther.Controls <- Controls.screen[!duplicated(Controls.screen),]
Validations.Esther.Controls.genesOnly <- Validations.Esther.Controls %>% dplyr::select(id,ensembl_id)

#write.table(Validation.list.Esther.noDupl.genesOnly, "Validation.list.Esther.noDupl.genesOnly.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
Validation.list.Esther.noDupl.genesOnly <- read.table("/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/Validation.list.Esther.noDupl.genesOnly.txt", header = FALSE)
## ----------  DESIGN VALIDATION LIST -------- 3 best performing sgRNAs------------------------------------------------------------

#invivo.guide.sum <- mageck_invivo_guides_pseudo.standard %>% group_by(Gene, guide_number) %>% summarise(meanLFC=mean(LFC))
#invivo.guide.sum <- invivo.guide.sum[order( invivo.guide.sum[,1], invivo.guide.sum[,3] ),]

#setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis')
#invivo.guide.best.performing <- invivo.guide.sum %>% group_by(Gene) %>%
# slice(1:3)
#write.table(invivo.guide.best.performing, "invivo.guide.best.performing.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#check <- invivo.guide.best.performing %>% group_by(Gene) %>% count()

all.pools.extraINFO <- all.pools
names(all.pools.extraINFO)[names(all.pools.extraINFO) == 'ensembl_gene_id'] <- 'ensembl_id'

Validation.list.Esther.noDupl.oldGuides <- merge(Validation.list.Esther.noDupl.genesOnly, all.pools.extraINFO, by='ensembl_id')
names(Validation.list.Esther.noDupl.oldGuides)[names(Validation.list.Esther.noDupl.oldGuides) == 'gene'] <- 'Gene'
names(Validation.list.Esther.noDupl.oldGuides)[names(Validation.list.Esther.noDupl.oldGuides) == 'Auto.pick.top.sgRNAs'] <- 'guide_number'
#Validation.list.Esther.noDupl.bestGuides <- merge(Validation.list.Esther.noDupl.oldGuides,invivo.guide.best.performing, by= c('Gene', 'guide_number'))

#check <- Validation.list.Esther.noDupl.bestGuides %>% group_by(Gene) %>% count()
#genes1guide<- check %>% filter(n == 1)
#GuidesToAdd1<- Validation.list.Esther.noDupl.NotUsedguides %>% filter(Gene %in% genes1guide$Gene) %>% group_by(Gene) %>%slice(1:2)
#genes2guides<- check %>% filter(n == 2)
#GuidesToAdd2<- Validation.list.Esther.noDupl.NotUsedguides %>% filter(Gene %in% genes2guides$Gene) %>% group_by(Gene) %>%slice(1:1)

#Validation.list.Esther.noDupl.3Guides <- bind_rows(Validation.list.Esther.noDupl.bestGuides, GuidesToAdd1, GuidesToAdd2)
#Validation.list.Esther.noDupl.3Guides <- subset(Validation.list.Esther.noDupl.3Guides, select = -c(7))
#Validation.list.Esther.noDupl.3Guides <- Validation.list.Esther.noDupl.3Guides %>% filter(!ensembl_id == 'ENSMUSG00000116933')
#Validation.list.Esther.noDupl.3Guides <- Validation.list.Esther.noDupl.3Guides %>% filter(!ensembl_id == 'ENSMUSG00000118560')
#write.table(Validation.list.Esther.noDupl.3Guides, "Validation.list.Esther.noDupl.3Guides.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(Validation.list.Esther.noDupl.oldGuides, "Validation.list.Esther.noDupl.oldGuides.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#abc <- Validation.list.Esther.noDupl.3Guides %>% group_by(Gene) %>% count()

all.pools.extraINFO$StickyEnd.Fw <- 'AGGACCATTACGACCTGTTCTCGTCTCACACC' 
all.pools.extraINFO$StickyEnd.Rv <- 'GTTTAGAGACGTGATGAAGAGCCTAGCACCTC'
all.pools.extraINFO$OrderingOligo <- paste(all.pools.extraINFO$StickyEnd.Fw, all.pools.extraINFO$cloning.FW_oligo, all.pools.extraINFO$StickyEnd.Rv, sep='')

Validation.list.Esther.noDupl.oldGuides$Gene_guideNr <- paste(Validation.list.Esther.noDupl.oldGuides$Gene, Validation.list.Esther.noDupl.oldGuides$guide_number, sep = '_')
Validation.list.Esther.noDupl.oldGuides <- Validation.list.Esther.noDupl.oldGuides %>% filter(!ensembl_id == 'ENSMUSG00000116933')
Validation.list.Esther.noDupl.oldGuides <- Validation.list.Esther.noDupl.oldGuides %>% filter(!ensembl_id == 'ENSMUSG00000118560')

#invivo.guide.best.performing$Gene_guideNr <- paste(invivo.guide.best.performing$Gene, invivo.guide.best.performing$guide_number, sep = '_')
all.pools.extraINFO$Gene_guideNr <- paste(all.pools.extraINFO$gene, all.pools.extraINFO$Auto.pick.top.sgRNAs, sep = '_')

#Validation.list.Esther.noDupl.NotUsedguides <- all.pools.extraINFO %>% filter(!Gene_guideNr %in% invivo.guide.best.performing$Gene_guideNr) 
#names(Validation.list.Esther.noDupl.NotUsedguides)[names(Validation.list.Esther.noDupl.NotUsedguides) == 'gene'] <- 'Gene'
#names(Validation.list.Esther.noDupl.NotUsedguides)[names(Validation.list.Esther.noDupl.NotUsedguides) == 'Auto.pick.top.sgRNAs'] <- 'guide_number'

Validations.Esther.woEnriched.guides <- merge(Validations.Esther.woEnriched.genesOnly, all.pools.extraINFO, by='ensembl_id')
Validations.Esther.woEnriched.guides <- Validations.Esther.woEnriched.guides %>% filter(!ensembl_id == 'ENSMUSG00000116933')

Validations.Esther.vivoEnriched.guides <- merge(Validations.Esther.vivoEnriched.genesOnly, all.pools.extraINFO, by='ensembl_id')
Validations.Esther.Controls.guides <- merge(Validations.Esther.Controls.genesOnly, all.pools.extraINFO, by='ensembl_id')
Validations.Esther.Controls.guides <- Validations.Esther.Controls.guides %>% filter(!ensembl_id == 'ENSMUSG00000118560')

Validations.Esther.woEnriched.guides.onlyContr <- Validations.Esther.Controls.guides[1:365,]
Validations.Esther.woEnriched.guides.onlyContr <- Validations.Esther.woEnriched.guides.onlyContr[!duplicated(Validations.Esther.woEnriched.guides.onlyContr),]
Validations.Esther.woEnriched.guides.Contr <- rbind(Validations.Esther.woEnriched.guides, Validations.Esther.Controls.guides[1:365,])
Validations.Esther.woEnriched.guides.Contr <- Validations.Esther.woEnriched.guides.Contr[!duplicated(Validations.Esther.woEnriched.guides.Contr),]
#write.table(Validations.Esther.woEnriched.guides.onlyContr, "Validations.Esther.woEnriched.guides.onlyContr.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

Validations.Esther.vivoEnriched.guides.onlyContr <- Validations.Esther.Controls.guides[366:485,]
Validations.Esther.vivoEnriched.guides.onlyContr <- Validations.Esther.vivoEnriched.guides.onlyContr[!duplicated(Validations.Esther.vivoEnriched.guides.onlyContr),]
Validations.Esther.vivoEnriched.guides.Contr <- rbind(Validations.Esther.vivoEnriched.guides, Validations.Esther.Controls.guides[366:485,])
Validations.Esther.vivoEnriched.guides.Contr <- Validations.Esther.vivoEnriched.guides.Contr[!duplicated(Validations.Esther.vivoEnriched.guides.Contr),]
#write.table(Validations.Esther.vivoEnriched.guides.onlyContr, "Validations.Esther.vivoEnriched.guides.onlyContr.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

Validations.Esther.vivoEnriched.guides.onlyContr

# abc <- Validations.Esther.vivoEnriched.guides.Contr %>% group_by(gene) %>% count()
n_distinct(Validations.Esther.vivoEnriched.guides.Contr$id)

Validations.Esther.woEnriched.guides.Contr$StickyEnd.Fw <- 'AGGACCATTACGACCTGTTCTCGTCTCACACC' #stickyEnd primer pool D mouse
Validations.Esther.woEnriched.guides.Contr$StickyEnd.Rv <- 'GTTTAGAGACGTGATGAAGAGCCTAGCACCTC'  #stickyEnd primer pool D mouse
Validations.Esther.woEnriched.guides.Contr$OrderingOligo <- paste(Validations.Esther.woEnriched.guides.Contr$StickyEnd.Fw, Validations.Esther.woEnriched.guides.Contr$cloning.FW_oligo, Validations.Esther.woEnriched.guides.Contr$StickyEnd.Rv, sep='')
#write.table(Validations.Esther.woEnriched.guides.Contr, "Validations.Esther.woEnriched.guides.Contr.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
setwd('/Volumes/elling/Esther/NGS_results/Screen_results/20210816_Yumm_allApools_invivovitro_screen/Analysis/')
Validations.Esther.woEnriched.guides.Contr <- read.table('Validations.Esther.woEnriched.guides.Contr.txt', header = T, sep = "\t", quote = "")

Validations.Esther.vivoEnriched.guides.Contr$StickyEnd.Fw <- 'CCCCAGCCTTTACTCTTCCTACGTCTCACACC' #stickyEnd primer pool C mouse
Validations.Esther.vivoEnriched.guides.Contr$StickyEnd.Rv <- 'GTTTAGAGACGGACCTGGAACAATTCGCCTGG'  #stickyEnd primer pool C mouse
Validations.Esther.vivoEnriched.guides.Contr$OrderingOligo <- paste(Validations.Esther.vivoEnriched.guides.Contr$StickyEnd.Fw, Validations.Esther.vivoEnriched.guides.Contr$cloning.FW_oligo, Validations.Esther.vivoEnriched.guides.Contr$StickyEnd.Rv, sep='')
#write.table(Validations.Esther.vivoEnriched.guides.Contr, "Validations.Esther.vivoEnriched.guides.Contr.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
Validations.Esther.vivoEnriched.guides.Contr <- read.table('Validations.Esther.vivoEnriched.guides.Contr.txt', header = T, sep = "\t", quote = "")
