# Analysis for Sim et al 2026

# Load Libraries ----------------------------------------------------------

library(dplyr)
library(Biostrings)
library(stringr)
library(seqinr)
library(ggplot2)
library(gggenes)
library(DESeq2)
library(ggrepel)
library(ggfortify)
library(tidyr)
library(tibble)
library(microeco)
library(cowplot)
library(forcats)
library(Hmisc)
library(corrplot)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(scales)
library(ggnewscale)
library(ggh4x)
# GENERAL Color Sets for Plots ----------------------------------------------------
rnf_gene_colors <- list()
rnf_gene_colors <- c("#86f2d3","#86f2f2","#8aa0a0","#e8ef8f","#efa470","#6d5c4a")
names(rnf_gene_colors) <- c("RnfC","RnfD","RnfG","RnfE","RnfA","RnfB")

succinate_gene_colors <- list()
succinate_gene_colors <- c("#FEA1E0","#FEA1E0","#FBAF41","#FBAF41","#27A3DE","#27A3DE","#27A3DE","#27A3DE",
                           "#1162C1","#1162C1","#1162C1","#967bb6","#967bb6")
names(succinate_gene_colors) <- c("fumB","ttda_fumarase","frdA","frdB","flxA","flxB","flxC","flxD","hdrA","hdrB","hdrC","ltrA","erm")

RNA_seq_colors_hex <- c(
  "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", "#FF9896",
  "#9467BD", "#C5B0D5", "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7",
  "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"
)

WLP_HYD_colors <- list()
WLP_HYD_colors <- c("#6d5c4a","#6d5c4a","#6d5c4a","#6d5c4a","#6d5c4a", "#6d5c4a","#6d5c4a", "#6d5c4a","#6d5c4a","#6d5c4a",
                    "#FF9896", "#FF9896","#FF9896", "#FF9896",
                    "#FEA1E0","#FEA1E0","#FBAF41","#FBAF41","#27A3DE","#27A3DE","#27A3DE","#27A3DE","#1162C1","#1162C1","#1162C1","#967bb6","#967bb6",
                    "#e8ef8f","#BCBD22", "#DBDB8D","#e8ef8f","#BCBD22", "#DBDB8D")

names(WLP_HYD_colors) <- c("wlp3","wlp1","acsD","acsC","acsE","nqrF","fchA","folD","fhs","wlp2","hydB",
                           "hydM","hydA","hydC","fumB","ttda_fumarase","frdA","frdB","flxA","flxB","flxC","flxD","hdrA","hdrB","hdrC","ltrA","erm",
                           "RnfC","RnfD","RnfG","RnfE","RnfA","RnfB")

gene_color_set <- c(
  "GKGPCAPO_01438" = "#1162C2",
  "GKGPCAPO_01437" = "#1162C2",
  "GKGPCAPO_01436" = "#1162C2",
  "GKGPCAPO_01435" = "#27A3DF",
  "GKGPCAPO_01434" = "#27A3DF",
  "GKGPCAPO_01433" = "#27A3DF",
  "GKGPCAPO_01432" = "#27A3DF",
  "GKGPCAPO_01431" = "#FCAF41",
  "GKGPCAPO_01430" = "#FCAF41",
  "GKGPCAPO_01429" = "#FFA1E1",
  "GKGPCAPO_01428" = "#FFA1E1",
  "E.coli_K12_frdA" = "#529985",
  "E.coli_K12_frdB" = "#6E9F6D",
  "E.coli_K12_frdC" = "#81A665",
  "E.coli_K12_frdD" = "#ACB955",
  "B.theta frdC" = "#E7B04D",
  "B.theta frdA" = "#D08151",
  "B.theta frdB" = "#C26B51",
  "frdA prevotella" = "#ACA4E2",
  "frdB prevotella"  = "#55B8D0",
  "frdC prevotella"    = "#38BEB4",
  "Treponema succinifaciens DSM 2489, complete genome (frdA)"= "#AC7299",
  "Treponema succinifaciens DSM 2489, complete genome (frdB)" = "#D8B2C6",
  "Treponema succinifaciens DSM 2489, complete genome (frdC)" = "#C28AB1")

# GENERAL Final Plots -----------------------------------------------------

# FBEB-FR Clusters in Blautia:
Figure_1D
Figure_S1_Succinate

# RNF Clusters in Blautia
Figure_S1_RNF

# WLP HYD FBEB-FR Clusters in Blautia
Figure_1F

# Expression of FBEB-FR,WLP,HYD in log phase Blautia
Figure_3B

# Overall Transcriptome Changes with Loss of FBEB
Figure_4A
Figure_4B

# Changes in Specific Gene Clusters
Figure_4D

# Metagenomics Analysis
Figure_6f
Figure_6d

# GENERAL Load Genome Annotations & Other Metadata & Labels--------------------------------------------


    # Isolate genomes from Sorbara & Littmann Cell Host Microbe 2020 were downloaded from NCBI and annotated via Bakta and Prokka.  EBT Blautia strains were isolated
    # as described in Materials and Methods and annotated with Bakta and Prokka.
    
    # Set working direction then
    setwd("/Users/matthewsorbara/Dropbox/Faculty_Cloudtop/17_Lab_Papers/07_Succinate Paper/Manuscript 2026/R Analysis")
    
    # load bakta annotations
    load(file = "Sim_et_al_Isolate_Annotations.rdata")

    # Load Prokka Annotations across isolate biobanks
    load(file="Sim_et_al_Prokka_4910_Publications.rdata")
    
    # Load lists of Isolates for Main Text and Supplemental Figures
    figure_s1_working_group <- read.csv(file="blautia_working_group.csv",header=T,sep=",")
    figure_1_working_group <- read.csv(file="blautia_selected_group.csv",header=T,sep=",")
    
    # Load Labels for Schematic Heatmaps 
    labels <- read.csv(file="summary_operons.csv",header=T,sep=",")
    
# GENERAL: Output FASTA Used to Generate a Blast_Db from 303 Lachnospiraceae-----------------------------------------------

ms_lachno_genes <- total_annotations %>%
  filter(seq_id %in% total_working_list$seq_id)%>%
  select(locus_tag,seq_id)%>%
  left_join(total_sequences)%>%
  filter(!is.na(aa_sequence))%>%
  mutate(seq_index=paste("s",row_number(),sep = ""))

# Write the fasta file
for (i in 1:(nrow(ms_lachno_genes))){
  write.fasta(ms_lachno_genes[i,]$aa_sequence,ms_lachno_genes[i,]$seq_index,file.out="MS_lachno_genes.fasta",open="a")
}


# make blast db and run protein blast for succinate Gene Cluster, RNF Gene Cluster, WLP Gene Cluster and HYD Gene Cluster.

# First Make Blast DB:
    # makeblastdb -in MS_lachno_genes.fasta -title Lachno_DB -out Lachno.db -blastdb_version 5 -max_file_sz 4GB -dbtype prot
    
# Second Blast the different gene clusters
    # blastp -query succinate_gene_cluster.fasta -db Lachno.db -outfmt 6 -out Lachno_succinate_hits.txt -max_target_seqs 50000
    # blastp -query C_sporogene_RNF.faa -db Lachno.db -outfmt 6 -out Lachno_rnf_hits.txt -max_target_seqs 50000
    # blastp -query wlp_genes.fasta -db Lachno.db -outfmt 6 -out Lachno_wlp_hits.txt -max_target_seqs 50000
    # blastp -query hyd_genes.fasta -db Lachno.db -outfmt 6 -out Lachno_hyd_hits.txt -max_target_seqs 50000

# FIGURE 1D and S1 Succinate Operons and S1  ------------------------------------------------------

# Read in hits against the succinate gene cluster
      succinate_hits <- read.csv(file="Lachno_succinate_hits.txt",sep="\t",header=F)
      colnames(succinate_hits) <- c("qseqid","seq_index", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Identify the best hit per locus tag.
operon_targets <- ms_lachno_genes%>%
  select(seq_id,seq_index,locus_tag)%>%
  left_join(succinate_hits)%>%
  filter(!is.na(qseqid))%>%
  dplyr::arrange(desc(pident))%>%
  dplyr::group_by(locus_tag)%>%
  dplyr::slice(1)%>%
  dplyr::ungroup()%>%
  select(seq_id,locus_tag,gene_label=qseqid)

# Identify the best match to the frdB from the gene cluster,
# will use that to orient the operon plots 

      best_hit_frdB <- ms_lachno_genes %>%
        select(seq_id,seq_index,locus_tag)%>%
        left_join(succinate_hits)%>%
        filter(qseqid=="putative_frdB")%>%
        filter(!is.na(qseqid))%>%
        dplyr::arrange(desc(pident))%>%
        dplyr::group_by(locus_tag)%>%
        dplyr::slice(1)%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(pident))%>%
        dplyr::group_by(seq_id)%>%
        dplyr::slice(1)%>%
        dplyr::ungroup()%>%
        select(locus_tag,seq_id)
      
       
      flip <- total_annotations %>%
        filter(locus_tag %in% best_hit_frdB$locus_tag)%>%
        filter(strand=="+")
        
      flipped_orientations <- flip$seq_id
      rm(best_hit_frdB)

# Identify the largest cluster of Blast hits for the 11 genes
      
operon_df <- total_annotations %>%
  filter(locus_tag %in% operon_targets$locus_tag)%>%
  dplyr::arrange(start)%>%
  dplyr::group_by(seq_id,contig)%>%
  dplyr::mutate(difference=start-lag(start))%>% 
  dplyr::mutate(difference=ifelse(is.na(difference),0,difference))%>%
  dplyr::mutate(new_operon=ifelse(difference>10000,1,0))%>%
  dplyr::mutate(new_operon=ifelse(is.na(new_operon),0,new_operon))%>%
  dplyr::mutate(cumulative_new_operon=cumsum(new_operon))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id,contig,cumulative_new_operon)%>%
  dplyr::summarize(earliest_start = min(start)-5000,latest_end = max(stop)+5000,num_genes=n())%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id)%>%
  dplyr::arrange(desc(num_genes))%>%
  dplyr::mutate(ordered_operon=seq(1:n()))%>%
  left_join(total_annotations %>%
              select(seq_id,contig,start,stop,strand,product,gene,locus_tag))%>%
  filter(start>=earliest_start & stop <=latest_end)%>%
  dplyr::mutate(start_plot_x=ifelse(seq_id %in% flipped_orientations,start-earliest_start,latest_end-stop))%>%
  dplyr::mutate(end_plot_x=ifelse(seq_id %in% flipped_orientations,stop-earliest_start,latest_end-start))%>%
  dplyr::mutate(strand=ifelse(strand=="+",TRUE,FALSE))%>%
  dplyr::mutate(strand=ifelse(!seq_id %in% flipped_orientations,
                              ifelse(strand==TRUE,FALSE,TRUE),
                              strand))%>%
  filter(ordered_operon==1)

# Manual inspection to identify list of locus tags encoding HydABC/M genes that get identified as hits to Hdr-FLX subunits

list_to_correct <- c("BGDCJC_02300","BGDCJC_02301","HMOEIF_00768","HMOEIF_00769","AMCACH_01280",
                     "JMJOBN_00125","KLHFCO_03070","KLHFCO_03071","AIGJFA_03187","AIGJFA_03188",
                     "FOCAJE_03175","FOCAJE_03176","PHIFMM_03067","PHIFMM_03068",
                     "HOBNBL_02042","HOBNBL_02043","NJOFKI_03991","NJOFKI_03992",
                     "KKMGGJ_02073","KKMGGJ_02074","KKMGGJ_02075","KKMGGJ_02076",
                     "FOJDPM_01939","DGAOLD_01240",
                     "KCLKCA_01566","LKKDIN_02284","LKKDIN_02285","LKKDIN_02286","LKKDIN_02287","NGNIBL_00758")

Figure_S1_Succinate <- operon_df%>%
  left_join(operon_targets)%>%
  left_join(total_working_list)%>%
  mutate(gene_label=gsub("putative_","",gene_label))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  filter(grepl("Blautia",Species))%>%
  filter(msk_id %in%figure_s1_working_group$msk_id)%>%
  mutate(Species=ifelse(grepl("massiliensis",Species),"Blautia luti",Species))%>%
  mutate(gene_label=ifelse(locus_tag %in% list_to_correct,NA,gene_label))%>%
  ggplot(aes(xmin=start_plot_x,xmax=end_plot_x,y=msk_id))+
      geom_gene_arrow(aes(fill=gene_label,forward=strand))+

      facet_grid(Species~.,scales="free_y",space="free_y")+
   
      scale_fill_manual(values=succinate_gene_colors,na.value = "white")+

      xlab(label="Relative Position")+
      ylab(label="Isolates")+
      scale_x_continuous(expand=c(0,0))+
      theme_minimal()+
      theme(panel.border = element_rect(fill=NA,color="black"),
            legend.position = "none",
           
            strip.text.y = element_text(angle=0,hjust=0,vjust=0.5),
            strip.text = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            axis.title.y = element_blank())

# Output at 8 in x 6 inch to Match Sizing.

# PLOT FOR FIGURE 1D

p1 <- operon_df%>%
  left_join(operon_targets)%>%
  left_join(total_working_list)%>%
  mutate(gene_label=gsub("putative_","",gene_label))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  filter(grepl("Blautia",Species))%>%
  filter(msk_id %in%figure_1_working_group$msk_id)%>%
  mutate(Species=ifelse(grepl("massiliensis",Species),"Blautia luti",Species))%>%
  mutate(gene_label=ifelse(locus_tag %in% list_to_correct,NA,gene_label))%>%
  mutate(Species=factor(Species,levels=c("Blautia faecis",
                                         "Blautia glucerasea",
                                         "Blautia luti",
                                         "Blautia schinkii",
                                         "Blautia obeum",
                                         "Blautia wexlerae")))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  mutate(msk_id=factor(msk_id,
                       levels=c("MSK.11.45",
                                "EBT232",
                                "MSK.20.42",
                                "MSK.14.58",
                                "MSK.13.40",
                                "EBT147",
                                "MSK.15.25",
                                "EBT100",
                                "MSK.17.81",
                                "MSK.20.66",
                                "EBT223",
                                "EBT159",
                                "MSK.21.50",
                                "EBT193")))%>%
  ggplot(aes(xmin=start_plot_x,xmax=end_plot_x,y=msk_id))+
      geom_gene_arrow(aes(fill=gene_label,forward=strand))+
      facet_grid(Species~.,scales="free_y",space="free_y")+
      scale_fill_manual(values=succinate_gene_colors,na.value = "white")+
      xlab(label="Relative Position")+
      ylab(label="Isolates")+
      scale_x_continuous(expand=c(0,0))+
      theme_minimal()+
      theme(panel.border = element_rect(fill=NA,color="black"),
            legend.position = "none",
            axis.text.y = element_blank(),
            strip.text.y = element_text(angle=0,hjust=0,vjust=0.5),
            strip.text = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            axis.title.y = element_blank())

p2 <- operon_df%>%
  left_join(operon_targets)%>%
  left_join(total_working_list)%>%
  mutate(gene_label=gsub("putative_","",gene_label))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  filter(grepl("Blautia",Species))%>%
  filter(msk_id %in%figure_1_working_group$msk_id)%>%
  
  left_join(labels %>%
              select(-FBEB_Cluster)%>%
              gather(key="characteristic",value="value",-msk_id))%>%
  
  mutate(Species=ifelse(grepl("massiliensis",Species),"Blautia luti",Species))%>%
  mutate(Species=as.factor(Species))%>%
  mutate(Species=factor(Species,levels=c("Blautia faecis",
                                         "Blautia glucerasea",
                                         "Blautia luti",
                                         "Blautia schinkii",
                                         "Blautia obeum",
                                         "Blautia wexlerae")))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  mutate(msk_id=factor(msk_id,
                       levels=c("MSK.11.45",
                                "EBT232",
                                "MSK.20.42",
                                "MSK.14.58",
                                "MSK.13.40",
                                "EBT147",
                                "MSK.15.25",
                                "EBT100",
                                "MSK.17.81",
                                "MSK.20.66",
                                "EBT223",
                                "EBT159",
                                "MSK.21.50",
                                "EBT193")))%>%
  mutate(characteristic=gsub("\\_"," ",characteristic))%>%
  ggplot(aes(x=characteristic,y=msk_id,fill=as.factor(value)))+
      geom_tile(color="black")+
      facet_grid(Species~ordered_operon,scales="free_y",space="free_y")+
      scale_fill_manual(values=c("white","skyblue"),na.value = "white")+
      theme_minimal()+
      theme(panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            axis.text.x=element_text(angle=90,hjust=1),
            axis.title.x = element_blank(),
            strip.text.y = element_blank(),
            strip.text = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            axis.title.y = element_blank())


Figure_1D <- plot_grid(p2,p1,ncol=2,align = "h",axis="tb",rel_widths = c(0.7,4.5))

# Figure S1 RNF Gene Cluster  -----------------------------------------------

# Read in Blast hits
    rnf_hits <- read.csv(file="Lachno_rnf_hits.txt",sep="\t",header=F)
    colnames(rnf_hits) <- c("qseqid","seq_index", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Identify the best hit per locus tag.
operon_targets <- ms_lachno_genes%>%
  select(seq_id,seq_index,locus_tag)%>%
  left_join(rnf_hits)%>%
  filter(!is.na(qseqid))%>%
  dplyr::arrange(desc(pident))%>%
  dplyr::group_by(locus_tag)%>%
  dplyr::slice(1)%>%
  dplyr::ungroup()%>%
  select(seq_id,locus_tag,gene_label=qseqid)

# Identify the best match to the rnfD from the gene cluster,
# will use that to orient the operon plots 
      best_hit_rnfD <- ms_lachno_genes %>%
        select(seq_id,seq_index,locus_tag)%>%
        left_join(rnf_hits)%>%
        filter(qseqid=="RnfD")%>%
        filter(!is.na(qseqid))%>%
        dplyr::arrange(desc(pident))%>%
        dplyr::group_by(locus_tag)%>%
        dplyr::slice(1)%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(pident))%>%
        dplyr::group_by(seq_id)%>%
        dplyr::slice(1)%>%
        dplyr::ungroup()%>%
        select(locus_tag,seq_id)
      
      flip <- total_annotations %>%
        filter(locus_tag %in% best_hit_rnfD$locus_tag)%>%
        filter(strand=="+")
      
      flipped_orientations <- flip$seq_id
      rm(best_hit_rnfD)

# Identify the largest cluster of Blast hits for the RNF Complex
operon_df <- total_annotations %>%
  filter(locus_tag %in% operon_targets$locus_tag)%>%
  dplyr::arrange(start)%>%
  dplyr::group_by(seq_id,contig)%>%
  dplyr::mutate(difference=start-lag(start))%>% 
  dplyr::mutate(difference=ifelse(is.na(difference),0,difference))%>%
  dplyr::mutate(new_operon=ifelse(difference>10000,1,0))%>%
  dplyr::mutate(new_operon=ifelse(is.na(new_operon),0,new_operon))%>%
  dplyr::mutate(cumulative_new_operon=cumsum(new_operon))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id,contig,cumulative_new_operon)%>%
  dplyr::summarize(earliest_start = min(start)-5000,latest_end = max(stop)+5000,num_genes=n())%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id)%>%
  dplyr::arrange(desc(num_genes))%>%
  dplyr::mutate(ordered_operon=seq(1:n()))%>%
  left_join(total_annotations %>%
              select(seq_id,contig,start,stop,strand,product,gene,locus_tag))%>%
  filter(start>=earliest_start & stop <=latest_end)%>%
  dplyr::mutate(start_plot_x=ifelse(seq_id %in% flipped_orientations,start-earliest_start,latest_end-stop))%>%
  dplyr::mutate(end_plot_x=ifelse(seq_id %in% flipped_orientations,stop-earliest_start,latest_end-start))%>%
  dplyr::mutate(strand=ifelse(strand=="+",TRUE,FALSE))%>%
  dplyr::mutate(strand=ifelse(!seq_id %in% flipped_orientations,
                              ifelse(strand==TRUE,FALSE,TRUE),
                              strand))%>%
  filter(ordered_operon==1)


Figure_S1_RNF <- operon_df%>%
  left_join(operon_targets)%>%
  left_join(total_working_list)%>%
  mutate(gene_label=gsub("putative_","",gene_label))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  filter(grepl("Blautia",Species))%>%
  filter(msk_id %in%figure_s1_working_group$msk_id)%>%
  mutate(Species=ifelse(grepl("massiliensis",Species),"Blautia luti",Species))%>%
  ggplot(aes(xmin=start_plot_x,xmax=end_plot_x,y=msk_id))+
  geom_gene_arrow(aes(fill=gene_label,forward=strand))+
  facet_grid(Species~.,scales="free_y",space="free_y")+
  scale_fill_manual(values=rnf_gene_colors,na.value = "white")+
  xlab(label="Relative Position")+
  ylab(label="Isolates")+
  scale_x_continuous(expand=c(0,0))+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA,color="black"),
        legend.position = "none",
        strip.text.y = element_text(angle=0,hjust=0,vjust=0.5),
        strip.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        axis.title.y = element_blank())



# Figure 1F WLP and HYD Genes ----------------------------------------------

# Read in Blast hits
    WLP_hits <- read.csv(file="Lachno_WLP_hits.txt",sep="\t",header=F)
    colnames(WLP_hits) <- c("qseqid","seq_index", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

    HYD_hits <- read.csv(file="Lachno_HYD_hits.txt",sep="\t",header=F)
    colnames(HYD_hits) <- c("qseqid","seq_index", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Identify the best hit per locus tag across HYD, WLP, FBEB-FR, RNF Blast results
    operon_targets <- ms_lachno_genes%>%
      select(seq_id,seq_index,locus_tag)%>%
      left_join(WLP_hits)%>%
      filter(!is.na(qseqid))%>%
      dplyr::arrange(desc(pident))%>%
      dplyr::group_by(locus_tag)%>%
      dplyr::slice(1)%>%
      dplyr::ungroup()%>%
      select(pident,seq_id,locus_tag,gene_label=qseqid)%>%
      rbind(ms_lachno_genes%>%
              select(seq_id,seq_index,locus_tag)%>%
              left_join(HYD_hits)%>%
              filter(!is.na(qseqid))%>%
              dplyr::arrange(desc(pident))%>%
              dplyr::group_by(locus_tag)%>%
              dplyr::slice(1)%>%
              dplyr::ungroup()%>%
              select(pident,seq_id,locus_tag,gene_label=qseqid))%>%
      rbind(ms_lachno_genes%>%
              select(seq_id,seq_index,locus_tag)%>%
              left_join(succinate_hits)%>%
              filter(!is.na(qseqid))%>%
              dplyr::arrange(desc(pident))%>%
              dplyr::group_by(locus_tag)%>%
              dplyr::slice(1)%>%
              dplyr::ungroup()%>%
              select(pident,seq_id,locus_tag,gene_label=qseqid))%>%
      rbind(ms_lachno_genes%>%
              select(seq_id,seq_index,locus_tag)%>%
              left_join(rnf_hits)%>%
              filter(!is.na(qseqid))%>%
              dplyr::arrange(desc(pident))%>%
              dplyr::group_by(locus_tag)%>%
              dplyr::slice(1)%>%
              dplyr::ungroup()%>%
              select(pident,seq_id,locus_tag,gene_label=qseqid))%>%
      dplyr::arrange(desc(pident))%>%
      dplyr::group_by(locus_tag)%>%
      dplyr::slice(1)
      

# Identify the best match to the ascC,
# will use that to orient the operon plots 

      best_hit_acsC <- ms_lachno_genes %>%
        select(seq_id,seq_index,locus_tag)%>%
        left_join(WLP_hits)%>%
        filter(qseqid=="acsC")%>%
        filter(!is.na(qseqid))%>%
        dplyr::arrange(desc(pident))%>%
        dplyr::group_by(locus_tag)%>%
        dplyr::slice(1)%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(pident))%>%
        dplyr::group_by(seq_id)%>%
        dplyr::slice(1)%>%
        dplyr::ungroup()%>%
        select(locus_tag,seq_id)
      
      flip <- total_annotations %>%
        filter(locus_tag %in% best_hit_acsC$locus_tag)%>%
        filter(strand=="+")
      
      flipped_orientations <- flip$seq_id
      rm(best_hit_acsC)
      
      # Add two Blautia wexlerae that do not have ascC beside HydABC
      flipped_orientations <- c(flipped_orientations,"GCF_013304405.1")
      flipped_orientations <- c(flipped_orientations,"EBT193")

# Identify the largest cluster of Blast hits
operon_df <- total_annotations %>%
  filter(locus_tag %in% operon_targets$locus_tag)%>%
  dplyr::arrange(start)%>%
  dplyr::group_by(seq_id,contig)%>%
  dplyr::mutate(difference=start-lag(start))%>% 
  dplyr::mutate(difference=ifelse(is.na(difference),0,difference))%>%
  dplyr::mutate(new_operon=ifelse(difference>10000,1,0))%>%
  dplyr::mutate(new_operon=ifelse(is.na(new_operon),0,new_operon))%>%
  dplyr::mutate(cumulative_new_operon=cumsum(new_operon))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id,contig,cumulative_new_operon)%>%
  dplyr::summarize(earliest_start = min(start)-5000,latest_end = max(stop)+5000,num_genes=n())%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id)%>%
  dplyr::arrange(desc(num_genes))%>%
  dplyr::mutate(ordered_operon=seq(1:n()))%>%
  left_join(total_annotations %>%
              select(seq_id,contig,start,stop,strand,product,gene,locus_tag))%>%
  filter(start>=earliest_start & stop <=latest_end)%>%
  dplyr::mutate(start_plot_x=ifelse(seq_id %in% flipped_orientations,start-earliest_start,latest_end-stop))%>%
  dplyr::mutate(end_plot_x=ifelse(seq_id %in% flipped_orientations,stop-earliest_start,latest_end-start))%>%
  dplyr::mutate(strand=ifelse(strand=="+",TRUE,FALSE))%>%
  dplyr::mutate(strand=ifelse(!seq_id %in% flipped_orientations,
                              ifelse(strand==TRUE,FALSE,TRUE),
                              strand))%>%
  filter(ordered_operon<2)

# Generate Plots

p1 <- operon_df%>%
  left_join(operon_targets)%>%
  left_join(total_working_list)%>%
  mutate(gene_label=gsub("putative_","",gene_label))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  filter(grepl("Blautia",Species))%>%
  filter(msk_id %in%figure_1_working_group$msk_id)%>%
  mutate(Species=ifelse(grepl("massiliensis",Species),"Blautia luti",Species))%>%
  mutate(Species=as.factor(Species))%>%
  mutate(Species=factor(Species,levels=c("Blautia faecis",
                                         "Blautia glucerasea",
                                         "Blautia luti",
                                         "Blautia schinkii",
                                         "Blautia obeum",
                                         "Blautia wexlerae")))%>%
  ggplot(aes(xmin=start_plot_x,xmax=end_plot_x,y=msk_id))+
      geom_gene_arrow(aes(fill=gene_label,forward=strand))+
      facet_grid(Species~ordered_operon,scales="free_y",space="free_y")+
      scale_fill_manual(values=WLP_HYD_colors,na.value = "white")+
      xlab(label="Relative Position")+
      ylab(label="Isolates")+
      scale_x_continuous(expand=c(0,0))+
      theme_minimal()+
      theme(panel.border = element_rect(fill=NA,color="black"),
            legend.position = "none",
            axis.text.y = element_blank(),
            strip.text.y = element_text(angle=0,hjust=0,vjust=0.5),
            strip.text = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            axis.title.y = element_blank())
p2 <- operon_df%>%
  left_join(operon_targets)%>%
  left_join(total_working_list)%>%
  mutate(gene_label=gsub("putative_","",gene_label))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  filter(grepl("Blautia",Species))%>%
  filter(msk_id %in%figure_1_working_group$msk_id)%>%
  left_join(labels %>%
              gather(key="characteristic",value="value",-msk_id))%>%
  mutate(Species=ifelse(grepl("massiliensis",Species),"Blautia luti",Species))%>%
  mutate(Species=as.factor(Species))%>%
  mutate(Species=factor(Species,levels=c("Blautia faecis",
                                         "Blautia glucerasea",
                                         "Blautia luti",
                                         "Blautia schinkii",
                                         "Blautia obeum",
                                         "Blautia wexlerae")))%>%
  mutate(msk_id=as.factor(msk_id))%>%
  mutate(msk_id=factor(msk_id,
                       levels=c("MSK.11.45",
                                "EBT232",
                                "MSK.20.42",
                                "MSK.14.58",
                                "MSK.13.40",
                                "EBT147",
                                "MSK.15.25",
                                "EBT100",
                                "MSK.17.81",
                                "MSK.20.66",
                                "EBT223",
                                "EBT159",
                                "MSK.21.50",
                                "EBT193")))%>%
  mutate(characteristic=gsub("\\_"," ",characteristic))%>%
  ggplot(aes(x=characteristic,y=msk_id,fill=as.factor(value)))+
  geom_tile(color="black")+
  facet_grid(Species~ordered_operon,scales="free_y",space="free_y")+
  scale_fill_manual(values=c("white","skyblue"),na.value = "white")+
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x=element_text(angle=90,hjust=1),
        axis.title.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        axis.title.y = element_blank())

Figure_1F <- plot_grid(p2,p1,ncol=2,align = "h",axis="tb",rel_widths = c(0.75,5))
Figure_1F
# Figure 3B and Figure 4 RNA Seq Analysis ------------------------------------------------

# Read in outputs and Process
    rna_counts <- read.csv("16s_v_hdrb.csv") %>%
      select(!"X")
    
    rna_colData <- as.data.frame(matrix(c("16s", "16s", "16s", "hdrB", "hdrB", "hdrB"), nrow=6,ncol=1))
    
    row.names(rna_colData) <- c("m1_16S", "m2_16s", "m3_16s", "m4_hdrB","m5_hdrB", "m6_hdrB")
    colnames(rna_colData) <- c('condition')
    
    dds <- DESeqDataSetFromMatrix(countData = rna_counts, colData = rna_colData, design = ~condition, tidy = TRUE)

    dds <- dds[rowSums(counts(dds)) > 10, ]
    
    dds <- DESeq(dds, betaPrior = TRUE, fitType = 'parametric')                             
    dds_table_counts_normalized <- counts(dds, normalized = TRUE)
    
    res <- results(dds)
    head(results(dds, tidy = TRUE))
    
    rld <- rlog(dds, blind = FALSE)
    
    norm_counts <- as.data.frame(assay(rld))
    
# Overall Transcriptome Changes
    
Figure_4A <- plotPCA(rld, intgroup = c("condition")) +
  theme_bw(base_size = 12) +
  geom_point(pch=21,size = 4,color="black") +
  scale_color_manual(values=c("gray","#967bb6"))+
  theme(aspect.ratio = 1,
        legend.position = "none") + 
  ggtitle('Transcriptional Response')
Figure_4A

# Test for Significance
      
      LRT_dds <- DESeq(dds, test = "LRT", reduced = ~1)
      
      res_LRT <- as.data.frame(results(LRT_dds))
      
      res_LRT$gene <- rownames(res_LRT)
      
      res_sign <- filter(res_LRT, padj < .05 & abs(log2FoldChange) > 1)
      
      res_sign <- res_sign$gene
      
      norm_counts_LRT <- norm_counts[which(rownames(norm_counts) %in% res_sign), ]
      LRT_prcomp <- prcomp(x = t(norm_counts_LRT), center = TRUE, scale = TRUE)
      
      
      t_norm_counts_LRT <- as.data.frame(t(norm_counts_LRT)) %>%
        mutate(Strain = str_extract(rownames(LRT_prcomp$x), pattern = "16s|16S|sdhE")) %>%
        mutate(Strain=ifelse(Strain=="16s","16S",Strain))
  
      
      df_16_vs_hdrB <- as.data.frame(results(dds, contrast = c("condition", "16s", "hdrB"))) %>%
        rownames_to_column('Gene') %>%
        mutate(Comparison = "df_16_vs_hdrB")%>%
        mutate(Significant = case_when(
          padj > 0.05 | (abs(log2FoldChange) < 1) ~ 'Not Sign',
          padj < 0.05 & log2FoldChange > 1 ~ 'Sign',
          padj < 0.05 & log2FoldChange < -1 ~ 'Sign'
        )) %>%
        mutate(locus_tag=gsub("n[0-9]*","",Gene))



# Make a volcano plot with the genes in the top X clusters labelled.  Note; colors being manually assigned, so would need to adjust/add colors to label additional"

volcano_annotations <- read.csv(file="Sim_et_al_Volcano_Annotations.csv",header=T,sep=",")

transcriptome_df <- df_16_vs_hdrB  %>%
  left_join(volcano_annotations%>%
              select(locus_tag,Label.on.volcano.plot.,Cluster.))

FBEB <- c( "GKGPCAPO_01438", "GKGPCAPO_01437", "GKGPCAPO_01436", "GKGPCAPO_01435",
           "GKGPCAPO_01434", "GKGPCAPO_01433", "GKGPCAPO_01432", "GKGPCAPO_01431",
           "GKGPCAPO_01430", "GKGPCAPO_01429", "GKGPCAPO_01428")

RNF <- c("GKGPCAPO_00263",
         "GKGPCAPO_00264",
         "GKGPCAPO_00265",
         "GKGPCAPO_00266",
         "GKGPCAPO_00267",
         "GKGPCAPO_00268")

HYD <- c("GKGPCAPO_02059",
         "GKGPCAPO_02060",
         "GKGPCAPO_02061",
         "GKGPCAPO_02062")

WLP <- c("GKGPCAPO_02048","GKGPCAPO_02049","GKGPCAPO_02050","GKGPCAPO_02051","GKGPCAPO_02052","GKGPCAPO_02053","GKGPCAPO_02054",
         "GKGPCAPO_02055","GKGPCAPO_02056","GKGPCAPO_02057","GKGPCAPO_02916")

label_colors_ms <- c("yes (other-yellow)","yes (orange)","no","yes (pink)","yes (green)",
                     "yes (red)","yes (purple)","yes (black)","yes (light blue)","yes (brown)",
                     "yes (dark blue)","WLP","RNF","HYD")
volcano_labels <- c("#FFFF9B","#F6C6AD","white",
                     "#E59EDD","#B4E6A2","#C1584B",
                     "#9B81AF","black","#53C9E9",
                    "#6D5C49","#4150A4","#6D5C49","#A0EFD5","#FF9896")

names(volcano_labels) <- label_colors_ms

labelled_genes <- transcriptome_df %>%
  mutate(Label.on.volcano.plot.=ifelse(Label.on.volcano.plot.=="","no",Label.on.volcano.plot.))%>%
  mutate(Label.on.volcano.plot.=ifelse(locus_tag %in% FBEB,"yes (purple)",Label.on.volcano.plot.))%>%
  mutate(Label.on.volcano.plot.=ifelse(locus_tag %in% WLP,"WLP",Label.on.volcano.plot.))%>%
  mutate(Label.on.volcano.plot.=ifelse(locus_tag %in% RNF,"RNF",Label.on.volcano.plot.))%>%
  mutate(Label.on.volcano.plot.=ifelse(locus_tag %in% HYD,"HYD",Label.on.volcano.plot.))%>%
  filter(!is.na(Label.on.volcano.plot.))%>%
  filter(Label.on.volcano.plot.!="no")

Figure_4B <- transcriptome_df %>%
  ggplot(aes(x=log2FoldChange,y=-log10(padj),shape=Significant))+
    geom_point(size=2,fill="white")+
    geom_point(data=labelled_genes,aes(fill=Label.on.volcano.plot.),size=3)+
  geom_vline(xintercept = -1, colour = "#383838", linetype = "dashed") +
  geom_vline(xintercept = 1, colour = "#383838", linetype = "dashed")+
  scale_fill_manual(values=volcano_labels,na.value = "#e6e6e6")+
  scale_shape_manual(values=c(22,21))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

Figure_4B 

# Plot Transcript Abundance for Genes in WLP,HYD,FBEB-FR in Log Phase Blautia

rna_counts <- read.csv("16s_v_hdrb.csv") %>%
  select(!"X")%>%
  mutate(locus_tag=gsub("n[0-9]*","",Locus))%>%
  select(-Locus)%>%
  filter(grepl("GKGPCAPO",locus_tag))

t <- publication_annotations %>%
  filter(locus_tag %in% rna_counts$locus_tag)%>%
  select(locus_tag,length_bp)

total_counts <- rna_counts %>%
  gather(key="condition",value="reads",-locus_tag)%>%
  dplyr::group_by(condition)%>%
  dplyr::summarise(total_counts=sum(reads))

rna_table <- rna_counts %>%
  gather(key="condition",value="reads",-locus_tag)%>%
  left_join(t)%>%
  left_join(total_counts)

rna_table <- rna_table %>%
  mutate(rpkm = reads / ((length_bp /1000) * (total_counts/1000000)))

rna_table <- rna_table %>%
  filter(grepl("16",condition))%>%
  group_by(locus_tag)%>%
  dplyr::summarise(average_rpkm=mean(rpkm))%>%
  ungroup()%>%
  dplyr::arrange(average_rpkm)%>%
  mutate(rank=row_number())%>%
  mutate(percentile=(100*rank/3508))
  
labels <- rna_table %>%
  mutate(label=ifelse(locus_tag %in% WLP,"WLP",
                      ifelse(locus_tag %in% HYD,"HYD",
                             ifelse(locus_tag %in% FBEB,"FBEB",
                                    NA))))%>%
  filter(!is.na(label))

Figure_3B <- rna_table %>%
    ggplot(aes(x=percentile,y=log(average_rpkm,10)))+
          geom_point(size=1,pch=21,fill="white",color="black",stroke=0.2)+
          geom_point(data=labels,aes(fill=label,color=label),pch=21,size=3)+
          
          scale_fill_manual(values=c("#9B81AF","#FF9896","#53C9E9"))+
          scale_color_manual(values=c("#9B81AF","#FF9896","#53C9E9"))+
          ylab(label="RPKM (log 10)")+
          xlab(label="Percentile")+
          theme_bw()
Figure_3B
# Calculating average percentile by cluster:

labels %>%
  dplyr::group_by(label)%>%
  dplyr::summarise(average_percentile=mean(percentile))

# Figure 4D RNA Seq Cluster Plots ---------------------------------------------------

volcano_annotations <- read.csv(file="Sim_et_al_Volcano_Annotations.csv",header=T,sep=",")

target_RNA_cluster <- volcano_annotations %>%
  filter(Cluster.=="yes")

# Create list of targets
operon_targets <- target_RNA_cluster %>%
  select(locus_tag)%>%
  left_join(publication_annotations)

# Now map operons

flipped_orientations <- "TM475"

operon_df <- operon_targets %>%
  select(locus_tag)%>%
  left_join(publicationlocations)%>%
  select(-source,-note,-phase,-filename,-score)%>%
  
  dplyr::select(seq_id,seqid)%>%
  unique()%>%
  dplyr::group_by(seq_id)%>%
  dplyr::mutate(contig_num=seq(1:n()))%>%
  dplyr::left_join(publicationlocations %>%
                     select(-source,-note,-phase,-filename,-score)%>%
                     filter(locus_tag %in% operon_targets$locus_tag))%>%
  dplyr::ungroup()%>%
  dplyr::arrange(start)%>%
  dplyr::group_by(seq_id,seqid,contig_num)%>%
  dplyr::mutate(difference=start-lag(start))%>% 
  dplyr::mutate(difference=ifelse(is.na(difference),0,difference))%>%
  dplyr::mutate(new_operon=ifelse(difference>10000,1,0))%>%
  dplyr::mutate(new_operon=ifelse(is.na(new_operon),0,new_operon))%>%
  dplyr::mutate(cumulative_new_operon=cumsum(new_operon))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id,seqid,contig_num,cumulative_new_operon)%>%
  dplyr::summarize(earliest_start = min(start)-5000,latest_end = max(end)+5000,num_genes=n())%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id)%>%
  dplyr::arrange(desc(num_genes))%>%
  dplyr::mutate(ordered_operon=seq(1:n()))%>%
  left_join(publicationlocations %>%
              select(seq_id,seqid,start,end,strand,product,gene,eC_number,locus_tag))%>%
  filter(start>=earliest_start & end <=latest_end)%>%
  dplyr::mutate(start_plot_x=ifelse(seq_id %in% flipped_orientations,start-earliest_start,latest_end-end))%>%
  dplyr::mutate(end_plot_x=ifelse(seq_id %in% flipped_orientations,end-earliest_start,latest_end-start))%>%
  dplyr::mutate(strand=ifelse(strand=="+",TRUE,FALSE))%>%
  dplyr::mutate(strand=ifelse(!seq_id %in% flipped_orientations,
                              ifelse(strand==TRUE,FALSE,TRUE),
                              strand))%>%
  mutate(ordered_operon=as.factor(ordered_operon))


p2 <- operon_df%>%
  left_join(df_16_vs_sdh%>%
              select(locus_tag,log2FoldChange))%>%
  mutate(Significant=ifelse(locus_tag%in%target_RNA_cluster$locus_tag,"Yes","No"))%>%
  mutate(ordered_operon=as.factor(ordered_operon))%>%
  mutate(ordered_operon=factor(ordered_operon,levels=c("1","3","4","5","2","6","7","8")))%>%
  filter(ordered_operon %in% c("1","3","4","5"))%>%
  ggplot(aes(xmin=start_plot_x,xmax=end_plot_x,y=fct_rev(ordered_operon)))+
  geom_gene_arrow(aes(fill=log2FoldChange,forward=strand,color=Significant,size=Significant))+
  # geom_gene_label(aes(label=gene),size=5)+
  xlab(label="Relative Position")+
  ylab(label="Isolates")+
  scale_fill_gradient2(low = "darkred", mid = "white", high = "skyblue1", midpoint = 0)+
  scale_x_continuous(expand=c(0,0))+
  scale_size_manual(values=c(0.2,0.6))+
  scale_color_manual(values=c("black","forestgreen"),na.value="#e6e6e6")+
  theme_minimal()+
  theme(
    panel.border = element_rect(fill=NA,color="black"),
    panel.grid.major.y = element_line(color="black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
 
    plot.margin = unit(c(0, 0, 0, 0), "lines"),
    axis.title.y = element_blank(),
    legend.position = "bottom")+
  theme(axis.line = element_line(colour="black"), axis.ticks = element_line(colour="black"))+
  theme(legend.position = "bottom",axis.text.x = element_text(angle=0),
        plot.title = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank())
p1 <- operon_df%>%
  left_join(df_16_vs_hdrB%>%
              select(locus_tag,log2FoldChange))%>%
  mutate(Significant=ifelse(locus_tag%in%target_RNA_cluster$locus_tag,"Yes","No"))%>%
  left_join(volcano_annotations%>%
              select(locus_tag,Label.on.volcano.plot.))%>%
  filter(!is.na(Label.on.volcano.plot.))%>%
  mutate(ordered_operon=as.factor(ordered_operon))%>%
  mutate(ordered_operon=factor(ordered_operon,levels=c("1","3","4","5","2","6","7","8")))%>%
  filter(ordered_operon %in% c("1","3","4","5"))%>%
  ggplot(aes(x=1,y=fct_rev(ordered_operon),fill=Label.on.volcano.plot.))+
        geom_tile(color="black")+
          scale_fill_manual(values=volcano_labels)+
        theme_minimal()+
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                legend.position = "none")

Figure_4D <- cowplot::plot_grid(p1,p2,ncol=2,align="h",axis="tb",rel_widths = c(1,13))
  
Figure_4D
# Figure 6BC Clusters Across Isolates  ---------------------------------------

all_protein_sequences <-  publication_isolates %>%
  select(msk_id,seq_id)%>%
  left_join(publication_annotations%>%
              filter(ftype=="CDS")%>%
              select(seq_id,locus_tag))%>%
  left_join(publication_sequences%>%
              select(locus_tag,prot_sequence))

all_protein_sequences<- all_protein_sequences%>%
  mutate(sequence_num=paste("s",row_number(),sep=""))

aa_set <- AAStringSet(setNames(all_protein_sequences$prot_sequence,all_protein_sequences$sequence_num))

writeXStringSet(aa_set, "publication_isolate_biobanks.fasta", format="fasta")
# 
# makeblastdb -in publication_isolate_biobanks.fasta -title PubBiobankDB -out PubBiobank.db -max_file_sz 4GB -dbtype prot
# blastp -query succinate_gene_cluster.fasta -db PubBiobank.db -outfmt 6 -out Publication_succinate_hits.txt -max_target_seqs 50000


succinate_blast_results_all <- read.csv(file="Publication_succinate_hits.txt",sep="\t",header = F)
colnames(succinate_blast_results_all) <- c("qseqid","sequence_num", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# HAVE BLAST HITS FOR ENTIRE 11 GENE CLUSTER; GROUP THEM BY LOCUS_TAG AND TAKE THE HIT WITH THE HIGHEST %ID

all_operon_targets <- all_protein_sequences %>%
  select(msk_id,sequence_num,locus_tag)%>%
  left_join(succinate_blast_results_all)%>%
  filter(!is.na(qseqid))%>%
  dplyr::arrange(desc(pident))%>%
  dplyr::group_by(locus_tag)%>%
  dplyr::slice(1)%>%
  dplyr::ungroup()%>%
  select(msk_id,locus_tag,gene_label=qseqid,)%>%
  mutate(gene_label=gsub("putative\\_","",gene_label))%>%
  left_join(publication_isolates%>%
              select(msk_id,seq_id,ms_id)%>%
              left_join(publication_annotations%>%
                          select(seq_id,locus_tag)))

#  GROUP THEM INTO POTENTIAL GENE CLUSTERS.  THOSE ARE  DEFINED AS TWO HITS WITHIN 10 KB OF EACH OTHER ON A CONTIG
# TAKE THE NUMBER OF HITS IN EACH CLUSTER IN EACH ISOLATE, AND THEN TAKE ONLY THE CLUSTER THAT IS BIGGEST 

all_operon_df <- all_operon_targets %>%
  select(ms_id,seq_id)%>%
  left_join(publication_annotations)%>%
  select(ms_id,seq_id,locus_tag)%>%
  left_join(publication_locations)%>%
  select(-source,-note,-phase,-filename,-score)%>%
  filter(locus_tag %in% all_operon_targets$locus_tag)%>%
  dplyr::select(seq_id,seqid,ms_id)%>%
  unique()%>%
  dplyr::group_by(seq_id)%>%
  dplyr::mutate(contig_num=seq(1:n()))%>%
  dplyr::left_join(publication_locations %>%
                     select(-source,-note,-phase,-filename,-score)%>%
                     filter(locus_tag %in% all_operon_targets$locus_tag))%>%
  dplyr::ungroup()%>%
  dplyr::arrange(start)%>%
  dplyr::group_by(seq_id,seqid,ms_id,contig_num)%>%
  dplyr::mutate(difference=start-lag(start))%>% 
  dplyr::mutate(difference=ifelse(is.na(difference),0,difference))%>%
  dplyr::mutate(new_operon=ifelse(difference>10000,1,0))%>%
  dplyr::mutate(new_operon=ifelse(is.na(new_operon),0,new_operon))%>%
  dplyr::mutate(cumulative_new_operon=cumsum(new_operon))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id,seqid,ms_id,contig_num,cumulative_new_operon)%>%
  dplyr::summarize(earliest_start = min(start)-5000,latest_end = max(end)+5000,num_genes=n())%>%
  dplyr::ungroup()%>%
  dplyr::group_by(seq_id,ms_id)%>%
  dplyr::arrange(desc(num_genes))%>%
  dplyr::mutate(ordered_operon=seq(1:n()))%>%
  left_join(publication_locations %>%
              select(seq_id,seqid,start,end,strand,product,gene,eC_number,locus_tag))%>%
  filter(start>=earliest_start & end <=latest_end)%>%
  dplyr::mutate(start_plot_x=start-earliest_start)%>%
  dplyr::mutate(end_plot_x=end-earliest_start)%>%
  dplyr::mutate(strand=ifelse(strand=="+",TRUE,FALSE))%>%
  filter(ordered_operon==1)

all_hit_percentage <- all_protein_sequences %>%
  select(msk_id,sequence_num,locus_tag)%>%
  left_join(succinate_blast_results_all)%>%
  filter(!is.na(qseqid))%>%
  dplyr::arrange(desc(pident))%>%
  dplyr::group_by(locus_tag)%>%
  dplyr::slice(1)%>%
  dplyr::ungroup()%>%
  select(msk_id,locus_tag,gene_label=qseqid,pident)%>%
  mutate(gene_label=gsub("putative\\_","",gene_label))%>%
  filter(locus_tag %in% all_operon_df$locus_tag)%>%
  left_join(publication_isolates%>%
              select(msk_id,Genus,seq_id,Family,Phylum)%>%
              left_join(publication_annotations%>%
                          select(locus_tag,seq_id)))%>%
  ungroup()%>%
  group_by(msk_id)%>%
  summarise(average_pident=mean(pident))


publication_isolates%>%
  mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum))%>%
  group_by(Phylum,Family)%>%
  summarise(num_in_Family=n())%>%
  filter(Phylum=="Firmicutes")%>%
  arrange(desc(num_in_Family))%>%
  mutate(rank=row_number())


p1 <- all_operon_df %>%
  select(seq_id,ms_id,num_genes)%>%
  unique()%>%
  right_join(publication_isolates%>%
               select(seq_id,msk_id,key_species,Genus,Family,Phylum)%>%
               mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum)))%>%
  mutate(num_genes=ifelse(is.na(num_genes),0,num_genes))%>%
  left_join(all_hit_percentage)%>%
  group_by(num_genes,Phylum,Family)%>%
  summarise(num_genes_count=n(),pident_by_cluster_size=mean(average_pident))%>%
  left_join(publication_isolates%>%
              mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum))%>%
              group_by(Phylum,Family)%>%
              summarise(num_in_Family=n())%>%
              filter(Phylum=="Firmicutes")%>%
              arrange(desc(num_in_Family))%>%
              mutate(rank=row_number()))%>%
  mutate(percent_of_isolates=100*num_genes_count/num_in_Family)%>%
  mutate(Family_label=paste(Family,"\n(Number in Family:",num_in_Family,")",sep=" "))%>%
  filter(!is.na(Phylum))%>%
  filter(Phylum=="Firmicutes")%>%
  ungroup()%>%
  mutate(Family=as.factor(Family))%>%
  mutate(Family=reorder(Family,rank))%>%
  filter(rank<9)%>%
  ggplot(aes(x=num_genes,y=1,fill=pident_by_cluster_size,size=percent_of_isolates))+
  geom_point(pch=21)+
  # geom_bar(stat="identity",color="black")+
  facet_grid(Family~.)+
  scale_fill_gradient(low="#383838",high = "#FEA1E0",limits=c(0,100))+
  # scale_fill_viridis(option="magma")+
  labs(fill="Average % Identity\nof hits in Cluster",size="% of Isolates")+
  scale_x_continuous(breaks = seq(0,16,2))+
  ylab(label="% of Isolates From Family")+
  xlab(label="# of Genes in the Largest Cluster\nof BLAST hits on the Genome")+
  theme_bw()+
  theme(strip.background = element_blank(),
        # strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(angle=0,hjust=0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.justification = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"))

        legend.position="none")

p1

p2 <- all_operon_df %>%
  select(seq_id,ms_id,num_genes)%>%
  unique()%>%
  right_join(publication_isolates%>%
               select(seq_id,msk_id,key_species,Genus,Family,Phylum)%>%
               mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum)))%>%
  mutate(num_genes=ifelse(is.na(num_genes),0,num_genes))%>%
  left_join(all_hit_percentage)%>%
  group_by(num_genes,Phylum,Family)%>%
  summarise(num_genes_count=n(),pident_by_cluster_size=mean(average_pident))%>%
  left_join(publication_isolates%>%
              mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum))%>%
              group_by(Phylum,Family)%>%
              summarise(num_in_Family=n())%>%
              filter(Phylum=="Firmicutes")%>%
              arrange(desc(num_in_Family))%>%
              mutate(rank=row_number()))%>%
  mutate(percent_of_isolates=100*num_genes_count/num_in_Family)%>%
  mutate(Family_label=paste(Family,"\n(Number in Family:",num_in_Family,")",sep=" "))%>%
  filter(!is.na(Phylum))%>%
  filter(Phylum=="Firmicutes")%>%
  ungroup()%>%
  mutate(Family=as.factor(Family))%>%
  mutate(Family=reorder(Family,rank))%>%
  filter(rank<9)%>%
  ggplot(aes(x=1,y=1,fill=num_in_Family))+
  geom_tile(color="black")+
  # geom_bar(stat="identity",color="black")+
  facet_grid(Family~.)+
  scale_fill_gradient(low="#e3e3e3",high = "#27A3DE")+
  # scale_fill_viridis(option="magma")+
  labs(fill="Number of isolates\nin Family",size="% of Isolates")+
  scale_x_continuous(breaks = seq(0,16,2))+
  ylab(label="% of Isolates From Family")+
  xlab(label="# of Genes in the Largest Cluster\nof BLAST hits on the Genome")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.justification = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        # legend.position="none",
        panel.border = element_blank())

p2

p3 <- all_operon_df %>%
  select(seq_id,ms_id,num_genes)%>%
  unique()%>%
  right_join(publication_isolates%>%
               # filter(grepl("Lachno",Family))%>%
               select(seq_id,msk_id,key_species,Genus,Family,Phylum)%>%
               mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum)))%>%
  mutate(num_genes=ifelse(is.na(num_genes),0,num_genes))%>%
  left_join(all_hit_percentage)%>%
  group_by(num_genes,Phylum)%>%
  summarise(num_genes_count=n(),pident_by_cluster_size=mean(average_pident))%>%
  left_join(publication_isolates%>%
              mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum))%>%
              # filter(grepl("Lachno",Family))%>%
              group_by(Phylum)%>%
              summarise(num_in_Phylum=n()))%>%
  mutate(percent_of_isolates=100*num_genes_count/num_in_Phylum)%>%
  mutate(Phylum_label=paste(Phylum,"\n(Number in Phylum:",num_in_Phylum,")",sep=" "))%>%
  filter(!is.na(Phylum))%>%
  ggplot(aes(x=num_genes,y=1,fill=pident_by_cluster_size,size=percent_of_isolates))+
  geom_point(pch=21)+
  # geom_bar(stat="identity",color="black")+
  facet_grid(Phylum~.)+
  scale_fill_gradient(low="#383838",high = "#FEA1E0",limits=c(0,100))+
  # scale_fill_viridis(option="magma")+
  labs(fill="Average % ID\nOf BLAST Hits\nIn Cluster",size="% of Isolates")+
  scale_x_continuous(breaks = seq(0,16,2))+
  ylab(label="% of Isolates From Family")+
  xlab(label="# of Genes in the Largest Cluster\nof BLAST hits on the Genome")+
  theme_bw()+
  theme(strip.background = element_blank(),
        # strip.text.y = element_blank(),
        strip.text.y = element_text(angle=0,hjust=0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.justification = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        legend.position="none")
p3


p4 <- all_operon_df %>%
  select(seq_id,ms_id,num_genes)%>%
  unique()%>%
  right_join(publication_isolates%>%
               # filter(grepl("Lachno",Family))%>%
               select(seq_id,msk_id,key_species,Genus,Family,Phylum)%>%
               mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum)))%>%
  mutate(num_genes=ifelse(is.na(num_genes),0,num_genes))%>%
  left_join(all_hit_percentage)%>%
  group_by(num_genes,Phylum)%>%
  summarise(num_genes_count=n(),pident_by_cluster_size=mean(average_pident))%>%
  left_join(publication_isolates%>%
              mutate(Phylum=ifelse(grepl("Bacillota",Phylum),"Firmicutes",Phylum))%>%
              group_by(Phylum)%>%
              summarise(num_in_Phylum=n()))%>%
  mutate(percent_of_isolates=100*num_genes_count/num_in_Phylum)%>%
  mutate(Phylum_label=paste(Phylum,"\n(Number in Phylum:",num_in_Phylum,")",sep=" "))%>%
  filter(!is.na(Phylum))%>%
  ggplot(aes(x=1,y=1,fill=num_in_Phylum))+
  geom_tile(color="black")+
  # geom_bar(stat="identity",color="black")+
  facet_grid(Phylum~.)+
  scale_fill_gradient(low="#e3e3e3",high = "#27A3DE")+
  # scale_fill_viridis(option="magma")+
  labs(fill="number of Isolate\nin Phylum",size="% of Isolates")+
  scale_x_continuous(breaks = seq(0,16,2))+
  ylab(label="% of Isolates From Family")+
  xlab(label="# of Genes in the Largest Cluster\nof BLAST hits on the Genome")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.title.x = element_blank(),
        legend.justification = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        # legend.position="none",
        panel.border = element_blank())
p4

Figure_6B <- cowplot::plot_grid(p4,p3,p2,p1,align="h",axis="tb",rel_widths = c(0.5,10,.5,10),nrow = 1)



Figure_6B
# Figure 7A Metagenomic Analysis ------------------------------------------
#Load gene length information
gene_lengths <- read.csv("Query_gene_length.csv")|>
  select(!X) 

#Total reads found in each sample
total_reads <- read.csv("MS_ERR_sample_reads.csv")|>
  select(!X) |>
  rename(sample = SAMPLE)

# Read in and Process Data
fbeb_frd_all <- read.csv("Bowtie_Industrial_and_NonIndustrial_June12th_2026.csv")|>
  select(!X)  

# Set analysis options 
mapq_min <- 30

filtered_reads <- fbeb_frd_all |>
  mutate(
    MAPQ = as.numeric(MAPQ)) |>
  filter(
    !is.na(MAPQ),
    MAPQ >= mapq_min,
    str_detect(SAMPLE, "^ERR"))

counts_by_sample <- filtered_reads |>
  dplyr::count(SAMPLE, RNAME, Diet, name = "reads") |>
  dplyr::arrange(SAMPLE, desc(reads))

#Define gene order and colors 
gene_colors <- c(
  "GKGPCAPO_01438" = "#1162C2",
  "GKGPCAPO_01437" = "#1162C2",
  "GKGPCAPO_01436" = "#1162C2",
  "GKGPCAPO_01435" = "#27A3DF",
  "GKGPCAPO_01434" = "#27A3DF",
  "GKGPCAPO_01433" = "#27A3DF",
  "GKGPCAPO_01432" = "#27A3DF",
  "GKGPCAPO_01431" = "#FCAF41",
  "GKGPCAPO_01430" = "#FCAF41",
  "GKGPCAPO_01429" = "#FFA1E1",
  "GKGPCAPO_01428" = "#FFA1E1",
  "E.coli_K12_frdA" = "#529985",
  "E.coli_K12_frdB" = "#6E9F6D",
  "E.coli_K12_frdC" = "#81A665",
  "E.coli_K12_frdD" = "#ACB955",
  "B.theta frdC" = "#E7B04D",
  "B.theta frdA" = "#D08151",
  "B.theta frdB" = "#C26B51",
  "frdA prevotella" = "#ACA4E2",
  "frdB prevotella"  = "#55B8D0",
  "frdC prevotella"    = "#38BEB4",
  "Treponema succinifaciens DSM 2489, complete genome (frdA)"= "#AC7299",
  "Treponema succinifaciens DSM 2489, complete genome (frdB)" = "#D8B2C6",
  "Treponema succinifaciens DSM 2489, complete genome (frdC)" = "#C28AB1")

all_genes <- (names(gene_colors))

counts_by_sample_Industrial <- counts_by_sample |>
  filter(Diet == "Industrial")|>
  complete(
    SAMPLE,
    RNAME = all_genes,
    fill = list(reads = 0)) |>
  select(!Diet) |>
  pivot_wider(
    names_from = SAMPLE,
    values_from = reads,
    values_fill = 0) |>
  column_to_rownames("RNAME")

counts_by_sample_Non_industrial <- counts_by_sample  |>
  filter(Diet == "Non-industrial")|>
  complete(
    SAMPLE,
    RNAME = all_genes,
    fill = list(reads = 0)) |>
  select(!Diet) |>
  pivot_wider(
    names_from = SAMPLE,
    values_from = reads,
    values_fill = 0) |>
  column_to_rownames("RNAME")


#Calculate RPKM 
#RPKM standardizes read counts by gene length and sequencing depth.
counts_rpkm_Industrial <- counts_by_sample_Industrial |>
  rownames_to_column("locus_tag") |>
  left_join(gene_lengths, by = "locus_tag") |>
  pivot_longer(
    cols = !c(locus_tag, length_bp),
    names_to = "sample",
    values_to = "counts")|>
  left_join(total_reads, by = "sample") |>
  mutate(
    counts = as.numeric(counts),
    length_bp = as.numeric(length_bp),
    total_reads = as.numeric(total_reads))|>
  mutate(
    rpkm = counts / ((length_bp / 1000) * (total_reads / 1e6))) |>
  select(locus_tag, sample, rpkm) |>
  pivot_wider(names_from = sample, values_from = rpkm)

counts_rpkm_Non_industrial <- counts_by_sample_Non_industrial |>
  rownames_to_column("locus_tag") |>
  left_join(gene_lengths, by = "locus_tag") |>
  pivot_longer(
    cols = !c(locus_tag, length_bp),
    names_to = "sample",
    values_to = "counts")|>
  left_join(total_reads, by = "sample") |>
  mutate(
    counts = as.numeric(counts),
    length_bp = as.numeric(length_bp),
    total_reads = as.numeric(total_reads))|>
  mutate(
    rpkm = counts / ((length_bp / 1000) * (total_reads / 1e6))) |>
  select(locus_tag, sample, rpkm) |>
  pivot_wider(names_from = sample, values_from = rpkm)

counts_rpkm_matrix_Industrial <- counts_rpkm_Industrial |>
  column_to_rownames("locus_tag")

counts_rpkm_matrix_Non_industrial <- counts_rpkm_Non_industrial |>
  column_to_rownames("locus_tag")

gene_order <- names(gene_colors)

sample_order_I <- colnames(counts_rpkm_matrix_Industrial) |>
  str_remove("ERR") |>
  as.numeric() |>
  order()
sample_order_NI <- colnames(counts_rpkm_matrix_Non_industrial) |>
  str_remove("ERR") |>
  as.numeric() |>
  order()


counts_rpkm_matrix_Industrial <- counts_rpkm_matrix_Industrial[gene_order, sample_order_I]
counts_rpkm_matrix_Non_industrial <- counts_rpkm_matrix_Non_industrial[gene_order, sample_order_NI]

#
mean_rpkm_long_Industrial <- counts_rpkm_matrix_Industrial |>
  as.data.frame() |>
  rownames_to_column("locus_tag") |>
  pivot_longer(
    cols = -locus_tag,
    names_to = "sample",
    values_to = "rpkm" ) |>
  mutate(locus_tag = factor(locus_tag, levels = gene_order))|>
  mutate(Diet = "Industrial")

mean_rpkm_long_Non_industrial <- counts_rpkm_matrix_Non_industrial |>
  as.data.frame() |>
  rownames_to_column("locus_tag") |>
  pivot_longer(
    cols = -locus_tag,
    names_to = "sample",
    values_to = "rpkm" ) |>
  mutate(locus_tag = factor(locus_tag, levels = gene_order))|>
  mutate(Diet = "Non_industrial")

mean_rpkm_long_Total <- mean_rpkm_long_Industrial |>
  rbind(mean_rpkm_long_Non_industrial)|>
  mutate(
    locus_tag = factor(
      locus_tag,
      levels = gene_order))




# Combine Western and Hazda into a single heatmap
# Log-transform
mat_log_I <- log2(counts_rpkm_matrix_Industrial + 1)
mat_log_NI <- log2(counts_rpkm_matrix_Non_industrial + 1)

# Preserve your existing sample ordering
mat_log_I <- mat_log_I[, sample_order_I]
mat_log_NI <- mat_log_NI[, sample_order_NI]

# Combine matrices
mat_log_combined <- cbind(
  mat_log_I,
  mat_log_NI)

# Single color scale across both cohorts
breaks_fixed <- seq(
  min(mat_log_combined, na.rm = TRUE),
  max(mat_log_combined, na.rm = TRUE),
  length.out = 101
)

# Row annotation (genes)
gene_groups <- data.frame(
  gene = rownames(mat_log_combined))

rownames(gene_groups) <- gene_groups$gene

# Use whichever color vector contains ALL genes
gene_colors <- gene_colors

# Column annotation (diet)

sample_annotation <- data.frame(
  Diet = c(
    rep("Industrial",
        ncol(mat_log_I)),
    rep("Non_industrial",
        ncol(mat_log_NI))))

rownames(sample_annotation) <- colnames(mat_log_combined)

# Annotation colors

annotation_colors <- list(
  gene = gene_colors,
  Diet = c(
    Industrial = "#3B6FB6",
    Non_industrial = "#D95F02"))

mat_log_combined <- as.matrix(mat_log_combined)

storage.mode(mat_log_combined) <- "numeric"

# Heatmap
heatmap_combined <- pheatmap(mat_log_combined,
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             annotation_row = gene_groups,
                             annotation_col = sample_annotation,
                             annotation_colors = annotation_colors,
                             border_color = NA,
                             annotation_legend = TRUE,
                             color = colorRampPalette(c(
                               "#26456E",
                               "#4A6FE3",
                               "white",
                               "#D33F6A",
                               "#9C0824"))(100),
                             breaks = breaks_fixed,
                             main = "",
                             silent = TRUE)$gtable

Figure_7A <- ggplotify::as.ggplot(heatmap_combined)

Figure_7A

ggsave(
  filename = "Heatmap_June12th_2026_Marissa.pdf",
  plot = Figure_7A,
  width = 20,
  height = 10,
  units = "in",
  bg = "white"
)

# Figure 7B Correlation Analysis --------------------------------------------------
#Convert long data to sample × gene matrix
mat <- mean_rpkm_long_Total |>
  select(sample, locus_tag, rpkm) |>
  pivot_wider(
    id_cols = sample,
    names_from = locus_tag,
    values_from = rpkm,
    values_fill = 0) |>
  column_to_rownames("sample") |>
  as.matrix()

# Remove genes with zero variance
mat_clean <- mat[,apply(mat, 2, function(x) {sd(x, na.rm = TRUE) > 0})]

# Correlation matrices
cor_mat <- cor(
  mat_clean,
  method = "spearman",
  use = "pairwise.complete.obs")

cor_pmat <- ggcorrplot::cor_pmat(
  mat_clean,
  method = "spearman")

# Desired plotting order

gene_order <- c(
  "GKGPCAPO_01438",
  "GKGPCAPO_01437",
  "GKGPCAPO_01436",
  "GKGPCAPO_01435",
  "GKGPCAPO_01434",
  "GKGPCAPO_01433",
  "GKGPCAPO_01432",
  "GKGPCAPO_01431",
  "GKGPCAPO_01430",
  "GKGPCAPO_01429",
  "GKGPCAPO_01428",
  "E.coli_K12_frdA",
  "E.coli_K12_frdB",
  "E.coli_K12_frdC",
  "E.coli_K12_frdD",
  "B.theta frdC",
  "B.theta frdA",
  "B.theta frdB",
  "frdA prevotella",
  "frdB prevotella",
  "frdC prevotella",
  "Treponema succinifaciens DSM 2489, complete genome (frdA)",
  "Treponema succinifaciens DSM 2489, complete genome (frdB)",
  "Treponema succinifaciens DSM 2489, complete genome (frdC)")

# Keep only genes that survived filtering

gene_order <- gene_order[
  gene_order %in% rownames(cor_mat)]

#debugging
missing_genes <- setdiff(
  c(
    "GKGPCAPO_01438",
    "GKGPCAPO_01437",
    "GKGPCAPO_01436",
    "GKGPCAPO_01435",
    "GKGPCAPO_01434",
    "GKGPCAPO_01433",
    "GKGPCAPO_01432",
    "GKGPCAPO_01431",
    "GKGPCAPO_01430",
    "GKGPCAPO_01429",
    "GKGPCAPO_01428",
    "E.coli_K12_frdA",
    "E.coli_K12_frdB",
    "E.coli_K12_frdC",
    "E.coli_K12_frdD",
    "B.theta frdC",
    "B.theta frdA",
    "B.theta frdB",
    "frdA prevotella",
    "frdB prevotella",
    "frdC prevotella",
    "Treponema succinifaciens DSM 2489, complete genome (frdA)",
    "Treponema succinifaciens DSM 2489, complete genome (frdB)",
    "Treponema succinifaciens DSM 2489, complete genome (frdC)"),
  
  rownames(cor_mat))
print(missing_genes)

#Reorder matrices

cor_mat <- cor_mat[
  gene_order,
  gene_order]

cor_pmat <- cor_pmat[
  gene_order,
  gene_order]

#Plot
Figure_7B <- as.ggplot(~corrplot::corrplot(
  cor_mat,
  method = "circle",
  type = "lower",
  order = "original",
  col = colorRampPalette(
    c("#D33F6A", "white", "#4A6FE3"))(100),
  tl.col = "black",
  tl.srt = 90,
  p.mat = cor_pmat,
  sig.level = c(0.001, 0.01, 0.05),
  insig = "label_sig",
  pch.cex = 0.7,
  pch.col = "black"))

Figure_7B

ggsave("Correlation_plot_June12th_2026.pdf", 
       Figure_7B,
       width = 14,
       height = 14)


# Figure 7C ---------------------------------------------------------------

t <- read.csv(file="Mean_rpkm_June_12th_2026 2.csv",header=TRUE,sep=",")%>%
  select(-X)%>%
  mutate(group=ifelse(grepl("GKGPCAP",locus_tag),"1. FBEB-FR",
                      ifelse(grepl("prevotella",locus_tag),"4. Prev.",
                             ifelse(grepl("theta",locus_tag),"2. Bact.",
                                     ifelse(grepl("coli",locus_tag),"3. E coli","5. Trep.")))))%>%
  dplyr::mutate(Diet=ifelse(Diet=="Industrial","2. Industrial","1. Non-industrial"))

library(ggh4x)
library(ggnewscale)
t %>%
  dplyr::group_by(sample,Diet,group)%>%
  dplyr::summarise(mean_rpkm=mean(rpkm))%>%
  write.csv(file="mean_rpkm.csv")

p1 <- t %>%

  dplyr::group_by(sample,Diet,group)%>%
  dplyr::summarise(mean_rpkm=mean(rpkm))%>%
  dplyr::ungroup()%>%
  filter(group=="1. FBEB-FR")%>%
    ggplot(aes(x=Diet,y=mean_rpkm,fill=Diet))+
    facet_grid(~group,scale="free_x",space = "free_x")+
    geom_boxplot(position = "dodge",outlier.shape = NA)+
    scale_fill_manual(values=c("#D95F02","#3B6FB6"))+
    geom_jitter(aes(x=Diet,y=mean_rpkm),size=0.5,position = position_jitterdodge(jitter.width = 0.5))+
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 2),
      breaks = 2^(0:10))+
    scale_x_discrete(expand = c(0,0))+
    theme_bw()+
    theme(legend.position = "none",
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  p1
p2 <- t %>%
  dplyr::group_by(sample,Diet,group)%>%
  dplyr::summarise(mean_rpkm=mean(rpkm))%>%
  dplyr::ungroup()%>%
  filter(group=="2. Bact.")%>%
  ggplot(aes(x=Diet,y=mean_rpkm,fill=Diet))+
  facet_grid(~group,scale="free_x",space = "free_x")+
  geom_boxplot(position = "dodge",outlier.shape = NA)+
  scale_fill_manual(values=c("#D95F02","#3B6FB6"))+
  geom_jitter(aes(x=Diet,y=mean_rpkm),size=0.5,position = position_jitterdodge(jitter.width = 0.5))+
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 2),
    breaks = 2^(0:10))+
  scale_x_discrete(expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

p3 <- t %>%
  dplyr::group_by(sample,Diet,group)%>%
  dplyr::summarise(mean_rpkm=mean(rpkm))%>%
  dplyr::ungroup()%>%
  filter(group=="3. E coli")%>%
  ggplot(aes(x=Diet,y=mean_rpkm,fill=Diet))+
  facet_grid(~group,scale="free_x",space = "free_x")+
  geom_boxplot(position = "dodge",outlier.shape = NA)+
  scale_fill_manual(values=c("#D95F02","#3B6FB6"))+
  geom_jitter(aes(x=Diet,y=mean_rpkm),size=0.5,position = position_jitterdodge(jitter.width = 0.5))+
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 2),
    breaks = 2^(0:10))+
  scale_x_discrete(expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
p3
p4 <- t %>%
  dplyr::group_by(sample,Diet,group)%>%
  dplyr::summarise(mean_rpkm=mean(rpkm))%>%
  dplyr::ungroup()%>%
  filter(group=="4. Prev.")%>%
  ggplot(aes(x=Diet,y=mean_rpkm,fill=Diet))+
  facet_grid(~group,scale="free_x",space = "free_x")+
  geom_boxplot(position = "dodge",outlier.shape = NA)+
  scale_fill_manual(values=c("#D95F02","#3B6FB6"))+
  geom_jitter(aes(x=Diet,y=mean_rpkm),size=0.5,position = position_jitterdodge(jitter.width = 0.5))+
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 2),
    breaks = 2^(0:10))+
  scale_x_discrete(expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
p4
p5 <- t %>%
  dplyr::group_by(sample,Diet,group)%>%
  dplyr::summarise(mean_rpkm=mean(rpkm))%>%
  dplyr::ungroup()%>%
  filter(group=="5. Trep.")%>%
  
  ggplot(aes(x=Diet,y=mean_rpkm,fill=Diet))+
  facet_grid(~group,scale="free_x",space = "free_x")+
  geom_boxplot(position = "dodge",outlier.shape = NA)+
  scale_fill_manual(values=c("#D95F02","#3B6FB6"))+
  geom_jitter(aes(x=Diet,y=mean_rpkm),size=0.5,position = position_jitterdodge(jitter.width = 0.5))+
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 2),
    breaks = 2^(0:10))+
  scale_x_discrete(expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
p5

cowplot::plot_grid(p1,p2,p3,p4,p5,ncol = 5,align="h",axis="tb")

