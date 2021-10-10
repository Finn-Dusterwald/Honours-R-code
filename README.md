
##################################################################################################################################

#SpoT Over-expression - RNA-Seq Analysis
#Differentially expressed genes after spoT overexpression.

##################################################################################################################################

#Remove all objects from workspace/environment.
rm(list = ls())


#Set working directory.
setwd("~/Documents/University Work/Research Project/R Code ")


#Load packages required for analysis.
library(tidyverse)
library(EnhancedVolcano)
library(GOfuncR)
library(ggplot2)
library(ggpubr)
library(gage)
library(pathview)
library(DESeq2)


#Read in dataset containing RNA-seq data on all differentially expressed genes in SpoT-overexpression .
Unfiltered_spoT <- read.csv("de_unfiltered_spoT_INvsUN_20210304.csv", 
                       header=TRUE)


#Filter list to remove genes with an FC (Fold change) ± 1.5 and an adjusted p-value equal to or less than 0.05.
Filtered_spoT <- filter(Unfiltered_spoT, 
                   FC >= 1.5 & padj <= 0.05 | FC <= -1.5 & padj <= 0.05)


#Split the upregulated and downregulated genes which are significant
UpregulatedGenes <- filter(Filtered_spoT,
                           FC >= 1.5)

DownregulatedGenes <- filter(Filtered_spoT,
                             FC <= -1.5)


#Read in Pseudomonas aeruginosa annotation file - Downloaded from www.psuedomonas.com.
PA_Annotation <- read.csv(file = 'Pseudomonas_aeruginosa_PAO1_107.tsv', 
                       header = TRUE, sep = '\t', 
                       fill = TRUE)


#Merge the filtered gene list with the annotation file based on the locus tag.
MergedData <- left_join(Filtered_spoT, PA_Annotation, by="Locus.Tag")


#Return the number of genes that are either down-regulated or up-regulated in the filtered dataset. Based on a fold change of 1.5.
sum(MergedData$FC < -1.5)
sum(MergedData$FC > 1.5)



##################################################################################################################################

#Volcano plot 

##################################################################################################################################

#Remove SpoT from dataset in order to better represent data on volcano plot.
Vplot <- Unfiltered_spoT[Unfiltered_spoT$name != "spoT", ]

#Color genes based on fold change(log2) in expression.
keyvals <- ifelse(
  Vplot$padj > 0.05, 'grey',
  ifelse(Vplot$log2FoldChange < (-0.58), 'red3',
         ifelse(Vplot$log2FoldChange > 0.58, 'green4',
                'grey')))

#Set colors to correspond to specific labels
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'green4'] <- 'Upregulated'
names(keyvals)[keyvals == 'grey'] <- 'Not Significant'
names(keyvals)[keyvals == 'red3'] <- 'Downregulated'

#Volcano Plot 
EnhancedVolcano(Vplot,
                lab = as.character(Vplot$name),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Dysregulated Genes in spoT Overexpression',
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 4.0,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'bottom',
                legendLabSize = 15)



##################################################################################################################################

#Gene Ontology Enrichment Analysis

##################################################################################################################################

#Read in GO annotation file for Pseudomonas aeruginosa PAO1 - Downloaded from www.pseudomonas.com.
GO_Annotation <- read.csv("gene_ontology_csv.csv", 
                       header=TRUE)


#Create custom annotation package for GO enrichment and merge GO annotation with Locus Tag ID.
custom_GO <- data.frame(Locus.Tag = GO_Annotation$Locus.Tag,
                        Accession = GO_Annotation$Accession)

#Count the number of each GO category in the PAO1 annotation file for gene ratio calculation
Annotation_Count <- rle(sort(GO_Annotation$GO.Term))


#Create a data frame containing the PAO1 annotation GO categories and counts of each category.
PAO1_annotation_count <- data.frame(GO.category = Annotation_Count$values,
                                    Count = Annotation_Count$lengths)



#############################################################################################

#GO analysis for up regulated genes

#############################################################################################
#Create Gene List from filtered dataset of upregulated genes.
gene_ids_up = c(UpregulatedGenes$Locus.Tag)


#Mandatory input for Go_enrich. Removes NA gene_IDs in the dataset.
input_hyper_up = data.frame(na.omit(gene_ids_up), 
                         is_candidate=1)

#Run GO enrichment using the upregulated genes and the custom GO annotation package
res_hyper_up = go_enrich(input_hyper_up,
                         annotations = custom_GO)


#Annotate all given genes to GO categories.
anno_up = get_anno_categories(c(Locus.Tag = UpregulatedGenes$Locus.Tag), 
                           annotations = custom_GO)


#Create a data frame containing the key GO enrichment results.
GOgraph_up <- data.frame(GO.category = res_hyper_up$results$node_name,
                      GO.number = res_hyper_up$results$node_id,
                      p.value = res_hyper_up$results$FWER_overrep,
                      ontology = res_hyper_up$results$ontology)


#Aggregate the number of times each GO category appears in the anno_up data frame.  
GO_count_up <- rle(sort(anno_up$name))


#Create a data frame containing the count for each GO category in the upregulated enriched data set.
Count_list_up <- data.frame(GO.category = GO_count_up$values,
                         Count = GO_count_up$lengths)


#Merge the counts of enriched GO categories with the total counts in the annotation file.
Combined_count_up <- data.frame(left_join(Count_list_up, PAO1_annotation_count, by = "GO.category"))


#Create a data frame containing each of the GO enriched genes, the count of each GO category and the gene ratio.
GO_count_up2 <- data.frame(GO.category = GO_count_up$values,
                        Count = GO_count_up$lengths,
                        GeneRatio = (Combined_count_up$Count.x / Combined_count_up$Count.y))


#Merge Go_count2 data frame with GOgraph
Dotplot_up <- data.frame(left_join(GO_count_up2, GOgraph_up, by = "GO.category"))


#Create a column that combines the GO number with the name.
Dotplot_up$merged <- paste(Dotplot_up$GO.category,Dotplot_up$GO.number, sep = '\n' )


###########################################################################

#Unregulated biological processes

###########################################################################

#Subset upregulated gene ontology categories into biological processes only
BioProGO_up <- subset(Dotplot_up, ontology == "biological_process")

#Filter based on a adjusted p-value of 0.1 and 5 more genes annotated
BioProGOfilter_up <- filter(BioProGO_up, Count >= 5 & p.value <= 0.05)

#Graph upregulated biological processes as dotplot and gene ratios
DotUpBP <- ggplot(BioProGOfilter_up, aes(x = GeneRatio, 
                                         y = fct_reorder(merged, GeneRatio), 
                                          color = ontology, size = Count)) +
  geom_point(color = "green4") +
  theme(axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 19),
        legend.title = element_text(color = "black", size = 24),
        legend.text = element_text(color = "black", size = 19),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'gray94', linetype = 'solid')) +
  xlab("Gene Ratio") 


###########################################################################

#Upregulated Molecular functions

###########################################################################

#Subset upregulated gene ontology categories into molecular function only
Molecular_functionGO_up <- subset(Dotplot_up, ontology == "molecular_function")

#Filter based on a adjusted p-value of 0.1 and 5 more genes annotated
Molecular_functionGOfilter_up <- filter(Molecular_functionGO_up, Count >= 5 & p.value <= 0.05)

#Graph upregulated molecular functions as dotplot with gene ratios
DotUpMF <- ggplot(Molecular_functionGOfilter_up, aes(x = GeneRatio, y = fct_reorder(merged, GeneRatio), 
                                                   color = ontology, size = Count)) +
  geom_point(color = "green4") +
  theme(axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 19),
        legend.title = element_text(color = "black", size = 24),
        legend.text = element_text(color = "black", size = 19),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'gray94', linetype = 'solid')) +
  xlab("Gene Ratio") 


###########################################################################

#Upregulated Cellular Components

###########################################################################

#Subset upregulated gene ontology categories into cellular component only
Cellular_componentGO_up <- subset(Dotplot_up, ontology == "cellular_component")

#Filter based on a adjusted p-value of 0.1 and 5 more genes annotated
Cellular_componentGOfilter_up <- filter(Cellular_componentGO_up, Count >= 5 & p.value <= 0.05)

#Graph upregulated cellular components as dotplot with gene ratios
DotUpCC <- ggplot(Cellular_componentGOfilter_up, aes(x = GeneRatio, y = fct_reorder(merged, GeneRatio), 
                                                     color = ontology, size = Count)) +
  geom_point(color = "green4") +
  theme(axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 19),
        legend.title = element_text(color = "black", size = 24),
        legend.text = element_text(color = "black", size = 19),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'gray94', linetype = 'solid')) +
  xlab("Gene Ratio") 




#############################################################################################

#GO analysis for down regulated genes

#############################################################################################
#Create Gene List from filtered dataset of upregulated genes.
gene_ids_down = c(DownregulatedGenes$Locus.Tag)


#Mandatory input for Go_enrich. Removes NA gene_IDs in the dataset.
input_hyper_down = data.frame(na.omit(gene_ids_down), 
                            is_candidate=1)

#Run GO enrichment using the downregulated genes and the custom GO annotation package
res_hyper_down = go_enrich(input_hyper_down,
                         annotations = custom_GO)


#Annotate all given genes to GO categories.
anno_down = get_anno_categories(c(Locus.Tag = DownregulatedGenes$Locus.Tag), 
                              annotations = custom_GO)


#Create a data frame containing the key GO enrichment results.
GOgraph_down <- data.frame(GO.category = res_hyper_down$results$node_name,
                         GO.number = res_hyper_down$results$node_id,
                         p.value = res_hyper_down$results$FWER_overrep,
                         ontology = res_hyper_down$results$ontology)


#Aggregate the number of times each GO category appears in the anno_up data frame.  
GO_count_down <- rle(sort(anno_down$name))


#Create a data frame containing the count for each GO category in the upregulated enriched data set.
Count_list_down <- data.frame(GO.category = GO_count_down$values,
                            Count = GO_count_down$lengths)


#Merge the counts of enriched GO categories with the total counts in the annotation file.
Combined_count_down <- data.frame(left_join(Count_list_down, PAO1_annotation_count, by = "GO.category"))


#Create a data frame containing each of the GO enriched genes, the count of each GO category and the gene ratio.
GO_count_down2 <- data.frame(GO.category = GO_count_down$values,
                           Count = GO_count_down$lengths,
                           GeneRatio = (Combined_count_down$Count.x / Combined_count_down$Count.y))


#Merge Go_count2 data frame with GOgraph
Dotplot_down <- data.frame(left_join(GO_count_down2, GOgraph_down, by = "GO.category"))


#Create a column that combines the GO number with the name.
Dotplot_down$merged <- paste(Dotplot_down$GO.category,Dotplot_down$GO.number, sep = "\n" )


###########################################################################

#Downregulated Biological Processes

###########################################################################

#Subset downregulated gene ontology categories into biological process only
Biological_processGO_down <- subset(Dotplot_down, ontology == "biological_process")

#Filter based on a adjusted p-value of 0.1 and 5 more genes annotated
Biological_processGOfilter_down <- filter(Biological_processGO_down, Count >= 5 & p.value <= 0.05)

#Graph downregulated biological processes with a dotplot and gene ratios
DotDownBP <- ggplot(Biological_processGOfilter_down, aes(x = GeneRatio, y = fct_reorder(merged, GeneRatio), 
                                          color = ontology, size = Count)) +
  geom_point(color = "red3") +
  theme(axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 19),
        legend.title = element_text(color = "black", size = 24),
        legend.text = element_text(color = "black", size = 19),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'gray94', linetype = 'solid')) +
  xlab("Gene Ratio") 

###########################################################################

#Downregulated Molecular function

###########################################################################

#Subset downregulated gene ontology categories into molecular function only
Molecular_functionGO_down <- subset(Dotplot_down, ontology == "molecular_function")

#Filter based on a adjusted p-value of 0.1 and 5 more genes annotated
Molecular_functionGOfilter_down <- filter(Molecular_functionGO_down, Count >= 5 & p.value <= 0.05)

#Graph downregulated molecular functions with dotplot with gene ratios
DotDownMF <- ggplot(MFF_down, aes(x = GeneRatio, y = fct_reorder(merged, GeneRatio), 
                                                         color = ontology, size = Count)) +
  geom_point(color = "red3") +
  theme(axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 19),
        legend.title = element_text(color = "black", size = 24),
        legend.text = element_text(color = "black", size = 19),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'gray94', linetype = 'solid')) +
  xlab("Gene Ratio")


###########################################################################

#Downregulated Cellular Component

###########################################################################

#Subset downregulated gene ontology categories into cellular component only
Cellular_componentGO_down <- subset(Dotplot_down, ontology == "cellular_component")

#Filter based on a adjusted p-value of 0.1 and 5 more genes annotated
Cellular_componentGOfilter_down <- filter(Cellular_componentGO_down, Count >= 5 & p.value <= 0.05)

#Graph biological processes with dotplot with gene ratios
DotDownCC <- ggplot(CCF_down, aes(x = GeneRatio, y = fct_reorder(merged, GeneRatio), 
                                                         color = ontology, size = Count)) +
  geom_point(color = "red3") +
  theme(axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        axis.title = element_text(color = "black", size = 19),
        legend.title = element_text(color = "black", size = 24),
        legend.text = element_text(color = "black", size = 19),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = 'gray94', linetype = 'solid')) +
  xlab("Gene Ratio")



######################################################################################################

#Biological processes up and down regulated merged into one figure
BPplot <- ggarrange(DotUpBP, DotDownBP,
                    ncol = 1, nrow = 2,
                    labels = c("A", "B"))

#Molecular function up and down regulated merged into one figure
MFplot <- ggarrange(DotUpMF, DotDownMF,
                    ncol = 1, nrow = 2,
                    labels = c("A", "B"))

#Cellular component up and down regulated merged into one figure
CCplot <- ggarrange (DotUpCC, DotDownCC,
                     ncol = 1, nrow = 2,
                     labels = c("A", "B"))


Mf##################################################################################################################################

#KEGG Pathway Analysis & pathview

##################################################################################################################################

#Create a value containing the fold change of the DE genes
deseq2.fc = Filtered_spoT$log2FoldChange


#Add locus tags to the fold change values in deseq2.fc
names(deseq2.fc) = Filtered_spoT$Locus.Tag


#Change KEGG species to Pseudomonas aeruginosa
kg.species <- kegg.gsets("pae")

kg.species$allpaths.idx <- 1:(length(kg.species$kg.sets))

kg.gs=kg.species$kg.sets[kg.species$allpaths.idx]


#Run the KEGG pathway enrichment analysis
fc.kegg.p <- gage(deseq2.fc, gsets = kg.gs, ref = NULL, samp = NULL)


#Return the number of significantly up- and down-regulated gene sets
fc.kegg.sig <- sigGeneSet(fc.kegg.p, outname="siggenes", cutoff=0.1, heatmap=TRUE)


#Return the top KEGG pathways which are up- and down-regulated 
lapply(fc.kegg.p, head)


#Create data frame containing up-regulated pathways from the KEGG enrichment
greater <- as.data.frame(fc.kegg.sig$greater)
setDT(greater, keep.rownames = T)
colnames(greater)[1] <- "Pathway"
greater <- lapply(greater, as.character)


#Create data frame containing downregulated pathways from the KEGG enrichment
less <- as.data.frame(fc.kegg.sig$less)
setDT(less, keep.rownames = T)
colnames(less)[1] <- "Pathway"
less <- lapply(less, as.character)


#Create a list of the up-regulated pathways
greater1 = substr(greater$Pathway, start=1, stop=8)
greater1


#Create a list of the down-regulated pathways
less1 = substr(less$Pathway, start = 1, stop = 8)
less1


#Run pathview on Greater pathways
tmp <- pathview(gene.data = deseq2.fc, pathway.id = greater1, species = "pae", gene.idtype = "kegg", limit = (gene = 2))


#Run pathview on Lesser pathways
tmp <- pathview(gene.data = deseq2.fc, pathway.id = less1, species = "pae", gene.idtype = "kegg", limit = (gene = 2))


