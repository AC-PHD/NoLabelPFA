# this will create subsets according to condition
# remove ribosomal and mitochondrial genes
# create PFA input tables which are labeled with A) the condition name and B) the a number for the condition (see # NAME variables)
# if you do not intend to use numbers as labels you can ignore ## Create Table using numbers as labels ---- 

# libraries
library(Seurat)
library(tidyverse)

# create directories
dir.create("Documentation_Plots")
dir.create("Output")
dir.create("PFA_Input")

# NAME variables
# e.g., for condition "age" with "young" and "old"
condition1_name <- c("Young")
condition2_name <- c("Old")
Cellinfo <- paste0(condition1_name,"_and_",condition2_name)
as_Number_condition1_name <- 0 # number for creating table -> 0_Young = 0 and 0_Old -> 1 for direct comparison  
as_Number_condition2_name <- 1

# # import data -----
# # .RDS format 
# we use the seurat object by Sole-Boldo et al. (2020) https://www.nature.com/articles/s42003-020-0922-4 
# the RDS-file is available via the GEO database: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130973 
SeuratObject <- readRDS("GSE130973_seurat_analysis_lyko.rds") # SeuratObject <- readRDS("path/to/inputfile")
str(SeuratObject)

# have a look at the metadata
View(SeuratObject@meta.data)
# see which "options" (cell types, conditions, etc.) are available (here "age")
colnames(SeuratObject@meta.data)

# select your criteria and your condition -> e.g. age in meta.data (YOUNG and OLD)
# and create your subsets (one subset for each condition)
criterion <- c("age") # row name of column of interest in Seurat object (e.g., interested in age, columname is age)
condition_1 <-  subset(SeuratObject, subset = age == "YOUNG")  # first option in column of interest in Seurat object (e.g., "YOUNG")
condition_2 <-  subset(SeuratObject, subset = age == "OLD")  # second option in column of interest in Seurat object (e.g., "OLD")

# if you want to plot UMAP with another column (e.g., cell types, or here cell types by age, with column name celltype.age) 
additional_criterion <- c("celltype.age")

# # visualize data -----
# all cells by criterion (e.g. age)
(dimPlot_start_criterion <- DimPlot(SeuratObject, 
                                    reduction = 'umap', group.by = criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",Cellinfo,"_",criterion,".png")
saveNamePDF <- paste0("Documentation_Plots/",Cellinfo,"_",criterion,".pdf")
png(file=saveNamePNG)
dimPlot_start_criterion
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_criterion
dev.off()

# cells by additional criterion (e.g. cell type and age)
(dimPlot_start_criterion2 <- DimPlot(SeuratObject, 
                                reduction = 'umap', group.by = additional_criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",Cellinfo,"_",additional_criterion,".png")
saveNamePDF <- paste0("Documentation_Plots/",Cellinfo,"_",additional_criterion,".pdf")
png(file=saveNamePNG)
dimPlot_start_criterion2
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_criterion2
dev.off()

# plots for condition 1 of your comparison (here comparison by age, all young cells)
(dimPlot_start_1_criterion <- DimPlot(condition_1, 
                                      reduction = 'umap', group.by = criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",condition1_name,"_",criterion,".png")
saveNamePDF <- paste0("Documentation_Plots/",condition1_name,"_",criterion,".pdf")
png(file=saveNamePNG)
dimPlot_start_1_criterion
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_1_criterion
dev.off()

# plots for condition 1 of your comparison using the additional criterion
(dimPlot_start_1_criterion2 <- DimPlot(condition_1, 
                                  reduction = 'umap', group.by = additional_criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",condition1_name,"_",additional_criterion,".png")
saveNamePDF <- paste0("Documentation_Plots/",condition1_name,"_",additional_criterion,".pdf")
png(file=saveNamePNG)
dimPlot_start_1_criterion2
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_1_criterion2
dev.off()

# plots for condition 2 of your comparison (here comparison by age, all old cells)
(dimPlot_start_2_criterion <- DimPlot(condition_2, 
                                      reduction = 'umap', group.by = criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",condition2_name,"_",criterion,".png")
saveNamePDF <- paste0("Documentation_Plots/",condition2_name,"_",criterion,".pdf")
png(file=saveNamePNG)
dimPlot_start_2_criterion
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_2_criterion
dev.off()

# plots for condition 2 of your comparison using the additional criterion
(dimPlot_start_2_criterion2 <- DimPlot(condition_2, 
                                  reduction = 'umap', group.by = additional_criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",condition2_name,"_",additional_criterion,".png")
saveNamePDF <- paste0("Documentation_Plots/",condition2_name,"_",additional_criterion,".pdf")
png(file=saveNamePNG)
dimPlot_start_2_criterion2
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_2_criterion2
dev.off()

## remove ribosomal and mitochondrial genes
# if you would like to save the resulting subsets that are ready for preparing the PFA input table as RDS objects, uncomment # saveRDS(...)

# Get the names of all mitochondrial and ribosomal genes in condition_1
genes_to_remove_condition_1 <- c(grep("^RP[LS]", rownames(condition_1@assays$RNA@counts), value = TRUE),
                                 grep("^MT-", rownames(condition_1@assays$RNA@counts), value = TRUE))
# Remove mitochondrial and ribosomal genes from condition_1
condition_1 <- subset(condition_1, features = setdiff(rownames(condition_1@assays$RNA@counts), genes_to_remove_condition_1))
# Get the names of all mitochondrial and ribosomal genes in condition_1
genes_to_remove_condition_1_after <- c(grep("^RP[LS]", rownames(condition_1@assays$RNA@counts), value = TRUE),
                                       grep("^MT-", rownames(condition_1@assays$RNA@counts), value = TRUE))

# # optional: save as RDS file
# saveNameRDS <- paste0("Output/subset_",condition1_name,".rds")
# saveRDS(condition_1, file = saveNameRDS)

# Get the names of all mitochondrial and ribosomal genes in condition_2
genes_to_remove_condition_2 <- c(grep("^RP[LS]", rownames(condition_2@assays$RNA@counts), value = TRUE),
                                 grep("^MT-", rownames(condition_2@assays$RNA@counts), value = TRUE))
# Remove mitochondrial and ribosomal genes from condition_2
condition_2 <- subset(condition_2, features = setdiff(rownames(condition_2@assays$RNA@counts), genes_to_remove_condition_2))
# Get the names of all mitochondrial and ribosomal genes in condition_2
genes_to_remove_condition_2_after <- c(grep("^RP[LS]", rownames(condition_2@assays$RNA@counts), value = TRUE),
                                       grep("^MT-", rownames(condition_2@assays$RNA@counts), value = TRUE))

# # optional: save as RDS file
# saveNameRDS <- paste0("Output/subset_",condition2_name,".rds")
# saveRDS(condition_2, file = saveNameRDS)

#plots after removing ribosomal and mitochondrial genes
# plots for condition 1 of your comparison (here comparison by age, all young cells)
(dimPlot_start_1_criterion_after <- DimPlot(condition_1, 
                                      reduction = 'umap', group.by = criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",condition1_name,"_",criterion,"_after.png")
saveNamePDF <- paste0("Documentation_Plots/",condition1_name,"_",criterion,"_after.pdf")
png(file=saveNamePNG)
dimPlot_start_1_criterion_after
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_1_criterion_after
dev.off()

# plots for condition 1 of your comparison using the additional criterion
(dimPlot_start_1_criterion2_after <- DimPlot(condition_1, 
                                       reduction = 'umap', group.by = additional_criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",condition1_name,"_",additional_criterion,"_after.png")
saveNamePDF <- paste0("Documentation_Plots/",condition1_name,"_",additional_criterion,"_after.pdf")
png(file=saveNamePNG)
dimPlot_start_1_criterion2_after
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_1_criterion2_after
dev.off()

# plots for condition 2 of your comparison (here comparison by age, all old cells)
(dimPlot_start_2_criterion_after <- DimPlot(condition_2, 
                                      reduction = 'umap', group.by = criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",condition2_name,"_",criterion,"_after.png")
saveNamePDF <- paste0("Documentation_Plots/",condition2_name,"_",criterion,"_after.pdf")
png(file=saveNamePNG)
dimPlot_start_2_criterion_after
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_2_criterion_after
dev.off()

# plots for condition 2 of your comparison using the additional criterion
(dimPlot_start_2_criterion2_after <- DimPlot(condition_2, 
                                       reduction = 'umap', group.by = additional_criterion))
# save
saveNamePNG <- paste0("Documentation_Plots/",condition2_name,"_",additional_criterion,"_after.png")
saveNamePDF <- paste0("Documentation_Plots/",condition2_name,"_",additional_criterion,"_after.pdf")
png(file=saveNamePNG)
dimPlot_start_2_criterion2_after
dev.off() 
pdf(file=saveNamePDF, width = 11.69, height = 8.27,  paper = "a4r" )
dimPlot_start_2_criterion2_after
dev.off()

# Create Table ----
# will create tables for each condition and a merged table using A) the name as label, B) the number as label
# if you prefer to use names as labels (recommended) you can ignore ## Create Table using numbers as labels ---- 

##condition_1
# get the counts of the raw RNA counts
Counts_DF_condition_1 <- as.data.frame(condition_1@assays$RNA@counts) %>%  rownames_to_column()
Counts_DF_names_condition_1 <- colnames(Counts_DF_condition_1)
# View(Counts_DF_condition_1)

# Check that sample names match in both files
checking_Counts_DF_condition_1 <- colnames(Counts_DF_condition_1)
checking_Counts_DF_condition_1 <- checking_Counts_DF_condition_1[-1]

# check
all(colnames(checking_Counts_DF_condition_1) == rownames(condition_1@meta.data))

## condition2
# get the counts of the raw RNA counts
Counts_DF_condition_2 <- as.data.frame(condition_2@assays$RNA@counts) %>%  rownames_to_column()
Counts_DF_names_condition_2 <- colnames(Counts_DF_condition_2)
# View(Counts_DF_condition_2)

# Check that sample names match in both files
checking_Counts_DF_condition_2 <- colnames(Counts_DF_condition_2)
checking_Counts_DF_condition_2 <- checking_Counts_DF_condition_2[-1]

# check
all(colnames(checking_Counts_DF_condition_2) == rownames(condition_2@meta.data))

## Create Table using names as labels ----
# create table with name for condition 1 -> variable as_Name
label_row_condition_1_name <- c(condition1_name)
sub_group_labeled_condition_1_name <- rbind(label_row_condition_1_name, Counts_DF_condition_1)
sub_group_labeled_condition_1_name[1,1] <- c("label")
# View(sub_group_labeled_condition_1_name)
# save labeled sub group in directory Output_Numbered_Tables
saveName <- paste0("Output/",condition1_name,"_as_",condition1_name,".csv")
write.table(sub_group_labeled_condition_1_name, file = saveName , sep=",",  col.names=FALSE, row.names = FALSE)

# create table with name for condition 2-> variable as_Name
label_row_condition_2_name <- c(condition2_name)
sub_group_labeled_condition_2_name <- rbind(label_row_condition_2_name, Counts_DF_condition_2)
sub_group_labeled_condition_2_name[1,1] <- c("label")
View(sub_group_labeled_condition_2_name)
# save labeled sub group in directory Output_Numbered_Tables
saveName <- paste0("Output/",condition2_name,"_as_",condition2_name,".csv")
write.table(sub_group_labeled_condition_2_name, file = saveName , sep=",",  col.names=FALSE, row.names = FALSE)

# merge table with name
Merged_table_name <- dplyr::left_join(sub_group_labeled_condition_1_name, sub_group_labeled_condition_2_name, 
                                      by = c('rowname' = 'rowname'))

# save the table => This table will be the input for the Randomizer script (the next step of the workflow)
TableName <- paste0("PFA_Input/",condition1_name,"_as_",condition1_name,"_and_",condition2_name,"_as_",condition2_name,".csv")
write.table(Merged_table_name, file = TableName , sep=",",  col.names=FALSE)

## Create Table using numbers as labels ----
# create table with number for condition 1 -> variable as_Number
label_row_condition_1 <- c(as_Number_condition1_name)
sub_group_labeled_condition_1 <- rbind(label_row_condition_1, Counts_DF_condition_1)
sub_group_labeled_condition_1[1,1] <- c("label")
# View(sub_group_labeled_condition_1)
# save labeled sub group in directory Output_Numbered_Tables
saveName <- paste0("Output/",condition1_name,"_as_",as_Number_condition1_name,".csv")
write.table(sub_group_labeled_condition_1, file = saveName , sep=",",  col.names=FALSE, row.names = FALSE)

# create table with number for condition 2-> variable as_Number
label_row_condition_2 <- c(as_Number_condition2_name)
sub_group_labeled_condition_2 <- rbind(label_row_condition_2, Counts_DF_condition_2)
sub_group_labeled_condition_2[1,1] <- c("label")
View(sub_group_labeled_condition_2)
# save labeled sub group in directory Output_Numbered_Tables
saveName <- paste0("Output/",condition2_name,"_as_",as_Number_condition2_name,".csv")
write.table(sub_group_labeled_condition_2, file = saveName , sep=",",  col.names=FALSE, row.names = FALSE)

# merge table with number
Merged_table <- dplyr::left_join(sub_group_labeled_condition_1, sub_group_labeled_condition_2, 
                                 by = c('rowname' = 'rowname'))

# save the table => This table will be the input for the Randomizer script (the next step of the workflow)
TableName <- paste0("PFA_Input/",condition1_name,"_as_",as_Number_condition1_name,"_and_",condition2_name,"_as_",as_Number_condition2_name,".csv")
write.table(Merged_table, file = TableName , sep=",",  col.names=FALSE)
