#### Install and load packages ####
pkgs_bioconductor = c("DESeq2", "GenomicFeatures")
pkgs.To.Install_bioconductor = 
	pkgs_bioconductor[! pkgs_bioconductor %in% installed.packages()]
pkgs = c("tidyverse", "devtools", "ggpubr")
pkgs.To.Install = pkgs[! pkgs %in% installed.packages()]

if (length(pkgs.To.Install) > 0) install.packages(pkgs.To.Install)
if (length(pkgs.To.Install_bioconductor) > 0) {
	BiocManager::install(pkgs.To.Install_bioconductor)
}

# CustomSelection has to be installed from github, so it cannot be grouped with the others yet
pkgs = c(pkgs, pkgs_bioconductor)

for (curPkg in pkgs) library(curPkg, character.only = T) 

##### VARIABLES #####
DAFS_ran = T

#### CustomSelection ####
Exp_table_files = paste0("Data/",
			 list.files(pattern = "kallisto_TPM_TMM_normalized_",
			 	   path = "Data")
			 )
Exp_table <- sapply(Exp_table_files, simplify = F, 
		    \(Exp) {
		    	read_csv(Exp, col_names = T) %>%
		    		dplyr::select(!group)
		    }) %>% list_rbind %>%
	as.data.frame()
rownames(Exp_table) = Exp_table$geneID
Exp_table = Exp_table[, -ncol(Exp_table)]

metadata <- read_delim("Data/samples_file.txt", 
		       col_names = T) |>
	dplyr::select(Sample_name, Replicate)

if (DAFS_ran) {
	cutv_table = read.delim("Data/cutv.tsv")
	cutv = cutv_table$thresholds
	names(cutv) = cutv_table$Replicates
} else {
	source("scripts/run_DAFS.r")
}

# We check, for each replicate, if the gene passed its expression threshold
passCutOff_by_rep = sapply(names(cutv), \(Replicate) {
	Exp_table[[Replicate]] > cutv[Replicate]^2
}) %>% cbind() %>% as.data.frame
rownames(passCutOff_by_rep) = rownames(Exp_table)

# Then, we check, for each sample, if the gene passed the threshold in all replicates
new_Exp = sapply(unique(metadata$Sample_name), simplify = F,
		       \(sampleName) {
	replicates = (metadata %>% filter(Sample_name == sampleName))$Replicate
	passedCutoff = apply(passCutOff_by_rep[, replicates], 1, all)
	sapply(replicates, \(replicate) {
		map2_dbl(passedCutoff, Exp_table[replicate],
			 \(x, y) {
			 	ifelse(x, y, 0)
			 })
	}) %>% cbind 
})

new_Exp_table = data.frame(cbind(
	new_Exp[[1]], new_Exp[[2]], new_Exp[[3]],
	new_Exp[[4]], new_Exp[[5]], new_Exp[[6]],
	new_Exp[[7]], new_Exp[[8]]
), row.names = rownames(Exp_table))


# we keep the genes that passed the threshold in at least one sample
Exp_filtered = new_Exp_table[rowSums(new_Exp_table) != 0,] 

Exp_filtered %>% 
	mutate(gene = rownames(Exp_filtered)) %>%
	write_csv("output_tables/Filtered_kallisto_TPM.csv",
		  quote = "none", append = F,
		  col_names = T)
	
