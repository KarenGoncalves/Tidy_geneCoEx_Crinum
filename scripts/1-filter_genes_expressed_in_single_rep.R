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

#### Gene selection by Ex95 ####
Exp_table_files = paste0("Data/",
			 list.files(pattern = "kallisto_TPM_TMM_normalized_",
			 	   path = "Data")
			 )
genes_Ex95_allSamples = read_delim("Data/genes_Ex95_allSamples.tsv", 
				   delim = "\t", col_names = "gene_id")

Exp_table <- sapply(Exp_table_files, simplify = F, 
		    \(Exp) {
		    	read_csv(Exp, col_names = T) %>%
		    		dplyr::select(!group)
		    }) %>% list_rbind %>%
	as.data.frame()
rownames(Exp_table) = Exp_table$geneID
Exp_table = Exp_table %>% 
	filter(geneID %in% genes_Ex95_allSamples$gene_id) %>%
	dplyr::select(!geneID)

metadata <- read_delim("Data/samples_file.txt", 
		       col_names = T) |>
	dplyr::select(Sample_name, Replicate)

# Then, we check, for each sample, if the gene passed the threshold in all replicates
new_Exp = sapply(unique(metadata$Sample_name), simplify = F,
		       \(sampleName) {
	replicates = (metadata %>% filter(Sample_name == sampleName))$Replicate
	passedCutoff = apply(Exp_table[, replicates], 1, \(x) all(x > 0))
	sapply(replicates, \(replicate) {
		sapply(1:length(passedCutoff), \(x) {
			ifelse(passedCutoff[x], Exp_table[x, replicate], 0)
		})
	}) %>% cbind %>% as.data.frame
}) 

new_Exp_table = new_Exp %>% do.call(what="cbind")
names(new_Exp_table) = gsub("(\\w+|Basal.plate)\\.\\1(_\\d)$", 
			    "\\1\\2", names(new_Exp_table) )

# we keep the genes that passed the threshold in at least one sample
Exp_filtered = new_Exp_table[rowSums(new_Exp_table) != 0,] 

Exp_filtered %>% 
	mutate(gene = rownames(Exp_filtered)) %>%
	write_csv("output_tables/Filtered_kallisto_TPM.csv",
		  quote = "none", append = F,
		  col_names = T)
	
