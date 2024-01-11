gc()
source("scripts/FUNCTIONS.R") # loads packages too,
metadata <- read_delim("Data/samples_file.txt", 
		       col_names = T, delim = "\t") 
#####  Reorder Sample_type and Sample_order columns, to set the modules in order in the plot ####
metadata$Sample_type = factor(metadata$Sample_type, 
			      levels = c("In vivo", "In vitro"))
metadata_sampleOrder = metadata %>%
	dplyr::select(Sample_name, Sample_type, Sample_order) %>%
	arrange(Sample_type, Sample_order) %>%
	unique
metadata$Sample_name = 
	factor(metadata$Sample_name, 
	       levels = metadata_sampleOrder$Sample_name)

Baits <- read_delim("Data/Genes_of_interest.txt")
funct_anno = read_delim("output_tables/Annotation_formatted.tsv")

r_cutoffs <- grep("MainAnalysis_Cor",
		  list.dirs("RDATA"), value = T) %>%
	gsub(pattern = "RDATA/MainAnalysis_Cor",
	     replacement = "")
resolutions <-  
	sapply(r_cutoffs, \(x) {
		list.files(path = paste0("RDATA/MainAnalysis_Cor", x),
			 pattern = "ModuleData") %>%
	gsub(pattern = ".+_res([0-9\\.]+).RData", 
	     replacement = "\\1") %>% unique
})

for (cur_rCutoff in r_cutoffs) {
	for (curResolution in resolutions) {
		RDATA = paste0("RDATA/MainAnalysis_Cor", cur_rCutoff, 
			       "/ModuleData_forPlot",
			       "_res", curResolution, ".RData")
		load(RDATA)
		annotationPath <- paste0("output_tables/Main_analysis_annotations_Cor",
					cur_rCutoff, "/Resolution",
					curResolution)
		dir.create(annotationPath, recursive = T)
		sapply(Expr_averaged_z_high_var_modules$module %>% unique,
		       \(curModule) {
		       	Expr_averaged_z_high_var_modules %>%
		       		dplyr::select(!all_of(c("Sample_name", "mean.logTPM", "z.score"))) %>%
		       		unique %>%
		       		filter(module == curModule) %>% 
		       		write_delim(
		       			file = paste0(annotationPath, "/",
		       				      "Module_", curModule, ".tsv"),
		       			delim = "\t", na = "", append = F,
		       			quote = "none",
		       			col_names = T
		       		)
		       })
		# Now the modules in the table module_peak_exp are a factor ordered by the 
		# Sample in which they are most expressed, so we join the tables:
		# metadata, module_mean_z and module_peak_exp
		module_peak_exp_mod <- 
			mutate(module_peak_exp, 
			       orderedSamples = sapply(Sample_name, \(x) {
			       	which(levels(metadata$Sample_name) == x)
			       })) %>%
			dplyr::select(module, Sample_name,
				      orderedSamples)
		
		heatmap_data = 
			full_join(modules_mean_z, module_peak_exp_mod[,-c(2, 3)], 
				  by = c('module')) %>%
			# Join with metadata, removing the column Replicate
			full_join(unique(metadata[, -2]), by = "Sample_name")
		
		# genes per module
		genes_per_module = Expr_averaged_z_high_var_modules %>% 
			group_by(module) %>% 
			summarize(gene_ID = unique(gene_ID)) %>% 
			count()
		
		heatmap_data = full_join(heatmap_data, genes_per_module,
					 by = "module") %>%
			arrange(desc(n)) %>%
			mutate(ordered_modules = factor(module,
							levels = unique(.$module))
			)
		
		### Heatmap of module mean expression z-score
		module_heatmap <- heatmap_data %>% 
			ggplot(aes(x = paste0(Sample_name, "\n", 
					      ifelse(is.na(Hormone), "", Hormone)),
				   y = ordered_modules)) +
			geom_tile(aes(fill = mean.z), color = "grey80") +
			scale_fill_gradient2(mid = "white",
					     high = "#67001F",
					     low = "#053061",
					     breaks = c(-1.5, 0, 1.5), 
					     labels = c("< -1.5", "0", "> 1.5")) +
			labs(x = NULL, y = "Module", fill = "z-score",
			     caption = paste0("r threshold = ", cur_rCutoff, 
			     		 "\nresolution = ", curResolution)
			) +
			facet_grid(cols = vars(Sample_type),
				   scales = "free") +
			heatmap_theme
		
		module_heatmap_nGenes = heatmap_data %>%
			select(ordered_modules, n) %>%
			unique %>%
			ggplot(aes(y = ordered_modules, x = "",
				   label = n)) +
			geom_tile(fill = "white", color = "white") +
			geom_text() + labs(x = "Genes\nin\nmodule") +
			annotation_theme
		
		wrap_plots(module_heatmap, 
			   module_heatmap_nGenes, 
			   ncol = 2, 
			   widths = c(1, 0.08))
		plotName = paste0("plots/MainAnalysis/Heatmap_Modules_",
				  "r", cur_rCutoff, "_res", curResolution, 
				  "highVarTop5pct.tiff")
		ggsave(plotName, height = 8.14, width = 8.86)
		
		annotation <- funct_anno %>%
			inner_join(my_network_modules,
				   by = join_by("geneID" == 'gene_ID')) %>%
			mutate(Module_name = curModule) %>%
			arrange(module) %>%
			write_delim(file = paste0("output_tables/Module", 
						  curModule,
						  "_annotations.tsv"),
				    delim = "\t", col_names = T, na = "", 
				    quote = "none")
		
	}
}
