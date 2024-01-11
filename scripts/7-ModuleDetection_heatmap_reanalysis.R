rm(list = ls())
gc()

source("scripts/FUNCTIONS.R") # loads packages too
rm("pkg", "pkgs", "pkgs.to.install")
# future::plan() - sets parallel session (using multisession())
plan(multisession, workers = detectCores())

curRes <- 2
r_cutoff <- 0.83
res_module <- 1
modulesToAnalyze <- c("2.1", "5.1")
dirs_to_create <- 
	paste0(c("output_tables/", "RDATA/", "plots/"),
	       paste0("Cor", r_cutoff, 
	              "_res", curRes, "_Subanalysis1/Module")
	)

dirs_to_create = sapply(modulesToAnalyze, \(x) {
	paste0(dirs_to_create, x)
}) %>% c() 
sapply(dirs_to_create, dir.create)

theme_set(theme_classic())
## Input files ##
load("RDATA/GeneSelection_objects.RData")
RDATA_plot <- paste0("RDATA/MainAnalysis_Cor", r_cutoff,
		     "/ModuleData_forPlot_res", curRes, ".RData") 
load(RDATA_plot)
rm(high_var_genes_5pct)

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

for (curSubmodule in modulesToAnalyze) {
	dirSubmodule = grep(value = T,
			    curSubmodule, dirs_to_create)
	names(dirSubmodule) = c("output_tables", "RDATA", "plots")
	curModule <- gsub("(\\d)\\.\\d",
			  "\\1", curSubmodule)
	genesModuleInterest <- 
		(read_delim(paste0(dirSubmodule["output_tables"],
				  "/../Module", curModule,
				  "_annotations.tsv")) %>%
		filter(Module_name == curSubmodule))$geneID
	Exp_table_long_averaged_z <- 
		Exp_table_long_averaged_z_high_var %>%
		filter(gene_ID %in% genesModuleInterest)
	
	z_score_wide <- Exp_table_long_averaged_z %>%
		dplyr::select(gene_ID, Sample_name, z.score) %>% 
		pivot_wider(names_from = Sample_name, 
			    values_from = z.score) %>% 
		as.data.frame()
	
	row.names(z_score_wide) <- z_score_wide$gene_ID
	nreps <- ncol(z_score_wide) - 1
	# Gene correlation matrix
	cor_matrix <- cor(t(z_score_wide[, -1]))
	
	##### Heatmap all genes of module #####
	heatmapGenes <- 
		heatmap_all_genes(cor_matrix, 
				  left_join(Exp_table_long_averaged_z, 
				  	  metadata_sampleOrder),
				  caption = NULL, 
				  title = paste("All genes in module", curSubmodule)) +
		facet_grid(cols = vars(Sample_type),
			   scales = "free") +
		heatmap_theme + 
		theme(axis.text.y = element_blank(),
		      axis.ticks.y = element_blank())
	ggsave(heatmapGenes, 
	       filename = paste0(dirSubmodule["plots"], 
	       		  "/Module", curSubmodule, "_heatmap_allGenes.png"),
	       height = 6.14, width = 6.86)
	
	# Transform the lower triangle of the matrix into NA, since it is a duplication of the upper part
	cor_matrix_upper_tri <- cor_matrix
	cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA
	
	
	# Edge selection: select only statistically significant correlations
	#' t-distribution approximation
	#' For each correlation coefficient r, you approximate a t statistics.
	#' The equation is t = r * ( (n-2) / (1 - r^2) )^0.5
	edge_table <- cor_matrix_upper_tri %>% 
		as.data.frame() %>% 
		mutate(from = row.names(cor_matrix)) %>% 
		pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% 
		filter(!is.na(r)) %>% # remove the lower triangle
		filter(from != to) %>% # remove self-to-self correlations
		mutate(t = get_T_value(r, nreps)) %>% 
		mutate(p.value = get_p_value(t, nreps)) %>% 
		mutate(FDR = p.adjust(p.value, method = "fdr"),
		       significant = ifelse(FDR < 0.05, T, F))
	
	#### DECIDE r_cutoff ####
	sub_r_cutoff = min((edge_table %>% 
			    	mutate(r_bins = round(r, digits = 3)) %>%
			    	group_by(r_bins) %>% 
			    	summarize(allSig = all(significant)) %>% 
			    	filter(allSig == T, 
			    	       r_bins > 0))$r_bins)
	
	ggplot(edge_table) +
		geom_histogram(aes(x = r,
				   fill = significant),
				   color = "white",
				   bins = 100) +
		geom_vline(xintercept = sub_r_cutoff, 
			   color = "black",
			   linewidth = 0.8) +
		# check at which r value, the number of correlations starts to drop fast
		labs(caption = paste0("r of gene correlations of module ", 
				      curSubmodule, "\nResolution ", curRes)) +
		labs(fill = "FDR < 0.05") +
		scale_x_continuous(breaks = seq(-.5, 1, .2)) +
		theme(text = element_text(size = 14),
		      axis.text = element_text(color = "black")
		)
	
	r_histo_name = paste0(dirSubmodule["plots"],
			      "/Module", curSubmodule)
	# ggsave(paste0(r_histo_name, ".svg"), height = 3.5, width = 5, bg = "white")
	ggsave(paste0(r_histo_name, ".png"), height = 3.5, width = 5, bg = "white")
	
	
	# co-expressed in the same direction
	split_at <- 90000
	dirPath = paste0(dirSubmodule["output_tables"],
			 "/Module", curSubmodule, "_edge_table"
	)
	if (dir.exists(dirPath)) {
		for (i in list.files(path = dirPath)) {
			paste0(dirPath, "/", i) %>%
				file.remove()
		}
	} else {
		dir.create(dirPath)
	}
	
	edge_table <- 
		edge_table %>% 
		filter(r > sub_r_cutoff)
	edge_table %>% 
		mutate(group = (row_number() - 1) %/% !! split_at) %>%
		group_split(group) %>%
		map(.f = ~{
			fileName = paste0(dirPath, "/edge_table_r",
					  sub_r_cutoff, "_",
					  unique(.x$group), ".tsv")
			write_delim(.x, fileName,
				    quote = "none", append = F,
				    col_names = T, delim = "\t")
		}) 
	
	node_table <- with(
		edge_table,
		data.frame(gene_ID = c(from, to) %>% unique()
		)) %>% 
		left_join(funct_anno, 
			  by = join_by("gene_ID" == "geneID")) %>%
		unique()

	cat("Current module:", curSubmodule, "-", nrow(node_table), "genes\n")
	
	# Create graph
	my_network <- graph_from_data_frame(
		edge_table,
		vertices = node_table,
		directed = F
	)
	
	optimization_results <- get_optimization_table(
		my_network, min_res = 0.25, max_res = 5, step = 0.25
	)
	
	## Plot optimization results ##
	transform_factor = 
		with(optimization_results,
		     max(num_contained_gene)/max(num_module)
		)
	
	optimization_results %>% 
		ggplot(aes(x = resolution)) +
		geom_line(aes(y = num_module),
			  linewidth = 1.1, alpha = 0.8, color = "dodgerblue2") +
		geom_point(aes(y = num_module),
			   size = 3, alpha = 0.7) +
		geom_line(aes(y = num_contained_gene/transform_factor),
			  linewidth = 1.1, alpha = 0.8, color = "violetred2") +
		geom_point(aes(y = num_contained_gene/transform_factor),
			   size = 3, alpha = 0.7) +
		scale_y_continuous(sec.axis = 
				   	sec_axis(~. * transform_factor,
				   		 name = "num. genes in\nmodules w/ >=3 genes"
				   	)) +
		labs(x = "resolution parameter",
		     y = "num. modules\nw/ >=3 genes") +
		theme_optimization
	
	ggsave(paste0(dirSubmodule["plots"], "/Module", curSubmodule,
		      "_Optimize_resolution.png"), 
	       height = 5, width = 5, bg = "white")
	
	#### Construct networks
	modules <- cluster_leiden_mod(
		my_network, r = res_module
	)
	
	my_network_modules <- data.frame(
		gene_ID = names(membership(modules)),
		module = as.vector(membership(modules)))
	
	# Continue only with modules with 3 or more genes
	module_3 <- my_network_modules %>% 
		group_by(module) %>% count() %>% 
		arrange(-n) %>% filter(n >= 3)
	
	my_network_modules <- my_network_modules %>% 
		filter(module %in% module_3$module)
	
	Expr_averaged_z_high_var_modules <- 
		Exp_table_long_averaged_z_high_var %>% 
		inner_join(my_network_modules, by = "gene_ID")
	
	modules_mean_z <- Expr_averaged_z_high_var_modules %>% 
		group_by(module, Sample_name) %>% 
		summarise(mean.z = mean(z.score)) %>% 
		ungroup()
	
	# This will create a table with the module number and the Sample in which the module is most expressed.
	module_peak_exp <- modules_mean_z %>% 
		group_by(module) %>% 
		slice_max(order_by = mean.z, n = 1) %>%
		mutate(orderedSamples = sapply(Sample_name, \(x) {
			which(levels(metadata$Sample_name) == x)
		})) %>%
		mutate(ordered_modules = 
		       	reorder(module, -orderedSamples)) %>%
		dplyr::select(module, Sample_name,
			      orderedSamples, ordered_modules)
	
	heatmap_data = 
		full_join(modules_mean_z, module_peak_exp[,-c(2, 3)], 
			  by = c('module')) %>%
		# Join with metadata, removing the column Replicate
		full_join(unique(metadata[, -2]), by = "Sample_name")
	
	# genes per module
	genes_per_module = Expr_averaged_z_high_var_modules %>% 
		group_by(module) %>% 
		summarize(gene_ID = unique(gene_ID)) %>% 
		count()
	
	heatmap_data = full_join(heatmap_data, genes_per_module,
				 by = "module")
	
	### Heatmap of module mean expression z-score
	module_heatmap <- heatmap_data %>% 
		ggplot(aes(x = Sample_name, y = ordered_modules)) +
		geom_tile(aes(fill = mean.z), color = "grey80") +
		scale_fill_gradient2(mid = "white",
				     high = "#67001F",
				     low = "#053061",
				     breaks = c(-1.5, 0, 1.5), 
				     labels = c("< -1.5", "0", "> 1.5")) +
		labs(x = NULL, y = "Module", fill = "z-score",
		     title = paste0("Nested analysis on module ", curSubmodule),
		     caption = paste0("r threshold = ", sub_r_cutoff, 
		     		 "\nresolution = ", res_module,
		     		 "\nTotal genes:",
		     		 sum(genes_per_module$n))
		) +
		facet_grid(cols = vars(Sample_type),
			   scales = "free") +
		heatmap_theme
	
	module_heatmap_nGenes = genes_per_module %>%
		ggplot(aes(y = as.factor(module), x = "",
			   label = n)) +
		geom_tile(fill = "white", color = "white") +
		geom_text() + 
		labs(x = "Genes\nin\nmodule") +
		annotation_theme
	
	wrap_plots(module_heatmap, module_heatmap_nGenes, 
		   ncol = 2, widths = c(1, 0.08))
	
	plotName = paste0(dirSubmodule["plots"],
			  "/Module", curSubmodule, 
			  "_heatmap_r", sub_r_cutoff, 
			  "_res", res_module, 
			  ".tiff")
	ggsave(plotName, height = 6.14, width = 6.86)
	
	annotation <- funct_anno %>%
		inner_join(my_network_modules,
			   by = join_by("geneID" == 'gene_ID')) %>%
		mutate(Module_name = paste0(curSubmodule, ".", module)) %>%
		arrange(module) %>%
		write_delim(file = paste0(dirSubmodule["output_tables"],
					  "/Module", curSubmodule,
					  "_annotations.tsv"),
			    delim = "\t", col_names = T, na = "", 
			    quote = "none")
}
