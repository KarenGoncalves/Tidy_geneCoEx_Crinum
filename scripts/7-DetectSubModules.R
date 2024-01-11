source("scripts/FUNCTIONS.R") # loads packages too
modulesToAnalyze = c(2, 3, 4)
r_cutoff = 0.8
# future::plan() - sets parallel session (using multisession())
plan(multisession, workers = detectCores())

## Input files ##
Exp_table_long_averaged_z <-
	read_delim("Data/Exp_table_long_average_z.tsv",
		   col_names = T)
metadata <- read_delim("Data/samples_file.txt", 
		       col_names = T, delim = "\t") 
Baits <- read_delim("Data/Genes_of_interest.txt",
		    col_names = T)
# Load gene annotation table
funct_anno = read_delim("output_tables/Annotation_formatted.tsv")

for (curModule in modulesToAnalyze) {
	edge_table = read_delim(paste0("output_tables/edgeTable_module_", 
				       curModule, ".tsv"))
	edge_table_select = edge_table %>% 
		filter(r >= r_cutoff) %>% unique
	
	dim(edge_table_select)
	
	
	# Merge annotation and edge tables
	node_table <- data.frame(
		gene_ID = c(edge_table_select$from, 
			    edge_table_select$to) %>% unique()
	) %>% 
		left_join(funct_anno, by = c("gene_ID" = "gene_ID")) %>%
		unique()
	
	head(node_table)
	dim(node_table)
	
	# Create graph
	my_network <- graph_from_data_frame(
		edge_table_select,
		vertices = node_table,
		directed = F
	)
	
	#' Too low resolution leads to forcing genes with different expression 
	#' patterns into the same module.
	#' Too high resolution leads to many genes not contained in any one module.
	
	optimization_results <- future_map(
		.x = seq(from = 0.25, to = 5, by = 0.25),
		\(x) optimize_resolution(network = my_network, 
					 resolution = x)) %>% 
		lapply(\(x) data.frame(V1 = x[1], V2 = x[2])) %>%
		list_rbind %>% 
		mutate(resolution = seq(from = 0.25, 
					to = 5, by = 0.25)) %>% 
		rename(num_module = V1,
		       num_contained_gene = V2)
	
	head(optimization_results)
	
	## Plot optimization results ##
	Optimize_num_module <- optimization_results %>% 
		ggplot(aes(x = resolution, y = num_module)) +
		geom_line(linewidth = 1.1, 
			  alpha = 0.8, 
			  color = "dodgerblue2") +
		geom_point(size = 3, alpha = 0.7) +
		geom_vline(xintercept = 2, 
			   linewidth = 1, linetype = 4) +
		labs(x = "resolution parameter",
		     y = "num. modules\nw/ >=5 genes") +
		theme_classic() +
		theme(text = element_text(size = 14),
		      axis.text = element_text(color = "black")
		)
	
	Optimize_num_gene <- optimization_results %>% 
		ggplot(aes(x = resolution, 
			   y = num_contained_gene)) +
		geom_line(linewidth = 1.1, 
			  alpha = 0.8, 
			  color = "violetred2") +
		geom_point(size = 3, alpha = 0.7) +
		geom_vline(xintercept = 2, 
			   linewidth = 1, 
			   linetype = 4) +
		labs(x = "resolution parameter",
		     y = "num. genes in\nmodules w/ >=5 genes") +
		theme_classic() +
		theme(text = element_text(size = 14),
		      axis.text = element_text(color = "black")
		)
	
	wrap_plots(Optimize_num_module, Optimize_num_gene, nrow = 2)
	
	optimization_results_plot = 
		paste0("plots/Optimize_resolution", 
		       curModule, ".png")
	# ggsave("plots/Optimize_resolution.svg", 
	#        height = 5, width = 3.2, bg = "white")
	ggsave(optimization_results_plot,
	       height = 5, width = 3.2, bg = "white")
	
	## Select resolution and detect modules
	resolution = 1.25# change for best one found above
	modules <- cluster_leiden(my_network, 
				  resolution_parameter = resolution, 
				  objective_function = "modularity")
	
	# Merge edge and node tables
	my_network_modules <- data.frame(
		gene_ID = names(membership(modules)),
		module = as.vector(membership(modules)) ) %>% 
		inner_join(node_table, by = "gene_ID")
	
	my_network_modules %>% 
		group_by(module) %>% 
		count() %>% 
		arrange(-n) %>% 
		filter(n >= 5)
	
	my_network_modules %>% 
		group_by(module) %>% 
		count() %>% 
		arrange(-n) %>% 
		filter(n >= 5) %>% 
		ungroup() %>% 
		summarise(sum = sum(n))
	
	# Continue only with modules with 5 or more genes
	module_5 <- my_network_modules %>% 
		group_by(module) %>% 
		count() %>% 
		arrange(-n) %>% 
		filter(n >= 5)
	
	my_network_modules <- my_network_modules %>% 
		filter(module %in% module_5$module)
	# Check how many bait genes are in the modules selected
	my_network_modules %>%
		filter(gene_ID %in% Baits$transcript) %>%
		nrow()
	
	## Peak expression of modules
	
	Exp_table_long_averaged_z_modules <- 
		Exp_table_long_averaged_z %>% 
		inner_join(my_network_modules, by = "gene_ID")
	
	write_delim(Exp_table_long_averaged_z_modules,
		    paste0("Data/Exp_table_averaged_module", curModule,
		           "_submodules.tsv"),
		    quote = "none", delim = "\t", na = '',
		    col_names = T)
	head(Exp_table_long_averaged_z_high_var_modules)
	
	modules_mean_z <- Exp_table_long_averaged_z_high_var_modules %>% 
		group_by(module, Sample_name) %>% 
		summarise(mean.z = mean(z.score)) %>% 
		ungroup()
	
	modules_mean_z %>%
		write_delim(
			paste0("Data/submodules_mean_z_module", 
			       curModule,".tsv"),
			quote = "none", delim = "\t", na = '',
			col_names = T)
	
	# This will create a table with the module number and the Sample in which the module is most expressed.
	module_peak_exp <- modules_mean_z %>% 
		group_by(module) %>% 
		slice_max(order_by = mean.z, n = 1)
	
	module_peak_exp %>%
		write_delim(
			paste0("Data/submodule_peak_exp_module", 
			       curModule,".tsv"),
			quote = "none", delim = "\t", na = '',
			col_names = T)
	
}