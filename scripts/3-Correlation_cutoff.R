source("scripts/FUNCTIONS.R") # loads packages too
Baits = read_delim("Data/Genes_of_interest.txt")
load("RDATA/GeneSelection_objects.RData")

# Wide z-score table for selected genes

z_score_wide <- Exp_table_long_averaged_z_high_var %>% 
	filter(gene_ID %in% high_var_genes_5pct$gene_ID) %>%
	select(gene_ID, Sample_name, z.score) %>% 
	pivot_wider(names_from = Sample_name, values_from = z.score) %>% 
	as.data.frame()

row.names(z_score_wide) <- z_score_wide$gene_ID
nreps = ncol(z_score_wide) - 1
# Gene correlation matrix
cor_matrix <- cor(t(z_score_wide[, -1]))

# Transform the lower triangle of the matrix into NA, since it is a duplication of the upper part
cor_matrix_upper_tri <- cor_matrix
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA

# Edge selection: select only statistically significant correlations
#' t-distribution approximation
#' For each correlation coefficient r, you approximate a t statistics.
#' The equation is t = r * ( (n-2) / (1 - r^2) )^0.5
# split cor table into 10 pieces
#splitSize = nrow(z_score_wide) %/% 10000
#colSplits = seq(splitSize, nrow(z_score_wide), splitSize)
colSplits = seq(100, nrow(z_score_wide), 100)

edge_table = matrix(nrow = 0, ncol = 6) %>%
	as.data.frame(col.names = 
	      	c("from", "to", "r", "t", 
	      	  "p.value", "FDR"))

for (i in seq_len(length(colSplits)) ) {
	start = ifelse(i == 1, 1, end + 1)
	end = ifelse(i == length(colSplits), nrow(z_score_wide), colSplits[i])

	tmp = cor_matrix_upper_tri[start:end,] %>% 
		as.data.frame() %>% 
		mutate(from = 
		       	row.names(cor_matrix)[start:end]) %>% 
		pivot_longer(cols = !from, 
			     names_to = "to", 
			     values_to = "r") %>% 
		# remove the lower triangle
		filter(!is.na(r)) %>% 
		# remove self-to-self correlations
		filter(from != to) %>% 
		mutate(t = r*sqrt( (nreps-2) / (1-r^2) ) ) %>% 
		mutate(p.value = case_when(
			t > 0 ~ pt(t, df = nreps-2, lower.tail = F),
			t <=0 ~ pt(t, df = nreps-2, lower.tail = T)
		)) %>% 
		mutate(FDR = p.adjust(p.value, method = "fdr"),
		       significant = ifelse(FDR < 0.01, T, F)) 

	edge_table = as.data.frame(
		rbind(edge_table, tmp)
	)
}

# edge_table <- cor_matrix_upper_tri %>% 
	# as.data.frame() %>% 
	# mutate(from = row.names(cor_matrix)) %>% 
	# pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% 
	# filter(!is.na(r)) %>% # remove the lower triangle
	# filter(from != to) %>% # remove self-to-self correlations
	# mutate(t = r*sqrt( (nreps-2) / (1-r^2) ) ) %>% 
	# mutate(p.value = case_when(
		# t > 0 ~ pt(t, df = nreps-2, lower.tail = F),
		# t <=0 ~ pt(t, df = nreps-2, lower.tail = T)
	# )) %>% 
	# mutate(FDR = p.adjust(p.value, method = "fdr")) 

# Plot the distribution of r values

#### DECIDE r_cutoff ####
r_cutoff = min((edge_table %>% 
			mutate(r_bins = round(r, digits = 3)) %>%
			group_by(r_bins) %>% 
			summarize(allSig = all(significant)) %>% 
			filter(allSig == T, 
			       r_bins > 0))$r_bins)
edge_table %>% 
	ggplot() +
	geom_histogram(aes(x = r,
			   fill = significant),
			   color = "white",
			   bins = 100) +
	geom_vline(xintercept = r_cutoff, 
		   color = "black",
		   linewidth = 0.5) +
	# check at which r value, the number of correlations starts to drop fast
	labs(fill = "FDR < 0.01", 
	     title = "Gene correlations") +
	scale_x_continuous(
		breaks = seq(-1, 1, .2)) +
	theme_classic() +
	theme(text = element_text(size = 14),
	      axis.text = element_text(color = "black")
	)

#ggsave("plots/r_histogram.svg", height = 3.5, width = 5, bg = "white")
ggsave("plots/MainAnalysis/r_histogram.png", height = 3.5, width = 5)

edge_table %>% 
	filter(from %in% Baits$transcript |
	       	to %in% Baits$transcript) %>% 
	ggplot() +
	geom_histogram(aes(x = r,
			   fill = significant),
			   color = "white",
			   bins = 100) +
	geom_vline(xintercept = r_cutoff, 
		   color = "black",
		   linewidth = 0.85) +
	# check at which r value, the number of correlations starts to drop fast
	labs(fill = "FDR < 0.01", 
	     title = "Gene correlations - Bait genes") +
	scale_x_continuous(
		breaks = seq(-.5, 1, .2)) +
	theme_classic() +
	theme(text = element_text(size = 14),
	      axis.text = element_text(color = "black")
	)
ggsave("plots/MainAnalysis/r_histogram_baitGenes.png", 
       height = 3.5, width = 5)

	
# co-expressed in the same direction
split_at <- 90000

if (dir.exists("output_tables/Main_analysis_edge_table")) {
 	for (i in list.files(path = "output_tables/Main_analysis_edge_table")) {
		paste0("output_tables/Main_analysis_edge_table/", i) %>%
	 	 	file.remove()
	}
} else {
	dir.create("output_tables/Main_analysis_edge_table")
}

edge_table %>% 
	filter(r > r_cutoff) %>%
	mutate(group = (row_number() - 1) %/% !! split_at) %>%
	group_split(group) %>%
	map(.f = ~{
		fileName = paste0("output_tables/Main_analysis_edge_table/edge_table_r",
				  r_cutoff, "_",
				  unique(.x$group), ".tsv")
		write_delim(.x, fileName,
			    quote = "none", append = F,
			    col_names = T, delim = "\t")
	}) 
