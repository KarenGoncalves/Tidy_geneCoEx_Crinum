source("scripts/FUNCTIONS.R") # loads packages too

## Input files ##
Exp_table <- read.csv("output_tables/Filtered_kallisto_TPM.csv",
		      header = T)
metadata <- read_delim("Data/samples_file.txt", 
		      col_names = T, delim = "\t") 

names(Exp_table)[ncol(Exp_table)] <- "gene_ID"

Baits <- read_delim("Data/Genes_of_interest.txt")

## Long exp_table ##

Exp_table_long <- Exp_table %>% 
	pivot_longer(cols = !gene_ID, 
		     names_to = "library", 
		     values_to = "tpm") %>% 
	mutate(logTPM = log10(tpm + 1)) 

Exp_table_log_wide <- Exp_table_long %>% 
	dplyr::select(gene_ID, library, logTPM) %>% 
	pivot_wider(names_from = library, 
		    values_from = logTPM, 
		    id_cols = gene_ID)

## PCA ##
my_pca <- prcomp(t(Exp_table_log_wide[, -1]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)

PCA_coord <- my_pca$x[, 1:10] %>% 
	as.data.frame() %>% 
	mutate(Replicate = row.names(.)) %>% 
	full_join(metadata, by = "Replicate")

axis_titles = sapply(1:2, \(x) {
	paste("PC", x, " (", 
	      pc_importance[x, 2] %>% signif(3)*100, 
	      "% of Variance)", sep = "")
})

PCA_by_sample <- PCA_coord %>% 
	ggplot(aes(x = PC1, y = PC2)) +
	geom_point(aes(fill = Sample_name), 
		   color = "grey20", 
		   shape = 21, 
		   size = 3, alpha = 0.8) +
	scale_fill_manual(values = brewer.pal(n = 8, "Accent")) +
	labs(x = axis_titles[1], 
	     y = axis_titles[2],
	     fill = NULL) +  
	theme_bw() +
	theme(text = element_text(size= 14),
	      axis.text = element_text(color = "black"),
	      legend.position = "bottom"
	)

PCA_by_sample

#ggsave("plots/PCA.svg", height = 5, width = 5.5, bg = "white")
ggsave("plots/MainAnalysis/PCA.png", height = 5, width = 5, bg = "white")

### Gene co-expression analysis ###
# Average up the reps #
Exp_table_long_averaged <- Exp_table_long %>% 
	full_join(metadata, 
		  by = join_by("library" == "Replicate")) %>% 
	group_by(gene_ID, Sample_name) %>% 
	summarise(mean.logTPM = mean(logTPM)) %>% 
	ungroup()  

head(Exp_table_long_averaged)

# Z-score #

Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
	group_by(gene_ID) %>% 
	mutate(z.score = zscore(mean.logTPM)) %>% 
	ungroup()
split_at <- 400000

Exp_table_long_averaged_z %>% 
	mutate(group = (row_number() - 1) %/% !! split_at) %>%
	group_split(group) %>%
	map(.f = ~{
		fileName = paste0("output_tables/Exp_table_long_averaged_z_",
				  unique(.x$group), ".tsv")
		write_delim(.x, fileName,
			    quote = "none", append = F,
			    col_names = T, delim = "\t")
	}) 

head(Exp_table_long_averaged_z)

# Gene selection #

all_var_and_ranks <- Exp_table_long_averaged_z %>% 
	group_by(gene_ID) %>% 
	summarise(var = var(mean.logTPM)) %>% 
	ungroup() %>% 
	mutate(rank = rank(var, ties.method = "average")) 
write_delim(all_var_and_ranks,
	    "output_tables/all_var_and_ranks.tsv",
	    col_names = T, quote = "none",
	    delim = "\t", na = '')

high_var_genes <- all_var_and_ranks %>%
	filter(var > quantile(var, 0.66))

high_var_genes_5pct <- high_var_genes %>% 
	slice_max(order_by = var, 
		  n = round(nrow(Exp_table)/20, digits = 0)
	)
write_delim(high_var_genes_5pct,
	    file = "output_tables/high_var_genes_5pct.tsv",
	    col_names = T, delim = "\t", quote = "none"
)

bait_var <- high_var_genes_5pct %>% 
	filter(gene_ID %in% Baits$transcript) %>% 
	group_by(gene_ID) %>% 
	slice_max(n = 1, order_by = var)

Exp_table_long_averaged_z_high_var <- 
	Exp_table_long_averaged_z %>% 
	filter(gene_ID %in% high_var_genes$gene_ID)

save(Exp_table_long_averaged_z_high_var,
     high_var_genes_5pct, 
     file = "RDATA/GeneSelection_objects.RData"
     )
# Check where the bait genes fall along the variance distribution
var_plot = all_var_and_ranks %>% 
	ggplot(aes(x = var, y = rank))  +
	geom_hline(
		data = bait_var, aes(yintercept = rank),
		color = "tomato1", linewidth = 0.1, alpha = 0.5
	) +
	geom_vline(
		data = bait_var, aes(xintercept = var), 
		color = "tomato1", linewidth = 0.1, alpha = 0.5
	) +
	geom_rect( 
		xmax = max(high_var_genes_5pct$var), 
		xmin = min(high_var_genes_5pct$var),
		ymax = nrow(all_var_and_ranks),
		ymin = nrow(all_var_and_ranks) - nrow(high_var_genes_5pct),
		fill = "dodgerblue2", alpha = 0.1
	) + 
	geom_line(linewidth = 0.8) +
	labs(y = "rank",
	     x = "var(log10(TPM))",
	     caption = "Blue box = top 5% high var genes.\nRed lines = bait genes.") +
	theme_classic() +
	theme(
		text = element_text(size = 11),
		axis.text = element_text(color = "black"),
		plot.caption = element_text(hjust = 0)
	) 

# ggsave(plot = var_plot,
#        filename = "plots/MainAnalysis/gene_var_distribution.svg", 
#        height = 5, width = 5)
ggsave(plot = var_plot, 
       filename = "plots/MainAnalysis/gene_var_distribution.png", 
       height = 5, width = 5)
