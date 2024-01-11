library(tidyverse)
sampleNames = c("Root", "Light3", "Light4",
		"Dark3", "Dark4", "Leaves", 
		"Bulb", "Basal.plate")
Ex95 = sapply(simplify = F, sampleNames, \(i) {
	ExN50_data = paste0("Data/ExN50.transcript.stats_", i) %>%
		read_delim() %>%
		filter(Ex <= 100)
	
	E80N50 = ExN50_data %>% filter(Ex == 80)
	E90N50 = ExN50_data %>% filter(Ex == 90)
	E95N50 = ExN50_data %>% filter(Ex == 95)
	
	ExN50_data[-1, ] %>%
		ggplot(aes(Ex, ExN50)) +
		geom_line(linewidth = 1.25) +
		geom_point(color = "red", size = 1.5) +
		geom_vline(aes(xintercept = 95),
			   color = "blue", linewidth = 1.5) +
		annotate(geom = "text",
			 y = 1250, x = 80,
			 label = paste0("ExN50 = ", E95N50$ExN50, "\n",
			 	       "N genes = ", E95N50$num_genes),
			 color = "blue"
		) +
		ggtitle(i) + 
		theme_bw() +
		theme(axis.text = element_text(size = 13),
		      axis.title = element_text(face = "bold", size = 14))
	
	
})

par(mfrow=c(4,4))	
pdf("plots/Ex95N50_perSample.pdf")
Ex95
dev.off()
