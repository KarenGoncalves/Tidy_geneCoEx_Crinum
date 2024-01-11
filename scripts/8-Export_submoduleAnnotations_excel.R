pkgs <- c("tidyverse", "openxlsx")
sapply(pkgs, \(x) library(x, character.only = T))

outFile = list.dirs(path = "output_tables",
	  recursive = F) %>%
	grep(pattern = "output_tables/Cor", value = T)

nodeTables = 
	list.files(path = outFile,
		   pattern = "_annotations") %>%
	sapply(simplify = F, \(x) {
		nodeTable = paste0(outFile, "/", x) %>%
		 read_delim()
	})
names(nodeTables) = gsub("_annotations.tsv", "", names(nodeTables))

write.xlsx(nodeTables, 
	   file = paste0(outFile, "_subModules.xlsx"))
