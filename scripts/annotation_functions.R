library(tidyverse)

blast_get_shortName = function(blastColumn) {
	str_split_i(blastColumn, pattern = "\\^", i = 1) 
}
blast_get_fullName = function(blastColumn) {
	gsub("^.+ Full=",
	     "", blastColumn) %>% 
		gsub(pattern = ";\\^.+$", replacement = "")
}
get_COG_category = function(EggNog_matchColumn) {
	gsub(".+COG Category: ([A-Z])\\^Des.+",
	     "\\1", EggNog_matchColumn)
}

get_COG_description = function(EggNog_matchColumn) {
	str_split_i(EggNog_matchColumn, "\\^", 5) %>%
		gsub(pattern = "Description: ",
		     replacement = "")
}

shorten_PFAM = function(PFAM_column) {
	sapply(PFAM_column, \(PFAM_data) {
		if (PFAM_data == ".") {
			"."
		} else {
			str_split(string = PFAM_data,
				  pattern = "`")[[1]] %>%
				sapply(\(x) {
					p1 = str_split_i(x, "\\^", 1)
					p2 = str_split_i(x, "\\^", 3)
					paste0(p1, "^", p2)
				}) %>% paste(collapse = '`', sep = "`")	
		}
	}) %>% unname
}

shorten_GO = function(GO_column) {
	sapply(GO_column, \(GO_data) {
		if (GO_data == ".") {
			"."
		} else {
			str_split(string = GO_data,
				  pattern = "`")[[1]] %>%
				sapply(\(x) {
					str_split_fixed(x, 
							pattern = "\\^",
							n = 1)
				}) %>% paste(collapse = '`', sep = "`")	
		}
	}) %>% unname
}
