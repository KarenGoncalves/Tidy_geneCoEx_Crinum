library(tidyverse)
source("scripts/annotation_functions.R")
functional_annotation_files = 
	paste0("Data/",
	       list.files(pattern = "Crinum_powellii_report_",
	       	   path = "Data"
	       )
	)

colNames = c("geneID", "transcriptID", "blastx", 
	     "protID", "prot_coord", "blastp", 
	     "PFAM", "empty", "TmHMM", "EggNOG_ID", 
	     "GO_blast", "GO_PFAM", "KEGG_blast", 
	     "EggNOG_match", "GO_emapper", "KEGG_emapper", 
	     "Emapper_domain", 
	     "protein_seq", "peptide_seqType", "protein_length")

colNames_keep = colNames[colNames != "empty"]

funct_anno <- sapply(functional_annotation_files,
		     simplify = F, \(curFile) {
		     	read_delim(curFile, 
		     		   delim = "\t", col_names = colNames) %>%
		     		dplyr::select(all_of(colNames_keep)) %>%
		     		filter(geneID != "X1" &
		     		       	geneID != "geneID")
		     }) %>% list_rbind %>%
	mutate(blastx_short = blast_get_shortName(blastx),
	       blastp_short = blast_get_shortName(blastp),
	       blastx_fullName = blast_get_fullName(blastx),
	       blastp_fullName = blast_get_fullName(blastp),
	       COG_category = get_COG_category(EggNOG_match),
	       COG_description = get_COG_description(EggNOG_match),
	       shortPFAM = shorten_PFAM(PFAM),
	       completeness = ifelse(peptide_seqType == "complete",
				     "complete", "incomplete") %>%
	       	factor(levels = c("complete", "incomplete")),
	       shortGO_blast = shorten_GO(GO_blast),
	       shortGO_PFAM = shorten_GO(GO_PFAM),
	       GO_emapper_n = gsub(",", "`", GO_emapper),
	       KEGG_emapper_n = gsub(",", "`", KEGG_emapper)
	) %>%
	dplyr::select(
		geneID, shortPFAM, TmHMM, EggNOG_ID, 
		shortGO_blast, shortGO_PFAM, KEGG_blast, 
		GO_emapper_n, KEGG_emapper_n, 
		Emapper_domain, peptide_seqType,
		blastx_short, blastp_short, 
		blastx_fullName, blastp_fullName,
		COG_category, COG_description, 
		protein_seq, completeness, protein_length
	) %>% unique 

funct_anno %>%
	group_by(geneID) %>% 
	arrange(completeness, protein_length) %>%
	summarize(blastx_short = blastx_short[1], 
		  blastp_short = blastp_short[1], 
		  blastx_fullName = blastx_fullName[1], 
		  blastp_fullName = blastp_fullName[1],
		  shortPFAM = shortPFAM[1], 
		  TmHMM = TmHMM[1], 
		  EggNOG_ID = EggNOG_ID[1], 
		  Emapper_domain = Emapper_domain[1], 
		  COG_category = COG_category[1], 
		  COG_description = COG_description[1], 
		  shortGO_blast = shortGO_blast[1], 
		  shortGO_PFAM = shortGO_PFAM[1], 
		  GO_emapper_n = GO_emapper_n[1], 
		  KEGG_blast = KEGG_blast[1], 
		  KEGG_emapper_n = KEGG_emapper_n[1], 
		  protein_seq = protein_seq[1], 
		  completeness = completeness[1], 
		  protein_length = protein_length[1]
	) %>%
	write_delim(
		file = "output_tables/Annotation_formatted.tsv", 
		col_names = T, quote = "none", 
		delim = "\t", append = F)
