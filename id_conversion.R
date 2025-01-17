#### HANDPICKED GENES ####

##### Get gene ids #####

# # Load gene list
# gene_list = fromJSON(file=paste0(paths$out, "/../gene_list.json"))
# 
# all_genes = do.call(c, gene_list)
# names(all_genes) = NULL
# all_genes = sort(unique(all_genes))
# 
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# human_to_mouse = getLDS(
#   attributes = c("hgnc_symbol"), 
#   filters = "hgnc_symbol", 
#   values = all_genes, 
#   mart = human, 
#   attributesL = c("mgi_symbol"), 
#   martL = mouse, 
#   uniqueRows=T) 
# human_to_mouse = human_to_mouse %>%
#   distinct(HGNC.symbol, .keep_all = T) %>%
#   distinct(MGI.symbol, .keep_all = T)
# 
# all_mouse_genes = getBM(attributes = "mgi_symbol", mart=mouse)[, 1]
# all_human_genes = getBM(attributes = "hgnc_symbol", mart=human)[, 1]
# 
# # Unconvered genes
# unconverted = setdiff(all_genes, human_to_mouse$HGNC.symbol)
# # Exclude genes that were already using the mouse IDs
# tmp = data.frame(
#   HGNC.symbol = intersect(unconverted, all_mouse_genes),
#   MGI.symbol = intersect(unconverted, all_mouse_genes)
# )
# human_to_mouse = distinct(rbind(human_to_mouse, tmp))
# unconverted = setdiff(unconverted, all_mouse_genes)
# # Some genes do not convert because mice don't have them
# unconverted = setdiff(unconverted, all_human_genes)
# cat("These genes were not converted")
# cat(paste(unconverted, collapse=", "))
# 
# rownames(human_to_mouse) = human_to_mouse$HGNC.symbol
# gene_list_conv = list()
# for (gset in names(gene_list)) {
#   gene_list_conv[[gset]] = setdiff(gene_list[[gset]], unconverted)
#   gene_list_conv[[gset]] = human_to_mouse[gene_list_conv[[gset]], "MGI.symbol" ]
#   gene_list_conv[[gset]] = sort(gene_list_conv[[gset]])
# }
# tmp = toJSON(gene_list_conv, indent=1)
# write_file(tmp, file="./gene_list_conv.json")
# 
# gene_list_conv = fromJSON(file=paste0(paths$out, "/../gene_list_conv.json"))