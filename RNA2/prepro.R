library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(patchwork)
library(edgeR)

setwd("/scratch/fmorandi/internal/Fathima/RNA2")

paths = list()
paths$data_aco = "./pipeline_out/05_counts/counts_aco"
paths$data_mus = "./pipeline_out/05_counts/counts_mus"
paths$meta = "./meta.txt"
paths$out = "./results"
paths$objects = "./results/objects"
paths$tables = "./results/tables"

dir.create(paths$out, showWarnings = F)
dir.create(paths$objects, showWarnings = F)
dir.create(paths$tables, showWarnings = F)

##### PLOTTING SETTINGS #####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

##### LOAD DATA #####

# Load metadata
meta = read.table(paths$meta, fill=T, header=T, sep="\t")

# Load counts
counts = list()
counts$aco = fread(paths$data_aco, data.table=F)
counts$mus = fread(paths$data_mus, data.table=F)

# Separate ginfo and counts, clean names
name_conv = meta$SampleName
names(name_conv) = meta$FileName
ginfo = list()
for (org in names(counts)) {
  # Split ginfo and counts
  ginfo[[org]] = counts[[org]][, 1:6]
  counts[[org]] = counts[[org]][, -c(1:6)]
  # Short gene id
  rownames(ginfo[[org]]) = ginfo[[org]]$Geneid
  rownames(counts[[org]]) = ginfo[[org]]$Geneid
  # Clean sample names
  colnames(counts[[org]]) = name_conv[str_extract(colnames(counts[[org]]), "/([^/]+).bam", group=1)]
}

# Get QC data
qc = data.frame()
# Aco
tmp = read.table(paste0(paths$data_aco, ".summary"), header = T) %>%
  column_to_rownames("Status")
colnames(tmp) = name_conv[str_extract(colnames(tmp), "\\.([^\\.]+).bam", group=1)]
qc = rbind(qc, t(tmp[rowSums(tmp) != 0, ]))
# Mus
tmp = read.table(paste0(paths$data_mus, ".summary"), header = T) %>%
  column_to_rownames("Status")
colnames(tmp) = name_conv[str_extract(colnames(tmp), "\\.([^\\.]+).bam", group=1)]
qc = rbind(qc, t(tmp[rowSums(tmp) != 0, ]))
# Sum total reads
qc$Total = rowSums(qc)
# Merge QC into meta
meta = meta %>%
  merge(., qc, by.x="SampleName", by.y=0) %>%
  column_to_rownames("SampleName")

rm(qc)

##### QC AND REMOVE OUTLIERS #####

# Number of reads
p1 = ggplot(meta, aes(Total))+
  geom_histogram(bins=50)+
  facet_wrap(~Species)
# Plot mapping rate vs dupication rate
p2 = ggplot(meta, aes(x=Assigned/Total, y=Unassigned_MultiMapping/Total))+
  geom_point()+
  facet_wrap(~Species)
p1/p2
ggsave(paste0(paths$out, "/qc.png"), width=w, height=0.8*h, units="mm")

# Define lowQ as the one sample with lower mapping rate
meta$PassesQC = meta$Assigned/meta$Total > 0.5
blacklist = rownames(meta)[!meta$PassesQC]

for (org in names(counts)) {
  counts[[org]] = counts[[org]][, !colnames(counts[[org]]) %in% blacklist]
}

write.table(meta, paste0(paths$tables, "/qc.tsv"), quote = F, sep="\t")
meta = subset(meta, PassesQC)

##### GENE SYBOLS #####

# Get gene symbol
ginfo$aco$GeneSymbol = ginfo$aco$Geneid
ginfo$mus$GeneSymbol = mapIds(org.Mm.eg.db, keys = ginfo$mus$Geneid, keytype = "ENSEMBL", column = "SYMBOL")

# Get Entrez IDs
ginfo$aco$Entrez = mapIds(org.Mm.eg.db, keys = ginfo$aco$GeneSymbol, keytype = "SYMBOL", column = "ENTREZID")
tmp = mapIds(org.Mm.eg.db, keys = ginfo$mus$GeneSymbol, keytype = "SYMBOL", column = "ENTREZID")
tmp = as.character(tmp)
tmp[tmp == "NULL"] = NA
ginfo$mus$Entrez = tmp

# Drop genes with no symbol
ginfo$mus = ginfo$mus %>%
  drop_na(GeneSymbol)
counts$mus = counts$mus[rownames(ginfo$mus), ]

##### NORMALIZED TABLES #####

norm = list()
for (org in names(counts)) {
  # Subset meta and define design
  this_meta = meta[colnames(counts[[org]]), ]
  design = model.matrix(~0+Treatment, data = this_meta)
  # Remove low expression and normalize
  dge = DGEList(counts[[org]], samples=this_meta)
  keep = filterByExpr(dge, design=design)
  dge = dge[keep, ]
  dge = calcNormFactors(dge)
  # Save outs
  norm[[org]] = t(cpm(dge, normalized.lib.sizes=T, log=T))
}

# Previously was only taking genes making it past filtering for each species indipendently
genes_intersection = intersect(ginfo$aco[colnames(norm$aco), "GeneSymbol"],
                               ginfo$mus[colnames(norm$mus), "GeneSymbol"])

tmp = list(
  ginfo$mus[, c("Geneid", "GeneSymbol")] %>%
    distinct(GeneSymbol, .keep_all = T) %>%
    rename("mus" = "Geneid"),
  ginfo$aco[, c("Geneid", "GeneSymbol")]%>%
    distinct(GeneSymbol, .keep_all = T) %>%
    rename("aco" = "Geneid"))
conv = Reduce(function(x, y) merge(x,y, by="GeneSymbol"), tmp)
any(duplicated(conv$GeneSymbol))
any(is.na(conv))

norm$union = data.frame()
for (org in names(counts)) {
  # Subset meta and define design
  this_meta = meta[colnames(counts[[org]]), ]
  design = model.matrix(~0+Treatment, data = this_meta)
  # Keep common genes, rename and normalize
  dge = DGEList(counts[[org]], samples=this_meta)
  dge = dge[conv[, org], ]
  rownames(dge) = conv$GeneSymbol
  dge = calcNormFactors(dge)
  # Save outs
  norm$union = rbind(norm$union, t(cpm(dge, normalized.lib.sizes=T, log=T)))
}
norm$union = norm$union[rownames(meta), ]
dim(norm$union)
norm$union[1:5, 1:5]
plot(norm$union[, "Col1a1"])
plot(norm$union[, "Tbp"])
norm$intersection = norm$union[, genes_intersection]
dim(norm$intersection)

##### COMMON GENE COUNT TABLE #####

# Across species
c_samples = rownames(meta[meta$Treatment == "Acetone", ])
counts$union = data.frame()
for (org in c("aco", "mus")) {
  tmp = counts[[org]]
  # Subset to common genes and control samples
  tmp = tmp[conv[,org], colnames(tmp) %in% c_samples]
  # Rename with mouse symbols
  rownames(tmp) = conv$GeneSymbol
  # Save
  counts$union = rbind(counts$union, t(tmp))
}
counts$union = counts$union[c_samples, ]
dim(counts$union)
counts$union[1:5, 1:5]
counts$intersection = counts$union[, genes_intersection]
dim(counts$intersection)

##### SAVE DATA #####

save(counts, ginfo, meta, norm, file=paste0(paths$objects, "/prepro.Rdata"))

