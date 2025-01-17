library(tidyverse)
library(edgeR)
library(data.table)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggVennDiagram)
library(biomaRt)
library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(ComplexHeatmap)
library(gtools)
library(colorspace)
library(BiocParallel)
library(gtools)
library(ggrepel)
library(GO.db)

setwd("/scratch/fmorandi/internal/Fathima/RNA1")

paths = list()
paths$meta = "./meta.txt"
paths$out = "./results_union"
paths$objects = "./results_union/objects"
paths$tables = "./results_union/tables"
paths$pathways = "./results_union/pathways"
paths$custom_gsets = "./results_union/custom_gsets"
dir.create(paths$pathways, showWarnings = F)
dir.create(paths$custom_gsets, showWarnings = F)

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### FUNCTIONS ####

my_gse = function(table, logFC_column, orgdb) {
  # Make gene list
  gene_list = table[[logFC_column]]
  names(gene_list) = table$Entrez
  gene_list = sort(gene_list, decreasing = TRUE)
  # Run GSEA
  set.seed(1337)
  res = gseGO(
    geneList=gene_list,
    ont ="BP",
    keyType = "ENTREZID",
    verbose = TRUE,
    OrgDb = orgdb,
    pvalueCutoff = 1.1,
    BPPARAM = SerialParam())
  return(res)
}

my_volcano = function(table, col, v1, v2, title=NULL) {
  tmp = table %>%
    group_by_at(col) %>%
    summarize(n=n()) %>%
    mutate(x = c(0, -3, 3))
  p = ggplot(table, aes(x=.data[[v1]], y=-log10(.data[[v2]]), color=.data[[col]]))+
    geom_point(size=0.1)+
    lims(x=c(-4, 4), y=c(0,50))+
    ggtitle(title)+
    geom_text(data=tmp, aes(x=x, y=40, label=n))
  if (!is.null(title)) p = p + ggtitle(title)
  return(p)
}

convert_ids_in_string = function(df, entrez, symbol) {
  conversions = symbol
  names(conversions) = entrez
  conversions = conversions[!is.na(names(conversions))]
  for (i in 1:nrow(df)) {
    this_str = df$core_enrichment[i]
    old_ids = unlist(strsplit(this_str, "/"))
    new_ids = conversions[old_ids]
    df[i, "core_enrichment"] = paste(new_ids, collapse="/")
  }
  return(df)
}

#### LOAD DATA ####

load(paste0(paths$objects, "/prepro.Rdata"))

#### REMOVE HSA ####

meta = subset(meta, Species != "HSA")
norm$hsa = NULL
counts$hsa = NULL
ginfo$hsa = NULL
counts$all = counts$all[rownames(counts$all) %in% rownames(meta), ]
norm$all = norm$all[rownames(norm$all) %in% rownames(meta), ]

#### PCA ####

pca = list()

tmp = prcomp(norm$all, center=T, scale.=T)
pca$all = merge(meta, tmp$x[,c("PC1", "PC2")], by=0)
tmp = prcomp(norm$aco, center=T, scale.=T)
pca$aco = merge(meta, tmp$x[,c("PC1", "PC2")], by=0)
tmp = prcomp(norm$mus, center=T, scale.=T)
pca$mus = merge(meta, tmp$x[,c("PC1", "PC2")], by=0)
tmp = prcomp(norm$musw, center=T, scale.=T)
pca$musw = merge(meta, tmp$x[,c("PC1", "PC2")], by=0)

# All samples PCA
p1 = ggplot(pca$all, aes(x=PC1, y=PC2, color=Species))+
  geom_point(size=0.2)
p2 = ggplot(pca$all, aes(x=PC1, y=PC2, color=Timepoint))+
  geom_point(size=0.2)
ps = align_patches(p1, p2)
pdf(paste0(paths$out, "/pca_all.pdf"), width=0.5*w_in, height=0.25*h_in)
ps
dev.off()

# Acomys only PCA
p1=ggplot(pca$aco, aes(x=PC1, y=PC2, color=Treatment))+
  geom_point()+
  ggtitle("Acomys")
p2=ggplot(pca$aco, aes(x=PC1, y=PC2, color=Timepoint))+
  geom_point()+
  ggtitle("Acomys")
p3=ggplot(pca$aco, aes(x=PC1, y=PC2, color=CellLine))+
  geom_point()+
  ggtitle("Acomys")
ps = align_patches(p1, p2, p3)
pdf(paste0(paths$out, "/pca_aco.pdf"), width=0.5*w_in, height=0.25*h_in)
ps
dev.off()

# Mouse only PCA
p1=ggplot(pca$mus, aes(x=PC1, y=PC2, color=Treatment))+
  geom_point()+
  ggtitle("Mouse")
p2=ggplot(pca$mus, aes(x=PC1, y=PC2, color=Timepoint))+
  geom_point()+
  ggtitle("Mouse")
p3=ggplot(pca$mus, aes(x=PC1, y=PC2, color=CellLine))+
  geom_point()+
  ggtitle("Mouse")
ps = align_patches(p1, p2, p3)
pdf(paste0(paths$out, "/pca_mus.pdf"), width=0.5*w_in, height=0.25*h_in)
ps
dev.off()

# Wild mouse only PCA
p1=ggplot(pca$musw, aes(x=PC1, y=PC2, color=Treatment))+
  geom_point()+
  ggtitle("Mouse")
p2=ggplot(pca$musw, aes(x=PC1, y=PC2, color=Timepoint))+
  geom_point()+
  ggtitle("Mouse")
p3=ggplot(pca$musw, aes(x=PC1, y=PC2, color=CellLine))+
  geom_point()+
  ggtitle("Mouse")
ps = align_patches(p1, p2, p3)
pdf(paste0(paths$out, "/pca_musw.pdf"), width=0.5*w_in, height=0.25*h_in)
ps
dev.off()

#### DIFFERENTIAL EXPRESSION ####

##### Treatment #####

des_treatment = list()

# Within species and timepoint, across treatments
for (org in c("aco", "mus", "musw")) {
  this_meta = subset(meta, tolower(Species) == org)
  for (t in unique(this_meta$Timepoint)) {
    # Subset meta and define design
    this_meta2 = subset(this_meta, Timepoint == t)
    design = model.matrix(~0+Treatment+CellLine, data = this_meta2)
    # Prepare dge
    dge = DGEList(counts[[org]][, rownames(this_meta2)], samples=this_meta2)
    # keep = filterByExpr(dge, design=design) # not when using union
    # dge = dge[keep, ]
    dge = calcNormFactors(dge)
    # Fit model
    dge = estimateDisp(dge, design)
    fit = glmFit(dge, design)
    # Get 10 vs C results
    conts = makeContrasts("Treatment10Gy - TreatmentCon", levels=design)
    res10vsC = as.data.frame(glmLRT(fit, contrast=conts))
    res10vsC$PAdj = p.adjust(res10vsC$PValue, method="BH")
    # Get 20 vs C results
    conts = makeContrasts("Treatment20Gy - TreatmentCon", levels=design)
    res20vsC = as.data.frame(glmLRT(fit, contrast=conts))
    res20vsC$PAdj = p.adjust(res20vsC$PValue, method="BH")
    # Save results with ginfo
    res = merge(res10vsC, res20vsC, by=0, suffixes = c("_10vsC","_20vsC")) %>%
      dplyr::select(-starts_with("LR"))
    des_treatment[[paste0(org, "_", t)]] = merge(ginfo[[org]], res, by.x=0, by.y="Row.names") %>%
      column_to_rownames("Row.names")
  }
}

rm(res10vsC, res20vsC)

##### Baseline #####

des_baseline = list()

# --- T1 ---
# Subset meta and define design
this_meta = meta[rownames(counts$all), ]
this_meta = subset(this_meta, Timepoint == "24h")
design = model.matrix(~0+Species, data = this_meta)
# Normalize
dge = DGEList(t(counts$all[rownames(this_meta), ]), samples=this_meta)
# keep = filterByExpr(dge, design=design) # not when using union
# dge = dge[keep, ]
dge = calcNormFactors(dge)
ggplot(dge$samples, aes(x=Species, y=norm.factors))+
  geom_violin()
# <!> Therefore it seems reasonabe to normalize with TMM, most are close to 1
# Fit model
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
# Get Aco vs Mus results
conts = makeContrasts("SpeciesACO - SpeciesMUS", levels=design)
res_aco_vs_mus = as.data.frame(glmLRT(fit, contrast=conts))
res_aco_vs_mus$PAdj = p.adjust(res_aco_vs_mus$PValue, method="BH")
# Get Aco vs MusW results
conts = makeContrasts("SpeciesACO - SpeciesMUSW", levels=design)
res_aco_vs_musw = as.data.frame(glmLRT(fit, contrast=conts))
res_aco_vs_musw$PAdj = p.adjust(res_aco_vs_musw$PValue, method="BH")
# Combine T1 results
des_baseline$T1 = merge(res_aco_vs_mus, res_aco_vs_musw, by=0, suffixes=c("_aco_vs_mus", "_aco_vs_musw")) %>%
  dplyr::select(-starts_with("LR")) %>%
  column_to_rownames("Row.names")
colSums(des_baseline$T1[c("PAdj_aco_vs_mus", "PAdj_aco_vs_musw")] < 0.05)

# --- T2 ---
# Subset meta and define design
this_meta = meta[rownames(counts$all), ]
this_meta = subset(this_meta, Timepoint == "12d")
design = model.matrix(~0+Species, data = this_meta)
# Normalize
dge = DGEList(t(counts$all[rownames(this_meta), ]), samples=this_meta)
# keep = filterByExpr(dge, design=design) # not when using union
# dge = dge[keep, ]
dge = calcNormFactors(dge)
ggplot(dge$samples, aes(x=Species, y=norm.factors))+
  geom_violin()
# <!> Therefore it seems reasonabe to normalize with TSS, most are close to 1
# Fit model
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
# Get Aco vs Mus results
conts = makeContrasts("SpeciesACO - SpeciesMUS", levels=design)
res_aco_vs_mus = as.data.frame(glmLRT(fit, contrast=conts))
res_aco_vs_mus$PAdj = p.adjust(res_aco_vs_mus$PValue, method="BH")
# Get Aco vs MusW results
conts = makeContrasts("SpeciesACO - SpeciesMUSW", levels=design)
res_aco_vs_musw = as.data.frame(glmLRT(fit, contrast=conts))
res_aco_vs_musw$PAdj = p.adjust(res_aco_vs_musw$PValue, method="BH")
# Combine T1 results
des_baseline$T2 = merge(res_aco_vs_mus, res_aco_vs_musw, by=0, suffixes=c("_aco_vs_mus", "_aco_vs_musw")) %>%
  dplyr::select(-starts_with("LR"))%>%
  column_to_rownames("Row.names")
colSums(des_baseline$T2[c("PAdj_aco_vs_mus", "PAdj_aco_vs_musw")] < 0.05)

rm(res_aco_vs_mus, res_aco_vs_musw)

##### Define DEGs #####

th_sig_gene = 0.05
th_logfc = 1

# --- Treatment ---
for (org in names(des_treatment)) {
  des_treatment[[org]] = des_treatment[[org]] %>%
    mutate(sig_10vsC = abs(logFC_10vsC) > th_logfc & PAdj_10vsC < th_sig_gene) %>%
    mutate(sig_20vsC = abs(logFC_20vsC) > th_logfc & PAdj_20vsC < th_sig_gene) %>%
    mutate(up_10vsC = sig_10vsC & logFC_10vsC > 0) %>%
    mutate(up_20vsC = sig_20vsC & logFC_20vsC > 0) %>%
    mutate(down_10vsC = sig_10vsC & logFC_10vsC < 0) %>%
    mutate(down_20vsC = sig_20vsC & logFC_20vsC < 0) %>%
    mutate(sig_col10 = interaction(sig_10vsC, up_10vsC)) %>%
    mutate(sig_col20 = interaction(sig_20vsC, up_20vsC)) %>%
    mutate(sig_col10 = fct_recode(sig_col10, 
                                  "Up" = "TRUE.TRUE",
                                  "Down" = "TRUE.FALSE",
                                  "NotSig" = "FALSE.FALSE")) %>%
    mutate(sig_col20 = fct_recode(sig_col20, 
                                  "Up" = "TRUE.TRUE",
                                  "Down" = "TRUE.FALSE",
                                  "NotSig" = "FALSE.FALSE"))
}

# --- Baseline ---
des_baseline$T1 = des_baseline$T1 %>%
  mutate(sig_aco_vs_mus = abs(logFC_aco_vs_mus) > th_logfc & PAdj_aco_vs_mus < th_sig_gene) %>%
  mutate(sig_aco_vs_musw = abs(logFC_aco_vs_musw) > th_logfc & PAdj_aco_vs_musw < th_sig_gene) %>%
  mutate(up_aco_vs_mus = sig_aco_vs_mus & logFC_aco_vs_mus > 0) %>%
  mutate(up_aco_vs_musw = sig_aco_vs_musw & logFC_aco_vs_musw > 0) %>%
  mutate(down_aco_vs_mus = sig_aco_vs_mus & logFC_aco_vs_mus < 0) %>%
  mutate(down_aco_vs_musw = sig_aco_vs_musw & logFC_aco_vs_musw < 0) %>%
  mutate(sig_aco_vs_mus = interaction(sig_aco_vs_mus, up_aco_vs_mus)) %>%
  mutate(sig_aco_vs_musw = interaction(sig_aco_vs_musw, up_aco_vs_musw)) %>%
  mutate(sig_aco_vs_mus = fct_recode(sig_aco_vs_mus, "Up" = "TRUE.TRUE", "Down" = "TRUE.FALSE", "NotSig" = "FALSE.FALSE")) %>%
  mutate(sig_aco_vs_musw = fct_recode(sig_aco_vs_musw, "Up" = "TRUE.TRUE", "Down" = "TRUE.FALSE", "NotSig" = "FALSE.FALSE"))
des_baseline$T2 = des_baseline$T2 %>%
  mutate(sig_aco_vs_mus = abs(logFC_aco_vs_mus) > th_logfc & PAdj_aco_vs_mus < th_sig_gene) %>%
  mutate(sig_aco_vs_musw = abs(logFC_aco_vs_musw) > th_logfc & PAdj_aco_vs_musw < th_sig_gene) %>%
  mutate(up_aco_vs_mus = sig_aco_vs_mus & logFC_aco_vs_mus > 0) %>%
  mutate(up_aco_vs_musw = sig_aco_vs_musw & logFC_aco_vs_musw > 0) %>%
  mutate(down_aco_vs_mus = sig_aco_vs_mus & logFC_aco_vs_mus < 0) %>%
  mutate(down_aco_vs_musw = sig_aco_vs_musw & logFC_aco_vs_musw < 0) %>%
  mutate(sig_aco_vs_mus = interaction(sig_aco_vs_mus, up_aco_vs_mus)) %>%
  mutate(sig_aco_vs_musw = interaction(sig_aco_vs_musw, up_aco_vs_musw)) %>%
  mutate(sig_aco_vs_mus = fct_recode(sig_aco_vs_mus, "Up" = "TRUE.TRUE", "Down" = "TRUE.FALSE", "NotSig" = "FALSE.FALSE")) %>%
  mutate(sig_aco_vs_musw = fct_recode(sig_aco_vs_musw, "Up" = "TRUE.TRUE", "Down" = "TRUE.FALSE", "NotSig" = "FALSE.FALSE"))

##### Plot volcanos #####

# --- Treatment: T1 ---
p1 = my_volcano(des_treatment$aco_24h, col="sig_col10", v1="logFC_10vsC", v2="PValue_10vsC", "Acomys")
p2 = my_volcano(des_treatment$mus_24h, col="sig_col10", v1="logFC_10vsC", v2="PValue_10vsC", "Mouse")
p3 = my_volcano(des_treatment$musw_24h, col="sig_col10", v1="logFC_10vsC", v2="PValue_10vsC", "Wild Mouse")

p4 = my_volcano(des_treatment$aco_24h, col="sig_col20", v1="logFC_20vsC", v2="PValue_20vsC")
p5 = my_volcano(des_treatment$mus_24h, col="sig_col20", v1="logFC_20vsC", v2="PValue_20vsC")
p6 = my_volcano(des_treatment$musw_24h, col="sig_col20", v1="logFC_20vsC", v2="PValue_20vsC")

p=p1+p2+p3+p4+p5+p6+
  plot_layout(guides = "collect", nrow=2) &
  scale_color_manual(values=c("#999999", "#5555cc", "#cc5555")) &
  guides(color="none")
ggsave(paste0(paths$out, "/volcanos_treatment_24h.png"), plot=p, width=w, height=0.6*h, units="mm")

# --- Treatment: T2 ---
p1 = my_volcano(des_treatment$aco_12d, col="sig_col10", v1="logFC_10vsC", v2="PValue_10vsC", "Acomys")
p2 = my_volcano(des_treatment$mus_12d, col="sig_col10", v1="logFC_10vsC", v2="PValue_10vsC", "Mouse")
p3 = my_volcano(des_treatment$musw_12d, col="sig_col10", v1="logFC_10vsC", v2="PValue_10vsC", "Wild Mouse")

p4 = my_volcano(des_treatment$aco_12d, col="sig_col20", v1="logFC_20vsC", v2="PValue_20vsC")
p5 = my_volcano(des_treatment$mus_12d, col="sig_col20", v1="logFC_20vsC", v2="PValue_20vsC")
p6 = my_volcano(des_treatment$musw_12d, col="sig_col20", v1="logFC_20vsC", v2="PValue_20vsC")

p=p1+p2+p3+p4+p5+p6+
  plot_layout(guides = "collect", nrow=2) &
  scale_color_manual(values=c("#999999", "#5555cc", "#cc5555")) &
  guides(color="none")
ggsave(paste0(paths$out, "/volcanos_treatment_12d.png"), plot=p, width=w, height=0.6*h, units="mm")

# --- Baseline: T1 ---
p1 = my_volcano(des_baseline$T1, col="sig_aco_vs_mus", v1="logFC_aco_vs_mus", v2="PValue_aco_vs_mus", "Acomys vs Mouse")
p2 = my_volcano(des_baseline$T1, col="sig_aco_vs_musw", v1="logFC_aco_vs_musw", v2="PValue_aco_vs_musw", "Acomys vs Wild Mouse")

p=p1+p2+
  plot_layout(guides = "collect", nrow=1) &
  scale_color_manual(values=c("#999999", "#5555cc", "#cc5555")) &
  guides(color="none")
ggsave(paste0(paths$out, "/volcanos_baseline_24h.png"), plot=p, width=w, height=0.4*h, units="mm")

# --- Baseline: T2 ---
p1 = my_volcano(des_baseline$T2, col="sig_aco_vs_mus", v1="logFC_aco_vs_mus", v2="PValue_aco_vs_mus", "Acomys vs Mouse")
p2 = my_volcano(des_baseline$T2, col="sig_aco_vs_musw", v1="logFC_aco_vs_musw", v2="PValue_aco_vs_musw", "Acomys vs Wild Mouse")

p=p1+p2+
  plot_layout(guides = "collect", nrow=1) &
  scale_color_manual(values=c("#999999", "#5555cc", "#cc5555")) &
  guides(color="none")
ggsave(paste0(paths$out, "/volcanos_baseline_12d.png"), plot=p, width=w, height=0.4*h, units="mm")

##### Combined table #####

des_combined = list()

# --- T1 ---
# Bring all DE results side by side on common genes
res = data.frame()
res = des_treatment$aco_24h %>%
  dplyr::select(GeneSymbol, logFC_10vsC:PAdj_20vsC) %>%
  mutate(Species = "aco") %>%
  rbind(res, .)
res = des_treatment$mus_24h %>%
  dplyr::select(GeneSymbol, logFC_10vsC:PAdj_20vsC) %>%
  mutate(Species = "mus") %>%
  rbind(res)
res = des_treatment$musw_24h %>%
  dplyr::select(GeneSymbol, logFC_10vsC:PAdj_20vsC) %>%
  mutate(Species = "musw") %>%
  rbind(res, .)
rownames(res) = NULL

# Split into 2 tables for 10gy and 20gy
des_combined$res24h_10vsC = res %>%
  distinct(GeneSymbol, Species, .keep_all = T) %>%
  dplyr::select(-ends_with("20vsC")) %>%
  pivot_wider(names_from=Species, values_from=ends_with("10vsC")) %>%
  mutate(Entrez = mapIds(org.Mm.eg.db, keys = GeneSymbol, keytype = "SYMBOL", column = "ENTREZID")) %>%
  column_to_rownames("GeneSymbol")
des_combined$res24h_20vsC = res %>%
  distinct(GeneSymbol, Species, .keep_all = T) %>%
  dplyr::select(-ends_with("10vsC")) %>%
  pivot_wider(names_from=Species, values_from=ends_with("20vsC")) %>%
  mutate(Entrez = mapIds(org.Mm.eg.db, keys = GeneSymbol, keytype = "SYMBOL", column = "ENTREZID")) %>%
  column_to_rownames("GeneSymbol")

# Add baseline results
des_combined$res24h_10vsC = des_baseline$T1 %>%
  dplyr::select(logFC_aco_vs_mus:PAdj_aco_vs_musw) %>%
  cbind(des_combined$res24h_10vsC[rownames(des_baseline$T1), ], .)
des_combined$res24h_20vsC = des_baseline$T1 %>%
  dplyr::select(logFC_aco_vs_mus:PAdj_aco_vs_musw) %>%
  cbind(des_combined$res24h_20vsC[rownames(des_baseline$T1), ], .)

# --- T2 ---
# Bring all DE results side by side on common genes
res = data.frame()
res = des_treatment$aco_12d %>%
  dplyr::select(GeneSymbol, logFC_10vsC:PAdj_20vsC) %>%
  mutate(Species = "aco") %>%
  rbind(res, .)
res = des_treatment$mus_12d %>%
  dplyr::select(GeneSymbol, logFC_10vsC:PAdj_20vsC) %>%
  mutate(Species = "mus") %>%
  rbind(res)
res = des_treatment$musw_12d %>%
  dplyr::select(GeneSymbol, logFC_10vsC:PAdj_20vsC) %>%
  mutate(Species = "musw") %>%
  rbind(res, .)
rownames(res) = NULL

# Split into 2 tables for 10gy and 20gy
des_combined$res12d_10vsC = res %>%
  distinct(GeneSymbol, Species, .keep_all = T) %>%
  dplyr::select(-ends_with("20vsC")) %>%
  pivot_wider(names_from=Species, values_from=ends_with("10vsC")) %>%
  mutate(Entrez = mapIds(org.Mm.eg.db, keys = GeneSymbol, keytype = "SYMBOL", column = "ENTREZID")) %>%
  column_to_rownames("GeneSymbol")
des_combined$res12d_20vsC = res %>%
  distinct(GeneSymbol, Species, .keep_all = T) %>%
  dplyr::select(-ends_with("10vsC")) %>%
  pivot_wider(names_from=Species, values_from=ends_with("20vsC")) %>%
  mutate(Entrez = mapIds(org.Mm.eg.db, keys = GeneSymbol, keytype = "SYMBOL", column = "ENTREZID")) %>%
  column_to_rownames("GeneSymbol")

# Add baseline results
colnames(des_baseline$T2)
des_combined$res12d_10vsC = des_baseline$T2 %>%
  dplyr::select(logFC_aco_vs_mus:PAdj_aco_vs_musw) %>%
  cbind(des_combined$res12d_10vsC[rownames(des_baseline$T2), ], .)
des_combined$res12d_20vsC = des_baseline$T2 %>%
  dplyr::select(logFC_aco_vs_mus:PAdj_aco_vs_musw) %>%
  cbind(des_combined$res12d_20vsC[rownames(des_baseline$T2), ], .)

#### CHECKPOINT 1 ####

# save.image(file=paste0(paths$objects, "/checkpoint1.Rdata"))
load(paste0(paths$objects, "/checkpoint1.Rdata"))

#### GSEA: TREATMENT ####

##### Run #####

# gse_treatment = list()
# 
# # --- Mouse ---
# gse_treatment$mus24h_10gy = my_gse(des_treatment[["mus_24h"]], "logFC_10vsC", orgdb = org.Mm.eg.db)
# gse_treatment$mus24h_20gy = my_gse(des_treatment[["mus_24h"]], "logFC_20vsC", orgdb = org.Mm.eg.db)
# gse_treatment$mus12d_10gy = my_gse(des_treatment[["mus_12d"]], "logFC_10vsC", orgdb = org.Mm.eg.db)
# gse_treatment$mus12d_20gy = my_gse(des_treatment[["mus_12d"]], "logFC_20vsC", orgdb = org.Mm.eg.db)
# 
# # --- Wild Mouse ---
# gse_treatment$musw24h_10gy = my_gse(des_treatment[["musw_24h"]], "logFC_10vsC", orgdb = org.Mm.eg.db)
# gse_treatment$musw24h_20gy = my_gse(des_treatment[["musw_24h"]], "logFC_20vsC", orgdb = org.Mm.eg.db)
# gse_treatment$musw12d_10gy = my_gse(des_treatment[["musw_12d"]], "logFC_10vsC", orgdb = org.Mm.eg.db)
# gse_treatment$musw12d_20gy = my_gse(des_treatment[["musw_12d"]], "logFC_20vsC", orgdb = org.Mm.eg.db)
# 
# # --- Acomys ---
# length(intersect(ginfo$aco$GeneSymbol, ginfo$mus$GeneSymbol))
# length(setdiff(ginfo$aco$GeneSymbol, ginfo$mus$GeneSymbol)) # Most genes share names with mice
# 
# gse_treatment$aco24h_10gy = my_gse(des_treatment[["aco_24h"]], "logFC_10vsC", orgdb = org.Mm.eg.db)
# gse_treatment$aco24h_20gy = my_gse(des_treatment[["aco_24h"]], "logFC_20vsC", orgdb = org.Mm.eg.db)
# gse_treatment$aco12d_10gy = my_gse(des_treatment[["aco_12d"]], "logFC_10vsC", orgdb = org.Mm.eg.db)
# gse_treatment$aco12d_20gy = my_gse(des_treatment[["aco_12d"]], "logFC_20vsC", orgdb = org.Mm.eg.db)
# 
# # Save results
# save(gse_treatment, file=paste0(paths$objects, "/gsea_treatment.Rdata"))

load(paste0(paths$objects, "/gsea_treatment.Rdata"))

##### Dotplots #####

# --- T1 ---
ps = list()
ps[[1]]=gse_treatment$aco24h_10gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys 10gy vs C")
ps[[2]]=gse_treatment$aco24h_20gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys 20gy vs C")
ps[[3]]=gse_treatment$mus24h_10gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Mouse 10gy vs C")
ps[[4]]=gse_treatment$mus24h_20gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Mouse 20gy vs C")
ps[[5]]=gse_treatment$musw24h_10gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Wild Mouse 10gy vs C")
ps[[6]]=gse_treatment$musw24h_20gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Wild Mouse 20gy vs C")
ps = align_patches(ps)
pdf(paste0(paths$out, "/gsea_treatment24h.pdf"), width=w_in*2, height=h_in)
ps
dev.off()

# --- T2 ---
ps = list()
ps[[1]]=gse_treatment$aco12d_10gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys 10gy vs C")
ps[[2]]=gse_treatment$aco12d_20gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys 20gy vs C")
ps[[3]]=gse_treatment$mus12d_10gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Mouse 10gy vs C")
ps[[4]]=gse_treatment$mus12d_20gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Mouse 20gy vs C")
ps[[5]]=gse_treatment$musw12d_10gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Wild Mouse 10gy vs C")
ps[[6]]=gse_treatment$musw12d_20gy %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Wild Mouse 20gy vs C")
ps = align_patches(ps)
pdf(paste0(paths$out, "/gsea_treatment12d.pdf"), width=w_in*2, height=h_in)
ps
dev.off()

##### Comparison #####

gse_combined = list()
gse_any_sig = list()

# --- T1 ---
gse_combined$T1 = data.frame()
for (comp in names(gse_treatment)[grepl("24h", names(gse_treatment))]) {
  tmp = data.frame(gse_treatment[[comp]]) %>%
    dplyr::select("ID", "Description", "NES", "p.adjust") %>%
    mutate(Comparison = comp)
  gse_combined$T1 = rbind(gse_combined$T1, tmp)
}
rownames(gse_combined$T1) = NULL
gse_combined$T1 = pivot_wider(gse_combined$T1, names_from=Comparison, values_from=c(NES, p.adjust))
# LogFC similarity across species and treatments
data.frame(cor(gse_combined$T1[, 3:8], use = "pairwise.complete.obs")) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols=-Var1, names_to = "Var2") %>%
  ggplot(., aes(x=Var1, y=Var2, fill=value, label=round(value,2)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_viridis_c()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(fill="Correlation")
ggsave(paste0(paths$out, "/gsea_response_similarity24h.png"), width=w, height=0.6*h, units="mm")

# Make min pval columns
colnames(gse_combined$T1) = str_replace_all(colnames(gse_combined$T1), "p.adjust", "padj")
gse_combined$T1$padj_mus = pmin(gse_combined$T1$padj_mus24h_10gy, gse_combined$T1$padj_mus24h_20gy, na.rm = T)
gse_combined$T1$padj_musw = pmin(gse_combined$T1$padj_musw24h_10gy, gse_combined$T1$padj_musw24h_20gy, na.rm = T)
gse_combined$T1$padj_aco = pmin(gse_combined$T1$padj_aco24h_10gy, gse_combined$T1$padj_aco24h_20gy, na.rm = T)
gse_combined$T1$padj_min = pmin(gse_combined$T1$padj_mus, gse_combined$T1$padj_musw, gse_combined$T1$padj_aco, na.rm=T)

gse_any_sig$T1 = subset(gse_combined$T1, padj_min < 0.01)
# gse_any_sig$T1 = arrange(gse_any_sig$T1, ID) # Previous ordering
gse_any_sig$T1 = arrange(gse_any_sig$T1, rowMeans(gse_any_sig$T1[, grepl("NES", colnames(gse_any_sig$T1))], na.rm = T))

n_per_page = 40
nes_range = c(
  min(gse_any_sig$T1[, grepl("NES", colnames(gse_any_sig$T1))], na.rm=T),
  max(gse_any_sig$T1[, grepl("NES", colnames(gse_any_sig$T1))], na.rm=T))
ps = list()
for (p in 1:ceiling(nrow(gse_any_sig$T1)/n_per_page)) {
  page_first = (p-1)*n_per_page+1
  page_ps = list()
  tmp = gse_any_sig$T1[page_first:min(page_first+n_per_page-1, nrow(gse_any_sig$T1)), 1:14] %>%
    mutate(ID = factor(ID, levels=ID)) %>%
    pivot_longer(
      names_to = c("Variable", "Species", "Timepoint", "Treatment"),
      names_pattern = "(.*)_(hsa|mus|musw|aco)(24h)_(10gy|20gy)",
      cols = -c(ID, Description)) %>%
    mutate(Comparison = paste(Species, Treatment, sep = "_")) %>%
    pivot_wider(names_from=Variable) %>%
    mutate(Sig = stars.pval(padj)) %>%
    mutate(Title = substr(Description, 1, 80)) %>%
    mutate(Title = fct_reorder(Title, as.numeric(ID)))
  ps[[p]] = ggplot(tmp, aes(x=Comparison, y=Title, fill=NES, label=Sig))+
    geom_tile()+
    geom_text()+
    scale_fill_gradient2(low="blue", high="red", limits=nes_range)+
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle=45, hjust=1))
}

ps = align_patches(ps)
pdf(paste0(paths$out, "/gsea_comparison24h.pdf"), width=w_in*1.2, height=h_in)
ps
dev.off()

# --- T2 ---
gse_combined$T2 = data.frame()
for (comp in names(gse_treatment)[grepl("12d", names(gse_treatment))]) {
  tmp = data.frame(gse_treatment[[comp]]) %>%
    dplyr::select("ID", "Description", "NES", "p.adjust") %>%
    mutate(Comparison = comp)
  gse_combined$T2 = rbind(gse_combined$T2, tmp)
}
rownames(gse_combined$T2) = NULL
gse_combined$T2 = pivot_wider(gse_combined$T2, names_from=Comparison, values_from=c(NES, p.adjust))
# LogFC similarity across species and treatments
data.frame(cor(gse_combined$T2[, 3:8], use = "pairwise.complete.obs")) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols=-Var1, names_to = "Var2") %>%
  ggplot(., aes(x=Var1, y=Var2, fill=value, label=round(value,2)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_viridis_c()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(fill="Correlation")
ggsave(paste0(paths$out, "/gsea_response_similarity12d.png"), width=w, height=0.6*h, units="mm")

# Make min pval columns
colnames(gse_combined$T2) = str_replace_all(colnames(gse_combined$T2), "p.adjust", "padj")
gse_combined$T2$padj_mus = pmin(gse_combined$T2$padj_mus12d_10gy, gse_combined$T2$padj_mus12d_20gy, na.rm = T)
gse_combined$T2$padj_musw = pmin(gse_combined$T2$padj_musw12d_10gy, gse_combined$T2$padj_musw12d_20gy, na.rm = T)
gse_combined$T2$padj_aco = pmin(gse_combined$T2$padj_aco12d_10gy, gse_combined$T2$padj_aco12d_20gy, na.rm = T)
gse_combined$T2$padj_min = pmin(gse_combined$T2$padj_mus, gse_combined$T2$padj_musw, gse_combined$T2$padj_aco, na.rm=T)

gse_any_sig$T2 = subset(gse_combined$T2, padj_min < 0.01)
gse_any_sig$T2 = arrange(gse_any_sig$T2, rowMeans(gse_any_sig$T2[, grepl("NES", colnames(gse_any_sig$T2))], na.rm = T))

n_per_page = 40
nes_range = c(
  min(gse_any_sig$T2[, grepl("NES", colnames(gse_any_sig$T2))], na.rm=T),
  max(gse_any_sig$T2[, grepl("NES", colnames(gse_any_sig$T2))], na.rm=T))
ps = list()
for (p in 1:ceiling(nrow(gse_any_sig$T2)/n_per_page)) {
  page_first = (p-1)*n_per_page+1
  page_ps = list()
  tmp = gse_any_sig$T2[page_first:min(page_first+n_per_page-1, nrow(gse_any_sig$T2)), 1:14] %>%
    mutate(ID = factor(ID, levels=ID)) %>%
    pivot_longer(
      names_to = c("Variable", "Species", "Timepoint", "Treatment"),
      names_pattern = "(.*)_(mus|musw|aco)(12d)_(10gy|20gy)",
      cols = -c(ID, Description)) %>%
    mutate(Comparison = paste(Species, Treatment, sep = "_")) %>%
    pivot_wider(names_from=Variable) %>%
    mutate(Sig = stars.pval(padj)) %>%
    mutate(Title = substr(Description, 1, 80)) %>%
    mutate(Title = fct_reorder(Title, as.numeric(ID)))
  ps[[p]] = ggplot(tmp, aes(x=Comparison, y=Title, fill=NES, label=Sig))+
    geom_tile()+
    geom_text()+
    scale_fill_gradient2(low="blue", high="red", limits=nes_range)+
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle=45, hjust=1))
}

ps = align_patches(ps)
pdf(paste0(paths$out, "/gsea_comparison12d.pdf"), width=w_in*1.2, height=h_in)
ps
dev.off()

#### GSEA: BASELINE ####

##### Run #####

# Run gsea on pairwise logFC baseline differences
des_baseline$T1 = des_baseline$T1 %>%
  rownames_to_column("Row.names") %>%
  mutate(Entrez = mapIds(org.Mm.eg.db, keys = Row.names, keytype = "SYMBOL", column = "ENTREZID"))%>%
  column_to_rownames("Row.names")
des_baseline$T2 = des_baseline$T2 %>%
  rownames_to_column("Row.names") %>%
  mutate(Entrez = mapIds(org.Mm.eg.db, keys = Row.names, keytype = "SYMBOL", column = "ENTREZID"))%>%
  column_to_rownames("Row.names")

# gse_baseline = list()
# 
# # --- Aco vs Mus
# gse_baseline$aco_vs_mus_24h = my_gse(des_baseline$T1, "logFC_aco_vs_mus", org.Mm.eg.db)
# gse_baseline$aco_vs_mus_12d = my_gse(des_baseline$T2, "logFC_aco_vs_mus", org.Mm.eg.db)
# 
# # --- Aco vs MusW
# gse_baseline$aco_vs_musw_24h = my_gse(des_baseline$T1, "logFC_aco_vs_musw", org.Mm.eg.db)
# gse_baseline$aco_vs_musw_12d = my_gse(des_baseline$T2, "logFC_aco_vs_musw", org.Mm.eg.db)
# 
# save(gse_baseline, file=paste0(paths$objects, "/gsea_baseline.Rdata"))
load(paste0(paths$objects, "/gsea_baseline.Rdata"))

##### Dotplots #####

facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Mouse"
)
p1=gse_baseline$aco_vs_mus_24h %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Mouse Baseline")

facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Wild Mouse"
)
p2=gse_baseline$aco_vs_musw_24h %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Wild Mouse Baseline")

ps = align_patches(p1, p2)
pdf(paste0(paths$out, "/gsea_baseline24h.pdf"), width=w_in*2, height=h_in)
ps
dev.off()

# --- T2 ---
facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Mouse"
)
p1=gse_baseline$aco_vs_mus_12d %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Mouse Baseline")

facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Wild Mouse"
)
p2=gse_baseline$aco_vs_musw_12d %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Wild Mouse Baseline")

ps = align_patches(p1, p2)
pdf(paste0(paths$out, "/gsea_baseline12d.pdf"), width=w_in*2, height=h_in)
ps
dev.off()

#### CHECKPOINT 2 ####

# save.image(file=paste0(paths$objects, "/checkpoint2.Rdata"))
load(paste0(paths$objects, "/checkpoint2.Rdata"))

#### GSEA: DIFFERENCE IN RESPONSE ####

##### Run #####
# 
# # Run gsea on pairwise logFC differences
# gse_diff_response = list()
# 
# # --- Aco vs Mus ---
# gse_diff_response$aco_vs_mus10_24h = des_combined$res24h_10vsC %>%
#   mutate(logFCdiff = logFC_10vsC_aco - logFC_10vsC_mus) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# gse_diff_response$aco_vs_mus20_24h = des_combined$res24h_20vsC %>%
#   mutate(logFCdiff = logFC_20vsC_aco - logFC_20vsC_mus) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# gse_diff_response$aco_vs_mus10_12d = des_combined$res12d_10vsC %>%
#   mutate(logFCdiff = logFC_10vsC_aco - logFC_10vsC_mus) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# gse_diff_response$aco_vs_mus20_12d = des_combined$res12d_20vsC %>%
#   mutate(logFCdiff = logFC_20vsC_aco - logFC_20vsC_mus) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# 
# # --- Aco vs MusW ---
# gse_diff_response$aco_vs_musw10_24h = des_combined$res24h_10vsC %>%
#   mutate(logFCdiff = logFC_10vsC_aco - logFC_10vsC_musw) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# gse_diff_response$aco_vs_musw20_24h = des_combined$res24h_20vsC %>%
#   mutate(logFCdiff = logFC_20vsC_aco - logFC_20vsC_musw) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# gse_diff_response$aco_vs_musw10_12d = des_combined$res12d_10vsC %>%
#   mutate(logFCdiff = logFC_10vsC_aco - logFC_10vsC_musw) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# gse_diff_response$aco_vs_musw20_12d = des_combined$res12d_20vsC %>%
#   mutate(logFCdiff = logFC_20vsC_aco - logFC_20vsC_musw) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# 
# save(gse_diff_response, file=paste0(paths$objects, "/gsea_diff_response.Rdata"))
load(paste0(paths$objects, "/gsea_diff_response.Rdata"))

##### Dotplots #####

# --- T1 ---
facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Mouse"
)
p1=gse_diff_response$aco_vs_mus10_24h %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Mouse 10gy")
p2=gse_diff_response$aco_vs_mus20_24h %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Mouse 20gy")

facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Wild Mouse"
)
p3=gse_diff_response$aco_vs_musw10_24h %>%
  dplyr::filter(p.adjust < 0.2) %>%#<!>
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Wild Mouse 10gy")
p4=gse_diff_response$aco_vs_musw20_24h %>%
  dplyr::filter(p.adjust < 0.2) %>%#<!>
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Wild Mouse 20gy")

ps = align_patches(p1, p2, p3, p4)
pdf(paste0(paths$out, "/gsea_diff_response24h.pdf"), width=w_in*2, height=h_in)
ps
dev.off()

# --- T2 ---
facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Mouse"
)
p1=gse_diff_response$aco_vs_mus10_12d %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Mouse 10gy")
p2=gse_diff_response$aco_vs_mus20_12d %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Mouse 20gy")

facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Wild Mouse"
)
p3=gse_diff_response$aco_vs_musw10_12d %>%
  dplyr::filter(p.adjust < 0.3) %>% #<!>
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Wild Mouse 10gy")
p4=gse_diff_response$aco_vs_musw20_12d %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Wild Mouse 20gy")

ps = align_patches(p1, p2, p3, p4)
pdf(paste0(paths$out, "/gsea_diff_response12d.pdf"), width=w_in*2, height=h_in)
ps
dev.off()

#### CHECKPOINT 3 ####

# save.image(file=paste0(paths$objects, "/checkpoint3.Rdata"))
load(paste0(paths$objects, "/checkpoint3.Rdata"))

#### PATHWAYS OF INTEREST ####

##### Prepare #####

# This is based on any sig, which required padj < 0.01
# Previously, I required padj < 0.05 so its a bit more stringent now
gse_diff_sig = list()
gse_diff_sig$T1 = gse_any_sig$T1[, 1:14] %>% 
  pivot_longer(
    names_to = c("Variable", "Species", "Timepoint", "Treatment"),
    names_pattern = "(.*)_(hsa|mus|musw|aco)(24h)_(10gy|20gy)",
    cols = -c(ID, Description)) %>%
  mutate(Comparison = paste(Species, Treatment, sep = "_")) %>%
  pivot_wider(names_from=Variable) %>%
  mutate(Title = str_replace_all(Description, "regulation", "reg")) %>%
  mutate(Title = str_wrap(Title, 30)) %>%
  mutate(Sig = stars.pval(padj))
gse_diff_sig$T2 = gse_any_sig$T2[, 1:14] %>% 
  pivot_longer(
    names_to = c("Variable", "Species", "Timepoint", "Treatment"),
    names_pattern = "(.*)_(mus|musw|aco)(12d)_(10gy|20gy)",
    cols = -c(ID, Description)) %>%
  mutate(Comparison = paste(Species, Treatment, sep = "_")) %>%
  pivot_wider(names_from=Variable) %>%
  mutate(Title = str_replace_all(Description, "regulation", "reg")) %>%
  mutate(Title = str_wrap(Title, 30)) %>%
  mutate(Sig = stars.pval(padj))

##### Select #####

pathways_of_interest = list(
  "repair" = c( # Redundant / hard to interpret terms not included
    "GO:0006281" = "DNA repair",
    "GO:0006284" = "base-excision repair",
    "GO:0006289" = "nucleotide-excision repair",
    "GO:0006298" = "mismatch repair",
    "GO:0006301" = "postreplication repair",
    "GO:0036297" = "interstrand cross-link repair",
    "GO:0006302" = "double-strand break repair",
    "GO:0006303" = "double-strand break repair via nonhomologous end joining",
    "GO:0000724" = "double-strand break repair via homologous recombination",
    "GO:0000727" = "double-strand break repair via break-induced replication"
  ),
  "damage" = c( # All which appear in T1 or T2
    "GO:0000077" = "DNA damage checkpoint signaling",
    "GO:0007095" = "mitotic G2 DNA damage checkpoint signaling",
    "GO:0030330" = "DNA damage response, signal transduction by p53 class mediator",
    "GO:0042770" = "signal transduction in response to DNA damage",
    "GO:0044773" = "mitotic DNA damage checkpoint signaling",
    "GO:2000001" = "regulation of DNA damage checkpoint",
    "GO:2001020" = "regulation of response to DNA damage stimulus",
    "GO:2001021" = "negative regulation of response to DNA damage stimulus",
    "GO:2001022" = "positive regulation of response to DNA damage stimulus"
  ),
  "telomeres" = c(
    "GO:0032205" = "negative regulation of telomere maintenance",
    "GO:0032206" = "positive regulation of telomere maintenance",
    "GO:0032200" = "telomere organization",
    "GO:0000723" = "telomere maintenance",
    "GO:0043247" = "telomere maintenance in response to DNA damage",
    "GO:0007004" = "telomere maintenance via telomerase",
    "GO:0010833" = "telomere maintenance via telomere lengthening"
  ), 
  "mapk_erk" = c(
    "GO:0070374" = "positive reg. of ERK1 and ERK2 cascade",
    "GO:0070372" = "reg. of ERK1 and ERK2 cascade",
    "GO:0070371" = "ERK1 and ERK2 cascade",
    "GO:0043410" = "positive reg. of MAPK cascade"
  )
  # "p53" = c(), # Used to have some but increased stringency removed them
  # "p21" = c() # Never had any sig
)

# Function to help search terms
gse_diff_sig$T1 %>%
  dplyr::filter(grepl("telo", Description, ignore.case = T)) %>%
  dplyr::select(ID, Description) %>%
  distinct() %>%
  View()

##### Plot #####

cols = brewer.pal(n = 6, name = "Paired")

# --- T1 ---
for (concept in names(pathways_of_interest)) {
  ps = list()
  for (term in names(pathways_of_interest[[concept]])) {
    tmp = gse_diff_sig$T1 %>%
      dplyr::filter(ID == term)
    if (nrow(tmp) == 0) next
    p=ggplot(tmp, aes(x=Comparison, y=NES, fill=Comparison, label=Sig))+
      geom_bar(stat="identity")+
      geom_text(aes(y=NES/2), color="white")+
      scale_fill_manual(values=cols)+
      ggtitle(str_wrap(pathways_of_interest[[concept]][term], 25))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    ps[[length(ps)+1]] = p
  }
  plt = wrap_plots(ps, guides = "collect", ncol=3)&
    theme(axis.title=element_text(size=8),
          legend.position = "bottom",
          legend.box = "vertical")
  page_height = ceiling(length(ps) / 3) / 4 # 1 if 4 rows (i.e. 0.25 per row)
  ggsave(plot=plt, paste0(paths$pathways, "/", concept, "_24h.png"), width=w, height=page_height*h, units="mm")
}

# --- T2 ---
for (concept in names(pathways_of_interest)) {
  ps = list()
  for (term in names(pathways_of_interest[[concept]])) {
    tmp = gse_diff_sig$T2 %>%
      dplyr::filter(ID == term)
    if (nrow(tmp) == 0) next
    p=ggplot(tmp, aes(x=Comparison, y=NES, fill=Comparison, label=Sig))+
      geom_bar(stat="identity")+
      geom_text(aes(y=NES/2), color="white")+
      scale_fill_manual(values=cols)+
      ggtitle(str_wrap(pathways_of_interest[[concept]][term], 25))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    ps[[length(ps)+1]] = p
  }
  plt = wrap_plots(ps, guides = "collect", ncol=3)&
    theme(axis.title=element_text(size=8),
          legend.position = "bottom",
          legend.box = "vertical")
  page_height = ceiling(length(ps) / 3) / 4 # 1 if 4 rows (i.e. 0.25 per row)
  ggsave(plot=plt, paste0(paths$pathways, "/", concept, "_12d.png"), width=w, height=page_height*h, units="mm")
}

#### CUSTOM GENE SETS ####

##### Load sets #####

my_read_gmt = function(file) {
  lines = readLines(file)
  lines = strsplit(lines, "\t")
  ids = vapply(lines, function(y) y[1], character(1))
  descs = vapply(lines, function(y) y[2], character(1))
  genes = lapply(lines, "[", -c(1:2))
  names(genes) = ids
  names(descs) = ids
  gmt = stack(genes)
  gmt$desc = descs[gmt$ind]
  colnames(gmt) = c("gene", "term", "desc")
  return(gmt[, c("term", "gene", "desc")])
}

gmts = list()
gmts$sen_mayo = my_read_gmt(paste0(paths$out, "/../../extra_gene_sets/SAUL_SEN_MAYO.v2023.2.Mm.edited.gmt"))
gmts$fathima_custom = my_read_gmt(paste0(paths$out, "/../../extra_gene_sets/fathima_custom_sets.gmt"))

##### Functions #####

my_gse_custom = function(table, logFC_column, sym_col, gmt) { # Have to specify symbol col because humans 
  options(warn = 1)
  # Make gene list
  gene_list = table[[logFC_column]]
  names(gene_list) = table[[sym_col]]
  gene_list = sort(gene_list, decreasing = TRUE)
  # Run GSEA
  set.seed(1337)
  res = GSEA(
    geneList=gene_list,
    TERM2GENE = gmt,
    verbose = TRUE,
    pvalueCutoff = 1.1,
    minGSSize = 5, # Important or some small custom lists were not shown
    maxGSSize = 1400,
    BPPARAM = SerialParam())
  return(res)
}
get_geneset_zscores = function(df, meta, gmt) {
  all_genes = unique(gmt$gene)
  all_genes = intersect(all_genes, colnames(df))
  zscores = df[rownames(meta), all_genes]
  zscores = apply(zscores, 2, scale)
  return(zscores)
}

##### Run GSEA and get Zscores #####

res_custom_gsets = list()
for (gmt in names(gmts)) {
  res_custom_gsets[[gmt]] = list()
  for (comp in names(des_treatment)) {
    geneids = ifelse(grepl("hsa", comp), "MouseSymbol", "GeneSymbol")
    res_custom_gsets[[gmt]]$gsea = rbind(
      res_custom_gsets[[gmt]]$gsea,
      my_gse_custom(des_treatment[[comp]], "logFC_10vsC", geneids, gmts[[gmt]]) %>%
        as.data.frame() %>%
        mutate(Sig = stars.pval(p.adjust)) %>%
        mutate(Context = comp, Comparison = "10vsC") %>%
        mutate(Condition = paste(Context, Comparison, sep="_")),
      my_gse_custom(des_treatment[[comp]], "logFC_20vsC", geneids, gmts[[gmt]]) %>%
        as.data.frame() %>%
        mutate(Sig = stars.pval(p.adjust)) %>%
        mutate(Context = comp, Comparison = "20vsC") %>%
        mutate(Condition = paste(Context, Comparison, sep="_")))
  }
  res_custom_gsets[[gmt]]$zscores = get_geneset_zscores(norm$union, meta, gmts[[gmt]])
}

##### Plot #####

for (gmt in names(gmts)) {
  for (tp in c("24h", "12d")) { # Separate pdfs for each timepoint
    pdf(paste0(paths$custom_gsets, "/", gmt, "_", tp, ".pdf"), height = h_in, width = w_in)
    for (gset in unique(gmts[[gmt]]$term)) {
      gset_desc = gmts[[gmt]][gmts[[gmt]]$term == gset, "desc"]
      gset_genes = gmts[[gmt]] %>%
        dplyr::filter(term == gset) %>%
        pull(gene)
      gset_zscores = res_custom_gsets[[gmt]]$zscores
      gset_zscores = gset_zscores[, colnames(gset_zscores) %in% gset_genes]
      cols = brewer.pal(n = 8, name = "Paired")
      if (tp == "12d") {cols = cols[-c(3:4)]}
      # GSEA NES bars
      p1 = res_custom_gsets[[gmt]]$gsea %>%
        dplyr::filter(grepl(tp, Context)) %>% # Filter to this timepoint
        dplyr::filter(ID == gset) %>%
        ggplot(., aes(x=Condition, y=NES, fill=Condition, label=Sig))+
        geom_bar(stat="identity")+
        geom_text(aes(y=NES/2), color="white")+
        scale_fill_manual(values=cols)+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        ggtitle("GSEA")
      # Average pathway zscore boxes
      p2 = meta %>% 
        mutate(setZscore = rowMeans(gset_zscores)) %>%
        dplyr::filter(Timepoint == tp) %>%
        dplyr::filter(Treatment != "ConP") %>% # We're not use the ConP condition
        mutate(Treatment = factor(Treatment, levels = c("Con", "10Gy", "20Gy"))) %>%
        ggplot(., aes(x=Species, y=setZscore, fill=Treatment))+
        geom_boxplot()+
        ggtitle("Gene set Zscore")
      # Zscore heatmap
      tmp = data.frame(gset_zscores) %>%
        mutate(Species = meta$Species, Treatment = meta$Treatment, Timepoint = meta$Timepoint) %>%
        dplyr::filter(Timepoint == tp) %>%
        dplyr::filter(Treatment != "ConP") %>% # We're not use the ConP condition
        dplyr::select(-Timepoint) %>%
        mutate(Species = factor(Species, levels = unique(Species))) %>% # Lot of mess to order columns nicely
        mutate(Treatment = factor(Treatment, levels = c("Con", "10Gy", "20Gy"))) %>%
        mutate(Condition = paste(Species, Treatment, sep="\n")) %>%
        mutate(Condition = factor(Condition, levels=unique(Condition[order(Species, Treatment)]))) %>%
        pivot_longer(cols=-c(Species, Treatment, Condition), names_to = "gene")
      tmp2 = tmp %>%
        group_by(Condition, gene) %>%
        summarize(meanExpr = mean(value)) %>%
        mutate(Condition = str_replace(Condition, " ", "_")) %>%
        pivot_wider(names_from = Condition, values_from = meanExpr) %>%
        column_to_rownames("gene")
      clust = hclust(dist(tmp2))
      p3 = tmp %>%
        group_by(Species, Treatment, Condition, gene) %>%
        summarize(meanExpr = mean(value)) %>% 
        mutate(gene = factor(gene, levels =rownames(tmp2)[clust$order])) %>%
        # mutate(gene = factor(gene, levels = rownames(stats))) %>%
        ggplot(., aes(Condition, gene, fill=meanExpr))+
        geom_tile()+
        scale_fill_gradient2(low="blue", high="red")+
        ggtitle("Gene heatmap")
      ptop = p1+p2+plot_layout(widths=c(1,2))
      print(ptop/p3+plot_layout(heights=c(1,4))+plot_annotation(title=gset_desc))
    }
    dev.off()
  }
}

#### CHECKPOINT 4 ####

# save.image(file=paste0(paths$objects, "/checkpoint4.Rdata"))
load(paste0(paths$objects, "/checkpoint4.Rdata"))

#### ACOMYS SPECIFIC GENES ####

##### Prep #####

comp24h = rbind(
  des_treatment$aco_24h %>%
    dplyr::select(GeneSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    mutate(Species = "ACO", Timepoint = "24h"),
  des_treatment$mus_24h %>%
    dplyr::select(GeneSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    mutate(Species = "MUS", Timepoint = "24h"),
  des_treatment$musw_24h %>%
    dplyr::select(GeneSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    mutate(Species = "MUSW", Timepoint = "24h"),
  des_treatment$hsa_24h %>%
    dplyr::select(MouseSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    dplyr::rename(GeneSymbol = "MouseSymbol") %>%
    mutate(Species = "HSA", Timepoint = "24h") %>%
    drop_na(GeneSymbol)
)

comp12d = rbind(
  des_treatment$aco_12d %>%
    dplyr::select(GeneSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    mutate(Species = "ACO", Timepoint = "12d"),
  des_treatment$mus_12d %>%
    dplyr::select(GeneSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    mutate(Species = "MUS", Timepoint = "12d"),
  des_treatment$musw_12d %>%
    dplyr::select(GeneSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    mutate(Species = "MUSW", Timepoint = "12d")
)

comp24h_wide = comp24h %>%
  group_by(GeneSymbol, Species) %>%
  mutate(row = row_number()) %>%
  dplyr::filter(row == 1) %>%
  pivot_wider(names_from = "Species", values_from = c("logFC_10vsC", "logFC_20vsC", "PAdj_10vsC", "PAdj_20vsC"), values_fill=NA) %>%
  dplyr::select(-row)
comp12d_wide = comp12d %>%
  group_by(GeneSymbol, Species) %>%
  mutate(row = row_number()) %>%
  dplyr::filter(row == 1) %>%
  pivot_wider(names_from = "Species", values_from = c("logFC_10vsC", "logFC_20vsC", "PAdj_10vsC", "PAdj_20vsC"), values_fill=NA) %>%
  dplyr::select(-row)

##### Aco specific response #####

th_high = 1
th_low = 0.5
th_pval = 0.1

regions = data.frame(
  xmin = c(-5, th_high),
  ymin = c(-th_low, -5),
  xmax = c(-th_high, 5),
  ymax = c(5, th_low),
  color = c("blue", "red")
)

ps = list()
aco_de24h = comp24h_wide %>%
  mutate(Passes = abs(logFC_10vsC_ACO) > th_high & PAdj_10vsC_ACO < th_pval)
ps[[1]] = ggplot(aco_de24h, aes(logFC_10vsC_ACO, logFC_10vsC_MUS, color=Passes))+
  geom_point(size=0.5)+
  geom_vline(xintercept = th_high)+
  geom_vline(xintercept = -th_high)+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de24h = subset(aco_de24h, Passes) %>%
  mutate(Passes = abs(logFC_20vsC_ACO) > th_high & PAdj_20vsC_ACO < th_pval)
ps[[2]] = ggplot(aco_de24h, aes(logFC_20vsC_ACO, logFC_20vsC_MUS, color=Passes))+
  geom_point(size=0.5)+
  geom_vline(xintercept = th_high)+
  geom_vline(xintercept = -th_high)+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de24h = subset(aco_de24h, Passes) %>%
  mutate(Passes = is.na(logFC_10vsC_MUS) | (abs(logFC_10vsC_MUS) < th_low & PAdj_10vsC_MUS > th_pval)| logFC_10vsC_MUS * logFC_10vsC_ACO < 0)
ps[[3]] = ggplot()+
  geom_point(data = aco_de24h, aes(logFC_10vsC_ACO, logFC_10vsC_MUS, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de24h = subset(aco_de24h, Passes) %>%
  mutate(Passes = is.na(logFC_20vsC_MUS) | (abs(logFC_20vsC_MUS) < th_low & PAdj_20vsC_MUS > th_pval)| logFC_20vsC_MUS * logFC_20vsC_ACO < 0)
ps[[4]] = ggplot()+
  geom_point(data = aco_de24h, aes(logFC_20vsC_ACO, logFC_20vsC_MUS, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de24h = subset(aco_de24h, Passes) %>%
  mutate(Passes = is.na(logFC_10vsC_MUSW) | (abs(logFC_10vsC_MUSW) < th_low & PAdj_10vsC_MUSW > th_pval)| logFC_10vsC_MUSW * logFC_10vsC_ACO < 0)
ps[[5]] = ggplot()+
  geom_point(data = aco_de24h, aes(logFC_10vsC_ACO, logFC_10vsC_MUSW, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de24h = subset(aco_de24h, Passes) %>%
  mutate(Passes = is.na(logFC_20vsC_MUSW) | (abs(logFC_20vsC_MUSW) < th_low & PAdj_20vsC_MUSW > th_pval)| logFC_20vsC_MUSW * logFC_20vsC_ACO < 0)
ps[[6]] = ggplot()+
  geom_point(data = aco_de24h, aes(logFC_20vsC_ACO, logFC_20vsC_MUSW, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de24h = subset(aco_de24h, Passes) %>%
  mutate(Passes = is.na(logFC_10vsC_HSA) | (abs(logFC_10vsC_HSA) < th_low & PAdj_10vsC_HSA > th_pval)| logFC_10vsC_HSA * logFC_10vsC_ACO < 0)
ps[[7]] = ggplot()+
  geom_point(data = aco_de24h, aes(logFC_10vsC_ACO, logFC_10vsC_HSA, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de24h = subset(aco_de24h, Passes) %>%
  mutate(Passes = is.na(logFC_20vsC_HSA) | (abs(logFC_20vsC_HSA) < th_low & PAdj_20vsC_HSA > th_pval)| logFC_20vsC_HSA * logFC_20vsC_ACO < 0)
ps[[8]] = ggplot()+
  geom_point(data = aco_de24h, aes(logFC_20vsC_ACO, logFC_20vsC_HSA, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de24h = subset(aco_de24h, Passes)
wrap_plots(ps, ncol=2)+plot_layout(guides="collect")
ggsave(paste0(paths$out, "/aco_specific_genes_filtering24h.png"), width = w, height = h, units="mm")

aco_de24h_fewna = aco_de24h[rowSums(is.na(aco_de24h))<=4, ] # Can be NA in one species (2x logFC, 2x pvals)


## Too complicated to define a set of genes that only do not chage in ACO in a way that is meaningful
# aco_not_de24h = comp24h_wide %>%
#   dplyr::filter(abs(logFC_10vsC_ACO) < th_low & PAdj_10vsC_ACO > th_pval) %>% # Not sig in ACO 10vsC
#   dplyr::filter(abs(logFC_20vsC_ACO) < th_low  & PAdj_20vsC_ACO > th_pval) %>% # Not sig in ACO 20vsC
#   dplyr::filter(abs(logFC_10vsC_MUS) > th_high & PAdj_10vsC_MUS < th_pval)
#   dplyr::filter(abs(logFC_20vsC_MUS) > th_high & PAdj_20vsC_MUS < th_pval)
#   dplyr::filter(abs(logFC_10vsC_MUSW) > th_high & PAdj_10vsC_MUSW < th_pval)
#   dplyr::filter(abs(logFC_20vsC_MUSW) > th_high & PAdj_20vsC_MUSW < th_pval)
#   dplyr::filter(abs(logFC_10vsC_HSA) > th_high & PAdj_10vsC_HSA < th_pval)
#   dplyr::filter(abs(logFC_20vsC_HSA) > th_high & PAdj_20vsC_HSA < th_pval)

ps = list()
aco_de12d = comp12d_wide %>%
  mutate(Passes = abs(logFC_10vsC_ACO) > th_high & PAdj_10vsC_ACO < th_pval)
ps[[1]] = ggplot(aco_de12d, aes(logFC_10vsC_ACO, logFC_10vsC_MUS, color=Passes))+
  geom_point(size=0.5)+
  geom_vline(xintercept = th_high)+
  geom_vline(xintercept = -th_high)+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de12d = subset(aco_de12d, Passes) %>%
  mutate(Passes = abs(logFC_20vsC_ACO) > th_high & PAdj_20vsC_ACO < th_pval)
ps[[2]] = ggplot(aco_de12d, aes(logFC_20vsC_ACO, logFC_20vsC_MUS, color=Passes))+
  geom_point(size=0.5)+
  geom_vline(xintercept = th_high)+
  geom_vline(xintercept = -th_high)+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de12d = subset(aco_de12d, Passes) %>%
  mutate(Passes = is.na(logFC_10vsC_MUS) | (abs(logFC_10vsC_MUS) < th_low & PAdj_10vsC_MUS > th_pval)| logFC_10vsC_MUS * logFC_10vsC_ACO < 0)
ps[[3]] = ggplot()+
  geom_point(data = aco_de12d, aes(logFC_10vsC_ACO, logFC_10vsC_MUS, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de12d = subset(aco_de12d, Passes) %>%
  mutate(Passes = is.na(logFC_20vsC_MUS) | (abs(logFC_20vsC_MUS) < th_low & PAdj_20vsC_MUS > th_pval)| logFC_20vsC_MUS * logFC_20vsC_ACO < 0)
ps[[4]] = ggplot()+
  geom_point(data = aco_de12d, aes(logFC_20vsC_ACO, logFC_20vsC_MUS, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de12d = subset(aco_de12d, Passes) %>%
  mutate(Passes = is.na(logFC_10vsC_MUSW) | (abs(logFC_10vsC_MUSW) < th_low & PAdj_10vsC_MUSW > th_pval)| logFC_10vsC_MUSW * logFC_10vsC_ACO < 0)
ps[[5]] = ggplot()+
  geom_point(data = aco_de12d, aes(logFC_10vsC_ACO, logFC_10vsC_MUSW, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de12d = subset(aco_de12d, Passes) %>%
  mutate(Passes = is.na(logFC_20vsC_MUSW) | (abs(logFC_20vsC_MUSW) < th_low & PAdj_20vsC_MUSW > th_pval)| logFC_20vsC_MUSW * logFC_20vsC_ACO < 0)
ps[[6]] = ggplot()+
  geom_point(data = aco_de12d, aes(logFC_20vsC_ACO, logFC_20vsC_MUSW, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4))
aco_de12d = subset(aco_de12d, Passes)
wrap_plots(ps, ncol=2)+plot_layout(guides="collect")
ggsave(paste0(paths$out, "/aco_specific_genes_filtering12d.png"), width = w, height = 0.8*h, units="mm")

aco_de12d_fewna = aco_de12d[rowSums(is.na(aco_de12d))<=4, ] # Can be NA in one species (2x logFC, 2x pvals)

# Too complicated to define a set of genes that only do not chage in ACO in a way that is meaningful
# aco_not_de12d = comp12d_wide %>%
#   dplyr::filter(abs(logFC_10vsC_ACO) < th_low & abs(logFC_20vsC_ACO) < th_low & 
#                 PAdj_10vsC_ACO > th_pval & PAdj_20vsC_ACO > th_pval) %>%
#   dplyr::filter((abs(logFC_10vsC_MUS) > th_high & PAdj_10vsC_MUS < th_pval) & 
#                 (abs(logFC_20vsC_MUS) > th_high & PAdj_20vsC_MUS < th_pval)) %>%
#   dplyr::filter((abs(logFC_10vsC_MUSW) > th_high & PAdj_10vsC_MUSW < th_pval) & 
#                 (abs(logFC_20vsC_MUSW) > th_high & PAdj_20vsC_MUSW < th_pval))

##### Save lists for STRING #####

dir.create(paste0(paths$out, "/lists"), showWarnings = F)

aco_de24h_fewna %>%
  dplyr::filter(logFC_10vsC_ACO > 0 & logFC_20vsC_ACO > 0) %>%
  pull(GeneSymbol) %>%
  write_lines(file=paste0(paths$out, "/lists/aco_only_de_up24h.txt"))
aco_de24h_fewna %>%
  dplyr::filter(logFC_10vsC_ACO < 0 & logFC_20vsC_ACO < 0) %>%
  pull(GeneSymbol) %>%
  write_lines(file=paste0(paths$out, "/lists/aco_only_de_down24h.txt"))
# aco_not_de24h %>%
#   dplyr::filter(logFC_10vsC_ACO > 0 & logFC_20vsC_ACO > 0) %>%
#   pull(GeneSymbol) %>%
#   write_lines(file=paste0(paths$out, "/lists/aco_not_de24h.txt"))
aco_de12d_fewna %>%
  dplyr::filter(logFC_10vsC_ACO > 0 & logFC_20vsC_ACO > 0) %>%
  pull(GeneSymbol) %>%
  write_lines(file=paste0(paths$out, "/lists/aco_only_de_up12d.txt"))
aco_de12d_fewna %>%
  dplyr::filter(logFC_10vsC_ACO < 0 & logFC_20vsC_ACO < 0) %>%
  pull(GeneSymbol) %>%
  write_lines(file=paste0(paths$out, "/lists/aco_only_de_down12d.txt"))
# aco_not_de12d %>%
#   dplyr::filter(logFC_10vsC_ACO > 0 & logFC_20vsC_ACO > 0) %>%
#   pull(GeneSymbol) %>%
#   write_lines(file=paste0(paths$out, "/lists/aco_not_de12d.txt"))

##### Plot heatmap #####

aco_de24h_fewna %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_10vsC_ACO)) %>%
  pivot_longer(cols=-c(GeneSymbol, Timepoint), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus", "musw", "hsa"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_genes_heatmap24h.png"), width = w, height = 0.5*h, units="mm")

aco_de12d_fewna %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_10vsC_ACO)) %>%
  pivot_longer(cols=-c(GeneSymbol, Timepoint), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus", "musw", "hsa"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_genes_heatmap12d.png"), width = w, height = 2*h, units="mm")

##### Gene Ontology #####

my_go = function(fg, bg, org, key = "SYMBOL") {
  res = enrichGO(
    gene          = fg,
    universe      = bg,
    OrgDb         = org,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1.1,
    qvalueCutoff = 1.1,
    readable      = TRUE,
    keyType = key)
  return(res)
}

go_aco_specific = list()
bg24h = comp24h_wide$GeneSymbol
bg12d = comp12d_wide$GeneSymbol

go_aco_specific$up24h = aco_de24h_fewna %>%
  dplyr::filter(logFC_10vsC_ACO > 0 & logFC_20vsC_ACO > 0) %>%
  pull(GeneSymbol) %>%
  my_go(., bg24h, org.Mm.eg.db)
go_aco_specific$down24h = aco_de24h_fewna %>%
  dplyr::filter(logFC_10vsC_ACO < 0 & logFC_20vsC_ACO < 0) %>%
  pull(GeneSymbol) %>%
  my_go(., bg24h, org.Mm.eg.db)
go_aco_specific$up12d = aco_de12d_fewna %>%
  dplyr::filter(logFC_10vsC_ACO > 0 & logFC_20vsC_ACO > 0) %>%
  pull(GeneSymbol) %>%
  my_go(., bg12d, org.Mm.eg.db)
go_aco_specific$down12d = aco_de12d_fewna %>%
  dplyr::filter(logFC_10vsC_ACO < 0 & logFC_20vsC_ACO < 0) %>%
  pull(GeneSymbol) %>%
  my_go(., bg12d, org.Mm.eg.db)

# save(go_aco_specific, file=paste0(paths$objects, "/go_aco_specific.Rdata"))
load(paste0(paths$objects, "/go_aco_specific.Rdata"))

go_aco_specific = rbind(
  as.data.frame(go_aco_specific$up24h) %>%
    mutate(Direction = "Up", Timepoint = "24h"),
  as.data.frame(go_aco_specific$up12d) %>%
    mutate(Direction = "Up", Timepoint = "12d"),
  as.data.frame(go_aco_specific$down24h) %>%
    mutate(Direction = "Down", Timepoint = "24h"),
  as.data.frame(go_aco_specific$down12d) %>%
    mutate(Direction = "Down", Timepoint = "12d")
)

go_aco_specific2 = go_aco_specific %>%
  dplyr::filter(p.adjust < 0.05) 
# Nothing solid: all have GeneRatio 1/3. Previously it was nothing at all, unclear why.

#### ACOMYS SPECIFIC PATHWAYS ####

##### Aco specific response #####

# Exploration to figure out thresholds
names(gse_combined$T1)
ggplot(gse_combined$T1, aes(NES_aco24h_10gy, -log10(padj_aco24h_10gy), color=padj_aco24h_10gy<0.01))+
  geom_point()
sum(gse_combined$T1$padj_aco24h_10gy<0.01, na.rm=T)
nrow(gse_combined$T1)

min(abs(gse_combined$T1[gse_combined$T1$padj_aco24h_10gy<0.01, "NES_aco24h_10gy"]), na.rm=T)
summary(abs(gse_combined$T1[gse_combined$T1$padj_aco24h_10gy>0.01, "NES_aco24h_10gy"]), na.rm=T)
th_high_nes = 1.5
th_low_nes = 1
th_pval_nes = 0.05

ggplot(gse_combined$T1, aes(NES_aco24h_10gy, NES_mus24h_10gy, color=padj_aco24h_10gy<0.01 | padj_mus24h_10gy<0.01))+
  geom_point(size=0.2)

aco_nes24h = gse_combined$T1 %>%
  dplyr::filter(abs(NES_aco24h_10gy) > th_high_nes & abs(NES_aco24h_20gy) > th_high_nes &
                padj_aco24h_10gy < th_pval_nes & padj_aco24h_20gy < th_pval_nes) %>%
  dplyr::filter((is.na(NES_mus24h_10gy) | (abs(NES_mus24h_10gy) < th_low_nes & padj_mus24h_10gy > th_pval_nes) | NES_mus24h_10gy * NES_aco24h_10gy < 0) &
                (is.na(NES_mus24h_20gy) | (abs(NES_mus24h_20gy) < th_low_nes & padj_mus24h_20gy > th_pval_nes) | NES_mus24h_20gy * NES_aco24h_20gy < 0)) %>%
  dplyr::filter((is.na(NES_musw24h_10gy) | (abs(NES_musw24h_10gy) < th_low_nes & padj_musw24h_10gy > th_pval_nes) | NES_musw24h_10gy * NES_aco24h_10gy < 0) &
                (is.na(NES_musw24h_20gy) | (abs(NES_musw24h_20gy) < th_low_nes & padj_musw24h_20gy > th_pval_nes)| NES_musw24h_20gy * NES_aco24h_20gy < 0)) %>%
  dplyr::filter((is.na(NES_hsa24h_10gy) | (abs(NES_hsa24h_10gy) < th_low_nes & padj_hsa24h_10gy > th_pval_nes)| NES_hsa24h_10gy * NES_aco24h_10gy < 0) &
                (is.na(NES_hsa24h_20gy) | (abs(NES_hsa24h_20gy) < th_low_nes & padj_hsa24h_20gy > th_pval_nes)| NES_hsa24h_20gy * NES_aco24h_20gy < 0)) %>%
  dplyr::select(ID:padj_aco24h_20gy)
aco_nes24h_fewna = aco_nes24h[rowSums(is.na(aco_nes24h))<=4, ]

aco_nes12d = gse_combined$T2 %>%
  dplyr::filter(abs(NES_aco12d_10gy) > th_high_nes & abs(NES_aco12d_20gy) > th_high_nes &
                  padj_aco12d_10gy < th_pval_nes & padj_aco12d_20gy < th_pval_nes) %>%
  dplyr::filter((is.na(NES_mus12d_10gy) | (abs(NES_mus12d_10gy) < th_low_nes & padj_mus12d_10gy > th_pval_nes) | NES_mus12d_10gy * NES_aco12d_10gy < 0) &
                  (is.na(NES_mus12d_20gy) | (abs(NES_mus12d_20gy) < th_low_nes & padj_mus12d_20gy > th_pval_nes) | NES_mus12d_20gy * NES_aco12d_20gy < 0)) %>%
  dplyr::filter((is.na(NES_musw12d_10gy) | (abs(NES_musw12d_10gy) < th_low_nes & padj_musw12d_10gy > th_pval_nes) | NES_musw12d_10gy * NES_aco12d_10gy < 0) &
                  (is.na(NES_musw12d_20gy) | (abs(NES_musw12d_20gy) < th_low_nes & padj_musw12d_20gy > th_pval_nes)| NES_musw12d_20gy * NES_aco12d_20gy < 0)) %>%
  dplyr::select(ID:padj_aco12d_20gy)
aco_nes12d_fewna = aco_nes12d[rowSums(is.na(aco_nes12d))<=4, ]

##### Plot examples #####

regions_nes = data.frame(
  xmin = c(-5, th_high_nes),
  ymin = c(-th_low_nes, -5),
  xmax = c(-th_high_nes, 5),
  ymax = c(5, th_low_nes),
  color = c("blue", "red")
)

p1 = ggplot()+
  geom_point(data = gse_combined$T1, aes(NES_aco24h_10gy, NES_mus24h_10gy), size=0.2)+
  geom_rect(data = regions_nes, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, color=color), fill=NA)+
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3))+
  scale_color_manual(values=c("blue", "red"), guide="none")
p2 = ggplot()+
  geom_text_repel(data = aco_nes24h_fewna, aes(NES_aco24h_10gy, NES_mus24h_10gy, label=str_wrap(Description, 30)), size=2)+
  geom_point(data = aco_nes24h_fewna, aes(NES_aco24h_10gy, NES_mus24h_10gy), size=0.2)+
  geom_rect(data = regions_nes, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, color=color), fill=NA)+
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3))+
  scale_color_manual(values=c("blue", "red"), guide="none")
p=p1+p2
ggsave(plot=p, paste0(paths$out, "/aco_specific_gsea_filtering_ex24h.png"), width = w, height = 0.4*h, units="mm")

p1 = ggplot()+
  geom_point(data = gse_combined$T2, aes(NES_aco12d_10gy, NES_mus12d_10gy), size=0.2)+
  geom_rect(data = regions_nes, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, color=color), fill=NA)+
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3))+
  scale_color_manual(values=c("blue", "red"), guide="none")
p2 = ggplot()+
  geom_text_repel(data = aco_nes12d_fewna, aes(NES_aco12d_10gy, NES_mus12d_10gy, label=str_wrap(Description, 30)), size=2)+
  geom_point(data = aco_nes12d_fewna, aes(NES_aco12d_10gy, NES_mus12d_10gy), size=0.2)+
  geom_rect(data = regions_nes, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, color=color), fill=NA)+
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3))+
  scale_color_manual(values=c("blue", "red"), guide="none")
p=p1+p2
ggsave(plot=p, paste0(paths$out, "/aco_specific_gsea_filtering_ex12d.png"), width = w, height = 0.4*h, units="mm")

##### Plot heatmap #####

aco_nes24h_fewna %>%
  mutate(Description = fct_reorder(Description, NES_aco24h_10gy)) %>%
  pivot_longer(cols=-c(ID, Description), names_pattern = "(.*)_(.*)24h_(.*)", names_to = c("Var", "Species", "Dose")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(Species, levels=c("aco", "mus", "musw", "hsa"))) %>%
  mutate(Dose = factor(Dose, levels=c("10gy", "20gy"))) %>%
  arrange(Species, Dose) %>%
  mutate(Condition = paste(Species, Dose, "vsC", sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(padj)) %>%
  ggplot(., aes(Condition, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_gsea_heatmap24h.png"), width = w, height = 0.5*h, units="mm")

aco_nes12d_fewna %>%
  mutate(Description = fct_reorder(Description, NES_aco12d_10gy)) %>%
  pivot_longer(cols=-c(ID, Description), names_pattern = "(.*)_(.*)12d_(.*)", names_to = c("Var", "Species", "Dose")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(Species, levels=c("aco", "mus", "musw"))) %>%
  mutate(Dose = factor(Dose, levels=c("10gy", "20gy"))) %>%
  arrange(Species, Dose) %>%
  mutate(Condition = paste(Species, Dose, "vsC", sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(padj)) %>%
  ggplot(., aes(Condition, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_gsea_heatmap12d.png"), width = w, height = 0.5*h, units="mm")

#### CHECKPOINT 5 ####

# save.image(file=paste0(paths$objects, "/checkpoint5.Rdata"))
load(paste0(paths$objects, "/checkpoint5.Rdata"))

#### TUMOR SUPPRESSOR DETAIL ####

##### Get list #####

tsgenes = fread("../extra_gene_sets/TSGene.txt", data.table = F)
oncokb = fread("../extra_gene_sets/OncoKB.tsv", data.table = F)

table(oncokb$`Is Oncogene`, oncokb$`Is Tumor Suppressor Gene`)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human_to_mouse = getLDS(
  attributes = c("hgnc_symbol"), 
  filters = "hgnc_symbol", 
  values = union(tsgenes$GeneSymbol, oncokb$`Hugo Symbol`), 
  mart = human, 
  attributesL = c("mgi_symbol"), 
  martL = mouse, 
  uniqueRows=T) 
human_to_mouse = human_to_mouse %>%
  distinct(HGNC.symbol, .keep_all = T) %>%
  distinct(MGI.symbol, .keep_all = T) %>%
  column_to_rownames("HGNC.symbol")
tsgenes$MouseSymbol = human_to_mouse[tsgenes$GeneSymbol, "MGI.symbol"]
oncokb$MouseSymbol = human_to_mouse[oncokb$`Hugo Symbol`, "MGI.symbol"]

tsgenes = subset(tsgenes, !is.na(tsgenes$MouseSymbol))
oncokb = subset(oncokb, !is.na(oncokb$MouseSymbol))

##### Zscores #####

# --- TSG ---
tmp = norm$union[rownames(meta), intersect(tsgenes$MouseSymbol, colnames(norm$union))]
tmp = apply(tmp, 2, scale)
tmp = meta %>%
  mutate(zScore = rowMeans(tmp))

tmp %>%
  dplyr::filter(Timepoint == "24h") %>%
  mutate(Treatment = factor(Treatment, levels=c("Con", "10Gy", "20Gy"))) %>%
  ggplot(., aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of TSG expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_TSG_24h.png"), width = w, height = 0.4*h, units="mm")
tmp %>%
  dplyr::filter(Timepoint == "12d") %>%
  # dplyr::filter(Treatment != "ConP") %>%
  mutate(Treatment = factor(Treatment, levels=c("ConP", "Con", "10Gy", "20Gy"))) %>%
  ggplot(., aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of TSG expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_TSG_12d.png"), width = w, height = 0.4*h, units="mm")

# --- OncoKB: oncogenes ---

which_genes = oncokb[oncokb$`Is Oncogene` == "Yes", "MouseSymbol"]
tmp = norm$union[rownames(meta), intersect(which_genes, colnames(norm$union))]
tmp = apply(tmp, 2, scale)
tmp = meta %>%
  mutate(zScore = rowMeans(tmp))

tmp %>%
  dplyr::filter(Timepoint == "24h") %>%
  mutate(Treatment = factor(Treatment, levels=c("Con", "10Gy", "20Gy"))) %>%
  ggplot(., aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of oncogene (OncoKB) expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_oncoKB_oncogene_24h.png"), width = w, height = 0.4*h, units="mm")
tmp %>%
  dplyr::filter(Timepoint == "12d") %>%
  # dplyr::filter(Treatment != "ConP") %>%
  mutate(Treatment = factor(Treatment, levels=c("ConP", "Con", "10Gy", "20Gy"))) %>%
  ggplot(., aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of oncogene (OncoKB) expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_oncoKB_oncogene_12d.png"), width = w, height = 0.4*h, units="mm")

# --- OncoKB: tumor suppressors ---

which_genes = oncokb[oncokb$`Is Tumor Suppressor Gene` == "Yes", "MouseSymbol"]
tmp = norm$union[rownames(meta), intersect(which_genes, colnames(norm$union))]
tmp = apply(tmp, 2, scale)
tmp = meta %>%
  mutate(zScore = rowMeans(tmp))

tmp %>%
  dplyr::filter(Timepoint == "24h") %>%
  mutate(Treatment = factor(Treatment, levels=c("Con", "10Gy", "20Gy"))) %>%
  ggplot(., aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of tumor suppressors (OncoKB) expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_oncoKB_tsupp_24h.png"), width = w, height = 0.4*h, units="mm")
tmp %>%
  dplyr::filter(Timepoint == "12d") %>%
  # dplyr::filter(Treatment != "ConP") %>%
  mutate(Treatment = factor(Treatment, levels=c("ConP", "Con", "10Gy", "20Gy"))) %>%
  ggplot(., aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of tumor suppressors (OncoKB) expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_oncoKB_tsupp_12d.png"), width = w, height = 0.4*h, units="mm")

##### Aco specific ones #####

aco_de24h_fewna$isTS = aco_de24h_fewna$GeneSymbol %in% tsgenes$MouseSymbol
aco_de12d_fewna$isTS = aco_de12d_fewna$GeneSymbol %in% tsgenes$MouseSymbol
aco_de24h_fewna$isTS2 = aco_de24h_fewna$GeneSymbol %in% oncokb[oncokb$`Is Tumor Suppressor Gene` == "Yes", "MouseSymbol"]
aco_de12d_fewna$isTS2 = aco_de12d_fewna$GeneSymbol %in% oncokb[oncokb$`Is Tumor Suppressor Gene` == "Yes", "MouseSymbol"]
aco_de24h_fewna$isOncogene = aco_de24h_fewna$GeneSymbol %in% oncokb[oncokb$`Is Oncogene` == "Yes", "MouseSymbol"]
aco_de12d_fewna$isOncogene = aco_de12d_fewna$GeneSymbol %in% oncokb[oncokb$`Is Oncogene` == "Yes", "MouseSymbol"]

# --- TSG ---
aco_de24h_fewna %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isTS) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_10vsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus", "musw", "hsa"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank(),
        legend.position = "top")+
  ggtitle("Tumor suppressors (TSG) with unique Acomys response")
ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap24h_TSG.png"), width = w, height = 0.3*h, units="mm")

aco_de12d_fewna %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isTS) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_10vsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus", "musw", "hsa"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())+
  ggtitle("Tumor suppressors (TSG) with unique Acomys response")
ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap12d_TSG.png"), width = w, height = 0.5*h, units="mm")

# --- OncoKB: tsupp ---
# aco_de24h_fewna %>%
#   ungroup() %>%
#   dplyr::select(-Passes) %>%
#   dplyr::filter(isTS2) %>%
#   mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_10vsC_ACO)) %>%
#   pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
#   pivot_wider(names_from = Var, values_from = value) %>%
#   mutate(Species = factor(tolower(Species), levels=c("aco", "mus", "musw", "hsa"))) %>%
#   arrange(Species, Comparison) %>%
#   mutate(Condition = paste(Species, Comparison, sep="_")) %>%
#   mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
#   mutate(Sig = stars.pval(PAdj)) %>%
#   ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
#   geom_tile()+
#   scale_fill_gradient2(low="blue", high="red")+
#   theme(axis.text.x = element_text(angle=45, hjust=1),
#         axis.title = element_blank(),
#         legend.position = "top")+
#   ggtitle("Tumor suppressors (OncoKB) with unique Acomys response")
# # Nothing
# ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap24h_OncoKB_tsupp.png"), width = w, height = 0.3*h, units="mm")

aco_de12d_fewna %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isTS2) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_10vsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus", "musw", "hsa"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())+
  ggtitle("Tumor suppressors (OncoKB) with unique Acomys response")
ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap12d_OncoKB_tsupp.png"), width = w, height = 0.25*h, units="mm")

# --- OncoKB: oncogene ---
aco_de24h_fewna %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isOncogene) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_10vsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus", "musw", "hsa"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank(),
        legend.position = "top")+
  ggtitle("Oncogenes (OncoKB) with unique Acomys response")
# Nothing
ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap24h_OncoKB_oncogene.png"), width = w, height = 0.3*h, units="mm")

aco_de12d_fewna %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isOncogene) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_10vsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus", "musw", "hsa"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())+
  ggtitle("Oncogenes (OncoKB) with unique Acomys response")
ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap12d_OncoKB_oncogene.png"), width = w, height = 0.5*h, units="mm")

#### MAPK-ERK detail ####

source("../utils.R")
pathwy_mapk = read_kegg_map("../pathway_maps/mmu04010.xml")

# 24h, 10gy
logFCs = des_combined$res24h_10vsC %>%
  dplyr::select(starts_with("logFC")) %>%
  dplyr::select(contains("vsC"))
pathwy_mapk = make_scores(pathwy_mapk, logFCs, "logFC_10vsC_aco")
plot_kegg_map(pathwy_mapk, "score")
ggsave(paste0(paths$out, "/detail_mapk/pathway_map_24h_10gy.png"), width = 2*w, height = 1*h, units="mm")

# 24h, 20gy
logFCs = des_combined$res24h_20vsC %>%
  dplyr::select(starts_with("logFC")) %>%
  dplyr::select(contains("vsC"))
pathwy_mapk = make_scores(pathwy_mapk, logFCs, "logFC_20vsC_aco")
plot_kegg_map(pathwy_mapk, "score")
ggsave(paste0(paths$out, "/detail_mapk/pathway_map_24h_20gy.png"), width = 2*w, height = 1*h, units="mm")

# 12d, 10gy
logFCs = des_combined$res12d_10vsC %>%
  dplyr::select(starts_with("logFC")) %>%
  dplyr::select(contains("vsC"))
pathwy_mapk = make_scores(pathwy_mapk, logFCs, "logFC_10vsC_aco")
plot_kegg_map(pathwy_mapk, "score")
ggsave(paste0(paths$out, "/detail_mapk/pathway_map_12d_10gy.png"), width = 2*w, height = 1*h, units="mm")

# 12d, 20gy
logFCs = des_combined$res12d_20vsC %>%
  dplyr::select(starts_with("logFC")) %>%
  dplyr::select(contains("vsC"))
pathwy_mapk = make_scores(pathwy_mapk, logFCs, "logFC_20vsC_aco")
plot_kegg_map(pathwy_mapk, "score")
ggsave(paste0(paths$out, "/detail_mapk/pathway_map_12d_20gy.png"), width = 2*w, height = 1*h, units="mm")



investigator = function(lst, des_wide, all_terms, terms_oi) {
  out = list()
  expressed = des_wide$Entrez #note that some are NA within comparison
  for (ptwy in names(terms_oi)) {
    out[[ptwy]] = list()
    all_ptwy_genes = all_terms[[terms_oi[ptwy]]]
    expr_ptwy_genes = intersect(all_ptwy_genes, expressed)
    # all_map_genes = lst$nodes$EntrezId1[!is.na(lst$nodes$EntrezId1)]
    # EntrezId1 only contains the first gene of the node (can be multiples) now i correct this:
    all_map_genes = str_extract_all(paste(lst$nodes$keggId, collapse=" "), "mmu:(\\d+)", simplify = T) # Paste all together and then split
    all_map_genes = unique(str_remove(all_map_genes, "mmu:"))
    expr_map_genes = intersect(all_map_genes, expressed)
    out[[ptwy]]$expr_ptwy_genes = expr_ptwy_genes
    out[[ptwy]]$expr_map_genes = expr_map_genes
    # Summary
    print(sprintf("========= %s =========", ptwy))
    print(sprintf("Gene set contains %i genes, of which %i are expressed", length(all_ptwy_genes), length(expr_ptwy_genes)))
    print(sprintf("Map contains %i genes, of which %i are expressed", length(all_map_genes), length(expr_map_genes)))
    print(sprintf("%.2f %% of map genes are in the gene set", 100*length(intersect(expr_ptwy_genes, expr_map_genes)) / length(expr_map_genes)))
    # GSEA on subsets
    gmtA = data.frame(
      "term" = rep("Map genes", length(expr_map_genes)), # Expr or not doesnt matter since implicit in GSEA (val empirically)
      "gene" = expr_map_genes)
    gmtB = data.frame(
      "term" = rep("Gset not in map", length(setdiff(expr_ptwy_genes, expr_map_genes))),
      "gene" = setdiff(expr_ptwy_genes, expr_map_genes))
    gmtC = data.frame(
      "term" = rep("Gset all", length(expr_ptwy_genes)),
      "gene" = expr_ptwy_genes)
    gmtD = data.frame(
      "term" = rep("Gset and map", length(intersect(expr_ptwy_genes, expr_map_genes))),
      "gene" = intersect(expr_ptwy_genes, expr_map_genes))
    out[[ptwy]]$subset_gsea = my_gse_custom(des_wide, "logFCdiff", "Entrez", rbind(gmtA, gmtB, gmtC, gmtD)) %>%
      as.data.frame()
    # Nicer gene symbols
    conversions = mapIds(org.Mm.eg.db, key = union(gmtA$gene, gmtC$gene), column = "SYMBOL", keytype = "ENTREZID")
    for (i in 1:nrow(out[[ptwy]]$subset_gsea)) {
      this_str = out[[ptwy]]$subset_gsea$core_enrichment[i]
      old_ids = unlist(strsplit(this_str, "/"))
      new_ids = conversions[old_ids]
      out[[ptwy]]$subset_gsea[i, "core_enrichment"] = paste(new_ids, collapse="/")
    }
    # So what is causing enrichment in "Gset not in map" ?
    # Gset not in map is a subset of gset, so gset is the bg
    out[[ptwy]]$go_diff = my_go(gmtB$gene, gmtC$gene, org.Mm.eg.db, "ENTREZID") %>%
      as.data.frame()%>%
      dplyr::filter(pvalue < 0.05)
    # And what is causing enrichment in "Gset and map" ?
    # Gset and map is a subset of the union of gset and map, so use the union as bg
    # Actually that doesnt work, it just gives enrichemnt of the obvious
    # Maybe just using gset?
    out[[ptwy]]$go_common = my_go(gmtD$gene, gmtC$gene, org.Mm.eg.db, "ENTREZID") %>%
      as.data.frame()%>%
      dplyr::filter(pvalue < 0.05)
  }
  return(out)
}

all_terms = as.list(org.Mm.egGO2ALLEGS)
termsoi = c(
  "mapk" = "GO:0000165",
  "pos_reg_mapk" = "GO:0043410",
  "reg_mapk" = "GO:0043408")[2] # Let's only go for pos_reg for interpretability
for (tp in c("24h", "12d")) {
  for (spe in c("mus", "musw", "hsa")) {
    if (spe == "hsa" & tp == "12d") next
    for (comp in c("10vsC", "20vsC")) {
      tmp = des_combined[[sprintf("res%s_%s", tp, comp)]]
      tmp$logFCdiff = tmp[, sprintf("logFC_%s_%s", comp, "aco")] - tmp[, sprintf("logFC_%s_%s", comp, spe)]
      test = investigator(pathwy_mapk, tmp, all_terms, termsoi)
      write.csv(test$pos_reg_mapk$subset_gsea, sprintf("%s/detail_mapk/subset_gsea_%s_%s_acoVs%s.csv", paths$out, tp, comp, spe))
      # Note: these are always the same because defined based on gene sets and not expression
      write.csv(test$pos_reg_mapk$go_diff, sprintf("%s/detail_mapk/go_on_gset-map.csv", paths$out))
    }
  }
}

View(test$pos_reg_mapk$subset_gsea)
View(test$pos_reg_mapk$go_diff)
View(test$pos_reg_mapk$go_common)

# pathwy_il17 = read_kegg_map("../pathway_maps/mmu04657.xml")
# logFCs = des_combined$res24h_10vsC %>%
#   dplyr::select(starts_with("logFC")) %>%
#   dplyr::select(contains("vsC"))
# pathwy_il17 = make_scores(pathwy_il17, logFCs, "logFC_10vsC_aco")
# plot_kegg_map(pathwy_il17, "score")
# termsoi = c(
#   "il17" = "GO:0032620")
# test = des_combined$res24h_10vsC %>%
#   mutate(logFCdiff = logFC_10vsC_aco - logFC_10vsC_hsa) %>%
#   investigator(pathwy_il17, ., all_terms, termsoi)
# View(test$il17$subset_gsea)
# View(test$il17$go_diff)
# View(gse_combined$T1[gse_combined$T1$ID == "GO:0032620", ])
# View(gse_diff_response$aco_vs_hsa10_24h[gse_diff_response$aco_vs_hsa10_24h$ID == "GO:0032620", ])
# # Doesnt show up too strongly because:
# # il17 appeared remarkable in aco_specific, which is based on gse_combined (from gse_treatment put side by side)
# # Here, however, i am working off of des_combined logFCdiff aka dLogFC
# # Still the question remains: why so little overlap?
# 
# ptwy_only = mapIds(org.Mm.eg.db, setdiff(test$il17$expr_ptwy_genes, test$il17$expr_map_genes), column="SYMBOL", keytype="ENTREZID")
# map_only = mapIds(org.Mm.eg.db, setdiff(test$il17$expr_map_genes, test$il17$expr_ptwy_genes), column="SYMBOL", keytype="ENTREZID")
# common = mapIds(org.Mm.eg.db, intersect(test$il17$expr_map_genes, test$il17$expr_ptwy_genes), column="SYMBOL", keytype="ENTREZID")
# 
# ptwy_only
# map_only
# common
# 
# # WNt 
# pathwy_wnt = read_kegg_map("../pathway_maps/mmu04310.xml")
# logFCs = des_combined$res24h_10vsC %>%
#   dplyr::select(starts_with("logFC")) %>%
#   dplyr::select(contains("vsC"))
# pathwy_il17 = make_scores(pathwy_il17, logFCs, "logFC_10vsC_aco")
# plot_kegg_map(pathwy_il17, "score")
# termsoi = c(
#   "il17" = "GO:0032620")
# test = des_combined$res24h_10vsC %>%
#   mutate(logFCdiff = logFC_10vsC_aco - logFC_10vsC_hsa) %>%
#   investigator(pathwy_il17, ., all_terms, termsoi)

#### TELOMERASE DETAIL ####

##### Load #####

telo_regulators = read.table("../extra_gene_sets/telo_regulators.tsv", header = T)
setdiff(telo_regulators$GeneId, ginfo$mus$GeneSymbol)

##### Tert #####

meta[rownames(counts$union), ] %>%
  mutate(Tert = counts$union[, "Tert"]) %>%
  ggplot(., aes(Species, Tert))+
  geom_boxplot()

meta[rownames(counts$union), ] %>%
  mutate(Tert = norm$union[rownames(counts$union), "Tert"]) %>%
  ggplot(., aes(Species, Tert))+
  geom_boxplot()
ggsave(paste0(paths$out, "/detail_telo/tert.png"))

this_meta = meta[rownames(counts$union), ]
this_meta = subset(this_meta, Timepoint == "24h")
design = model.matrix(~0+Species, data = this_meta)
dge = DGEList(t(counts$union[rownames(this_meta), ]), samples=this_meta)
# keep = filterByExpr(dge, design=design) # not when using union
# dge = dge[keep, ]
dge = calcNormFactors(dge)
this_meta$Tert = cpm(dge, normalized.lib.sizes=T, log=T)["Tert", ]
ggplot(this_meta, aes(Species, Tert))+
  geom_boxplot()

##### GO terms #####

go_sets = as.data.frame(GOTERM) %>%
  dplyr::select(-1) %>%
  dplyr::filter(grepl("telom", Term)) %>%
  distinct(go_id, .keep_all = T)

go_members = select(org.Mm.eg.db, keys=keys(org.Mm.eg.db), columns=c("GO", "ONTOLOGY", "SYMBOL")) %>%
  dplyr::select(-EVIDENCE) %>%
  dplyr::filter(GO %in% go_sets$go_id) %>%
  distinct(GO, SYMBOL, .keep_all = T) %>%
  mutate(gene = SYMBOL)
length(unique(go_members$ENTREZID))
setdiff(telo_regulators$GeneId, go_members$SYMBOL)
intersect(telo_regulators$GeneId, go_members$SYMBOL)

tmp = go_members %>%
  group_by(GO) %>%
  summarize(n = n()) %>%
  dplyr::filter(n >= 10, n <= 500)
go_members = go_members[go_members$GO %in% tmp$GO, ]

this_meta = meta %>%
  dplyr::filter(rownames(.) %in% rownames(counts$union))%>%
  dplyr::filter(Timepoint == "24h")
telo_zscores = get_geneset_zscores(norm$union, this_meta, go_members)
rownames(telo_zscores) = rownames(this_meta)

ps = list()
for (goid in unique(go_members$GO)) {
  godesc = go_sets[go_sets$go_id == goid, "Term"]
  genes_oi = go_members %>%
    dplyr::filter(GO == goid) %>%
    dplyr::filter(SYMBOL %in% colnames(telo_zscores)) %>%
    pull(SYMBOL)
  ps[[goid]] = this_meta %>%
    mutate(meanZscore = rowMeans(telo_zscores[, genes_oi])) %>%
    ggplot(., aes(Species, meanZscore, fill=Species)) +
    geom_boxplot()+
    ggtitle(godesc)
}

plt = wrap_plots(ps, guides = "collect", ncol=3)&
  theme(axis.title=element_text(size=8),
        legend.position = "bottom",
        legend.box = "vertical")
ggsave(plot=plt, paste0(paths$out, "/detail_telo/go_terms.png"), width=1.2*w, height=1.5*h, units="mm")

# Telomere capping similar in human and acomys
View(go_members[go_members$GO == "GO:0016233", ])


go_sets2 = as.data.frame(GOTERM) %>%
  dplyr::select(-1) %>%
  dplyr::filter(grepl("shelterin", Term)) %>%
  distinct(go_id, .keep_all = T)
intersect(go_sets$go_id, go_sets2$go_id)

go_members2 = select(org.Mm.eg.db, keys=keys(org.Mm.eg.db), columns=c("GO", "ONTOLOGY", "SYMBOL")) %>%
  dplyr::select(-EVIDENCE) %>%
  dplyr::filter(GO %in% go_sets2$go_id) %>%
  distinct(GO, SYMBOL, .keep_all = T) %>%
  mutate(gene = SYMBOL)

tmp = go_members2 %>%
  group_by(GO) %>%
  summarize(n = n()) %>%
  dplyr::filter(n >= 5, n <= 500)
go_members2 = go_members2[go_members2$GO %in% tmp$GO, ]

shelt_zscores = get_geneset_zscores(norm$union, this_meta, go_members2)
rownames(shelt_zscores) = rownames(this_meta)

ps = list()
for (goid in unique(go_members2$GO)) {
  godesc = go_sets2[go_sets2$go_id == goid, "Term"]
  genes_oi = go_members2 %>%
    dplyr::filter(GO == goid) %>%
    dplyr::filter(SYMBOL %in% colnames(shelt_zscores)) %>%
    pull(SYMBOL)
  ps[[goid]] = this_meta %>%
    mutate(meanZscore = rowMeans(shelt_zscores[, genes_oi])) %>%
    ggplot(., aes(Species, meanZscore, fill=Species)) +
    geom_boxplot()+
    ggtitle(godesc)
}

ps[[1]]
ggsave(plot=plt, paste0(paths$out, "/detail_telo/go_terms.png"), width=1.2*w, height=1.5*h, units="mm")

##### Regulators #####

this_meta = meta %>%
  dplyr::filter(rownames(.) %in% rownames(counts$union))%>%
  dplyr::filter(Timepoint == "24h")
telo_regulators$gene = telo_regulators$GeneId
teloreg_zscores = get_geneset_zscores(norm$union, this_meta, telo_regulators)
rownames(teloreg_zscores) = rownames(this_meta)
this_meta %>%
  mutate(meanZscore = rowMeans(teloreg_zscores)) %>%
  ggplot(., aes(Species, meanZscore, fill=Species)) +
  geom_boxplot()

design = model.matrix(~0+Species, data = this_meta)
dge = DGEList(t(counts$union[rownames(this_meta), ]), samples=this_meta)
# keep = filterByExpr(dge, design=design) # not when using union
# dge = dge[keep, ]
dge = calcNormFactors(dge)
tmp = t(cpm(dge, normalized.lib.sizes=T, log=T))
teloreg_zscores = get_geneset_zscores(tmp, this_meta, telo_regulators)
rownames(teloreg_zscores) = rownames(this_meta)
this_meta %>%
  mutate(meanZscore = rowMeans(teloreg_zscores)) %>%
  ggplot(., aes(Species, meanZscore, fill=Species)) +
  geom_boxplot()

#### PAPER FIGURES ####

# NES heatmap comparing ACO and MUS, MUSW response to irradiation
# Focus on topics like apoptosis, cell cycle, senescence, inflam, repair

# Help select terms based on gsea_diff_response and gsea_comparison.pdf
common_diff = intersect(
  gse_diff_response$aco_vs_mus10_24h %>% dplyr::filter(p.adjust < 0.05) %>% pull(ID),
  gse_diff_response$aco_vs_mus20_24h %>% dplyr::filter(p.adjust < 0.05) %>% pull(ID))
common_diff = intersect(common_diff,
  gse_diff_response$aco_vs_musw10_24h %>% dplyr::filter(p.adjust < 0.05) %>% pull(ID))
common_diff = intersect(common_diff,
  gse_diff_response$aco_vs_musw20_24h %>% dplyr::filter(p.adjust < 0.05) %>% pull(ID))

common_diff_strong = intersect(common_diff,
  gse_diff_response$aco_vs_mus10_24h %>% dplyr::filter(abs(NES) >= 2.1) %>% pull(ID))
common_diff_strong = intersect(common_diff_strong,
  gse_diff_response$aco_vs_mus20_24h %>% dplyr::filter(abs(NES) >= 2.1) %>% pull(ID))
common_diff_strong = intersect(common_diff_strong,
  gse_diff_response$aco_vs_musw10_24h %>% dplyr::filter(abs(NES) >= 2.1) %>% pull(ID))
common_diff_strong = intersect(common_diff_strong,
  gse_diff_response$aco_vs_musw20_24h %>% dplyr::filter(abs(NES) >= 2.1) %>% pull(ID))

tmp = gse_diff_sig$T1 %>%
  dplyr::filter(Species != "hsa") %>%
  dplyr::filter(ID %in% common_diff) %>%
  mutate(Description = str_sub(Description, 1, 57)) %>%
  mutate(Species = factor(Species, levels=c("aco", "musw", "mus")))
tmp2 = tmp %>%
  dplyr::filter(Species == "aco", Treatment == "10gy", NES > 0) %>%
  pull(Description)
tmp %>%
  mutate(NESpos = Description %in% tmp2) %>%
  ggplot(., aes(Species, Description, fill=NES, label=Sig))+
  geom_tile()+
  facet_grid(cols=vars(Treatment), rows=vars(NESpos), space="free", scales="free")+
  scale_fill_gradient2(low="blue", high="red")
ggsave(paste0(paths$out, "/paper_figs/helper1.pdf"), width =w, height = 2*h, units="mm")
tmp %>%
  dplyr::filter(ID %in% common_diff_strong) %>%
  mutate(NESpos = Description %in% tmp2) %>%
  ggplot(., aes(Species, Description, fill=NES, label=Sig))+
  geom_tile()+
  facet_grid(cols=vars(Treatment), rows=vars(NESpos), space="free", scales="free")+
  scale_fill_gradient2(low="blue", high="red")
ggsave(paste0(paths$out, "/paper_figs/helper2.pdf"), width =w, height = 1*h, units="mm")

selected_terms = c(
  "positive regulation of cell cycle phase transition", # From gsea_diff_response
  "DNA-templated DNA replication", # From gsea_diff_response
  "spindle organization", # From gsea_diff_response
  "sister chromatid segregation", # From gsea_diff_response
  "positive regulation of double-strand break repair", # From gsea_comparison
  "base-excision repair", # From gsea_comparison
  "nucleotide-excision repair", # From gsea_comparison
  "regulation of epithelial cell apoptotic process", # From gsea_comparison
  "regulation of endothelial cell apoptotic process", # From gsea_comparison
  "positive regulation of leukocyte chemotaxis", # From gsea_diff_response
  "innate immune response", # From gsea_diff_response
  "extracellular structure organization", # From gsea_diff_response
  "extracellular matrix organization") # From gsea_diff_response
gse_diff_sig$T1 %>%
  dplyr::filter(Species != "hsa") %>%
  mutate(Description = str_sub(Description, 1, 57)) %>%
  mutate(Species = factor(Species, levels=c("aco", "musw", "mus")))%>%
  dplyr::filter(Description %in% selected_terms) %>%
  mutate(Description = factor(Description, levels = selected_terms)) %>%
  ggplot(., aes(Species, Description, fill=NES))+
  geom_tile()+
  facet_grid(cols=vars(Treatment))+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme(axis.title = element_blank())
ggsave(paste0(paths$out, "/paper_figs/fig2d.pdf"), width = 0.8*w, height=0.3*h, units="mm")
gse_diff_sig$T1 %>%
  dplyr::filter(Species != "hsa") %>%
  mutate(Description = str_sub(Description, 1, 57)) %>%
  mutate(Species = factor(Species, levels=c("aco", "musw", "mus")))%>%
  dplyr::filter(Description %in% selected_terms) %>%
  mutate(Description = factor(Description, levels = selected_terms)) %>%
  ggplot(., aes(Species, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  facet_grid(cols=vars(Treatment))+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme(axis.title = element_blank())
ggsave(paste0(paths$out, "/paper_figs/fig2d_stars.pdf"), width = 0.8*w, height=0.3*h, units="mm")

go_sets 




# Checking if anything related to these topics has extreme NES differences
th_high_nes = 1.5 # 1.5
th_low_nes = 1 # 1
th_pval_nes = 0.05 # 0.05
aco_nes24h_nohsa = gse_combined$T1 %>%
  dplyr::filter(abs(NES_aco24h_10gy) > th_high_nes & abs(NES_aco24h_20gy) > th_high_nes &
                  padj_aco24h_10gy < th_pval_nes & padj_aco24h_20gy < th_pval_nes) %>%
  dplyr::filter((is.na(NES_mus24h_10gy) | (abs(NES_mus24h_10gy) < th_low_nes & padj_mus24h_10gy > th_pval_nes) | NES_mus24h_10gy * NES_aco24h_10gy < 0) &
                  (is.na(NES_mus24h_20gy) | (abs(NES_mus24h_20gy) < th_low_nes & padj_mus24h_20gy > th_pval_nes) | NES_mus24h_20gy * NES_aco24h_20gy < 0)) %>%
  dplyr::filter((is.na(NES_musw24h_10gy) | (abs(NES_musw24h_10gy) < th_low_nes & padj_musw24h_10gy > th_pval_nes) | NES_musw24h_10gy * NES_aco24h_10gy < 0) &
                  (is.na(NES_musw24h_20gy) | (abs(NES_musw24h_20gy) < th_low_nes & padj_musw24h_20gy > th_pval_nes)| NES_musw24h_20gy * NES_aco24h_20gy < 0)) %>%
  dplyr::select(ID:padj_aco24h_20gy)
aco_nes24h_nohsa = aco_nes24h_nohsa[rowSums(is.na(aco_nes24h_nohsa))<=2, ]

th_high_nes = 0 # 1.5
th_pval_nes = 0.05 # 0.05
aco_nes24h_nohsa = gse_combined$T1 %>%
  dplyr::filter(abs(NES_aco24h_10gy) > th_high_nes & abs(NES_aco24h_20gy) > th_high_nes &
                  padj_aco24h_10gy < th_pval_nes & padj_aco24h_20gy < th_pval_nes) %>%
  dplyr::filter((is.na(NES_mus24h_10gy) | NES_mus24h_10gy * NES_aco24h_10gy < 0) &
                  (is.na(NES_mus24h_20gy) | NES_mus24h_20gy * NES_aco24h_20gy < 0)) %>%
  dplyr::filter((is.na(NES_musw24h_10gy) | NES_musw24h_10gy * NES_aco24h_10gy < 0) &
                  (is.na(NES_musw24h_20gy) | NES_musw24h_20gy * NES_aco24h_20gy < 0)) %>%
  dplyr::select(ID:padj_aco24h_20gy)
aco_nes24h_nohsa = aco_nes24h_nohsa[rowSums(is.na(aco_nes24h_nohsa))<=2, ]

#### dLogFC>GSEA vs dGSEA ####

# gse_treatment contains GSEA results performed on logFCs
# gse_baseline contains GSEA on logFCs between Ctrl treated species
# gse_combined has them merged side by side
# gse_diff_response has GSEA performed on dLogFCs
# gse_diff_combined needs to be made

# dLogFC > GSEA
gse_diff_combined = rbind(
  data.frame(gse_diff_response$aco_vs_hsa10_24h) %>%
    mutate(Species1 = "aco", Species2 = "hsa", Dose = "10gy", Timepoint = "24h"),
  data.frame(gse_diff_response$aco_vs_hsa20_24h) %>%
    mutate(Species1 = "aco", Species2 = "hsa", Dose = "20gy", Timepoint = "24h"),
  data.frame(gse_diff_response$aco_vs_mus10_24h) %>%
    mutate(Species1 = "aco", Species2 = "mus", Dose = "10gy", Timepoint = "24h"),
  data.frame(gse_diff_response$aco_vs_mus20_24h) %>%
    mutate(Species1 = "aco", Species2 = "mus", Dose = "20gy", Timepoint = "24h"),
  data.frame(gse_diff_response$aco_vs_musw10_24h) %>%
    mutate(Species1 = "aco", Species2 = "musw", Dose = "10gy", Timepoint = "24h"),
  data.frame(gse_diff_response$aco_vs_musw20_24h) %>%
    mutate(Species1 = "aco", Species2 = "musw", Dose = "20gy", Timepoint = "24h"),
  data.frame(gse_diff_response$aco_vs_mus10_12d) %>%
    mutate(Species1 = "aco", Species2 = "mus", Dose = "10gy", Timepoint = "12d"),
  data.frame(gse_diff_response$aco_vs_mus20_12d) %>%
    mutate(Species1 = "aco", Species2 = "mus", Dose = "20gy", Timepoint = "12d"),
  data.frame(gse_diff_response$aco_vs_musw10_12d) %>%
    mutate(Species1 = "aco", Species2 = "musw", Dose = "10gy", Timepoint = "12d"),
  data.frame(gse_diff_response$aco_vs_musw20_12d) %>%
    mutate(Species1 = "aco", Species2 = "musw", Dose = "20gy", Timepoint = "12d"))

# GSEA > dNES
tmp1 = gse_combined$T1 %>%
  dplyr::select(ID, Description, starts_with("NES")) %>%
  pivot_longer(cols=-c(ID, Description), names_pattern = "NES_(.*)(24h|12d)_(.*)", 
               names_to = c("Species", "Timepoint", "Dose"), values_to = "NES") %>%
  pivot_wider(names_from = Species, values_from = NES) %>%
  mutate(dNES_aco_vs_mus = aco - mus, dNES_aco_vs_musw = aco - musw, dNES_aco_vs_hsa = aco - hsa) %>%
  dplyr::select(ID:Dose, starts_with("dNES"))
tmp2 = gse_combined$T2 %>%
  dplyr::select(ID, Description, starts_with("NES")) %>%
  pivot_longer(cols=-c(ID, Description), names_pattern = "NES_(.*)(24h|12d)_(.*)", 
               names_to = c("Species", "Timepoint", "Dose"), values_to = "NES")%>%
  pivot_wider(names_from = Species, values_from = NES) %>%
  mutate(dNES_aco_vs_mus = aco - mus, dNES_aco_vs_musw = aco - musw)%>%
  dplyr::select(ID:Dose, starts_with("dNES"))
dNES = plyr::rbind.fill(tmp1, tmp2)

NES_dLogFC = gse_diff_combined %>%
  dplyr::select(ID, Description, Species1:Timepoint, NES) %>%
  pivot_wider(names_from = c(Species1, Species2), names_glue = "NES_dLogFC_{Species1}_vs_{Species2}", values_from = NES)

tmp3 = merge(dNES, NES_dLogFC)

cor(tmp3$dNES_aco_vs_mus, tmp3$NES_dLogFC_aco_vs_mus, use = "pairwise.complete.obs")
cor(tmp3$dNES_aco_vs_musw, tmp3$NES_dLogFC_aco_vs_musw, use = "pairwise.complete.obs")
cor(tmp3$dNES_aco_vs_hsa, tmp3$NES_dLogFC_aco_vs_hsa, use = "pairwise.complete.obs")

#### SAVE OBJECTS ####

saveRDS(des_baseline, file = paste0(paths$objects, "/des_baseline.rds"))
saveRDS(des_treatment, file = paste0(paths$objects, "/des_treatment.rds"))
saveRDS(des_combined, file = paste0(paths$objects, "/des_combined.rds"))

#### SAVE TABLES ####

# gse_treatment content hand't been converted to table yet
for (comp in names(gse_treatment)) {
  gse_treatment[[comp]] = data.frame(gse_treatment[[comp]])
}

# same for gse_baseline
for (comp in names(gse_baseline)) {
  gse_baseline[[comp]] = data.frame(gse_baseline[[comp]])
}

# same for gse_diff_response
for (comp in names(gse_diff_response)) {
  gse_diff_response[[comp]] = data.frame(gse_diff_response[[comp]])
}

##### Response #####

# Acomys conversion
gse_treatment$aco24h_10gy = convert_ids_in_string(gse_treatment$aco24h_10gy, ginfo$aco$Entrez, ginfo$aco$GeneSymbol)
gse_treatment$aco24h_20gy = convert_ids_in_string(gse_treatment$aco24h_20gy, ginfo$aco$Entrez, ginfo$aco$GeneSymbol)
gse_treatment$aco12d_10gy = convert_ids_in_string(gse_treatment$aco12d_10gy, ginfo$aco$Entrez, ginfo$aco$GeneSymbol)
gse_treatment$aco12d_20gy = convert_ids_in_string(gse_treatment$aco12d_20gy, ginfo$aco$Entrez, ginfo$aco$GeneSymbol)
# Human conversion
gse_treatment$hsa24h_10gy = convert_ids_in_string(gse_treatment$hsa24h_10gy, ginfo$hsa$Entrez, ginfo$hsa$GeneSymbol)
gse_treatment$hsa24h_20gy = convert_ids_in_string(gse_treatment$hsa24h_20gy, ginfo$hsa$Entrez, ginfo$hsa$GeneSymbol)
# Mouse conversion
gse_treatment$mus24h_10gy = convert_ids_in_string(gse_treatment$mus24h_10gy, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_treatment$mus24h_20gy = convert_ids_in_string(gse_treatment$mus24h_20gy, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_treatment$mus12d_10gy = convert_ids_in_string(gse_treatment$mus12d_10gy, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_treatment$mus12d_20gy = convert_ids_in_string(gse_treatment$mus12d_20gy, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
# Wild mouse conversion
gse_treatment$musw24h_10gy = convert_ids_in_string(gse_treatment$musw24h_10gy, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_treatment$musw24h_20gy = convert_ids_in_string(gse_treatment$musw24h_20gy, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_treatment$musw12d_10gy = convert_ids_in_string(gse_treatment$musw12d_10gy, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_treatment$musw12d_20gy = convert_ids_in_string(gse_treatment$musw12d_20gy, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)

# Aco
write.table(subset(gse_treatment$aco24h_10gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_aco24h_10vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$aco24h_20gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_aco24h_20vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$aco12d_10gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_aco12d_10vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$aco12d_20gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_aco12d_20vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
# Mus
write.table(subset(gse_treatment$mus24h_10gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_mus24h_10vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$mus24h_20gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_mus24h_20vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$mus12d_10gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_mus12d_10vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$mus12d_20gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_mus12d_20vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
# MusW
write.table(subset(gse_treatment$musw24h_10gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_musw24h_10vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$musw24h_20gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_musw24h_20vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$musw12d_10gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_musw12d_10vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$musw12d_20gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_musw12d_20vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
# Hsa
write.table(subset(gse_treatment$hsa24h_10gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_hsa24h_10vsC.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_treatment$hsa24h_20gy, p.adjust < 0.01), paste0(paths$tables, "/gsea_treatment_hsa24h_20vsC.tsv"), 
            row.names = F, sep="\t", quote=F)

# Combined
write.table(gse_combined$T1, paste0(paths$tables, "/gsea_treatment_combined24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(gse_combined$T2, paste0(paths$tables, "/gsea_treatment_combined12d.tsv"), 
            row.names = F, sep="\t", quote=F)

##### Difference in baseline #####

# All genes use the mouse ids
gse_baseline$aco_vs_hsa_24h = convert_ids_in_string(gse_baseline$aco_vs_hsa_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_baseline$aco_vs_mus_24h = convert_ids_in_string(gse_baseline$aco_vs_mus_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_baseline$aco_vs_mus_12d = convert_ids_in_string(gse_baseline$aco_vs_mus_12d, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_baseline$aco_vs_musw_24h = convert_ids_in_string(gse_baseline$aco_vs_musw_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_baseline$aco_vs_musw_12d = convert_ids_in_string(gse_baseline$aco_vs_musw_12d, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)

write.table(subset(gse_baseline$aco_vs_hsa_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_baseline_aco_vs_hsa24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_baseline$aco_vs_mus_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_baseline_aco_vs_mus24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_baseline$aco_vs_mus_12d, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_baseline_aco_vs_mus12d.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_baseline$aco_vs_musw_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_baseline_aco_vs_musw24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_baseline$aco_vs_musw_12d, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_baseline_aco_vs_musw12d.tsv"), 
            row.names = F, sep="\t", quote=F)

##### Difference in response #####

# All genes use the mouse ids
gse_diff_response$aco_vs_hsa10_24h = convert_ids_in_string(gse_diff_response$aco_vs_hsa10_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_hsa20_24h = convert_ids_in_string(gse_diff_response$aco_vs_hsa20_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_mus10_24h = convert_ids_in_string(gse_diff_response$aco_vs_mus10_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_mus20_24h = convert_ids_in_string(gse_diff_response$aco_vs_mus20_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_mus10_12d = convert_ids_in_string(gse_diff_response$aco_vs_mus10_12d, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_mus20_12d = convert_ids_in_string(gse_diff_response$aco_vs_mus20_12d, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_musw10_24h = convert_ids_in_string(gse_diff_response$aco_vs_musw10_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_musw20_24h = convert_ids_in_string(gse_diff_response$aco_vs_musw20_24h, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_musw10_12d = convert_ids_in_string(gse_diff_response$aco_vs_musw10_12d, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
gse_diff_response$aco_vs_musw20_12d = convert_ids_in_string(gse_diff_response$aco_vs_musw20_12d, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)

write.table(subset(gse_diff_response$aco_vs_hsa10_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_hsa10_24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_hsa20_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_hsa20_24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_mus10_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_mus10_24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_mus20_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_mus20_24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_mus10_12d, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_mus10_12d.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_mus20_12d, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_mus20_12d.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_musw10_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_musw10_24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_musw20_24h, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_musw20_24h.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_musw10_12d, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_musw10_12d.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(subset(gse_diff_response$aco_vs_musw20_12d, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_musw20_12d.tsv"), 
            row.names = F, sep="\t", quote=F)
write.table(gse_diff_combined, paste0(paths$tables, "/gsea_diff_combined.tsv"), 
            row.names = F, sep="\t", quote=F)

##### Count tables and DEs #####

write.csv(norm$union, paste0(paths$tables, "/cpm_combined.csv"))

write.table(des_treatment$aco_24h[,-c(2:6)], paste0(paths$tables, "/des_aco24h.tsv"), row.names = F, sep="\t", quote=F)
write.table(des_treatment$aco_12d[,-c(2:6)], paste0(paths$tables, "/des_aco12d.tsv"), row.names = F, sep="\t", quote=F)
write.table(des_treatment$mus_24h[,-c(2:6)], paste0(paths$tables, "/des_mus24h.tsv"), row.names = F, sep="\t", quote=F)
write.table(des_treatment$mus_12d[,-c(2:6)], paste0(paths$tables, "/des_mus12d.tsv"), row.names = F, sep="\t", quote=F)
write.table(des_treatment$musw_24h[,-c(2:6)], paste0(paths$tables, "/des_musw24h.tsv"), row.names = F, sep="\t", quote=F)
write.table(des_treatment$musw_12d[,-c(2:6)], paste0(paths$tables, "/des_musw12d.tsv"), row.names = F, sep="\t", quote=F)
write.table(des_treatment$hsa_24h[,-c(2:6)], paste0(paths$tables, "/des_hsa24h.tsv"), row.names = F, sep="\t", quote=F)

#### OBJECT SUMMARY ####

# aco_de(12d/24h)(/_fewna) = acomys specific genes
# aco_nes(12d/24h)(/_fewna) = acomys specific pathways
# comp(12d/24h)(/_wide) = aggregation of des_treatment used to generate aco_de
# counts = raw counts
# des_baseline = DE results comparing species baseline expression
# des_combined = DE results side by side, used to generate gse_diff_response
# des_treatment = DE results comparing treated and untreated within species
# dNES = pairwise NES differences between species, used to valiate the dLogFC > GSEA strat
# ginfo = gene info
# go_aco_specific(/2) = GO results for aco specific lists
# gse_any_sig = from gse_combined s.t. sig for at least one comparison, used for gsea_comparion figures and to make gse_diff_sig
# gse_baseline = GSEA on des_baseline
# gse_combined = used for gsea_response_similarity figure, making gse_any_sig, determing aco_nes
# gse_diff_combined = aggregation of gse_diff_response used for dLogFC>GSEA vs dGSEA
# gse_diff_response = GSEA on des_combined dLogFC
# gse_diff_sig = from gse_any_sig, used for pathways of interest
# gse_treatment = GSEA on des_treatment
# meta = meta
# NES_dLogFC = used for validation
# norm = normalized counts
# pathways_of_interest = list of GO terms of interest
# pca = pca
# rest is auxilliary

