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

setwd("/scratch/fmorandi/internal/Fathima/RNA2")

paths = list()
paths$meta = "./meta.txt"
paths$out = "./results"
paths$objects = "./results/objects"
paths$tables = "./results/tables"
paths$pathways = "./results/pathways"
paths$custom_gsets = "./results/custom_gsets"
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

# common_genes = colnames(norm$all)
genes_intersection = colnames(norm$intersection)
genes_union = colnames(norm$union)

#### PCA ####

pca = list()
tmp = prcomp(norm$intersection, center=T, scale.=T)
pca$all = merge(meta, tmp$x[,c("PC1", "PC2")], by=0)

p1 = ggplot(pca$all, aes(x=PC1, y=PC2, color=Species))+
  geom_point(size=1)
p2 = ggplot(pca$all, aes(x=PC1, y=PC2, color=Treatment))+
  geom_point(size=1)
ps = align_patches(p1, p2)
pdf(paste0(paths$out, "/pca_all.pdf"), width=0.5*w_in, height=0.25*h_in)
ps
dev.off()

#### DIFFERENTIAL EXPRESSION ####

##### Treatment #####

des_treatment = list()

# Within species and timepoint, across treatments
for (org in c("aco", "mus")) {
  # Subset meta and define design
  this_meta = subset(meta, tolower(Species) == org)
  design = model.matrix(~Treatment, data = this_meta)
  # Prepare dge
  dge = DGEList(counts[[org]][, rownames(this_meta)], samples=this_meta)
  keep = filterByExpr(dge, design=design)
  dge = dge[keep, ]
  dge = calcNormFactors(dge)
  # Fit model
  dge = estimateDisp(dge, design)
  fit = glmFit(dge, design)
  # Get results
  conts = makeContrasts("TreatmentDMBA_TPA", levels=design)
  res = as.data.frame(glmLRT(fit, contrast=conts))
  res$PAdj = p.adjust(res$PValue, method="BH")
  colnames(res) = paste0(colnames(res), "_DTvsC")
  # Save results with ginfo
  res = res %>%
    dplyr::select(-starts_with("LR")) %>%
    rownames_to_column("Row.names")
  des_treatment[[org]] = merge(ginfo[[org]], res, by.x=0, by.y="Row.names") %>%
    column_to_rownames("Row.names")
}

##### Baseline #####

des_baseline = list()

# Subset meta and define design
this_meta = meta[rownames(counts$union), ]
design = model.matrix(~0+Species, data = this_meta)
# Normalize
dge = DGEList(t(counts$union[rownames(this_meta), ]), samples=this_meta)
dge = calcNormFactors(dge)
ggplot(dge$samples, aes(x=Species, y=norm.factors))+
  geom_violin()
# <!> Therefore it seems reasonabe to normalize with TMM, most are close to 1
# Fit model
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
# Get Aco vs Mus results
conts = makeContrasts("SpeciesAco - SpeciesMus", levels=design)
res_aco_vs_mus = as.data.frame(glmLRT(fit, contrast=conts))
res_aco_vs_mus$PAdj = p.adjust(res_aco_vs_mus$PValue, method="BH")
# Combine T1 results
des_baseline = res_aco_vs_mus %>%
  dplyr::select(-starts_with("LR"))
colnames(des_baseline) = paste0(colnames(des_baseline), "_aco_vs_mus")
sum(des_baseline$PAdj_aco_vs_mus < 0.05)

rm(res, res_aco_vs_mus)

##### Define DEGs #####

th_sig_gene = 0.05
th_logfc = 1

# --- Treatment ---
for (org in names(des_treatment)) {
  des_treatment[[org]] = des_treatment[[org]] %>%
    mutate(sig = abs(logFC_DTvsC) > th_logfc & PAdj_DTvsC < th_sig_gene) %>%
    mutate(up = sig & logFC_DTvsC > 0) %>%
    mutate(down = sig & logFC_DTvsC < 0) %>%
    mutate(sig_col = interaction(sig, up)) %>%
    mutate(sig_col = fct_recode(sig_col, 
                                "Up" = "TRUE.TRUE",
                                "Down" = "TRUE.FALSE",
                                "NotSig" = "FALSE.FALSE"))
}

# --- Baseline ---
des_baseline = des_baseline %>%
  mutate(sig = abs(logFC_aco_vs_mus) > th_logfc & PAdj_aco_vs_mus < th_sig_gene) %>%
  mutate(up = sig & logFC_aco_vs_mus > 0) %>%
  mutate(down = sig & logFC_aco_vs_mus < 0) %>%
  mutate(sig_col = interaction(sig, up)) %>%
  mutate(sig_col = fct_recode(sig_col, 
                              "Up" = "TRUE.TRUE",
                              "Down" = "TRUE.FALSE",
                              "NotSig" = "FALSE.FALSE"))

##### Plot volcanos #####

# --- Treatment ---
p1 = my_volcano(des_treatment$aco, col="sig_col", v1="logFC_DTvsC", v2="PValue_DTvsC", "Acomys")
p2 = my_volcano(des_treatment$mus, col="sig_col", v1="logFC_DTvsC", v2="PValue_DTvsC", "Mouse")

p= p1 + p2 +
  plot_layout(guides = "collect", ncol=2) &
  scale_color_manual(values=c("#999999", "#5555cc", "#cc5555")) &
  guides(color="none")
ggsave(paste0(paths$out, "/volcanos_treament.png"), plot=p, width=w, height=0.4*h, units="mm")

# --- Baseline ---
p = my_volcano(des_baseline, col="sig_col", v1="logFC_aco_vs_mus", v2="PValue_aco_vs_mus", "Acomys vs Mouse")
p = p +
  scale_color_manual(values=c("#999999", "#5555cc", "#cc5555")) &
  guides(color="none")
ggsave(paste0(paths$out, "/volcanos_baseline.png"), plot=p, width=0.5*w, height=0.3*h, units="mm")

##### Combined table #####

# Bring all DE results side by side on common genes
res = data.frame()
res = des_treatment$aco %>%
  dplyr::select(GeneSymbol, logFC_DTvsC:PAdj_DTvsC) %>%
  mutate(Species = "aco") %>%
  rbind(res, .)
res = des_treatment$mus %>%
  dplyr::select(GeneSymbol, logFC_DTvsC:PAdj_DTvsC) %>%
  mutate(Species = "mus") %>%
  rbind(res)
rownames(res) = NULL
res = subset(res, GeneSymbol %in% genes_intersection)

des_combined = res %>%
  distinct(GeneSymbol, Species, .keep_all = T) %>%
  pivot_wider(names_from=Species, values_from=ends_with("DTvsC")) %>%
  mutate(Entrez = mapIds(org.Mm.eg.db, keys = GeneSymbol, keytype = "SYMBOL", column = "ENTREZID")) %>%
  column_to_rownames("GeneSymbol")
des_combined = des_baseline[genes_intersection, ] %>%
  dplyr::select(logFC_aco_vs_mus:PAdj_aco_vs_mus) %>%
  cbind(des_combined[genes_intersection, ], .)

#### CHECKPOINT 1 ####

# save.image(file=paste0(paths$objects, "/checkpoint1.Rdata"))
load(paste0(paths$objects, "/checkpoint1.Rdata"))

#### GSEA: TREATMENT ####

##### Run #####
# 
# gse_treatment = list()
# 
# # --- Mouse ---
# gse_treatment$musDTvsC = my_gse(des_treatment[["mus"]], "logFC_DTvsC", orgdb = org.Mm.eg.db)
# 
# # --- Acomys ---
# length(intersect(ginfo$aco$GeneSymbol, ginfo$mus$GeneSymbol))
# length(setdiff(ginfo$aco$GeneSymbol, ginfo$mus$GeneSymbol)) # Most genes share names with mice
# gse_treatment$acoDTvsC = my_gse(des_treatment[["aco"]], "logFC_DTvsC", orgdb = org.Mm.eg.db)
# 
# # Save results
# save(gse_treatment, file=paste0(paths$objects, "/gsea_treatment.Rdata"))
load(paste0(paths$objects, "/gsea_treatment.Rdata"))

##### Dotplots #####

ps = list()
ps[[1]]= gse_treatment$aco %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys DT vs C")
ps[[2]]=gse_treatment$mus %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign)+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Mouse DT vs C")
ps = align_patches(ps)
pdf(paste0(paths$out, "/gsea_treatment.pdf"), width=w_in*2, height=h_in)
ps
dev.off()

##### Comparison #####

gse_combined = data.frame()

for (comp in names(gse_treatment)) {
  tmp = data.frame(gse_treatment[[comp]]) %>%
    dplyr::select("ID", "Description", "NES", "p.adjust") %>%
    mutate(Comparison = comp)
  gse_combined = rbind(gse_combined, tmp)
}
rownames(gse_combined) = NULL
gse_combined = pivot_wider(gse_combined, names_from=Comparison, values_from=c(NES, p.adjust))
# LogFC similarity across species and treatments
data.frame(cor(gse_combined[, 3:4], use = "pairwise.complete.obs")) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols=-Var1, names_to = "Var2") %>%
  ggplot(., aes(x=Var1, y=Var2, fill=value, label=round(value,2)))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_viridis_c()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(fill="Correlation")
ggsave(paste0(paths$out, "/gsea_response_similarity.png"), width=0.6*w, height=0.3*h, units="mm")

# Make min pval columns
colnames(gse_combined) = str_replace_all(colnames(gse_combined), "p.adjust", "padj")
gse_combined$padj_min = pmin(gse_combined$padj_musDTvsC, gse_combined$padj_acoDTvsC, na.rm=T)

gse_any_sig = subset(gse_combined, padj_min < 0.01)
gse_any_sig = arrange(gse_any_sig, rowMeans(gse_any_sig[, c("NES_musDTvsC", "NES_acoDTvsC")], na.rm = T))

n_per_page = 40
nes_range = c(
  min(c(gse_any_sig$NES_musDTvsC, gse_any_sig$NES_acoDTvsC), na.rm = T),
  max(c(gse_any_sig$NES_musDTvsC, gse_any_sig$NES_acoDTvsC), na.rm = T))
ps = list()
for (p in 1:ceiling(nrow(gse_any_sig)/n_per_page)) {
  page_first = (p-1)*n_per_page+1
  page_ps = list()
  tmp = gse_any_sig[page_first:min(page_first+n_per_page-1, nrow(gse_any_sig)), 1:6] %>%
    mutate(ID = factor(ID, levels=ID)) %>%
    pivot_longer(
      names_to = c("Variable", "Species", "Treatment"),
      names_pattern = "(.*)_(mus|aco)(DTvsC)",
      cols = -c(ID, Description)) %>%
    mutate(Comparison = paste(Species, Treatment, sep = "_")) %>%
    pivot_wider(names_from=Variable) %>%
    mutate(Sig = stars.pval(padj)) %>%
    mutate(Title = substr(Description, 1, 100)) %>%
    mutate(Title = fct_reorder(Title, as.numeric(ID)))
  ps[[p]] = ggplot(tmp, aes(x=Comparison, y=Title, fill=NES, label=Sig))+
    geom_tile()+
    geom_text()+
    scale_fill_gradient2(low="blue", high="red", limits=nes_range)+
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle=45, hjust=1))
}

ps = align_patches(ps)
pdf(paste0(paths$out, "/gsea_comparison.pdf"), width=w_in*1.2, height=h_in)
ps
dev.off()

#### GSEA: BASELINE ####

##### Run #####

# Run gsea on pairwise logFC baseline differences
des_baseline = des_baseline %>%
  rownames_to_column("Row.names") %>%
  mutate(Entrez = mapIds(org.Mm.eg.db, keys = Row.names, keytype = "SYMBOL", column = "ENTREZID")) %>%
  column_to_rownames("Row.names")

# # --- Aco vs Mus
# gse_baseline = my_gse(des_baseline, "logFC_aco_vs_mus", org.Mm.eg.db)
# 
# save(gse_baseline, file=paste0(paths$objects, "/gsea_baseline.Rdata"))
load(paste0(paths$objects, "/gsea_baseline.Rdata"))

##### Dotplots #####

facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Mouse"
)
p1= gse_baseline %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Mouse Baseline")

pdf(paste0(paths$out, "/gsea_baseline.pdf"), width=w_in*2, height=1*h_in)
p1
dev.off()

#### CHECKPOINT 2 ####

# save.image(file=paste0(paths$objects, "/checkpoint2.Rdata"))
load(paste0(paths$objects, "/checkpoint2.Rdata"))

#### GSEA: DIFFERENCE IN RESPONSE ####

##### Run #####
# 
# # Run gsea on pairwise logFC differences
# gse_diff_response = des_combined %>%
#   mutate(logFCdiff = logFC_DTvsC_aco - logFC_DTvsC_mus) %>%
#   my_gse(., "logFCdiff", org.Mm.eg.db)
# 
# save(gse_diff_response, file=paste0(paths$objects, "/gsea_diff_response.Rdata"))
load(paste0(paths$objects, "/gsea_diff_response.Rdata"))

##### Dotplots #####

facet_titles = c(
  activated = "Higher in Acomys",
  suppressed = "Higher in Mouse"
)
p1=gse_diff_response %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dotplot(., split=".sign", showCategory = 20)+
  facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
  theme(axis.text.y = element_text(size = 10))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  ggtitle("Acomys vs Mouse")

pdf(paste0(paths$out, "/gsea_diff_response.pdf"), width=w_in*2, height=h_in)
p1
dev.off()

#### CHECKPOINT 3 ####

# save.image(file=paste0(paths$objects, "/checkpoint3.Rdata"))
load(paste0(paths$objects, "/checkpoint3.Rdata"))

#### PATHWAYS OF INTEREST ####

##### Prepare #####

gse_diff_sig = gse_any_sig[, 1:6] %>% 
  pivot_longer(
    names_to = c("Variable", "Species", "Treatment"),
    names_pattern = "(.*)_(mus|aco)(DTvsC)",
    cols = -c(ID, Description)) %>%
  mutate(Comparison = paste(Species, Treatment, sep = "_")) %>%
  pivot_wider(names_from=Variable) %>%
  mutate(Title = str_replace_all(Description, "regulation", "reg")) %>%
  mutate(Title = str_wrap(Title, 30)) %>%
  mutate(Sig = stars.pval(padj))

##### Select #####

pathways_of_interest = list(
  "repair" = c(
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
  "damage" = c(
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
    "GO:0000723" = "Telomere maintenance",
    "GO:0032200" = "Telomere organization"
  ),
  "mapk_erk" = c(
    "GO:0070374" = "positive reg. of ERK1 and ERK2 cascade",
    "GO:0070372" = "reg. of ERK1 and ERK2 cascade",
    "GO:0070371" = "ERK1 and ERK2 cascade",
    "GO:0043410" = "positive reg. of MAPK cascade"
  ),
  "wnt_signalling" = c(
    "GO:0060070" = "canonical Wnt signaling pathway",
    "GO:0060828" = "regulation of canonical Wnt signaling pathway",
    "GO:0030111" = "regulation of Wnt signaling pathway",
    "GO:0030178" = "negative regulation of Wnt signaling pathway",
    "GO:0198738" = "cell-cell signaling by wnt",
    "GO:0016055" = "Wnt signaling pathway"
  ),
  "bmp_signalling" = c(
    "GO:0030510" = "regulation of BMP signaling pathway",
    "GO:0071772" = "response to BMP",
    "GO:0071773" = "cellular response to BMP stimulus",
    "GO:0030509" = "BMP signaling pathway"
  )
  # "p53" = c(),
  # "p21" = c()
)

# Function to help search terms
gse_diff_sig %>%
  dplyr::filter(grepl("bmp", Description, ignore.case = T)) %>%
  dplyr::select(ID, Description) %>%
  distinct()

##### Plot #####

for (concept in names(pathways_of_interest)) {
  ps = list()
  for (term in names(pathways_of_interest[[concept]])) {
    tmp = gse_diff_sig %>%
      dplyr::filter(ID == term)
    if (nrow(tmp) == 0) next
    p=ggplot(tmp, aes(x=Comparison, y=NES, fill=Comparison, label=Sig))+
      geom_bar(stat="identity")+
      geom_text(aes(y=NES/2), color="white")+
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
  ggsave(plot=plt, paste0(paths$pathways, "/", concept, ".png"), width=w, height=0.5*h, units="mm")
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

my_gse_custom = function(table, logFC_column, gmt) {
  options(warn = 1)
  # Make gene list
  gene_list = table[[logFC_column]]
  names(gene_list) = table$GeneSymbol
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
    res_custom_gsets[[gmt]]$gsea = rbind(
      res_custom_gsets[[gmt]]$gsea,
      my_gse_custom(des_treatment[[comp]], "logFC_DTvsC", gmts[[gmt]]) %>%
        as.data.frame() %>%
        mutate(Sig = stars.pval(p.adjust)) %>%
        mutate(Context = comp, Comparison = "DTvsC") %>%
        mutate(Condition = paste(Context, Comparison, sep="_")))
  }
  res_custom_gsets[[gmt]]$zscores = get_geneset_zscores(norm$union, meta, gmts[[gmt]])
}

##### Plot #####

for (gmt in names(gmts)) {
  pdf(paste0(paths$custom_gsets, "/", gmt, ".pdf"), height = h_in, width = w_in)
  for (gset in unique(gmts[[gmt]]$term)) {
    gset_desc = gmts[[gmt]][gmts[[gmt]]$term == gset, "desc"]
    gset_genes = gmts[[gmt]] %>%
      dplyr::filter(term == gset) %>%
      pull(gene)
    gset_zscores = res_custom_gsets[[gmt]]$zscores
    gset_zscores = gset_zscores[, colnames(gset_zscores) %in% gset_genes]
    # GSEA NES bars
    p1 = res_custom_gsets[[gmt]]$gsea %>%
      dplyr::filter(ID == gset) %>%
      ggplot(., aes(x=Condition, y=NES, fill=Condition, label=Sig))+
      geom_bar(stat="identity")+
      geom_text(aes(y=NES/2), color="white")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      ggtitle("GSEA")
    # Average pathway zscore boxes
    p2 = meta %>% 
      mutate(setZscore = rowMeans(gset_zscores)) %>%
      # mutate(Treatment = factor(Treatment, levels = c("Con", "10Gy", "20Gy"))) %>%
      ggplot(., aes(x=Species, y=setZscore, fill=Treatment))+
      geom_boxplot()+
      ggtitle("Gene set Zscore")
    # Zscore heatmap
    tmp = data.frame(gset_zscores) %>%
      mutate(Species = meta$Species, Treatment = meta$Treatment) %>%
      mutate(Species = factor(Species, levels = unique(Species))) %>% # Lot of mess to order columns nicely
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

#### CHECKPOINT 4 ####

# save.image(file=paste0(paths$objects, "/checkpoint4.Rdata"))
load(paste0(paths$objects, "/checkpoint4.Rdata"))

#### ACOMYS SPECIFIC GENES ####

##### Prep #####

comp = rbind(
  des_treatment$aco %>%
    dplyr::select(GeneSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    mutate(Species = "ACO"),
  des_treatment$mus %>%
    dplyr::select(GeneSymbol, starts_with("logFC"), starts_with("PAdj")) %>%
    mutate(Species = "MUS")
)

comp_wide = comp %>%
  group_by(GeneSymbol, Species) %>%
  mutate(row = row_number()) %>%
  dplyr::filter(row == 1) %>%
  pivot_wider(names_from = "Species", values_from = c("logFC_DTvsC", "PAdj_DTvsC"), values_fill=NA) %>%
  dplyr::select(-row)

##### Aco specific response #####

th_high = 1
th_low = 0.5
th_pval = 0.1

regions = data.frame(
  xmin = c(-14, th_high),
  ymin = c(-th_low, -14),
  xmax = c(-th_high, 14),
  ymax = c(14, th_low),
  color = c("blue", "red")
)

ps = list()
aco_de = comp_wide %>%
  mutate(Passes = abs(logFC_DTvsC_ACO) > th_high & PAdj_DTvsC_ACO < th_pval)
ps[[1]] = ggplot(aco_de, aes(logFC_DTvsC_ACO, logFC_DTvsC_MUS, color=Passes))+
  geom_point(size=0.5)+
  geom_vline(xintercept = th_high)+
  geom_vline(xintercept = -th_high)+
  coord_cartesian(xlim=c(-8,8), ylim=c(-12,12))
aco_de = subset(aco_de, Passes) %>%
  mutate(Passes = (abs(logFC_DTvsC_MUS) < th_low & PAdj_DTvsC_MUS > th_pval)| logFC_DTvsC_MUS * logFC_DTvsC_ACO < 0)
ps[[2]] = ggplot()+
  geom_point(data = aco_de, aes(logFC_DTvsC_ACO, logFC_DTvsC_MUS, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-8,8), ylim=c(-12,12))
aco_de = subset(aco_de, Passes)
wrap_plots(ps, ncol=2)+plot_layout(guides="collect")
ggsave(paste0(paths$out, "/aco_specific_genes_filtering.png"), width = w, height = 0.35*h, units="mm")

comp_wide = comp_wide %>%
  mutate(SideCloud = !is.na(logFC_DTvsC_ACO+logFC_DTvsC_MUS) & logFC_DTvsC_ACO < -1 & logFC_DTvsC_MUS > 7)
ggplot(comp_wide, aes(logFC_DTvsC_ACO, logFC_DTvsC_MUS, color=SideCloud))+
  geom_point(size=0.2)+
  stat_cor()
ggsave(paste0(paths$out, "/side_cloud.png"), width = 0.8*w, height = 0.5*h, units="mm")

side_cloud = comp_wide %>%
  dplyr::filter(SideCloud) %>%
  pull(GeneSymbol)

##### Mus specific response #####

ps = list()
mus_de = comp_wide %>%
  mutate(Passes = abs(logFC_DTvsC_MUS) > th_high & PAdj_DTvsC_MUS < th_pval)
ps[[1]] = ggplot(mus_de, aes(logFC_DTvsC_MUS, logFC_DTvsC_ACO, color=Passes))+
  geom_point(size=0.5)+
  geom_vline(xintercept = th_high)+
  geom_vline(xintercept = -th_high)+
  coord_cartesian(xlim=c(-12,12), ylim=c(-8,8))
mus_de = subset(mus_de, Passes) %>%
  mutate(Passes = (abs(logFC_DTvsC_ACO) < th_low & PAdj_DTvsC_ACO > th_pval)| logFC_DTvsC_ACO * logFC_DTvsC_MUS < 0)
ps[[2]] = ggplot()+
  geom_point(data = mus_de, aes(logFC_DTvsC_MUS, logFC_DTvsC_ACO, color=Passes), size=0.5)+
  geom_rect(data = regions, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), fill=NA, color="black")+
  coord_cartesian(xlim=c(-12,12), ylim=c(-8,8))
mus_de = subset(mus_de, Passes)
wrap_plots(ps, ncol=2)+plot_layout(guides="collect")
ggsave(paste0(paths$out, "/mus_specific_genes_filtering.png"), width = w, height = 0.35*h, units="mm")

##### Save lists for STRING #####

dir.create(paste0(paths$out, "/lists"), showWarnings = F)

aco_de %>%
  dplyr::filter(logFC_DTvsC_ACO > 0) %>%
  pull(GeneSymbol) %>%
  write_lines(file=paste0(paths$out, "/lists/aco_only_de_up.txt"))
aco_de %>%
  dplyr::filter(logFC_DTvsC_ACO < 0) %>%
  pull(GeneSymbol) %>%
  write_lines(file=paste0(paths$out, "/lists/aco_only_de_down.txt"))
mus_de %>%
  dplyr::filter(logFC_DTvsC_MUS > 0) %>%
  pull(GeneSymbol) %>%
  write_lines(file=paste0(paths$out, "/lists/mus_only_de_up.txt"))
mus_de %>%
  dplyr::filter(logFC_DTvsC_MUS < 0) %>%
  pull(GeneSymbol) %>%
  write_lines(file=paste0(paths$out, "/lists/mus_only_de_down.txt"))
write_lines(side_cloud, file=paste0(paths$out, "/lists/side_cloud.txt"))

##### Plot heatmap #####

aco_de %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
  pivot_longer(cols=-GeneSymbol, names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_genes_heatmap.png"), width = w, height = 2*h, units="mm")
# ^ too big so gotta restrict further to top genes

aco_de %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  mutate(dLogFC = logFC_DTvsC_ACO - logFC_DTvsC_MUS) %>%
  mutate(dLogFC_sign = sign(dLogFC)) %>%
  group_by(dLogFC_sign) %>%
  slice_max(abs(dLogFC), n=30) %>%
  # mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, dLogFC_sign)) %>%
  pivot_longer(cols=-c(GeneSymbol, dLogFC, dLogFC_sign), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_genes_heatmap_topN.png"), width = w, height = 1*h, units="mm")

aco_de %>%
  ungroup() %>%
  dplyr::select(GeneSymbol:Passes, -Passes) %>%
  mutate(dLogFC = logFC_DTvsC_ACO - logFC_DTvsC_MUS) %>%
  mutate(dLogFC_sign = sign(dLogFC)) %>%
  dplyr::filter(!grepl("krt", GeneSymbol, ignore.case = T))%>%
  dplyr::filter(!grepl("Gm", GeneSymbol))%>%
  group_by(dLogFC_sign) %>%
  slice_max(abs(dLogFC), n=30) %>%
  # mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, dLogFC_sign)) %>%
  pivot_longer(cols=-c(GeneSymbol, dLogFC, dLogFC_sign), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_genes_heatmap_topN_filtered.png"), width = w, height = 1*h, units="mm")

##### Gene Ontology #####

my_go = function(fg, bg, org) {
  res = enrichGO(
    gene          = fg,
    universe      = bg,
    OrgDb         = org,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1.1,
    qvalueCutoff = 1.1,
    readable      = TRUE,
    keyType = "SYMBOL")
  return(res)
}

go_aco_specific = list()
bg = comp_wide$GeneSymbol

go_aco_specific$up = aco_de %>%
  dplyr::filter(logFC_DTvsC_ACO > 0) %>%
  pull(GeneSymbol) %>%
  my_go(., bg, org.Mm.eg.db)
go_aco_specific$down = aco_de %>%
  dplyr::filter(logFC_DTvsC_ACO < 0) %>%
  pull(GeneSymbol) %>%
  my_go(., bg, org.Mm.eg.db)
go_aco_specific$side_cloud = my_go(side_cloud, bg, org.Mm.eg.db)

# save(go_aco_specific, file=paste0(paths$objects, "/go_aco_specific.Rdata"))
load(paste0(paths$objects, "/go_aco_specific.Rdata"))

go_aco_specific = rbind(
  as.data.frame(go_aco_specific$up) %>%
    mutate(Direction = "Up"),
  as.data.frame(go_aco_specific$down) %>%
    mutate(Direction = "Down"),
  as.data.frame(go_aco_specific$side_cloud) %>%
    mutate(Direction = "SideCloud")
)

go_aco_specific2 = go_aco_specific %>%
  dplyr::filter(p.adjust < 0.05)
table(go_aco_specific2$Direction)
# <!> There is something to plot here!

go_aco_specific2 %>%
  dplyr::filter(Direction != "SideCloud") %>%
  group_by(Direction) %>%
  mutate(Direction2 = ifelse(Direction == "Up", 1, -1))%>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
  mutate(BgRatio = sapply(BgRatio, function(x) eval(parse(text=x)))) %>%
  mutate(enrichment = GeneRatio / BgRatio) %>%
  slice_min(p.adjust, n=30) %>%
  mutate(Description = if_else(str_length(Description) > 60, str_c(str_sub(Description, 1, 57), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust)*Direction2)) %>%
  ggplot(., aes(Description, -log10(pvalue)*Direction2, fill=enrichment*Direction2))+
  geom_bar(stat="identity")+
  coord_flip()+
  scale_fill_gradient2(low="blue", high="red")+
  labs(y = "-log10(p) * Direction", fill = "Enrichment * Direction")+
  theme(legend.position = "top")
ggsave(paste0(paths$out, "/aco_specific_genes_go.png"), width = w, height = 1*h, units="mm")

#### ACOMYS SPECIFIC PATHWAYS ####

##### Aco specific response #####

names(gse_combined)
ggplot(gse_combined, aes(NES_acoDTvsC, -log10(padj_acoDTvsC), color=padj_acoDTvsC<0.01))+
  geom_point()
sum(gse_combined$padj_acoDTvsC<0.01, na.rm=T)
nrow(gse_combined)

min(abs(gse_combined[gse_combined$padj_acoDTvsC<0.01, "NES_acoDTvsC"]), na.rm=T)
summary(abs(gse_combined[gse_combined$padj_acoDTvsC >0.01, "NES_acoDTvsC"]), na.rm=T)
th_high_nes = 1.5
th_low_nes = 0.0 # <!> changed from 1 because we need more stringency here
th_pval_nes = 0.05

ggplot(gse_combined, aes(NES_acoDTvsC, NES_musDTvsC, color=padj_acoDTvsC<0.01 | padj_musDTvsC<0.01))+
  geom_point(size=0.2)

aco_nes = gse_combined %>%
  dplyr::filter(abs(NES_acoDTvsC) > th_high_nes & padj_acoDTvsC < th_pval_nes) %>%
  dplyr::filter((abs(NES_musDTvsC) < th_low_nes & padj_musDTvsC > th_pval_nes) | NES_musDTvsC * NES_acoDTvsC < 0) %>%
  dplyr::select(ID:padj_acoDTvsC)

##### Plot examples #####

regions_nes = data.frame(
  xmin = c(-5, th_high_nes),
  ymin = c(-th_low_nes, -5),
  xmax = c(-th_high_nes, 5),
  ymax = c(5, th_low_nes),
  color = c("blue", "red")
)

p1 = ggplot()+
  geom_point(data = gse_combined, aes(NES_acoDTvsC, NES_musDTvsC), size=0.2)+
  geom_rect(data = regions_nes, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, color=color), fill=NA)+
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3))+
  scale_color_manual(values=c("blue", "red"), guide="none")
p2 = ggplot()+
  geom_text_repel(data = aco_nes, aes(NES_acoDTvsC, NES_musDTvsC, label=str_wrap(Description, 30)), size=2)+
  geom_point(data = aco_nes, aes(NES_acoDTvsC, NES_musDTvsC), size=0.2)+
  geom_rect(data = regions_nes, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, color=color), fill=NA)+
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3))+
  scale_color_manual(values=c("blue", "red"), guide="none")
p=p1+p2
ggsave(plot=p, paste0(paths$out, "/aco_specific_gsea_filtering.png"), width = w, height = 0.4*h, units="mm")

##### Plot heatmap #####

aco_nes %>%
  mutate(Description = fct_reorder(Description, NES_acoDTvsC)) %>%
  pivot_longer(cols=-c(ID, Description), names_pattern = "(.*)_(aco|mus)(.*)", names_to = c("Var", "Species", "Treatment")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(Species, levels=c("aco", "mus"))) %>%
  arrange(Species) %>%
  mutate(Condition = paste(Species, Treatment, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(padj)) %>%
  ggplot(., aes(Condition, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_gsea_heatmap.png"), width = 2*w, height = 5*h, units="mm")

aco_nes %>%
  mutate(dNES = NES_acoDTvsC - NES_musDTvsC) %>%
  mutate(dNES_sign = sign(dNES)) %>%
  group_by(dNES_sign) %>%
  slice_max(abs(dNES), n=30) %>%
  # mutate(Description = fct_reorder(Description, NES_acoDTvsC)) %>%
  mutate(Description = fct_reorder(Description, dNES)) %>%
  pivot_longer(cols=-c(ID, Description, dNES, dNES_sign), names_pattern = "(.*)_(aco|mus)(.*)", names_to = c("Var", "Species", "Treatment")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(Species, levels=c("aco", "mus"))) %>%
  arrange(Species) %>%
  mutate(Condition = paste(Species, Treatment, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(padj)) %>%
  ggplot(., aes(Condition, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/aco_specific_gsea_heatmap_topN.png"), width = 2*w, height = 1*h, units="mm")

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
# ! Id conversion error: nMyc was matched to sMyc ! fix manually
oncokb[oncokb$MouseSymbol == "Mycs", "MouseSymbol"] = "Mycn"

gmts$tsgenes = data.frame(
  "term" = "TSGene",
  "gene" = tsgenes$MouseSymbol)
gmts$oncokb_tsup = data.frame(
  "term" = "OncoKB tumor suppressor",
  "gene" = oncokb[oncokb$`Is Tumor Suppressor Gene` == "Yes", "MouseSymbol"])
gmts$oncokb_onco = data.frame(
  "term" = "OncoKB oncogene",
  "gene" = oncokb[oncokb$`Is Oncogene` == "Yes", "MouseSymbol"])
gmts$carcinoma = read.gmt("../extra_gene_sets/carcinoma.gmt")

##### Zscores #####

# --- TSG ---
tmp = norm$union[rownames(meta), intersect(tsgenes$MouseSymbol, colnames(norm$union))]
tmp = apply(tmp, 2, scale)
tmp = meta %>%
  mutate(zScore = rowMeans(tmp))

ggplot(tmp, aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of TSG expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_TSG.png"), width = 0.6*w, height = 0.3*h, units="mm")

# --- OncoKB: oncogenes ---

which_genes = oncokb[oncokb$`Is Oncogene` == "Yes", "MouseSymbol"]
tmp = norm$union[rownames(meta), intersect(which_genes, colnames(norm$union))]
tmp = apply(tmp, 2, scale)
tmp = meta %>%
  mutate(zScore = rowMeans(tmp))

ggplot(tmp, aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of oncogene (OncoKB) expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_oncoKB_oncogene.png"), width = 0.6*w, height = 0.3*h, units="mm")

# --- OncoKB: tumor suppressors ---

which_genes = oncokb[oncokb$`Is Tumor Suppressor Gene` == "Yes", "MouseSymbol"]
tmp = norm$union[rownames(meta), intersect(which_genes, colnames(norm$union))]
tmp = apply(tmp, 2, scale)
tmp = meta %>%
  mutate(zScore = rowMeans(tmp))

ggplot(tmp, aes(Species, zScore, fill=Treatment)) +
  geom_boxplot()+
  ggtitle("Mean z-score of tumor suppressors (OncoKB) expression")
ggsave(paste0(paths$out, "/detail_onco/zscore_oncoKB_tsupp.png"), width = w, height = 0.4*h, units="mm")

# --- Carcinoma ---
ps = list()
for (subterm in unique(gmts$carcinoma$term)) {
  gmt = gmts$carcinoma[gmts$carcinoma$term == subterm, ]
  which_genes = unique(gmt$gene)
  tmp = norm$union[rownames(meta), intersect(which_genes, colnames(norm$union))]
  tmp = apply(tmp, 2, scale)
  ps[[length(ps)+1]] = meta %>%
    mutate(zScore = rowMeans(tmp)) %>%
    ggplot(., aes(Species, zScore, fill=Treatment)) +
    geom_boxplot()+
    ggtitle(subterm)
}
wrap_plots(ps, ncol=2, guides = "collect")
ggsave(paste0(paths$out, "/detail_onco/zscore_carcinoma.png"), width = w, height = 0.6*h, units="mm")

##### Aco specific ones #####

aco_de$isTS = aco_de$GeneSymbol %in% tsgenes$MouseSymbol
aco_de$isTS2 = aco_de$GeneSymbol %in% oncokb[oncokb$`Is Tumor Suppressor Gene` == "Yes", "MouseSymbol"]
aco_de$isOncogene = aco_de$GeneSymbol %in% oncokb[oncokb$`Is Oncogene` == "Yes", "MouseSymbol"]

# --- TSG ---
aco_de %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isTS) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())+
  ggtitle("Tumor suppressors (TSG)\nwith unique Acomys response")
ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap_TSG.png"), width = 0.6*w, height = 1*h, units="mm")

# --- OncoKB: tsupp ---
aco_de %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isTS2) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())+
  ggtitle("Tumor suppressors (OncoKB)\nwith unique Acomys response")
ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap_OncoKB_tsupp.png"), width = 0.6*w, height = 0.25*h, units="mm")

# --- OncoKB: oncogene ---
aco_de %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isOncogene) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(tolower(Species), levels=c("aco", "mus"))) %>%
  arrange(Species, Comparison) %>%
  mutate(Condition = paste(Species, Comparison, sep="_")) %>%
  mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())+
  ggtitle("Oncogenes (OncoKB)\nwith unique Acomys response")
ggsave(paste0(paths$out, "/detail_onco/aco_specific_genes_heatmap_OncoKB_oncogene.png"), width = 0.6*w, height = 1*h, units="mm")

# # --- Carcinoma: increased skin carcinoma ---
# length(unique(gmts$carcinoma$gene))
# merge(aco_de, gmts$carcinoma, by.x="GeneSymbol", by.y="gene")%>%
#   ungroup() %>%
#   dplyr::select(-Passes) %>%
#   mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
#   pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
#   pivot_wider(names_from = Var, values_from = value) %>%
#   mutate(Species = factor(tolower(Species), levels=c("aco", "mus"))) %>%
#   arrange(Species, Comparison) %>%
#   mutate(Condition = paste(Species, Comparison, sep="_")) %>%
#   mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
#   mutate(Sig = stars.pval(PAdj)) %>%
#   ggplot(., aes(Condition, GeneSymbol, fill=logFC, label=Sig))+
#   geom_tile()+
#   scale_fill_gradient2(low="blue", high="red")+
#   theme(axis.text.x = element_text(angle=45, hjust=1),
#         axis.title = element_blank())+
#   facet_wrap(~term)
#   ggtitle("Oncogenes (OncoKB)\nwith unique Acomys response")

##### Gsea #####

tmp = des_combined %>%
  rownames_to_column("GeneSymbol")
res_custom_gsets$carcinoma = list()
res_custom_gsets$carcinoma$gsea = rbind(
  my_gse_custom(tmp, "logFC_DTvsC_aco", gmts$carcinoma) %>%
    as.data.frame() %>%
    mutate(Sig = stars.pval(p.adjust)) %>%
    mutate(Context = "aco", Comparison = "DTvsC") %>%
    mutate(Condition = "A. cah\nDTvsC"),
  my_gse_custom(tmp, "logFC_DTvsC_mus", gmts$carcinoma) %>%
    as.data.frame() %>%
    mutate(Sig = stars.pval(p.adjust)) %>%
    mutate(Context = "mus", Comparison = "DTvsC") %>%
    mutate(Condition = "M. mus\nDTvsC"))

View(res_custom_gsets$carcinoma$gsea)
ggplot(res_custom_gsets$carcinoma$gsea, aes(Condition, ID, fill=NES, label=Sig))+
  geom_tile()+
  geom_text()+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme(axis.title=element_blank())
ggsave(paste0(paths$out, "/detail_onco/gsea_carcinoma.pdf"), width = 0.9*w, height = 0.2*h, units="mm")

View(gmts$carcinoma)

gmts$all_onco = rbind(
  gmts$tsgenes,
  gmts$oncokb_tsup,
  gmts$oncokb_onco,
  gmts$carcinoma
)
res_custom_gsets$all_onco = list()
res_custom_gsets$all_onco$gsea = rbind(
  my_gse_custom(tmp, "logFC_DTvsC_aco", gmts$all_onco) %>%
    as.data.frame() %>%
    mutate(Sig = stars.pval(p.adjust)) %>%
    mutate(Context = "aco", Comparison = "DTvsC") %>%
    mutate(Condition = "A. cah\nDTvsC"),
  my_gse_custom(tmp, "logFC_DTvsC_mus", gmts$all_onco) %>%
    as.data.frame() %>%
    mutate(Sig = stars.pval(p.adjust)) %>%
    mutate(Context = "mus", Comparison = "DTvsC") %>%
    mutate(Condition = "M. mus\nDTvsC"))
ggplot(res_custom_gsets$all_onco$gsea, aes(Condition, ID, fill=NES, label=Sig))+
  geom_tile()+
  geom_text()+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme(axis.title=element_blank())
ggsave(paste0(paths$out, "/detail_onco/gsea_all_onco.pdf"), width = 0.9*w, height = 0.3*h, units="mm")

#### PATHWAY MAP DETAIL ####

source("../utils.R")

##### Ignoring NA #####

logFCs = des_combined %>%
  dplyr::select(starts_with("logFC")) %>%
  dplyr::select(contains("vsC"))

pathwy_mapk = read_kegg_map("../pathway_maps/mmu04010.xml")
pathwy_mapk = make_scores(pathwy_mapk, logFCs, "logFC_DTvsC_aco")
plot_kegg_map(pathwy_mapk, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/mapk.png"), width = 2*w, height = 1*h, units="mm")

pathwy_wnt = read_kegg_map("../pathway_maps/mmu04310.xml")
pathwy_wnt = make_scores(pathwy_wnt, logFCs, "logFC_DTvsC_aco")
plot_kegg_map(pathwy_wnt, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/wnt.png"), width = 2*w, height = 1*h, units="mm")

pathwy_tgfb = read_kegg_map("../pathway_maps/mmu04350.xml")
pathwy_tgfb = make_scores(pathwy_tgfb, logFCs, "logFC_DTvsC_aco")
plot_kegg_map(pathwy_tgfb, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/tgfb.png"), width = 2*w, height = 1*h, units="mm")

##### NA to 0 #####

logFCs = comp_wide %>%
  dplyr::select(starts_with("logFC")) %>%
  dplyr::select(contains("vsC"))%>%
  column_to_rownames("GeneSymbol")
logFCs[is.na(logFCs)] = 0

pathwy_mapk = read_kegg_map("../pathway_maps/mmu04010.xml")
pathwy_mapk = make_scores(pathwy_mapk, logFCs, "logFC_DTvsC_ACO")
plot_kegg_map(pathwy_mapk, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/mapk_na0.png"), width = 2*w, height = 1*h, units="mm")

pathwy_wnt = read_kegg_map("../pathway_maps/mmu04310.xml")
pathwy_wnt = make_scores(pathwy_wnt, logFCs, "logFC_DTvsC_ACO")
plot_kegg_map(pathwy_wnt, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/wnt_na0.png"), width = 2*w, height = 1*h, units="mm")

pathwy_tgfb = read_kegg_map("../pathway_maps/mmu04350.xml")
pathwy_tgfb = make_scores(pathwy_tgfb, logFCs, "logFC_DTvsC_ACO")
plot_kegg_map(pathwy_tgfb, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/tgfb_na0.png"), width = 2*w, height = 1*h, units="mm")

##### Baseline #####

make_scores2 = function(lst, logFCs) {
  lst$nodes$EntrezId1 = str_extract(lst$nodes$keggId, "^mmu:(\\d+)", group=1)
  lst$nodes$Symbol1 = entrez_to_symbol(lst$nodes$EntrezId1)
  valid = !is.na(lst$nodes$Symbol1) & lst$nodes$Symbol1 %in% names(logFCs)
  which_genes = lst$nodes$Symbol1[valid]
  tmp = logFCs[which_genes]
  lst$nodes[valid, "score"] = tmp
  return(lst)
}

logFCs = des_baseline[, "logFC_aco_vs_mus"]
names(logFCs) = rownames(des_baseline)

pathwy_mapk = read_kegg_map("../pathway_maps/mmu04010.xml")
pathwy_mapk = make_scores2(pathwy_mapk, logFCs)
plot_kegg_map(pathwy_mapk, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/mapk_baseline.png"), width = 2*w, height = 1*h, units="mm")

pathwy_wnt = read_kegg_map("../pathway_maps/mmu04310.xml")
pathwy_wnt = make_scores2(pathwy_wnt, logFCs)
plot_kegg_map(pathwy_wnt, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/wnt_baseline.png"), width = 2*w, height = 1*h, units="mm")

pathwy_tgfb = read_kegg_map("../pathway_maps/mmu04350.xml")
pathwy_tgfb = make_scores2(pathwy_tgfb, logFCs)
plot_kegg_map(pathwy_tgfb, "score")
ggsave(paste0(paths$out, "/detail_pathway_maps/tgfb_baseline.png"), width = 2*w, height = 1*h, units="mm")

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

meta %>%
  mutate(Tert = norm$union[rownames(meta), "Tert"]) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "Aco", "M. mus" = "Mus")) %>%
  ggplot(., aes(Species, Tert, fill=Treatment))+
  geom_boxplot()+
  labs(y = "Tert expression [log(CPM)]")+
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$out, "/detail_telo/tert.pdf"), width = 0.6*w, height=0.3*h, units="mm")

##### GO terms zcores #####

go_sets = as.data.frame(GOTERM) %>%
  dplyr::select(-1) %>%
  dplyr::filter(grepl("telom", Term) | grepl("shelterin", Term)) %>%
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
  dplyr::filter(n >= 5, n <= 500)
go_members = go_members[go_members$GO %in% tmp$GO, ]

all.equal(rownames(norm$union), rownames(meta))
telo_zscores = get_geneset_zscores(norm$union, meta, go_members)
rownames(telo_zscores) = rownames(meta)

ps = list()
for (goid in unique(go_members$GO)) {
  godesc = go_sets[go_sets$go_id == goid, "Term"]
  genes_oi = go_members %>%
    dplyr::filter(GO == goid) %>%
    dplyr::filter(SYMBOL %in% colnames(telo_zscores)) %>%
    pull(SYMBOL)
  ps[[goid]] = meta %>%
    mutate(meanZscore = rowMeans(telo_zscores[, genes_oi])) %>%
    ggplot(., aes(Species, meanZscore, fill=Treatment)) +
    geom_boxplot()+
    ggtitle(str_wrap(godesc, width=35))
}

plt = wrap_plots(ps, guides = "collect", ncol=4)&
  theme(axis.title=element_text(size=8),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "vertical")
ggsave(plot=plt, paste0(paths$out, "/detail_telo/go_terms_zcores.pdf"), width=1.2*w, height=2.2*h, units="mm")

##### GO terms GSEA #####

gse_combined %>%
  dplyr::filter(ID %in% go_sets$go_id) %>%
  drop_na()%>%
  dplyr::select(-padj_min) %>%
  mutate(Description = fct_reorder(Description, NES_musDTvsC))%>%
  pivot_longer(cols=-c(ID, Description)) %>%
  mutate(Species = str_extract(name, "(mus|aco)DTvsC", group=1)) %>%
  mutate(name = str_remove(name, "_.*")) %>%
  pivot_wider(names_from = "name") %>%
  mutate(sig = stars.pval(padj))%>%
  ggplot(., aes(Species, Description, fill=NES, label=sig))+
  geom_tile()+
  geom_text()+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme(axis.title = element_blank())
ggsave(paste0(paths$out, "/detail_telo/go_terms_gsea.pdf"), width=0.9*w, height=0.8*h, units="mm")

##### Regulators #####

this_meta = meta %>%
  dplyr::filter(rownames(.) %in% rownames(counts$union))
telo_regulators$gene = telo_regulators$GeneId
teloreg_zscores = get_geneset_zscores(norm$union, this_meta, telo_regulators)
rownames(teloreg_zscores) = rownames(this_meta)
this_meta %>%
  mutate(meanZscore = rowMeans(teloreg_zscores)) %>%
  ggplot(., aes(Species, meanZscore, fill=Species)) +
  geom_boxplot()

design = model.matrix(~0+Species, data = this_meta)
dge = DGEList(t(counts$union[rownames(this_meta), ]), samples=this_meta)
keep = filterByExpr(dge, design=design)
dge = dge[keep, ]
dge = calcNormFactors(dge)
tmp = t(cpm(dge, normalized.lib.sizes=T, log=T))
teloreg_zscores = get_geneset_zscores(tmp, this_meta, telo_regulators)
rownames(teloreg_zscores) = rownames(this_meta)
this_meta %>%
  mutate(meanZscore = rowMeans(teloreg_zscores)) %>%
  ggplot(., aes(Species, meanZscore, fill=Species)) +
  geom_boxplot()

#### WNT DETAIL ####

##### Load lists #####

# Note: best to focus on bcat depenent in general
bcat_dependent = read.csv("../extra_gene_sets/bcat_dependent.csv")
tcf_dependent = read.csv("../extra_gene_sets/tcf_dependent.csv")
tcf_independent = read.csv("../extra_gene_sets/tcf_independent.csv")
wnt_targets = read.csv("../extra_gene_sets/wnt_targets.csv", header = F)
colnames(wnt_targets) = c("Symbol", "Desc", "Uniprot")

# Convert human gene symbols to mouse
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human_to_mouse = getLDS(
  attributes = c("hgnc_symbol"), 
  filters = "hgnc_symbol", 
  values = bcat_dependent$symbol, 
  mart = human, 
  attributesL = c("mgi_symbol"), 
  martL = mouse, 
  uniqueRows=T) 
human_to_mouse = human_to_mouse %>%
  distinct(HGNC.symbol, .keep_all = T) %>%
  distinct(MGI.symbol, .keep_all = T) %>%
  column_to_rownames("HGNC.symbol")
bcat_dependent$MouseSymbol = human_to_mouse[bcat_dependent$symbol, "MGI.symbol"]
tcf_dependent$MouseSymbol = human_to_mouse[tcf_dependent$symbol, "MGI.symbol"]
tcf_independent$MouseSymbol = human_to_mouse[tcf_independent$symbol, "MGI.symbol"]

intersect(bcat_dependent$MouseSymbol, wnt_tagets$Symbol)

# --- Properly with GSEA ---
gmts$wnt_target = read.csv("../extra_gene_sets/bcat_dependent.csv")
gmts$wnt_target$gene = human_to_mouse[gmts$wnt_target$symbol, "MGI.symbol"]
gmts$wnt_target[, "term"] = "Upregulated by Wnt"
gmts$wnt_target[gmts$wnt_target$logFC.WT < 0, "term"] = "Downregulated by Wnt"
gmts$wnt_target = drop_na(gmts$wnt_target)
gmts$wnt_target = gmts$wnt_target[, c("term", "gene")]

gmts$wnt_target_tcfdep = read.csv("../extra_gene_sets/tcf_dependent.csv")
gmts$wnt_target_tcfdep$gene = human_to_mouse[gmts$wnt_target_tcfdep$symbol, "MGI.symbol"]
gmts$wnt_target_tcfdep[, "term"] = "Upregulated by Wnt"
gmts$wnt_target_tcfdep[gmts$wnt_target_tcfdep$logFC.WT < 0, "term"] = "Downregulated by Wnt"
gmts$wnt_target_tcfdep = drop_na(gmts$wnt_target_tcfdep)
gmts$wnt_target_tcfdep = gmts$wnt_target_tcfdep[, c("term", "gene")]

##### Plot #####

# Main TCF/LEF factors
des_combined_fcNp = des_combined %>% # DEs combines, FCs and pvalues
  dplyr::select(starts_with("logFC"), starts_with("PAdj")) %>%
  dplyr::select(contains("vsC")) %>%
  rownames_to_column("Symbol") %>%
  pivot_longer(cols=-Symbol) %>%
  mutate(Species = str_extract(name, "_(aco|mus)", group=1))%>%
  mutate(name = str_remove(name, "_DTvsC_(aco|mus)")) %>%
  pivot_wider(names_from = name, values_from = value)%>%
  mutate(Sig = stars.pval(PAdj))

# All TCF/LEF factors
des_combined_fcNp %>%
  dplyr::filter(grepl("^(tcf|lef1)", Symbol, ignore.case = T)) %>%
  ggplot(.,aes(Species, Symbol, fill=logFC, label=Sig))+
  geom_tile()+
  geom_text()+
  scale_fill_gradient2(low="blue", high="red")

# Main TCF/LEF factors
des_combined_fcNp %>%
  dplyr::filter(Symbol %in% c("Tcf7", "Tcf7l1", "Tcf7l2", "Lef1")) %>%
  ggplot(.,aes(Species, Symbol, fill=logFC, label=Sig))+
  geom_tile()+
  geom_text()+
  scale_fill_gradient2(low="blue", high="red")

wnt_targets = merge(des_combined_fcNp, wnt_targets, by="Symbol")
bcat_dependent = bcat_dependent %>%
  dplyr::filter(!is.na(MouseSymbol))%>%
  dplyr::select(-c(EntrezID:gene)) %>%
  merge(des_combined_fcNp, ., by.x="Symbol", by.y="MouseSymbol")
tcf_dependent = tcf_dependent %>%
  dplyr::filter(!is.na(MouseSymbol))%>%
  dplyr::select(-c(EntrezID:gene)) %>%
  merge(des_combined_fcNp, ., by.x="Symbol", by.y="MouseSymbol")
tcf_independent = tcf_independent %>%
  dplyr::filter(!is.na(MouseSymbol))%>%
  dplyr::select(-c(EntrezID:gene)) %>%
  merge(des_combined_fcNp, ., by.x="Symbol", by.y="MouseSymbol")

ggplot(bcat_dependent, aes(logFC.WT, logFC))+
  geom_point()+
  facet_wrap(~Species)+
  stat_cor()
ggplot(tcf_dependent, aes(logFC.WT, logFC))+
  geom_point()+
  facet_wrap(~Species)+
  stat_cor()
ggplot(tcf_independent, aes(logFC.WT, logFC))+
  geom_point()+
  facet_wrap(~Species)+
  stat_cor()

bcat_dependent %>%
  mutate(signCHIR = ifelse(logFC.WT > 0, "Upregulated by Wnt", "Downregulated by Wnt")) %>%
  ggplot(., aes(Species, logFC))+
  geom_boxplot()+
  facet_wrap(~signCHIR)
tcf_dependent %>%
  mutate(signCHIR = ifelse(logFC.WT > 0, "Upregulated by Wnt", "Downregulated by Wnt")) %>%
  ggplot(., aes(Species, logFC))+
  geom_boxplot()+
  facet_wrap(~signCHIR)
tcf_independent %>%
  mutate(signCHIR = ifelse(logFC.WT > 0, "Upregulated by Wnt", "Downregulated by Wnt")) %>%
  ggplot(., aes(Species, logFC))+
  geom_boxplot()+
  facet_wrap(~signCHIR)

# bcat_dependent = bcat_dependent %>%
#   mutate(signCHIR = ifelse(logFC.WT > 0, "Upregulated by Wnt", "Downregulated by Wnt"))
# tmp1 = bcat_dependent %>%
#   group_by(signCHIR, Species) %>%
#   summarize(meanLogFC = mean(logFC), sdLogFC = sd(logFC))
# ggplot()+
#   geom_bar(data=tmp1, aes(Species, meanLogFC), stat="identity")+
#   geom_jitter(data=bcat_dependent, aes(Species, logFC), width = 0.2)+
#   facet_wrap(~signCHIR)

gene_order = bcat_dependent %>%
  dplyr::filter(Species == "mus") %>%
  group_by(sign(logFC)) %>%
  arrange(logFC) %>%
  pull(Symbol)
bcat_dependent %>%
  mutate(signCHIR = ifelse(logFC.WT > 0, "Up", "Down")) %>%
  mutate(Symbol = factor(Symbol, levels=gene_order)) %>%
  ggplot(., aes(Species, Symbol, fill=logFC, label=Sig))+
  geom_tile()+
  geom_text()+
  facet_wrap(~signCHIR, scales="free")+
  scale_fill_gradient2(low="blue", high="red")
tcf_dependent %>%
  mutate(signCHIR = ifelse(logFC.WT > 0, "Up", "Down")) %>%
  ggplot(., aes(Species, Symbol, fill=logFC))+
  geom_tile()+
  facet_wrap(~signCHIR, scales="free")+
  scale_fill_gradient2(low="blue", high="red")
tcf_independent %>%
  mutate(signCHIR = ifelse(logFC.WT > 0, "Up", "Down")) %>%
  ggplot(., aes(Species, Symbol, fill=logFC))+
  geom_tile()+
  facet_wrap(~signCHIR, scales="free")+
  scale_fill_gradient2(low="blue", high="red")

##### Gsea #####

View(gmts$wnt_target)

tmp = des_combined %>%
  rownames_to_column("GeneSymbol")
res_custom_gsets$wnt_target = list()
res_custom_gsets$wnt_target$gsea = rbind(
  my_gse_custom(tmp, "logFC_DTvsC_aco", gmts$wnt_target) %>%
    as.data.frame() %>%
    mutate(Sig = stars.pval(p.adjust)) %>%
    mutate(Context = "aco", Comparison = "DTvsC") %>%
    mutate(Condition = "A. cah\nDTvsC"),
  my_gse_custom(tmp, "logFC_DTvsC_mus", gmts$wnt_target) %>%
    as.data.frame() %>%
    mutate(Sig = stars.pval(p.adjust)) %>%
    mutate(Context = "mus", Comparison = "DTvsC") %>%
    mutate(Condition = "M. mus\nDTvsC"))
res_custom_gsets$wnt_target_tcfdep = list()
res_custom_gsets$wnt_target_tcfdep$gsea = rbind(
  my_gse_custom(tmp, "logFC_DTvsC_aco", gmts$wnt_target_tcfdep) %>%
    as.data.frame() %>%
    mutate(Sig = stars.pval(p.adjust)) %>%
    mutate(Context = "aco", Comparison = "DTvsC") %>%
    mutate(Condition = "A. cah\nDTvsC"),
  my_gse_custom(tmp, "logFC_DTvsC_mus", gmts$wnt_target_tcfdep) %>%
    as.data.frame() %>%
    mutate(Sig = stars.pval(p.adjust)) %>%
    mutate(Context = "mus", Comparison = "DTvsC") %>%
    mutate(Condition = "M. mus\nDTvsC"))

# View(res_custom_gsets$wnt_target$gsea)
ggplot(res_custom_gsets$wnt_target$gsea, aes(Condition, NES, fill=Condition, label=Sig))+
  geom_bar(stat="identity")+
  geom_text(aes(y=NES/2), color="white")+
  facet_wrap(~Description)
ggplot(res_custom_gsets$wnt_target2$gsea, aes(Condition, NES, fill=Condition, label=Sig))+
  geom_bar(stat="identity")+
  geom_text(aes(y=NES/2), color="white")+
  facet_wrap(~Description)

#### DIFFERENTIAL REGULATORS ####

##### Write seqs #####

# Extract sequence around TSS of aco only up and down genes (done same for mus for uniformity)
paths$seqs = paste0(paths$out, "/seqs")
dir.create(paths$seqs, showWarnings = F)
# script to extract tss from gff3 file in in Documents/scripts
# bedtools slop -l 500 -r 0 -s -i acahirinus_TSS.bed -g acahirinus.scaffold.fa.fai > acahirinus_TSS500.bed
# bedtools slop -l 500 -r 0 -s -i Mus_musculus.GRCm39.108_TSS.bed -g Mus_musculus.GRCm39.dna.primary_assembly.fa.fai > Mus_musculus.GRCm39.108_TSS500.bed
# bedtools getfasta -name -fi acahirinus.scaffold.fa -bed acahirinus_TSS500.bed > acahirinus_TSS500.fa
# bedtools getfasta -name -fi Mus_musculus.GRCm39.dna.primary_assembly.fa -bed Mus_musculus.GRCm39.108_TSS500.bed > Mus_musculus.GRCm39.108_TSS500.fa

# bedtools slop -l 2000 -r 0 -s -i acahirinus_TSS.bed -g acahirinus.scaffold.fa.fai > acahirinus_TSS2000.bed
# bedtools slop -l 2000 -r 0 -s -i Mus_musculus.GRCm39.108_TSS.bed -g Mus_musculus.GRCm39.dna.primary_assembly.fa.fai > Mus_musculus.GRCm39.108_TSS2000.bed
# bedtools getfasta -name -fi acahirinus.scaffold.fa -bed acahirinus_TSS2000.bed > acahirinus_TSS2000.fa
# bedtools getfasta -name -fi Mus_musculus.GRCm39.dna.primary_assembly.fa -bed Mus_musculus.GRCm39.108_TSS2000.bed > Mus_musculus.GRCm39.108_TSS2000.fa

# Aco specific
tss_seqs = readDNAStringSet("/scratch/fmorandi/external/references/AcoCah2/acahirinus_TSS500.fa")
tss_seq_genes = str_extract(names(tss_seqs), "([^:]+)::", group=1)
any(is.na(tss_seq_genes))
any(tss_seq_genes == "")
tss_seqs_up = tss_seqs[tss_seq_genes %in% unlist(aco_de[aco_de$logFC_DTvsC_ACO > 0, "GeneSymbol"])]
tss_seqs_down = tss_seqs[tss_seq_genes %in% unlist(aco_de[aco_de$logFC_DTvsC_ACO < 0, "GeneSymbol"])]
tss_seqs_bg = tss_seqs[tss_seq_genes %in% comp_wide$GeneSymbol & !tss_seq_genes %in% unlist(aco_de[, "GeneSymbol"])]
writeXStringSet(tss_seqs_up, filepath = paste0(paths$seqs, "/aco_only_de_up.fa"))
writeXStringSet(tss_seqs_down, filepath = paste0(paths$seqs, "/aco_only_de_down.fa"))
writeXStringSet(tss_seqs_bg, filepath = paste0(paths$seqs, "/aco_only_de_bg.fa"))

# Mus specific
tss_seqs = readDNAStringSet("/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.108_TSS500.fa")
tss_seq_genes = str_extract(names(tss_seqs), "([^:]+)::", group=1)
any(is.na(tss_seq_genes))
any(tss_seq_genes == "")
tss_seqs_up = tss_seqs[tss_seq_genes %in% unlist(mus_de[mus_de$logFC_DTvsC_MUS > 0, "GeneSymbol"])]
tss_seqs_down = tss_seqs[tss_seq_genes %in% unlist(mus_de[mus_de$logFC_DTvsC_MUS < 0, "GeneSymbol"])]
tss_seqs_bg = tss_seqs[tss_seq_genes %in% comp_wide$GeneSymbol &!tss_seq_genes %in% unlist(mus_de[, "GeneSymbol"])]
tss_seqs_bg = tss_seqs_bg[!duplicated(names(tss_seqs_bg))]
writeXStringSet(tss_seqs_up, filepath = paste0(paths$seqs, "/mus_only_de_up.fa"))
writeXStringSet(tss_seqs_down, filepath = paste0(paths$seqs, "/mus_only_de_down.fa"))
writeXStringSet(tss_seqs_bg, filepath = paste0(paths$seqs, "/mus_only_de_bg.fa"))

# Lef1 promoter analysis
tss_seqs = readDNAStringSet("/scratch/fmorandi/external/references/AcoCah2/acahirinus_TSS2000.fa")
writeXStringSet(tss_seqs[grepl("Lef1", names(tss_seqs))], filepath = paste0(paths$seqs, "/aco_lef1_2000bp.fa"))
tss_seqs = readDNAStringSet("/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.108_TSS2000.fa")
writeXStringSet(tss_seqs[grepl("Lef1:", names(tss_seqs))], filepath = paste0(paths$seqs, "/mus_lef1_2000bp.fa"))

##### Lef1 promoter alignment #####

lef1_pro = readDNAStringSet(paste0(paths$seqs, "/lef1_2000bp.txt"))
# Aco lef1 is on -strand
# Mus lef1 is on +strand
aln = pairwiseAlignment(reverseComplement(lef1_pro[1]), lef1_pro[2], type = "global", gapExtension=0)
aln

aln_comp = data.frame(
  "Aco" = as.vector(str_split(as.character(alignedPattern(aln)), "", simplify = T)),
  "Mus" = as.vector(str_split(as.character(alignedSubject(aln)), "", simplify = T))) %>%
  mutate("Position" = row_number()-nrow(.)) %>%
  mutate(Matched = Aco == Mus)
window = 25
for (i in 1:nrow(aln_comp)){
  si = max(1, i-window)
  se = min(nrow(aln_comp), i+window)
  nmatch = sum(aln_comp[si:se, "Matched"])
  aln_comp[i, "pMatched"] = 100*nmatch / (2*window)
}
posMus = posAco = 0
for (i in nrow(aln_comp):1){
  if (aln_comp[i, "Aco"] != "-") posAco = posAco -1
  if (aln_comp[i, "Mus"] != "-") posMus = posMus -1
  aln_comp[i, "posAco"] = posAco
  aln_comp[i, "posMus"] = posMus
}
sum(aln_comp$Aco != "-")
sum(aln_comp$Mus != "-")

aln_comp = aln_comp %>%
  mutate(AlnAco = ifelse(Matched, "Match", "Mismatch"))%>%
  mutate(AlnMus = ifelse(Matched, "Match", "Mismatch"))%>%
  mutate(AlnAco = ifelse(Mus == "-", "Insertion", AlnAco))%>%
  mutate(AlnAco = ifelse(Aco == "-", "Deletion", AlnAco))%>%
  mutate(AlnMus = ifelse(Mus == "-", "Deletion", AlnMus))%>%
  mutate(AlnMus = ifelse(Aco == "-", "Insertion", AlnMus))
  

p1=ggplot(aln_comp, aes(Position, pMatched))+
  geom_line()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p2=aln_comp %>%
  pivot_longer(cols=c(Aco, Mus))%>%
  mutate(Gap = value == "-") %>%
  ggplot(., aes(Position, name, fill=Gap))+
  geom_tile()
p1/p2

##### Fimo #####

# by default list is ranked by pval
fimo = read.table(paste0(paths$seqs, "/fimo_lef1_2000bp.tsv"), header = T) %>%
  mutate(sequence_name = str_remove(sequence_name, "_Lef1"))
fimo[2, ]
substr(lef1_pro[1], 29, 40)
substr(reverseComplement(lef1_pro[1]), 2001-40+1, 2001-29+1)
# fimo[fimo$sequence_name == "Aco", "start"] = 2001 - fimo[fimo$sequence_name == "Aco", "start"]
# fimo[fimo$sequence_name == "Aco", "stop"] = 2001 - fimo[fimo$sequence_name == "Aco", "stop"]
tmp = 2001 - fimo[fimo$sequence_name == "Aco", "stop"]+1
fimo[fimo$sequence_name == "Aco", "stop"] = 2001 - fimo[fimo$sequence_name == "Aco", "start"]+1
fimo[fimo$sequence_name == "Aco", "start"] = tmp
fimo$start = fimo$start - 2001
fimo$stop = fimo$stop - 2001
length(unique(fimo$motif_id))

ggplot(fimo, aes(start))+
  geom_histogram(bins=50)+
  facet_wrap(~sequence_name)
ggplot(fimo, aes(score))+
  geom_histogram()

fimo_wide = fimo %>%
  group_by(sequence_name, motif_id) %>%
  arrange(p.value) %>%
  do(head(., n = 1)) %>%
  pivot_wider(names_from="sequence_name", values_from=c(start:matched_sequence))

ggplot(fimo_wide, aes(start_Aco, start_Mus, color=-log10(p.value_Aco)))+
  geom_point()
ggplot(fimo_wide, aes(score_Aco, score_Mus))+
  geom_point()

fimo_unique = fimo_wide %>%
  dplyr::filter(is.na(score_Aco) | is.na(score_Mus)) %>%
  mutate(motif_alt_id = str_to_title(motif_alt_id))
fimo_moderate = fimo_unique %>%
  merge(., comp_wide, by.x="motif_alt_id", by.y="GeneSymbol") %>%
  mutate(uniquelyBinds = ifelse(is.na(score_Aco), "M. mus pro. only", "A. cah pro. only"))%>%
  dplyr::filter(logFC_DTvsC_ACO * logFC_DTvsC_MUS < 0)%>%
  dplyr::rename("logFC_DTvsC_Aco" = "logFC_DTvsC_ACO", "logFC_DTvsC_Mus"="logFC_DTvsC_MUS") %>%
  dplyr::rename("PAdj_DTvsC_Aco" = "PAdj_DTvsC_ACO", "PAdj_DTvsC_Mus"="PAdj_DTvsC_MUS")
fimo_extreme = fimo_unique %>%
  merge(., aco_de, by.x="motif_alt_id", by.y="GeneSymbol") %>%
  mutate(uniquelyBinds = ifelse(is.na(score_Aco), "M. mus pro. only", "A. cah pro. only"))%>%
  dplyr::rename("logFC_DTvsC_Aco" = "logFC_DTvsC_ACO", "logFC_DTvsC_Mus"="logFC_DTvsC_MUS") %>%
  dplyr::rename("PAdj_DTvsC_Aco" = "PAdj_DTvsC_ACO", "PAdj_DTvsC_Mus"="PAdj_DTvsC_MUS")



fimo_moderate %>%
  mutate(motif_alt_id = fct_reorder(motif_alt_id, logFC_DTvsC_Aco - logFC_DTvsC_Mus))%>%
  pivot_longer(cols=-c(motif_alt_id, motif_id, SideCloud, uniquelyBinds), names_to=c(".value", "spe"), names_pattern = "(.*)_(Aco|Mus)")%>%
  mutate(sig = stars.pval(PAdj_DTvsC))%>%
  ggplot(., aes(spe, motif_alt_id, fill=logFC_DTvsC, label=sig))+
  geom_tile()+
  geom_text(color="white")+
  facet_grid(rows=vars(uniquelyBinds), space = "free", scale="free")+
  scale_fill_gradient2(low="blue", high="red")
fimo_extreme %>%
  mutate(motif_alt_id = fct_reorder(motif_alt_id, logFC_DTvsC_Aco - logFC_DTvsC_Mus))%>%
  pivot_longer(cols=-c(motif_alt_id, motif_id, Passes, uniquelyBinds), names_to=c(".value", "spe"), names_pattern = "(.*)_(Aco|Mus)")%>%
  mutate(sig = stars.pval(PAdj_DTvsC))%>%
  ggplot(., aes(spe, motif_alt_id, fill=logFC_DTvsC, label=sig))+
    geom_tile()+
    geom_text(color="white")+
    facet_grid(rows=vars(uniquelyBinds), space = "free", scale="free")+
  scale_fill_gradient2(low="blue", high="red")
# Foxq1 and Foxe1 have similar seqs as usual
# Shh (sonic hedgehog is also in aco_de and mixed up with Gli1 in wnt and more)
# Mycn/mycs is also very different between aco and mus! oncokb missed it due to ambiguous ortho!!!
cat(paste(fimo_unique$motif_alt_id, collapse="\n"))

tmp = fimo_extreme %>%
  dplyr::select(motif_alt_id, starts_with("start"), starts_with("stop"))%>%
  pivot_longer(cols=-c(motif_alt_id), names_to = c(".value", "seq"), names_sep = "_")%>%
  drop_na()%>%
  mutate(position = (start + stop) / 2) %>%
  mutate(nudge_dir = ifelse(seq == "Mus", +1, -1))%>%
  mutate(seq = ifelse(seq == "Aco", "A. cah", "M. mus"))
# Convert sequence position to aln position

as.integer(tmp[tmp$seq == "A. cah", ]$position)
aln_comp[match(as.integer(tmp[tmp$seq == "A. cah", ]$position), aln_comp$posAco), "posAco"]
aln_comp[match(as.integer(tmp[tmp$seq == "A. cah", ]$position), aln_comp$posAco), "Position"]
match(as.integer(tmp[tmp$seq == "A. cah", ]$position), aln_comp$posAco)
substr(lef1_pro[1], 2001-426, 2001-415) # Correct
substr(reverseComplement(lef1_pro[1]), 2001-1585, 2001-1574) #correct
a = paste(aln_comp$Aco, collapse="")
substr(a, 2638-2187-5, 2638-2187+5)
grepl("GGCAATTATTT", a)
b = str_replace_all(a, "-", "")
grepl("GGCAATTATTT", b)


tmp[tmp$seq == "A. cah", "position2"] = aln_comp[match(as.integer(tmp[tmp$seq == "A. cah", ]$position), aln_comp$posAco), "Position"]
tmp[tmp$seq == "M. mus", "position2"] = aln_comp[match(as.integer(tmp[tmp$seq == "M. mus", ]$position), aln_comp$posMus), "Position"]

p1=ggplot(aln_comp, aes(Position, pMatched))+
  geom_line()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(y="% matches\n(50 bp window)")
p2=aln_comp %>%
  pivot_longer(cols=c(Aco, Mus))%>%
  mutate(Gap = value == "-") %>%
  mutate(name = ifelse(name == "Aco", "A. cah", "M. mus"))%>%
  ggplot(.)+
  geom_tile(aes(Position, name, fill=Gap), height=0.8)+
  geom_text_repel(data=tmp, aes(position, seq, label=motif_alt_id), nudge_y = 1*tmp$nudge_dir)+
  labs(x = "Consensus position (from Lef1 TSS)", y="Promoter sequence")+
  scale_fill_manual(values=c("#66cc66", "#888888"))
p3 = fimo_extreme %>%
  mutate(motif_alt_id = fct_reorder(motif_alt_id, logFC_DTvsC_Aco - logFC_DTvsC_Mus))%>%
  pivot_longer(cols=-c(motif_alt_id, motif_id, Passes, uniquelyBinds), names_to=c(".value", "spe"), names_pattern = "(.*)_(Aco|Mus)")%>%
  mutate(sig = stars.pval(PAdj_DTvsC))%>%
  mutate(spe = ifelse(spe == "Aco", "A. cah\nDTvsC", "M. mus\nDTvsC"))%>%
  ggplot(., aes(spe, motif_alt_id, fill=logFC_DTvsC, label=sig))+
  geom_tile()+
  geom_text(color="white")+
  facet_grid(rows=vars(uniquelyBinds), space = "free", scale="free")+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.title=element_blank())+
  labs(fill="logFC")

p2 = aln_comp %>%
  pivot_longer(cols=c(Aco, Mus, AlnAco, AlnMus))%>%
  mutate(spe = ifelse(grepl("Mus", name), "Mus", "Aco"))%>%
  mutate(type = ifelse(grepl("Aln", name), "Aln", "Seq"))%>%
  dplyr::select(-name)%>%
  pivot_wider(names_from = type)%>%
  mutate(spe = ifelse(spe == "Aco", "A. cah", "M. mus"))%>%
  mutate(Aln = factor(Aln, levels=c("Match", "Mismatch", "Insertion", "Deletion")))%>%
  ggplot(.)+
  geom_tile(aes(Position, spe, fill=Aln), height=0.8)+
  geom_text_repel(data=tmp, aes(position2, seq, label=motif_alt_id), nudge_y = 1*tmp$nudge_dir)+
  labs(x = "Alignment position (from Lef1 TSS)", y="Promoter sequence", fill="Alignment")+
  scale_fill_manual(values=c("#55bb55", "#dd5555", "#dddd66", "#666666")) #"#66cc66", "#cc6666", "#888888", "#444444"

design = "AAC
          BBC
          BBC
          BBC"
p=p1+p2+p3+plot_layout(design = design, guides = "collect")
p
ggsave(plot=p, paste0(paths$out, "/paper_figs/promoter_analysis.pdf"), width = 1*w, height = 0.5*h, units="mm")

#### PAPER FIGURES ####

##### Pca #####

pca = prcomp(norm$intersection, center=T, scale.=T)
pc_importance = summary(pca)$importance[2,1:2] * 100
pc_titles = sprintf(c("PC1 (%.1f%%)", "PC2 (%.1f%%)"), pc_importance)
pca = merge(meta, pca$x[,c("PC1", "PC2")], by=0)
ggplot(pca, aes(x=PC1, y=PC2, color=Treatment, shape=Species))+
  geom_point(size=1.5)+
  labs(x=pc_titles[1], y=pc_titles[2])
ggsave(paste0(paths$out, "/paper_figs/pca.pdf"), width = 0.5*w, height = 0.3*h, units="mm")

##### Volcanos #####

# --- Treatment ---
p1 = my_volcano(des_treatment$aco, col="sig_col", v1="logFC_DTvsC", v2="PValue_DTvsC", "A. cahirinus")
p2 = my_volcano(des_treatment$mus, col="sig_col", v1="logFC_DTvsC", v2="PValue_DTvsC", "M. musculus")

p= p1 + p2 +
  plot_layout(guides = "collect", ncol=2) &
  scale_color_manual(values=c("#999999", "#5555cc", "#cc5555")) &
  guides(color="none")&
  labs(x="logFC", y = "-log10(pvalue)")
ggsave(paste0(paths$out, "/paper_figs/volcanos_treament.pdf"), plot=p, width=w, height=0.4*h, units="mm")

##### Figure 4c (aco specific gsea) #####

# Heatmaps to Fathima can select terms easier

aco_nes %>%
  dplyr::filter(NES_musDTvsC > 0)%>%
  dplyr::filter(padj_musDTvsC < 0.05)

aco_nes %>%
  mutate(dNES = NES_acoDTvsC - NES_musDTvsC) %>%
  mutate(dNES_sign = sign(dNES)) %>%
  group_by(dNES_sign) %>%
  # slice_max(abs(dNES), n=30) %>%
  # mutate(Description = fct_reorder(Description, NES_acoDTvsC)) %>%
  mutate(Description = ifelse(nchar(Description) > 60, paste0(substr(Description, 1, 60), "..."), Description))%>%
  mutate(Description = fct_reorder(Description, dNES)) %>%
  pivot_longer(cols=-c(ID, Description, dNES, dNES_sign), names_pattern = "(.*)_(aco|mus)(.*)", names_to = c("Var", "Species", "Treatment")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(Species, levels=c("aco", "mus"))) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "aco", "M. mus" = "mus")) %>%
  mutate(Treatment = "DTvsC") %>%
  # arrange(Species) %>%
  # mutate(Condition = paste(Species, Treatment, sep="_")) %>%
  # mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(padj)) %>%
  ggplot(., aes(Species, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  facet_grid(cols=vars(Treatment))+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/paper_figs/aco_specificDT_gsea_heatmap.pdf"), width = 0.8*w, height = 7*h, units="mm", limitsize =F)

aco_nes %>%
  mutate(dNES = NES_acoDTvsC - NES_musDTvsC) %>%
  mutate(dNES_sign = sign(dNES)) %>%
  group_by(dNES_sign) %>%
  slice_max(abs(dNES), n=30) %>%
  mutate(Description = ifelse(nchar(Description) > 80, paste0(substr(Description, 1, 80), "..."), Description))%>%
  mutate(Description = fct_reorder(Description, dNES)) %>%
  pivot_longer(cols=-c(ID, Description, dNES, dNES_sign), names_pattern = "(.*)_(aco|mus)(.*)", names_to = c("Var", "Species", "Treatment")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(Species, levels=c("aco", "mus"))) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "aco", "M. mus" = "mus")) %>%
  mutate(Treatment = "DTvsC") %>%
  # arrange(Species) %>%
  # mutate(Condition = paste(Species, Treatment, sep="_")) %>%
  # mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(padj)) %>%
  ggplot(., aes(Species, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  facet_grid(cols=vars(Treatment))+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/paper_figs/aco_specificDT_gsea_heatmap_topN.pdf"), width = 0.75*w, height = 1*h, units="mm")

selected_terms = fromJSON(file=paste0(paths$out, "/paper_figs/selected_terms.txt"))
selected_terms = stack(selected_terms)
colnames(selected_terms) = c("Description", "Topic")
setdiff(selected_terms$Description, aco_nes$Description)
intersect(selected_terms$Description, aco_nes$Description)
length(intersect(selected_terms$Description, aco_nes$Description))
length(aco_nes$Description)

aco_nes %>%
  merge(., selected_terms, by="Description")%>%
  mutate(Description = ifelse(nchar(Description) > 80, paste0(substr(Description, 1, 80), "..."), Description))%>%
  mutate(Description = fct_reorder(Description, NES_acoDTvsC)) %>%
  pivot_longer(cols=-c(ID, Description, Topic), names_pattern = "(.*)_(.*)DTvsC", names_to = c("Var", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = factor(Species, levels=c("aco", "mus"))) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "aco", "M. mus" = "mus")) %>%
  arrange(Species) %>%
  # mutate(Condition = paste(Species, Dose, "vsC", sep="_")) %>%
  # mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
  mutate(Sig = stars.pval(padj)) %>%
  mutate(Comp = "DT vs C") %>%
  ggplot(., aes(Species, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  facet_grid(cols=vars(Comp), rows=vars(Topic), scale="free", space="free")+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/paper_figs/aco_specificDT_gsea_heatmap_byTopic.pdf"), width = 0.75*w, height = 1*h, units="mm")

##### Figure 5 (wnt) #####

des_combined_fcNp %>%
  mutate(Species = fct_recode(Species, "A. cah" = "aco", "M. mus" = "mus"))%>%
  dplyr::filter(Symbol %in% c("Tcf7", "Tcf7l1", "Tcf7l2", "Lef1")) %>%
  mutate(Symbol = factor(Symbol, c("Tcf7l1", "Tcf7l2", "Tcf7", "Lef1"))) %>%
  ggplot(.,aes(Species, Symbol, fill=logFC, label=Sig))+
  geom_tile()+
  geom_text()+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  theme(axis.title=element_blank())
ggsave(paste0(paths$out, "/paper_figs/tcf_lef.pdf"), width = 0.35*w, height = 0.2*h, units="mm")

bcat_dependent %>%
  mutate(Species = fct_recode(Species, "A. cah" = "aco", "M. mus" = "mus")) %>%
  ggplot(., aes(logFC.WT, logFC))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Species)+
  stat_cor()+
  ylim(-4, 7)+
  labs(x="logFC CHIRvsC", y="logFC DTvsC")
ggsave(paste0(paths$out, "/paper_figs/wnt_targets.pdf"), width = 0.6*w, height = 0.3*h, units="mm")

ggplot(res_custom_gsets$wnt_target$gsea, aes(Condition, NES, fill=Condition, label=Sig))+
  geom_bar(stat="identity")+
  geom_text(aes(y=NES/2), color="white")+
  facet_wrap(~Description)+
  guides(fill="none")+
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$out, "/paper_figs/wnt_targets_gsea1.pdf"), width = 0.5*w, height = 0.3*h, units="mm")

all.equal(rownames(norm$union), rownames(meta))
meta %>%
  mutate(Lef1 = norm$union[, "Lef1"]) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "Aco", "M. mus" = "Mus")) %>%
  ggplot(., aes(Species, Lef1, fill=Treatment))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())+
  labs(y = "Lef1 [log(CPM)]")
ggsave(paste0(paths$out, "/paper_figs/lef1_expression.pdf"), width = 0.6*w, height = 0.3*h, units="mm")

##### dLogFC > GSEA #####

# Dotplot column doesnt look good because almost all pvalues are extreme
# data.frame(gse_diff_response)%>%
#   dplyr::filter(p.adjust < 0.05) %>% 
#   mutate(signNES = ifelse(NES > 0, "Higher in A. cahirinus", "Higher in M.musculus")) %>%
#   group_by(signNES) %>%
#   top_n(n=20, -pvalue) %>%
#   arrange(sign(NES) * pvalue) %>%
#   mutate(Description = factor(Description, levels=unique(Description))) %>%
#   ggplot(., aes(1, Description, fill=NES, size=-log10(p.adjust)))+
#   geom_point(pch=21)

# Noticed that the "immune response to tumor cell is enriched, lets check genes:
tmp = data.frame(gse_diff_response) %>%
  dplyr::filter(grepl("immune response to tumor cell", Description)) %>%
  head(n=1) %>%
  pull(core_enrichment)
tmp = mapIds(org.Mm.eg.db, keys=str_split(tmp, "/", simplify = T), column = "SYMBOL", keytype = "ENTREZID")
intersect(aco_de$GeneSymbol, tmp)
# Basically some NK receptors

tmp = data.frame(gse_diff_response)%>%
    dplyr::filter(p.adjust < 0.05)%>% 
    mutate(signNES = ifelse(NES > 0, "Higher in A.cahirinus", "Higher in M.musculus")) 
tmp = tmp %>%
    group_by(signNES) %>%
    arrange(pvalue, -abs(NES)) %>%
    do(head(., n = 40)) # Top most significant, breaking ties by NES
selected = c(
  "intermediate filament organization",
  "epidermis development",
  "keratinocyte differentiation",
  "hair cycle",
  "BMP signaling pathway",
  "Wnt signaling pathway",
  
  "lymphocyte mediated immunity",
  "response to interferon-gamma",
  "response to interferon-beta",
  "cell killing",
  "positive regulation of cytokine production",
  "activation of immune response"
)
tmp %>%
  ungroup()%>%
  dplyr::filter(Description %in% selected) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  ggplot(., aes(Description, NES, fill = signNES))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme(axis.title.y=element_blank())
ggsave(paste0(paths$out, "/paper_figs/gsea_diff_response.pdf"), width = 1*w, height = 0.4*h, units="mm")
tmp %>%
  ungroup()%>%
  dplyr::filter(Description %in% selected) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  mutate(Comp = "DT vs C") %>%
  mutate(sig = stars.pval(p.adjust))%>%
  ggplot(., aes(1, Description, fill=NES, label = sig))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  facet_grid(rows=vars(signNES), cols = vars(Comp), space="free", scale="free")+
  theme(axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(paste0(paths$out, "/paper_figs/gsea_diff_response2.pdf"), width = 0.6*w, height = 0.4*h, units="mm")

# facet_titles = c(
#   activated = "Higher in Acomys",
#   suppressed = "Higher in Mouse"
# )
# gse_diff_response %>%
#   dplyr::filter(p.adjust < 0.05) %>%
#   dotplot(., split=".sign", showCategory = 20)+
#   facet_grid(.~.sign, labeller=as_labeller(facet_titles))+
#   theme(axis.text.y = element_text(size = 10))+
#   scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
#   ggtitle("Acomys vs Mouse")
# 
# pdf(paste0(paths$out, "/paper_figs/gsea_diff_response.pdf"), width=w_in*2, height=h_in)
# p1
# dev.off()

##### Senescence #####


# --- SenMayo ---
tmp1 = res_custom_gsets$sen_mayo$gsea %>%
  mutate(Source = "SenMayo") %>%
  mutate(Description = str_remove(Description, "SenMayo_")) %>%
  mutate(Description = str_replace_all(Description, "_", " ")) %>%
  mutate(Description = fct_rev(Description)) %>%
  mutate(Species = str_remove(Context, "_.+"))%>%
  mutate(Species = factor(Species, levels=c("aco", "mus"))) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "aco", "M. mus" = "mus")) %>%
  mutate(Comparison = "DTvsC")
# --- Fathima custom ---
tmp2 = res_custom_gsets$fathima_custom$gsea %>%
  merge(., distinct(gmts$fathima_custom[, c("term", "desc")]), by.x="ID", by.y="term") %>%
  mutate(Source = "Senescence Signatures") %>%
  mutate(Description = str_remove(desc, "SenMayo_")) %>%
  # mutate(Description = str_replace_all(Description, "_", " ")) %>%
  mutate(Species = str_remove(Context, "_.+"))%>%
  mutate(Species = factor(Species, levels=c("aco", "mus"))) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "aco", "M. mus" = "mus"))%>%
  mutate(Comparison = "DTvsC") %>%
  mutate(Description = factor(Description, levels=c(
    "NFkB regulated SASP",
    "Fibroblast SASP factors",
    "Fibroblast secretory SASP",
    "Core factors",
    "Core I",
    "Core II",
    "P53 regulated SASP",
    "p53 regulators",
    "Top10 secreted proteins exclusive to IR secretomes",
    "Top10 secreted proteins exclusive to RAS secretomes",
    "Modifiers and downstream"))) %>%
  mutate(Description = fct_rev(Description)) %>%
  dplyr::select(colnames(tmp1))
tmp = rbind(tmp1, tmp2)

tmp %>%
  # dplyr::filter(Description != "Protein modifying enzyme") %>%
  ggplot(., aes(Species, Description, fill=NES, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  facet_grid(cols=vars(Comparison), rows=vars(Source), scales = "free", space="free")+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        axis.title = element_blank())
ggsave(paste0(paths$out, "/paper_figs/sen_sigsDT.pdf"), width = 0.75*w, height = 0.5*h, units="mm")

colnames(norm$union)[grepl("cdkn", colnames(norm$union), ignore.case = T)]
colnames(norm$union)[grepl("p16", colnames(norm$union), ignore.case = T)]

all.equal(rownames(meta), rownames(norm$union))
meta %>%
  mutate(p21 = norm$union[rownames(meta), "Cdkn1a"]) %>%
  dplyr::filter(Treatment != "ConP")%>%
  dplyr::filter(Species != "HSA") %>%
  mutate(Species = fct_recode(Species, "A. cah" = "Aco", "M. mus" = "Mus")) %>%
  ggplot(., aes(Species, p21, fill=Treatment))+
  geom_boxplot()+
  labs(y = "p21 expression [log(CPM)]")+
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$out, "/paper_figs/p21.pdf"), width = 0.5*w, height=0.3*h, units="mm")

meta %>%
  mutate(p16 = norm$union[rownames(meta), "Cdkn2a"]) %>%
  dplyr::filter(Treatment != "ConP")%>%
  dplyr::filter(Species != "HSA") %>%
  mutate(Species = fct_recode(Species, "A. cah" = "Aco", "M. mus" = "Mus")) %>%
  ggplot(., aes(Species, p16, fill=Treatment))+
  geom_boxplot()+
  labs(y = "p16 expression [log(CPM)]")+
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$out, "/paper_figs/p16.pdf"), width = 0.5*w, height=0.3*h, units="mm")

##### Onco aco specific #####

# --- OncoKB: tsupp ---
aco_de %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isTS2) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "ACO", "M. mus" = "MUS")) %>%
  arrange(Species, Comparison) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Species, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        axis.title = element_blank())
  # ggtitle("Tumor suppressors (OncoKB)\nwith unique Acomys response")
ggsave(paste0(paths$out, "/paper_figs/aco_specific_OncoKB_tsupp.pdf"), width = 0.35*w, height = 0.25*h, units="mm")

# --- OncoKB: oncogene ---
aco_de %>%
  ungroup() %>%
  dplyr::select(-Passes) %>%
  dplyr::filter(isOncogene) %>%
  mutate(GeneSymbol = fct_reorder(GeneSymbol, logFC_DTvsC_ACO)) %>%
  pivot_longer(cols=c(starts_with("logFC"), starts_with("PAdj")), names_pattern = "(.*)_(.*)_(.*)", names_to = c("Var", "Comparison", "Species")) %>%
  pivot_wider(names_from = Var, values_from = value) %>%
  mutate(Species = fct_recode(Species, "A. cah" = "ACO", "M. mus" = "MUS")) %>%
  arrange(Species, Comparison) %>%
  mutate(Sig = stars.pval(PAdj)) %>%
  ggplot(., aes(Species, GeneSymbol, fill=logFC, label=Sig))+
  geom_tile()+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        axis.title = element_blank())
  # ggtitle("Oncogenes (OncoKB)\nwith unique Acomys response")
ggsave(paste0(paths$out, "/paper_figs/aco_specific_OncoKB_oncogene.pdf"), width = 0.35*w, height = 0.6*h, units="mm")

#### dLogFC>GSEA vs dGSEA ####

# gse_treatment contains GSEA results performed on logFCs
# gse_combined has them merged side by side
# gse_diff_response and gse_baseline have GSEA performed on dLogFCs
# gse_diff has them merged side by side

# View(gse_combined$T1)
# View(gse_diff)
# tmp1 = gse_combined$T1 %>%
#   dplyr::select(ID, Description, starts_with("NES")) %>%
#   pivot_longer(cols=-c(ID, Description), names_pattern = "NES_(.*)(24h|12d)_(.*)", 
#                names_to = c("Species", "Timepoint", "Dose"), values_to = "NES") %>%
#   pivot_wider(names_from = Species, values_from = NES) %>%
#   mutate(dNES_aco_vs_mus = aco - mus, dNES_aco_vs_musw = aco - musw, dNES_aco_vs_hsa = aco - hsa) %>%
#   dplyr::select(ID:Dose, starts_with("dNES"))
# tmp2 = gse_combined$T2 %>%
#   dplyr::select(ID, Description, starts_with("NES")) %>%
#   pivot_longer(cols=-c(ID, Description), names_pattern = "NES_(.*)(24h|12d)_(.*)", 
#                names_to = c("Species", "Timepoint", "Dose"), values_to = "NES")%>%
#   pivot_wider(names_from = Species, values_from = NES) %>%
#   mutate(dNES_aco_vs_mus = aco - mus, dNES_aco_vs_musw = aco - musw)%>%
#   dplyr::select(ID:Dose, starts_with("dNES"))
# tmp1 = plyr::rbind.fill(tmp1, tmp2)
# tmp2 = gse_diff %>%
#   dplyr::select(ID, Description, starts_with("NES")) %>% 
#   dplyr::select(ID, Description, contains("10"), contains("20")) %>% # Remove baseline NES
#   pivot_longer(cols=-c(ID, Description), names_pattern = "NES_(.*)_vs_(.*)(10|20)_(24h|12d)", 
#                names_to = c("Species1", "Species2", "Dose", "Timepoint"), values_to = "NES") %>%
#   mutate(Dose = paste0(Dose, "gy")) %>%
#   pivot_wider(names_from = c(Species1, Species2), names_glue = "dLogFC_{Species1}_vs_{Species2}", values_from = "NES")
# tmp3 = merge(tmp1, tmp2)
# cor(tmp3$dNES_aco_vs_mus, tmp3$dLogFC_aco_vs_mus, use = "pairwise.complete.obs")
# cor(tmp3$dNES_aco_vs_musw, tmp3$dLogFC_aco_vs_musw, use = "pairwise.complete.obs")
# cor(tmp3$dNES_aco_vs_hsa, tmp3$dLogFC_aco_vs_hsa, use = "pairwise.complete.obs")

#### SAVE OBJECTS ####

saveRDS(des_treatment, file = paste0(paths$objects, "/des_treatment.rds"))
saveRDS(des_baseline, file = paste0(paths$objects, "/des_baseline.rds"))
saveRDS(des_combined, file = paste0(paths$objects, "/des_combined.rds"))

#### SAVE TABLES ####

# gse_treatment content hand't been converted to table yet
for (comp in names(gse_treatment)) {
  gse_treatment[[comp]] = data.frame(gse_treatment[[comp]])
}

##### Response #####

# Acomys conversion
gse_treatment$acoDTvsC = convert_ids_in_string(gse_treatment$acoDTvsC, ginfo$aco$Entrez, ginfo$aco$GeneSymbol)
# Mouse conversion
gse_treatment$musDTvsC = convert_ids_in_string(gse_treatment$musDTvsC, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)

# Aco
write.table(subset(gse_treatment$acoDTvsC, p.adjust < 0.01), paste0(paths$tables, "/gsea_res_acoDTvsC.tsv"), 
            row.names = F, sep="\t", quote=F)
# Mus
write.table(subset(gse_treatment$musDTvsC, p.adjust < 0.01), paste0(paths$tables, "/gsea_res_musDTvsC.tsv"), 
            row.names = F, sep="\t", quote=F)

# Combined
write.table(gse_combined, paste0(paths$tables, "/gsea_res_combined.tsv"), 
            row.names = F, sep="\t", quote=F)

##### Difference in baseline #####

# All genes use the mouse ids
gse_baseline2 = data.frame(gse_baseline)
gse_baseline2 = convert_ids_in_string(gse_baseline2, ginfo$mus$Entrez, ginfo$mus$GeneSymbol)
write.table(subset(gse_baseline2, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_baseline.tsv"), 
            row.names = F, sep="\t", quote=F)

##### Difference in response #####

# All genes use the mouse ids
gse_diff_response2 = convert_ids_in_string(as.data.frame(gse_diff_response), ginfo$mus$Entrez, ginfo$mus$GeneSymbol)

write.table(subset(gse_diff_response2, p.adjust < 0.01), 
            paste0(paths$tables, "/gsea_diff_response_aco_vs_mus.tsv"), 
            row.names = F, sep="\t", quote=F)

##### Count tables and DEs #####

write.csv(counts$aco, paste0(paths$tables, "/counts_aco.csv"))
write.csv(counts$mus, paste0(paths$tables, "/counts_mus.csv"))
write.csv(norm$all, paste0(paths$tables, "/cpm_combined.csv"))

write.table(des_treatment$aco[,-c(2:6)], paste0(paths$tables, "/des_aco.tsv"), row.names = F, sep="\t", quote=F)
write.table(des_treatment$mus[,-c(2:6)], paste0(paths$tables, "/des_mus.tsv"), row.names = F, sep="\t", quote=F)
write.table(des_baseline[,-c(2:6)], paste0(paths$tables, "/des_baseline.tsv"), row.names = F, sep="\t", quote=F)