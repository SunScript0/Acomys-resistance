library(XML)
library(ggplot2)
library(ggforce)

read_kegg_map = function(path) {
  doc = xmlParse(path)
  doc = xmlToList(doc, addAttributes = TRUE, simplify = FALSE)
  entries = doc[names(doc) == "entry"]
  relations = doc[names(doc) == "relation"]
  # --- Node df ---
  mapdf = list()
  for (i in 1:length(entries)) {
    mapdf[[i]] = list(
      "Id" = entries[[i]]$.attrs["id"],
      "type" = entries[[i]]$.attrs["type"],
      "keggId" = entries[[i]]$.attrs["name"],
      "aliases" = entries[[i]]$graphics["name"],
      "fgcolor" = entries[[i]]$graphics["fgcolor"],
      "bgcolor" = entries[[i]]$graphics["bgcolor"],
      "shape" = entries[[i]]$graphics["type"],
      "x" = entries[[i]]$graphics["x"],
      "y" = entries[[i]]$graphics["y"],
      "w" = entries[[i]]$graphics["width"],
      "h" = entries[[i]]$graphics["height"]
    )
  }
  mapdf = Reduce(function(x,y) rbind(x, y), mapdf)
  mapdf = apply(mapdf, 2, trimws)
  mapdf = data.frame(mapdf)
  mapdf[, c("x", "y", "w", "h")] = sapply(mapdf[, c("x", "y", "w", "h")], as.numeric)
  mapdf$GeneName = str_extract(mapdf$aliases, "[[:alnum:]]+")
  # --- Node df ---
  reldf = list()
  for (i in 1:length(relations)) {
    reldf[[i]] = list(
      "Id1" = relations[[i]]$.attrs["entry1"],
      "Id2" = relations[[i]]$.attrs["entry2"],
      "type" = relations[[i]]$.attrs["type"],
      "subtype" = relations[[i]]$subtype["name"],
      "arrow" = relations[[i]]$subtype["value"]
    )
  }
  reldf = Reduce(function(x,y) rbind(x, y), reldf)
  reldf = apply(reldf, 2, trimws)
  reldf = data.frame(reldf)
  # Calculate better arrow positions
  rownames(mapdf) = mapdf$Id
  reldf$x1 = mapdf[reldf$Id1, "x"]
  reldf$y1 = mapdf[reldf$Id1, "y"]
  reldf$w1 = mapdf[reldf$Id1, "w"]
  reldf$h1 = mapdf[reldf$Id1, "h"]
  reldf$x2 = mapdf[reldf$Id2, "x"]
  reldf$y2 = mapdf[reldf$Id2, "y"]
  reldf$w2 = mapdf[reldf$Id2, "w"]
  reldf$h2 = mapdf[reldf$Id2, "h"]
  reldf$Gene1 = mapdf[reldf$Id1, "GeneName"]
  reldf$Gene2 = mapdf[reldf$Id2, "GeneName"]
  reldf = reldf %>%
    mutate(dx = x2-x1, dy = y2-y1) %>%
    mutate(xstart = ifelse(abs(dy) > sqrt(3) * abs(dx), x1, x1 + sign(x2-x1) * w1/2)) %>%
    mutate(ystart = ifelse(abs(dy) > sqrt(3) * abs(dx), y1 + sign(y2-y1) * h1/2, y1))
  for (i in 1:nrow(reldf)) {
    Xe = lineIntersectionOnRect2(
      reldf[i, "xstart"],
      reldf[i, "ystart"],
      reldf[i, "x2"],
      reldf[i, "y2"],
      reldf[i, "w2"],
      reldf[i, "h2"],
      margin=0.05)
    reldf[i, "xend"] = Xe[1]
    reldf[i, "yend"] = Xe[2]
  }
  return(list(
    "nodes" = mapdf,
    "edges" = reldf
  ))
}


lineIntersectionOnRect2 = function(xs, ys, xe, ye, w, h, margin=0) {
  # Horizontal lines
  t1 = (ye + (h/2) - ys) / (ye - ys)
  t2 = (ye - (h/2) - ys) / (ye - ys)
  # Vertical lines
  t3 = (xe + (w/2) - xs) / (xe - xs)
  t4 = (xe - (w/2) - xs) / (xe - xs)
  # Check if horz or vert intersection
  dx = abs(xe-xs)
  dy = abs(ye-ys)
  if (dy/dx > h/w) { # intersects top
    t = min(max(0, t1), max(0, t2))
  } else {
    t = min(max(0, t3), max(0, t4))
  }
  x = xs + (xe-xs)*t*(1-margin)
  y = ys + (ye-ys)*t*(1-margin)
  return(c(x,y))
}


plot_kegg_map = function(lst, gene_color=NULL) {
  nodes_genes = subset(lst$nodes, type == "gene")
  nodes_comps = subset(lst$nodes, type == "compound")
  nodes_map = subset(lst$nodes, type == "map")
  nodes_grp = subset(lst$nodes, type == "group")
  # table(lst$nodes$type, lst$nodes$shape)
  # table(lst$edges$subtype, lst$edges$arrow)
  # distinct(lst$edges[, c("Id1", "Id2")])
  edges_act = subset(lst$edges, arrow == "-->" | arrow == "+p")
  edges_repr = subset(lst$edges, arrow == "--|")
  edges_line = subset(lst$edges, arrow == "---")
  edges_indir = subset(lst$edges, arrow == "..>")
  p = ggplot()+
    geom_segment(data=edges_act, aes(x = xstart, y = ystart, xend = xend, yend = yend), arrow=arrow(type="closed", length = unit(.2,"cm")))+
    geom_segment(data=edges_repr, aes(x = xstart, y = ystart, xend = xend, yend = yend), arrow=arrow(type="closed", length = unit(.2,"cm"), angle = 90))+
    geom_segment(data=edges_line, aes(x = xstart, y = ystart, xend = xend, yend = yend))+
    geom_segment(data=edges_indir, aes(x = xstart, y = ystart, xend = xend, yend = yend), linetype="dashed")+
    geom_rect(data=nodes_grp, aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2), fill=NA, color="black")+
    geom_circle(data=nodes_comps, aes(x0=x, y0=y, r=w/2), fill="yellow")+
    geom_text(data=nodes_comps, aes(x, y, label = GeneName), size=3, vjust=1.5)+
    # geom_rect(data=nodes_map, aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2), fill="#bbccdd")+
    geom_label(data=nodes_map, aes(x, y, label = aliases), size=3, fill="#bbccdd")+
    scale_y_reverse()
  if (is.null(gene_color)) {
    p = p + 
      geom_rect(data=nodes_genes, aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2), fill="white")+
      geom_text(data=nodes_genes, aes(x, y, label = GeneName), size=3)
  } else {
    maxabs = max(abs(range(lst$nodes[[gene_color]], na.rm = T)))
    p = p + 
      geom_rect(data=nodes_genes, aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2, fill=.data[[gene_color]]))+
      geom_text(data=nodes_genes, aes(x, y, label = GeneName), size=3)+
      scale_fill_gradient2(low="blue", high="red", limits=c(-maxabs,maxabs))
  }
  
  return(p)
}

entrez_to_symbol = function(entrez) {
  # Use mapIds, but only for valid (non-NA, non-empty) gene IDs
  valid_indices = !is.na(entrez) & entrez != ""
  out = rep(NA, length(entrez))  # Pre-fill with NA
  out[valid_indices] = mapIds(
    org.Mm.eg.db,
    keys = entrez[valid_indices],
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first")
  return(out)
}

make_scores = function(lst, logFCs, mainCol, center=T, scale=F) {
  lst$nodes$EntrezId1 = str_extract(lst$nodes$keggId, "^mmu:(\\d+)", group=1)
  lst$nodes$Symbol1 = entrez_to_symbol(lst$nodes$EntrezId1)
  valid = !is.na(lst$nodes$Symbol1) & lst$nodes$Symbol1 %in% rownames(logFCs)
  which_genes = lst$nodes$Symbol1[valid]
  tmp = logFCs[which_genes, ] %>%
    dplyr::select(starts_with("logFC")) %>%
    dplyr::select(contains("vsC"))
  scr = data.frame(t(apply(tmp, 1, scale, center = center, scale = scale)))
  colnames(scr) = colnames(tmp)
  lst$nodes[valid, "score"] = scr[, mainCol]
  return(lst)
}


