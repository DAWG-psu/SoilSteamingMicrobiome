#Figure 4&5. Heat tree matrix (Bacteria, Fungi) and Table S5 & S6
table(tax_table(ps.bac)[, "Kingdom"], exclude = NULL) #For Domain/Kingdom
##Filter <NA> Kingdom - meaningless
ps.bac.filt <- subset_taxa(ps.bac, Kingdom == "Bacteria")

library(devtools)
library(microbiome)
library(ggplot2)
library(metagMisc)
library(metacoder)
sample_names(ps.fun) <- sample_names(ps.bac)
ps.fun <- phyloseq(otu_table(ps.fun), tax_table(ps.fun), sample_data(ps.bac))

metacoder_creating_tax_tree <- function(phyloseq_object, tax_level, Kingdom) {
  require(phyloseq)
  require(metagMisc)
  if (Kingdom == "bacteria"){
  ps.bac.filt <- subset_taxa(ps.bac, Kingdom == "Bacteria")}
  if (Kingdom == "fungi"){
  ps.bac.filt <- subset_taxa(ps.fun, Kingdom == "k__Fungi")
  tax_table(ps.bac.filt) = gsub(pattern = "k__", replace= "", x = tax_table(ps.bac.filt))
  tax_table(ps.bac.filt) = gsub(pattern = "g__", replace= "", x = tax_table(ps.bac.filt))
  tax_table(ps.bac.filt) = gsub(pattern = "c__", replace= "", x = tax_table(ps.bac.filt))
  tax_table(ps.bac.filt) = gsub(pattern = "f__", replace= "", x = tax_table(ps.bac.filt))
  tax_table(ps.bac.filt) = gsub(pattern = "g__", replace= "", x = tax_table(ps.bac.filt))
  tax_table(ps.bac.filt) = gsub(pattern = "p__", replace= "", x = tax_table(ps.bac.filt))
  tax_table(ps.bac.filt) = gsub(pattern = "o__", replace= "", x = tax_table(ps.bac.filt))
  }
  ps.bac.filt -> ps
  
  print("Loading Phyloseq object...")
  ps_filter = phyloseq_filter_prevalence(ps, threshold_condition = "OR",
                                         abund.trh = 10,abund.type = "total",
                                         prev.trh = 0.2)
  phyloseq_prevalence_plot(ps_filter,taxcolor = "Phylum", facet=T)
  ps_tax_level = tax_glom(ps_filter, taxrank = tax_level)
  return(ps_tax_level)
  }

ps <- metacoder_creating_tax_tree(ps.bac,tax_level = "Genus", Kingdom = "Bacteria")
ps.bac_metacoder <- parse_phyloseq(ps)
ps.bac_metacoder$data$otu_props <- calc_obs_props(ps.bac_metacoder, data = "otu_table", cols = ps.bac_metacoder$data$sample_data$sample_id)   ####Convert to proportions
ps.bac_metacoder$data$tax_table <- calc_taxon_abund(ps.bac_metacoder, data = "otu_table", cols = ps.bac_metacoder$data$sample_data$sample_id)   ####Convert to proportions
ps.bac_metacoder$data$sample_data <- ps.bac_metacoder$data$sample_data[-18,] #Remove negative control 
ps.bac_metacoder$data$diff_table <- calc_diff_abund_deseq2(ps.bac_metacoder, data = "tax_table", 
                                                           cols = ps.bac_metacoder$data$sample_data$sample_id, 
                                                           groups = ps.bac_metacoder$data$sample_data$Depth)
x <- ps.bac_metacoder
# Remove taxa with only small differences 
per_taxon_fold_changes <- obs(x, data = 'diff_table', value = 'log2FoldChange')
per_taxon_p_values <- obs(x, data = 'diff_table', value = 'pvalue')
per_taxon_max_change <- unlist(lapply(seq_along(per_taxon_fold_changes), function(i) {
  pvals <- per_taxon_p_values[[i]]
  pvals <- pvals[! is.na(pvals)]
  if (length(pvals) == 0 || all(pvals > 0.05)) {
    return(0)
  } else {
    return(max(abs(per_taxon_fold_changes[[i]][pvals <= 0.05])))
  }
}))
x <- filter_taxa(x, per_taxon_max_change > 2.5, supertaxa = TRUE, reassign_obs = c(diff_table = FALSE))

# Plot results (might take a few minutes)
heatmap(x,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = ifelse(is.na(padj) | padj > 0.05, 0, log2FoldChange),
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 fold change")

ps.bac_metacoder$data$otu_props <- calc_obs_props(ps.bac_metacoder, data = "otu_table", cols = ps.bac_metacoder$data$sample_data$sample_id)   ####Convert to proportions
ps.bac_metacoder$data$tax_table <- calc_taxon_abund(ps.bac_metacoder, data = "otu_table", cols = ps.bac_metacoder$data$sample_data$sample_id)   ####Convert to proportions
ps.bac_metacoder$data$sample_data <- ps.bac_metacoder$data$sample_data[-18,] #Remove negative control 
Comparison_group <- "Time"
group <- unlist(ps.bac_metacoder$data$sample_data[,Comparison_group])
ps.bac_metacoder$data$diff_table <- calc_diff_abund_deseq2(ps.bac_metacoder, data = "tax_table", 
                                                     cols = ps.bac_metacoder$data$sample_data$sample_id, 
                                                     groups = group)
range(ps.bac_metacoder$data$diff_table$padj, finite = TRUE)   ###Calculate differences between groups and multiple test correction.
range(ps.bac_metacoder$data$diff_table$log2FoldChange, finite = TRUE)
ps.bac_metacoder$data$diff_table$log2FoldChange[ps.bac_metacoder$data$diff_table$padj > 0.05] <- 0
ps.bac_metacoder$data$diff_table$log2FoldChange[is.na(ps.bac_metacoder$data$diff_table$padj) == TRUE] <- 0

Tree <- heat_tree_matrix(ps.bac_metacoder,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = ifelse(is.na(padj) | padj > 0.05, 0, log2FoldChange),
                 node_color_range = c("#238b45","grey","#6a51a3"),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 fold change")
Tree

ggsave(Tree, "tree.tiff", device = "tiff", dpi = 600)
ps.bac_metacoder$data$diff_table <- subset(ps.bac_metacoder$data$diff_table, 
                                           abs(log2FoldChange) >0.5)

