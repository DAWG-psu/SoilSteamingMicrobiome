##Figure 3 - Principal Coordinates analysis (PCoA)
PCoA_16s <- ordinate(ps.bac.rare, "PCoA", "bray")
PCoA_16s_plot <- plot_ordination(ps.bac.2, PCoA_16s, shape = "Depth", color = "Time") +
  geom_point(size=4) + theme_bw() +scale_color_manual(values = c("#238B45","#6A51A3", "#A1D99B","#BCBDDC"))
PCoA_16s_plot

PCoA_16s <- ordinate(ps.fun.rare, "PCoA", "bray")
PCoA_16s_plot <- plot_ordination(ps.fun.rare, PCoA_16s, shape = "Depth", color = "Time") +
  geom_point(size=4) + theme_bw() +scale_color_manual(values = c("#238B45","#6A51A3", "#A1D99B","#BCBDDC"))
PCoA_16s_plot

##Table S4 - PErmutational multivariate analysis of variance and beta-dispersion
ibrary(vegan)
metadata <- as(sample_data(ps.bac.rare), "data.frame")
adonis(phyloseq::distance(ps.bac.rare, method="bray") ~ SampleID,
       data = metadata)
permutest(betadisper(phyloseq::distance(ps.bac.rare, method="bray"), metadata$SampleID), permutations= 99)

metadata <- as(sample_data(ps.fun.rare), "data.frame")
adonis(phyloseq::distance(ps.fun.rare, method="bray") ~ SampleID,
       data = metadata)
permutest(betadisper(phyloseq::distance(ps.fun.rare, method="bray"), metadata$SampleID), permutations= 99)
