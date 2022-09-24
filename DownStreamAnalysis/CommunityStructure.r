#For bacterial community (Phylum level)
ps.bac <- subset_samples(ps.bac, SampleID != "Negative_control")
glom <- tax_glom(ps.bac, taxrank = "Phylum", NArm=TRUE)
phylum_name <- as.character(tax_table(glom)[,"Phylum"])
taxa_names(glom) <- phylum_name
phylum_10 = names(sort(taxa_sums(glom), FALSE)[1:(length(phylum_name)-10)])
glom_rel <- transform_sample_counts(glom, function(x) x/sum(x))
data_glom <- psmelt(glom_rel)
for (i in 1:length(phylum_10)){
  data_glom$OTU[data_glom$OTU == phylum_10[i]] <- "Other"
}

plot <- ggplot(data=data_glom, 
               aes(x=factor(Sample, 
                            level = c("1A.2cm.1D", 
                                      "2A.2cm.1D",
                                      "3A.2cm.1D", 
                                      "1B.15cm.1D",
                                      "2B.15cm.1D", 
                                      "3B.15cm.1D",
                                      "1A.2cm.1M", 
                                      "2A.2cm.1M",
                                      "3A.2cm.1M", 
                                      "1B.15cm.1M",
                                      "2B.15cm.1M", 
                                      "3B.15cm.1M",
                                      "1A.2cm.2M", 
                                      "2A.2cm.2M",
                                      "3A.2cm.2M", 
                                      "1B.15cm.2M",
                                      "2B.15cm.2M", 
                                      "3B.15cm.2M",
                                      "1A.2cm.5M", 
                                      "2A.2cm.5M",
                                      "3A.2cm.5M", 
                                      "1B.15cm.5M",
                                      "2B.15cm.5M", 
                                      "3B.15cm.5M")),
                                   y=Abundance, fill = forcats::fct_infreq(OTU))) + 
  geom_bar(stat="identity", position="stack") + 
  guides(fill=guide_legend(nrow=11)) + 
  scale_fill_manual(values=col)+
  theme(axis.text.x = element_text(angle=90), 
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        legend.title = element_blank()) +
  xlab("Samples") + 
  ylab("Relative abundance")
plot

                                    
#For bacterial community (Genus level)
glom <- tax_glom(ps.bac, taxrank = "Genus", NArm=TRUE)
genus_name <- as.character(tax_table(glom)[,"Genus"])
taxa_names(glom) <- genus_name
genus_10 = names(sort(taxa_sums(glom), FALSE)[1:(length(genus_name)-10)])
glom_rel <- transform_sample_counts(glom, function(x) x/sum(x))
data_glom <- psmelt(glom_rel)
for (i in 1:length(genus_10)){
  data_glom$Genus[data_glom$OTU == genus_10[i]] <- "Other"
}

plot <- ggplot(data=data_glom, 
               aes(x=factor(Sample, 
                            level = c("1A.2cm.1D", 
                                      "2A.2cm.1D",
                                      "3A.2cm.1D", 
                                      "1B.15cm.1D",
                                      "2B.15cm.1D", 
                                      "3B.15cm.1D",
                                      "1A.2cm.1M", 
                                      "2A.2cm.1M",
                                      "3A.2cm.1M", 
                                      "1B.15cm.1M",
                                      "2B.15cm.1M", 
                                      "3B.15cm.1M",
                                      "1A.2cm.2M", 
                                      "2A.2cm.2M",
                                      "3A.2cm.2M", 
                                      "1B.15cm.2M",
                                      "2B.15cm.2M", 
                                      "3B.15cm.2M",
                                      "1A.2cm.5M", 
                                      "2A.2cm.5M",
                                      "3A.2cm.5M", 
                                      "1B.15cm.5M",
                                      "2B.15cm.5M", 
                                      "3B.15cm.5M")),
                                   y=Abundance, fill = forcats::fct_infreq(Genus))) + 
  geom_bar(stat="identity", position="stack") + 
  guides(fill=guide_legend(nrow=11)) + 
  scale_fill_manual(values=col)+
  theme(axis.text.x = element_text(angle=90), 
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        legend.title = element_blank()) +
  xlab("Samples") + 
  ylab("Relative abundance")
plot                                    
                                    
                                    
                                    
                                    
#For fungal community
ps.fun <- subset_samples(ps.fun, SampleID != "Negative_control")
glom <- tax_glom(ps.fun, taxrank = "Phylum", NArm=TRUE)
phylum_name <- as.character(tax_table(glom)[,"Phylum"])
taxa_names(glom) <- phylum_name
glom_rel <- transform_sample_counts(glom, function(x) x/sum(x))
data_glom <- psmelt(glom_rel)
plot <- ggplot(data=data_glom, 
               aes(x=factor(Sample, 
                            level = c("1A.2cm.1D", 
                                      "2A.2cm.1D",
                                      "3A.2cm.1D", 
                                      "1B.15cm.1D",
                                      "2B.15cm.1D", 
                                      "3B.15cm.1D",
                                      "1A.2cm.1M", 
                                      "2A.2cm.1M",
                                      "3A.2cm.1M", 
                                      "1B.15cm.1M",
                                      "2B.15cm.1M", 
                                      "3B.15cm.1M",
                                      "1A.2cm.2M", 
                                      "2A.2cm.2M",
                                      "3A.2cm.2M", 
                                      "1B.15cm.2M",
                                      "2B.15cm.2M", 
                                      "3B.15cm.2M",
                                      "1A.2cm.5M", 
                                      "2A.2cm.5M",
                                      "3A.2cm.5M", 
                                      "1B.15cm.5M",
                                      "2B.15cm.5M", 
                                      "3B.15cm.5M")),
                                   y=Abundance, fill = forcats::fct_infreq(OTU))) + 
  geom_bar(stat="identity", position="stack") + 
  guides(fill=guide_legend(nrow=11)) + 
  scale_fill_manual(values=col)+
  theme(axis.text.x = element_text(angle=90), 
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        legend.title = element_blank()) +
  xlab("Samples") + 
  ylab("Relative abundance")
plot                             
                                    
#For fungal community (Genus level)
glom <- tax_glom(ps.fun, taxrank = "Genus", NArm=TRUE)
genus_name <- as.character(tax_table(glom)[,"Genus"])
taxa_names(glom) <- genus_name
genus_10 = names(sort(taxa_sums(glom), FALSE)[1:(length(genus_name)-10)])
glom_rel <- transform_sample_counts(glom, function(x) x/sum(x))
data_glom <- psmelt(glom_rel)
for (i in 1:length(genus_10)){
  data_glom$Genus[data_glom$OTU == genus_10[i]] <- "Other"
}

plot <- ggplot(data=data_glom, 
               aes(x=factor(Sample, 
                            level = c("1A.2cm.1D", 
                                      "2A.2cm.1D",
                                      "3A.2cm.1D", 
                                      "1B.15cm.1D",
                                      "2B.15cm.1D", 
                                      "3B.15cm.1D",
                                      "1A.2cm.1M", 
                                      "2A.2cm.1M",
                                      "3A.2cm.1M", 
                                      "1B.15cm.1M",
                                      "2B.15cm.1M", 
                                      "3B.15cm.1M",
                                      "1A.2cm.2M", 
                                      "2A.2cm.2M",
                                      "3A.2cm.2M", 
                                      "1B.15cm.2M",
                                      "2B.15cm.2M", 
                                      "3B.15cm.2M",
                                      "1A.2cm.5M", 
                                      "2A.2cm.5M",
                                      "3A.2cm.5M", 
                                      "1B.15cm.5M",
                                      "2B.15cm.5M", 
                                      "3B.15cm.5M")),
                                   y=Abundance, fill = forcats::fct_infreq(Genus))) + 
  geom_bar(stat="identity", position="stack") + 
  guides(fill=guide_legend(nrow=11)) + 
  scale_fill_manual(values=col)+
  theme(axis.text.x = element_text(angle=90), 
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        legend.title = element_blank()) +
  xlab("Samples") + 
  ylab("Relative abundance")
plot                                     
                                    
                                    
                                    
