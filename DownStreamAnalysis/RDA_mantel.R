# How to integrate environmental/climate variables with sequencing data

library(corrplot)
library(microbiome)
library(phyloseq)
library(vegan)

# How many samples do you have? How many variables do you have?
# if number of variables > number of samples, RDA/CCA are not appropriate
ps <- readRDS("ps.bac.rds")
ps.bac <- subset_taxa(ps, Kingdom == "Bacteria")

ps.rel  = transform_sample_counts(ps.bac, function(x) x / sum(x) )

# first remove samples that do not have the ancillary data
ps.chem <- subset_samples(ps.rel, Manganese_PPM != "NA")
# let's start with a PCoA ####

bac.bray = ordinate(ps.chem, "PCoA", "bray")
p <- plot_ordination(ps.chem, bac.bray, "samples", color="Depth", shape = "Time")

# assess whether any chemical variables are associated with the distribution along the x-axis
A1 <- p$data$Axis.1

chemistry <- data.frame(sample_data(ps.chem))
chemistry[1:8] <- NULL
chemistry$CEC.H_. <- NULL
chemistry$Limestone_. <- NULL
chemistry$Moisture_Analysis <- NULL
chemistry$A1 <- A1
chemistry$LimeReq_Tons.AF <- NULL
chemistry$GypReq_Calc_Tons.AF <- NULL
chemistry$CEC.Na_. <- NULL
# rda with stepwise selection
filt.asv <- data.frame(otu_table(ps.chem))
# check axis lengths (if > 4, transform via hellinger or use CCA)
decorana(filt.asv)
chemistry$A1 <- NULL
chemistry.log  <- log(chemistry+1)
rda.1 <- rda(decostand(filt.asv, 'hellinger') ~ ., chemistry.log, scaling = 'sites')
envfit(rda.1 ~., chemistry)

# let's keep significant variables
chemistry.log[,1:10] <- NULL
chemistry.log[,3:16] <- NULL

rda.2 <- rda(decostand(filt.asv, 'hellinger') ~ ., chemistry.log, scaling = 'sites')

safe_colorblind <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                     "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

metadata <- data.frame(sample_data(ps.chem))
pch2 <- c(15, 16)

plot(rda.2, type="n", display="sites",xlab = paste("PCA1"),
     ylab = paste("PCA2"), cex.lab = 1.5, cex.axis = 1.2)
points(rda.2, display="bp", lwd=2, col="black")
text(rda.2, display="bp", col="black", font=2, cex = 0.5)
points(rda.2, display="sites", col=safe_colorblind[metadata$Time],cex=3, font=1, pch = pch2[metadata$Depth])
legend("bottomright", legend = levels(metadata$Time), col = safe_colorblind[metadata$Time], pch = 16, cex = 1, title = "Time")
legend("bottomright", legend = levels(metadata$Depth), col = "black", pch = c(15,16), cex = 1, title = "Depth", inset = c(0.37, 0.00))


# Mantel Tests - Bacteria (Sulfate)####

x.dist<- vegdist(chemistry$Sulfate_meq.L, "euclidean")
y.dist <- vegdist(filt.asv, "bray")

mantel(x.dist, y.dist, method ="spearman", permutations = 999)

lm.all <- lm(y.dist~x.dist)
plot(y.dist~x.dist, ylab = "Bray-Curtis Dissimilarity (16S)", xlab = "Sulfate (meq)", pch = 19, frame = FALSE, cex = 1, font = 2)
abline(lm.all, col = "blue", lty =1, lwd=4 )

x.dist<- vegdist(chemistry, "euclidean")
y.dist <- vegdist(filt.asv, "bray")

# Mantel Tests - Bacteria (Zinc)####
x.dist <- vegdist(chemistry$Zinc_PPM, "euclidean")
y.dist <- vegdist(filt.asv, "bray")
mantel(x.dist, y.dist, method ="spearman", permutations = 999)
lm.all <- lm(y.dist~x.dist)
plot(y.dist~x.dist, ylab = "Bray-Curtis Dissimilarity (16S)", xlab = "Zn (ppn) (distance/dissimality)", pch = 19, frame = FALSE, cex = 1, font = 2)
abline(lm.all, col = "blue", lty =1, lwd=4 )




# fungal dataset ####
ps.fun.rel  = transform_sample_counts(ps.fun, function(x) x / sum(x) )

#  remove samples that do not have the ancillary data 
ps.fun.chem <- subset_samples(ps.fun.rel, Chemistry_data == "Y")
ps.fun.chem.1 = filter_taxa(ps.fun.chem, function(x) sum(x) > 0, TRUE)

# Ordinations: PCoA ASV-Level ####

fun.bray = ordinate(ps.fun.chem.1, "PCoA", "bray")
plot_ordination(ps.fun.chem.1, fun.bray, "biplot", color="Depth", shape = "Time")
p <- plot_ordination(ps.fun.chem.1, fun.bray, "sites", color="Depth", shape = "Time")
p

# assess whether any chemical variables are associated with the distribution along the x-axis
A1 <- p$data$Axis.1

fun.chemistry <- data.frame(sample_data(ps.fun.chem))
fun.chemistry[1:8] <- NULL
fun.chemistry$CEC.H_. <- NULL
fun.chemistry$Limestone_. <- NULL
fun.chemistry$Moisture_Analysis <- NULL
fun.chemistry$A1 <- A1
fun.chemistry$LimeReq_Tons.AF <- NULL
fun.chemistry$GypReq_Calc_Tons.AF <- NULL
fun.chemistry$CEC.Na_. <- NULL

# rda with stepwise selection
filt.fun.asv <- data.frame(otu_table(ps.fun.chem.1))
# check axis lengths (if > 4, transform via hellinger or use CCA)
decorana(filt.fun.asv)
fun.chemistry$A1 <- NULL
fun.chemistry.log  <- log(fun.chemistry+1)
rda.1 <- rda(decostand(filt.fun.asv, 'hellinger') ~ ., fun.chemistry.log, scale = TRUE)
envfit(rda.1 ~., chemistry)


fun.chemistry.log[1:14] <- NULL
fun.chemistry.log[2:10] <- NULL
fun.chemistry.log[3] <- NULL

rda.2 <- rda(decostand(filt.fun.asv, 'hellinger') ~ ., fun.chemistry.log)
fun.metadata <- data.frame(sample_data(ps.fun.chem))

plot(rda.2, type="n", display="sites",xlab = paste("PCA1"),
     ylab = paste("PCA2"), cex.lab = 1.5, cex.axis = 1.2)
points(rda.2, display="bp", lwd=2, col="black")
text(rda.2, display="bp", col="black", font=2, cex = 0.5)
points(rda.2, display="sites", col=safe_colorblind[fun.metadata$Time],cex=3, font=1, pch = pch2[fun.metadata$Depth])
legend("bottomright", legend = levels(fun.metadata$Time), col = safe_colorblind[fun.metadata$Time], pch = 16, cex = 0.7, title = "Time")
legend("bottomright", legend = levels(fun.metadata$Depth), col = "black", pch = c(15,16), cex = 0.7, title = "Depth", inset = c(0, 0.2))



