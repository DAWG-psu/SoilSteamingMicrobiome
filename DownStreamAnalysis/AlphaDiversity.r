##Alpha diversity (Figure S1)

ps.bac.2 <- subset_samples(ps.bac, SampleID != "Negative_control")
alpha <- plot_richness(ps.bac.2, measures= c("Chao1"), x = "SampleID", color = "SampleID",)
alpha_data <- alpha$data
ggplot(alpha_data, aes(x=factor(SampleID, level = c("2cm.1D", "15cm.1D", 
                                                    "2cm.1M", "15cm.1M",
                                                    "2cm.2M", "15cm.2M",
                                                    "2cm.5M", "15cm.5M")), y = value, fill=value)) +
  stat_summary(fun=mean, geom="bar", fill="skyblue") + 
  stat_summary(geom = "errorbar", width=0.4) + theme_bw() +
  xlab("") + ylab("Bacteria Chao1")

ps.fun.2 <- subset_samples(ps.fun, SampleID != "Negative_control")
alpha <- plot_richness(ps.fun.2, measures= c("Chao1"), x = "SampleID", color = "SampleID",)
alpha_data <- alpha$data
ggplot(alpha_data, aes(x=factor(SampleID, level = c("2cm.1D", "15cm.1D", 
                                                    "2cm.1M", "15cm.1M",
                                                    "2cm.2M", "15cm.2M",
                                                    "2cm.5M", "15cm.5M")), y = value, fill=value)) +
  stat_summary(fun=mean, geom="bar", fill="#84ac4f") + 
  stat_summary(geom = "errorbar", width=0.4) + theme_bw() +
  xlab("") + ylab("Fungi Chao1")

##Table S2, S3 - alpha diversity of each community
bac_alpha <- estimate_richness(ps.bac.2)
fun_alpha <- estimate_richness(ps.fun.2)
