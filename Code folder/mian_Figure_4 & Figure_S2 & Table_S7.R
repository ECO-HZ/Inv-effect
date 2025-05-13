library(ggplot2)
library(openxlsx)
library(picante)
library(ggExtra)
library(ggrepel)
library(caper)
library(dplyr)
library(gghalves)
library(psych)
library(corrplot)

mytheme = theme(
  legend.position = "none",
  panel.grid=element_blank(), 
  strip.background = element_rect(color="black", fill="white", size=0.5, linetype="solid"),
  strip.text.x = element_text(size = 9, color = "black"), # face = "bold.italic"
  legend.title = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 10),
  legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
  axis.ticks = element_line(color='black'),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(colour='black', size=12),
  axis.title.y = element_text(colour='black', size=12),
  axis.text = element_text(colour='black',size=11),
  plot.background = element_blank(), 
  plot.tag = element_text(size = 13, face = "bold")) 

# loading functional traits data
traits_data <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
colnames(traits_data)
traits_mean <- unique(traits_data[,c("Species","Origin","Hmax", "Chol", "LA", "SLA", "LDMC", "SRL", "FRR", "RMF")])
rownames(traits_mean) <- traits_mean$Species

# Loading phylogenetic files of study species
plant_tree <- read.tree("IQ_tree_plant_2025.NEWICK")
to_drop <- c ("Amborella_trichopoda","")
plant_tree <- drop.tip(as.phylo(plant_tree), to_drop) 
#plant_dist <- cophenetic.phylo(plant_tree)

traits_mean = traits_mean[plant_tree$tip.label,]
# Data transformation
traits_mean$SRL = log10(traits_mean$SRL)

################## Phylogenetic principal component analysis ###################
# functional traits
select_traits = c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF") 
phylpcaall6 <- phytools::phyl.pca(plant_tree, traits_mean[,select_traits], mode="corr", method="lambda")
#print(phylpcaall6)   #All relevant information
#plot(phylpcaall6)    #Eigenvalue
pca_loadings <- data.frame(Variables = rownames(phylpcaall6$L),phylpcaall6$L)

################################# Table S7 #####################################
## Cumulative Proportion
ephylall6 <- summary(phylpcaall6)$importance
#print(ephylall6)

## loading
print(phylpcaall6$L)

## created the database
df_pca <- as.data.frame(phylpcaall6$S[ ,1:4])
df_pca$Species <- rownames(df_pca)

PCA_fun <- df_pca %>% left_join(traits_mean[,c("Species", "Origin")], by = "Species") 
PCA_fun$Origin <- factor(PCA_fun$Origin, levels = c("Native","Exotic"))

pgls_data <- comparative.data(phy = plant_tree, data = PCA_fun, 
                              names.col = Species, vcv = TRUE, 
                              na.omit = FALSE, warn.dropped = TRUE)

pgls_PC1 <- pgls(PC1 ~ Origin, data = pgls_data, lambda = "ML")
anova(pgls_PC1)

pgls_PC2 <- pgls(PC2 ~ Origin, data = pgls_data, lambda = "ML")
anova(pgls_PC2)

#
ggplot(data=PCA_fun, aes(x = PC1, y = PC2))+
  geom_point(aes(shape = Origin, color = Origin, fill = Origin), size=2.2, show.legend = T)+
  geom_segment(data = pca_loadings, mapping = aes(x = 0, y = 0, xend = (PC1*1.5), yend = (PC2*1.5)), 
               color = "black", alpha=0.8) +
  annotate("text", x = (pca_loadings$PC1*2), y = (pca_loadings$PC2*2),label = pca_loadings$Variables, size = 3.5) +
  labs(x = paste("PCA1 (", sprintf("%.1f", ephylall6[2, 1] * 100), "%)", sep = ""),
       y = paste("PCA2 (", sprintf("%.1f", ephylall6[2, 2] * 100), "%)", sep = ""), tag = "b") +
  geom_hline(yintercept = 0,linetype = 2)+ 
  geom_vline(xintercept = 0,linetype = 2)+
  theme_classic() + mytheme + 
  theme(legend.position=c(0.15,0.85)) + 
  scale_shape_manual(values = c(16,15)) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  scale_color_manual(values = c("#70A7C3","#A67C2A")) + 
  scale_fill_manual(values = c("#70A7C3","#A67C2A")) -> mian_Fig_4b#; mian_Fig_4b

ggMarginal(mian_Fig_4b, type = "boxplot", groupColour = F, 
           groupFill = TRUE, alpha = 1, color = "black", size = 8) -> mian_Fig_4b; mian_Fig_4b

################################################################################
# Loading phylogenetic files of study species
phylo_cov <- vcv.phylo(plant_tree, corr = FALSE)
#if(any(eigen(phylo_cov)$values < 0)) stop("非正定矩阵需校正")
pca <- prcomp(phylo_cov, scale. = FALSE)
phylo_impor <- summary(pca)$importance
print(as.data.frame(summary(pca)$importance)[,c(1:6)])

df_pca <- as.data.frame(pca$x)[,c(1,2)] 
df_pca$Species <- rownames(df_pca)

PCA_phylo <- df_pca %>% left_join(traits_mean[,c("Species", "Origin")], by = "Species") 
PCA_phylo$Origin <- factor(PCA_phylo$Origin, levels = c("Native","Exotic"))

pgls_data <- comparative.data(phy = plant_tree, data = PCA_phylo, 
                              names.col = Species, vcv = TRUE, 
                              na.omit = FALSE, warn.dropped = TRUE)

pgls_PC1 <- pgls(PC1 ~ Origin, data = pgls_data, lambda = "ML")
anova(pgls_PC1)

pgls_PC2 <- pgls(PC2 ~ Origin, data = pgls_data, lambda = "ML")
anova(pgls_PC2)

ggplot(data = PCA_phylo, aes(x = PC1, y = PC2))+
  geom_point(aes(shape = Origin, color = Origin, fill = Origin), size=2.2, show.legend = T)+
  labs(x = paste("PCA1 (", sprintf("%.1f", phylo_impor[2, 1] * 100), "%)", sep = ""),
       y = paste("PCA2 (", sprintf("%.1f", phylo_impor[2, 2] * 100), "%)", sep = ""), tag = "c") +
  geom_hline(yintercept=0,linetype=2)+ 
  geom_vline(xintercept=0,linetype=2)+
  theme_classic() + mytheme + 
  scale_shape_manual(values = c(16,15)) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  scale_color_manual(values = c("#70A7C3","#A67C2A")) + 
  scale_fill_manual(values = c("#70A7C3","#A67C2A")) -> mian_Fig_4c#; mian_Fig_4c

ggMarginal(mian_Fig_4c, type = "boxplot", groupColour = F, 
           groupFill = TRUE, alpha = 1, color = "black", size = 8) -> mian_Fig_4c; mian_Fig_4c


################## Phylogenetic Generalized Least Squares ######################
PCA_phylo2 <- PCA_phylo
PCA_fun2 <- PCA_fun
colnames(PCA_fun2)[1:4] <- c("Funct_PC1", "Funct_PC2", "Funct_PC3", "Funct_PC4")
colnames(PCA_phylo2)[1:2] <- c("Phylo_PC1", "Phylo_PC2")
traits_mean_pca <- traits_mean %>% left_join(PCA_fun2) %>% left_join(PCA_phylo2)

vector <- c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF","Funct_PC1", "Funct_PC2","Phylo_PC1", "Phylo_PC2")

pgls_final <- NULL

for (i in vector) {
  traits_all1 <- traits_mean_pca[,c(i, "Origin", "Species")]
  pgls_data <- comparative.data(phy = plant_tree, data = traits_all1, 
                                names.col = Species, vcv = TRUE, 
                                na.omit = FALSE, warn.dropped = TRUE)
  colnames(pgls_data$data)[1] <- "Traits"
  model.pgls <- pgls(Traits ~ Origin, data = pgls_data, lambda = "ML")
  aov_data <- as.data.frame(anova(model.pgls))
  pgls_summry <- data.frame(Traits = i, df = "1,52", F = round(aov_data$`F value`[1], 2), p = round(aov_data$`Pr(>F)`[1], 3))
  pgls_final <- rbind(pgls_final, pgls_summry)
}

print(pgls_final)

traits_mean_long <- traits_mean_pca %>%
  tidyr::pivot_longer(cols = c(Chol, SLA, LDMC, SRL, FRR, RMF),
                      names_to = "Traits",   
                      values_to = "Traits_values")

traits_mean_long$Origin = factor(traits_mean_long$Origin, levels = c("Native", "Exotic"))
traits_mean_long$Traits = factor(traits_mean_long$Traits, levels = c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF"))

# Replace the labels
Traits_label = c(Chol = "Leaf chlorophyll (SPAD)", SLA = "Specific leaf area (cm2 g-1)", LDMC = "Leaf dry matter content (g g-1)",
                 SRL = "Specific root length (cm g-1, log10)", FRR = "Fine-to-total root mass (g g-1)", RMF = "Root mass fraction (g g-1)")

ggplot(traits_mean_long, aes(x = Origin,y = Traits_values, fill = Origin, shape = Origin))+
  geom_half_violin(position=position_nudge(x=0.15, y=0),side='R',adjust=1.2,trim=T,color=NA,alpha=0.8) +
  geom_boxplot(width = 0.2, alpha = 1, outliers = FALSE) + 
  geom_point(aes(x = Origin, y = Traits_values, fill = Origin), color = "black",  
             position = position_jitter(width = 0.08), size = 1.5,alpha = 0.6) + 
  labs(x = '', y = 'Traits value', tag = "a") + 
  theme_bw() + mytheme + 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 35)) + 
  scale_shape_manual(values = c(21, 22)) + 
  scale_color_manual(values = c("#70A7C3","#A67C2A")) + 
  scale_fill_manual(values = c("#70A7C3","#A67C2A")) + 
  facet_wrap(Traits ~., ncol = 2, nrow = 5, scales = "free_y",
             labeller = labeller(Traits = Traits_label)) -> mian_Fig_4a; mian_Fig_4a


################################### Fig.S2 #####################################
# Phylogenetic pearson’s correlation coefficients among six plant functional traits
rownames(traits_mean_pca) <- traits_mean_pca$Species
traits_mean_pca <- traits_mean_pca[plant_tree$tip.label,]
corr_plot <- traits_mean_pca[,c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF","Funct_PC1", "Funct_PC2","Phylo_PC1", "Phylo_PC2")]

colnames(corr_plot) <- c("Leaf chlorophyll", "Specific leaf area", "Leaf dry matter content", 
                         "Specific root length", "Fine-to-total root mass", "Root mass fraction",
                         "The first PCA axes of plant traits", "The second PCA axes of plant traits",
                         "The first PCA axes of plant phylogeny", "The second PCA axes of plant phylogeny")


obj <- phytools::phyl.vcv(as.matrix(corr_plot), vcv(plant_tree),1)
## correlation between x & y
r.xy<-cov2cor(obj$R)
diag(r.xy) = 0
r_corr_matrix <- r.xy

## t-statistic & P-value
t.xy<-r.xy*sqrt((Ntip(plant_tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(plant_tree)-2,lower.tail=F)
p_corr_matrix <- P.xy

#corr_matrix <- corr.test(as.matrix(scale(corr_plot)), use = "pairwise", method = "pearson",adjust="none", alpha=.05,ci=FALSE)
#r_corr_matrix <- as.matrix(corr_matrix$r)
#p_corr_matrix <- as.matrix(corr_matrix$p)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r_corr_matrix, p.mat =p_corr_matrix, sig.level = 0.05, insig = 'blank', method = 'number',type = 'upper',
         diag = FALSE, col=col(200), tl.cex = 0.8,tl.col = "black", number.cex = 0.6, order = "original",tl.srt = 45)

corrplot(r_corr_matrix, p.mat = p_corr_matrix, sig.level = 0.05, insig = 'blank', method = 'square',
         add = TRUE, type = 'upper', diag = FALSE, col=col(200), tl.pos = 'n', cl.pos = 'n',outline = F, order = "original",
         addCoef.col = "black", number.cex = 0.6)

# Notice that,
# For more picture details, we have further adjusted it in Adobe illustrator.








