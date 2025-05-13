library(funrar)
library(ggplot2)
library(openxlsx)
library(vegan)
library(phytools)
library(picante)
library(geiger)
library(GUniFrac)
library(patchwork)


###########################  Greenhouse experiment #############################
green_otu <- read.xlsx("Greenhouse_fungi_Flattening.xlsx", sheet = "green_flattening", colNames = T, rowNames = T)
green_otu[1:6,1:6]
green_otu <- green_otu[,-c(1)]

## load group data
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)
colnames(Green_group)
Green_group = Green_group[colnames(green_otu), ]
rownames(t(green_otu)) %in% rownames(Green_group)
Green_group$SRL <- log10(Green_group$SRL)

## Richness of overall fungi
green_richness <- as.data.frame(specnumber(t(green_otu)))
green_richness$Sample_ID = rownames(green_richness); colnames(green_richness)[1] = "SR"
Green_group = Green_group %>% left_join(green_richness)
rownames(Green_group) = Green_group$Sample_ID

############################## Part of Table S4 ################################
## Richness 
green_SR_mod = lm(SR ~ Origin/Species, data = Green_group)
#shapiro.test(residuals(green_SR_mod))
#anova(green_SR_mod)
Richness_lm_sum = as.data.frame(anova(green_SR_mod))
Richness_lm_sum$`p.adj`=p.adjust(Richness_lm_sum$`Pr(>F)`, method = "BH")
print(Richness_lm_sum)

## composition of rhizosphere overall fungi
green_hel_no = t(green_otu)/9690
rowSums(green_hel_no)
Bray_dist_green_no <- vegdist(green_hel_no, method = 'bray')

############################## Part of Table S4 ################################
## PERMANOVA
Green_group = Green_group[rownames(as.matrix(Bray_dist_green_no)),]

set.seed(1234)
green_perMAONA2 = GUniFrac::adonis3(green_hel_no ~ Origin/Species, method = "bray", by = "margin", 
                                    data = Green_group, permutations = 9999)
green_perMAONA2$aov.tab
permanova_green = as.data.frame(green_perMAONA2$aov.tab)
permanova_green$R2 = round(permanova_green$R2,3)
permanova_green$`F` = round(permanova_green$`F`,2)
permanova_green$`q-vaules`=p.adjust(permanova_green$`Pr(>F)`, method = "BH")
print(permanova_green)


############################## Part of Table S7 ################################
# Constructing species identity matrix
species_factor <- as.character(Green_group$Species)
sample_names <- Green_group$Sample_ID
species_matrix <- outer(species_factor, species_factor, FUN = function(x, y) ifelse(x == y, 1, 0))
rownames(species_matrix) <- sample_names
colnames(species_matrix) <- sample_names
species_dist <- as.dist(species_matrix)

rownames(as.matrix(Bray_dist_green_no)) %in% rownames(species_matrix)
rownames(as.matrix(Bray_dist_green_no)) %in% rownames(as.matrix(traits_dis))

# partial Mantel statistic
trait_names <- c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF",
                "Funct_pPC1", "Funct_pPC2", "Phylo_vPC1", "Phylo_vPC2")
Green_Bray_Curtis_mantel <- NULL

for (i in trait_names) {
  all_otu_TEST <- Green_group[,c(i, "Species")]
  all_otu_TEST$Species <- NULL
  #traits_dis <- compute_dist_matrix(all_otu_TEST, metric = "euclidean", scale = TRUE, center = TRUE) 
  traits_dis <- vegdist(all_otu_TEST, method = "euclidean") 
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel.partial(Bray_dist_green_no, (traits_dis), species_dist, method = 'spearman', permutations = 9999, na.rm = TRUE)
  mantel_Bray_Curtis_result <- data.frame(mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif,i,"Fungal composition","Traits",'Bray_Curtis')
  colnames(mantel_Bray_Curtis_result) <- c("Mantel_R","P_value","from","to","Group","dist_type")
  Green_Bray_Curtis_mantel <- rbind(Green_Bray_Curtis_mantel,mantel_Bray_Curtis_result)
}

Green_Bray_Curtis_mantel$p.adj <- round(p.adjust(Green_Bray_Curtis_mantel$P_value, method = "BH"), 3)
print(Green_Bray_Curtis_mantel)


##############################  Field experiment ###############################
# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Data Transformation
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Year <- as.factor(Field_group$Year)
Field_group$Site <- factor(Field_group$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Field_group$Origin <- factor(Field_group$Origin, levels = c("Native","Exotic"))

# notes: I have completed the above work, so I directly load the completed file
fungi_Flattening <- read.xlsx("Field_fungi_Flattening.xlsx", sheet = "field_flattening", rowNames = T, colNames = T)
fungi_Flattening <- fungi_Flattening[,-c(1)]
fungi_Flattening[1:6, 1:6]
colSums(fungi_Flattening)

# Z-scored standardized before analyses
Field_group <- Field_group[colnames(fungi_Flattening), ] 
colnames(Field_group)
select_var_scale <- c("Site_pool","Soil_ph", "Wcont","Soil_N","Tave","Prec","CV_Tave","CV_Prec",
                      "Chol","SLA","LDMC","SRL","FRR","RMF","Funct_pPC1","Funct_pPC2","Phylo_vPC1","Phylo_vPC2")

pd_attributes_variable <- attributes(scale(Field_group[select_var_scale]))

Field_group_scale = Field_group
colnames(Field_group_scale)

Field_group_scale[select_var_scale] = scale(Field_group_scale[select_var_scale])


############################## Part of Table S4 ################################
############################## Fungal richness #################################
field_richness <- as.data.frame(specnumber(t(fungi_Flattening)))
colnames(field_richness) <- "SR"
field_richness$Sample_ID <- rownames(field_richness)
Field_group <- Field_group %>% left_join(field_richness)
rownames(Field_group) <- Field_group$Sample_ID

Fungal_SR_mod <- lm(SR ~ Year + Site + Origin/Species + Year:Site + 
                     Origin:Year + Origin:Site + Origin:Year:Site + 
                     (Origin/Species):Year + (Origin/Species):Site, data = Field_group)
anova(Fungal_SR_mod)
SR_mod_anova <- as.data.frame(anova(Fungal_SR_mod))
SR_mod_anova$p_adj <- round(p.adjust(SR_mod_anova$`Pr(>F)`, method = "BH"), 3) 
SR_mod_anova[, 2:4] <- round(SR_mod_anova[, 2:4], digits = 2)
SR_mod_anova$`Pr(>F)` <- round(SR_mod_anova$`Pr(>F)`, 3)
SR_mod_anova <- SR_mod_anova[c(1:7,9,10,8),] # Reorder
print(SR_mod_anova) 

############################## Part of Table S4 ################################
# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
field_hel_no = t(fungi_Flattening)/9477
rowSums(field_hel_no)
Bray_dist_field_no <- vegdist(field_hel_no, method = 'bray')

set.seed(1234)
permanova_mod = GUniFrac::adonis3(field_hel_no ~ Year + Site + Origin/Species + Year:Site + 
                                    Origin:Year + Origin:Site + Origin:Year:Site + 
                                    (Origin/Species):Year + (Origin/Species):Site, method = "bray", by = "margin",
                                  data = Field_group_scale, permutations = 9999)

permanova_mod_sum = as.data.frame(permanova_mod$aov.tab)
permanova_mod_sum$R2 = round(permanova_mod_sum$R2,3)
permanova_mod_sum$df = paste0(permanova_mod_sum$Df, ",", permanova_mod_sum$Df[length(permanova_mod_sum$Df)-1])
permanova_mod_sum$`F` = round(permanova_mod_sum$F.Model,2)
permanova_mod_sum$`p.adj`= round(p.adjust(permanova_mod_sum$`Pr(>F)`, method = "BH"), 3)
print(permanova_mod_sum)

############################## Part of Table S7 ################################
# Constructing species identity matrix
species_factor <- as.character(Field_group$Species)
sample_names <- Field_group$Sample_ID
species_matrix <- outer(species_factor, species_factor, FUN = function(x, y) ifelse(x == y, 1, 0))
rownames(species_matrix) <- sample_names
colnames(species_matrix) <- sample_names
species_dist <- as.dist(species_matrix)

colnames(Field_group)
# partial Mantel statistic
trait_names <- c("Chol", "SLA", "LDMC", "SRL", "FRR", "RMF",
                 "Funct_pPC1", "Funct_pPC2", "Phylo_vPC1", "Phylo_vPC2")

env_vals <- c("Soil_ph", "Wcont", "Soil_N", "Tave", "Prec", "CV_Tave", "CV_Prec")

Field_Bray_Curtis_mantel <- NULL

#for (i in trait_names) {
for (i in env_vals) {
  all_otu_TEST <- Field_group[,c(i, "Species")]
  all_otu_TEST <- all_otu_TEST[rownames(as.matrix(Bray_dist_field_no)), ]
  all_otu_TEST$Species <- NULL
  #factors_dis <- compute_dist_matrix(all_otu_TEST, metric = "euclidean", scale = TRUE, center = TRUE) 
  factors_dis <- vegdist(all_otu_TEST, method = "euclidean") 
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel.partial((Bray_dist_field_no), (factors_dis), species_dist, method = 'spearman', permutations = 9999, na.rm = TRUE)
  mantel_Bray_Curtis_result <- data.frame(mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif,i,"Fungal composition","Factors",'Bray_Curtis')
  colnames(mantel_Bray_Curtis_result) <- c("Mantel_R","P_value","from","to","Group","dist_type")
  Field_Bray_Curtis_mantel <- rbind(Field_Bray_Curtis_mantel,mantel_Bray_Curtis_result)
}

Field_Bray_Curtis_mantel$p.adj <- round(p.adjust(Field_Bray_Curtis_mantel$P_value, method = "BH"), 3)
print(Field_Bray_Curtis_mantel)


env_vals <- c("Tave", "Prec", "Soil_ph", "Wcont", "Soil_N", "CV_Tave", "CV_Prec")
env_labs <- c("Average annual temperature (°C)", "Annual precipitation (mm)",
              "Soil pH", "Soil water content (g/g, sqrt)", "Soil total nitrogen content (%, sqrt)",
              "CV of average annual temperature", "CV of Annual precipitation")

plot_list <- list()

for (i in seq_along(env_vals)) {
  var <- env_vals[i]
  var_lab <- env_labs[i]
  
  all_otu_TEST <- Field_group[, c(var, "Species")]
  all_otu_TEST <- all_otu_TEST[rownames(as.matrix(Bray_dist_field_no)), ]
  all_otu_TEST$Species <- NULL
  
  #factors_dis <- compute_dist_matrix(all_otu_TEST, metric = "euclidean", scale = TRUE, center = TRUE) 
  factors_dis <- vegdist(all_otu_TEST, method = "euclidean") 
  #View(as.matrix(factors_dis))
  env_vect <- as.vector((factors_dis))
  fungal_vect <- as.vector(Bray_dist_field_no)
  ggplot_data <- data.frame(env_vect = env_vect, fungal_vect = fungal_vect)
  
  plot_list[[i]] <- ggplot(ggplot_data, aes(x = env_vect, y = fungal_vect)) + 
    geom_point(size = 1.5, pch = 21) + 
    geom_smooth(method = "lm", color = "#880000", se = TRUE) + 
    labs(
      x = "Environmental distance (euclidean distance)",
      y = "Fungal community dissimilarity\n(Bray–Curtis distance)",
      tag = letters[i], title = var_lab) +
    theme_bw() + mytheme + 
    theme(plot.title = element_textbox(
      size = 14, color = "black", fill = "grey90",
      box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
      halign = 0.5, width = grid::unit(1, "npc")))
  assign(paste0("p", i), plot_list[[i]])
}

(p1|p2|p3)/(p4|p5|p5)
