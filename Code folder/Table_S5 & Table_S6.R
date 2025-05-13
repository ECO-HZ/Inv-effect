library(ggplot2)
library(ggtext)
library(openxlsx)
library(adespatial)
library(vegan)
library(phytools)
library(picante)
library(geiger)
library(dplyr)

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
  axis.title.x = element_text(colour='black', size=13),
  axis.title.y = element_text(colour='black', size=13),
  axis.text = element_text(colour='black',size=11),
  plot.title = element_textbox(
    size = 14, color = "black", fill = "grey90",
    box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
    halign = 0.5, width = grid::unit(1, "npc")), #r = unit(3, "pt")     
  plot.tag = element_text(size = 14, face = "bold")) 

################################################################################
################################ Field survey ##################################
# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Data Transformation
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Year <- as.factor(Field_group$Year)
Field_group <- Field_group[colnames(fungi_Flattening), ] 
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

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
field_hel_no <- t(fungi_Flattening)/9477
#Bray_dist_field_no <- vegdist(field_hel_no, method = 'bray')

## Tax INFORMATION
tax_default <- read.xlsx("Field_fungi_Flattening.xlsx", sheet = "ASV_tax", rowNames = F, colNames = T)[,1:3]
head(tax_default)
rownames(tax_default) <- tax_default$ASV_name
tax_default <- tax_default[colnames(field_hel_no),]
OTU_tax <- tax_default %>% tidyr::separate(col = taxonomy_all, 
                                                    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                                    sep = ";\\s*", fill = "right")
head(OTU_tax[,1:6])
OTU_tax$GENUS <- sub("g__", "", OTU_tax$Genus)

# Loading FungalTraits database
FungalTraits <- read.xlsx("FungalTraits.xlsx", sheet = "FungalTraits", rowNames = F, colNames = T)
head(FungalTraits[,c(2:6,8)])
OTU_tax2 <- OTU_tax %>% left_join(FungalTraits[,c("GENUS","primary_lifestyle")], by = "GENUS")
head(OTU_tax2)

####  |Plant Pathogen|
Pathogen_id1 <- subset(OTU_tax2, primary_lifestyle == "plant_pathogen")$ASV_name
#### Fusarium
Pathogen_id2 <- subset(OTU_tax, Genus == "g__Fusarium")$ASV_name
### |Plant Pathogen| & Fusarium
Pathogens <- as.data.frame(unique(c(Pathogen_id1, Pathogen_id2)))
Pathogens$guild <- "Plant Pathogen"; colnames(Pathogens)[1] <- "ASV_name"

############################ Arbuscular Mycorrhizal ############################
####  Arbuscular Mycorrhizal
AMF_id1 <- subset(OTU_tax2, primary_lifestyle == "arbuscular_mycorrhizal")$ASV_name; length(AMF_id1)
#### Glomeromycota
AMF_id2 <- subset(OTU_tax, Phylum == "p__Glomeromycota")$ASV_name; length(AMF_id2)
### Arbuscular Mycorrhizal & Glomeromycota
AMF <- as.data.frame(unique(c(AMF_id1, AMF_id2)))
AMF$guild <- "Arbuscular Mycorrhizal"; colnames(AMF)[1] <- "ASV_name"

################################################################################
#### Plant Pathogen
Path_field_hel_no <- as.data.frame(field_hel_no)[ ,Pathogens$ASV_name]

## Richness
field_path_SR <- as.data.frame(specnumber(Path_field_hel_no))
colnames(field_path_SR) <- "Path_SR"
field_path_SR$Sample_ID <- rownames(field_path_SR)
Field_group <- Field_group %>% left_join(field_path_SR)
rownames(Field_group)  <- Field_group$Sample_ID

Path_field_SR_mod <- lm(Path_SR ~ Year + Site + Origin/Species + Year:Site + 
                         Origin:Year + Origin:Site + Origin:Year:Site + 
                         (Origin/Species):Year + (Origin/Species):Site, data = Field_group)
anova(Path_field_SR_mod)
shapiro.test(residuals(Path_field_SR_mod))
Path_SR_mod_anova <- as.data.frame(anova(Path_field_SR_mod))
Path_SR_mod_anova$p_adj <- round(p.adjust(Path_SR_mod_anova$`Pr(>F)`, method = "BH"), 3) 
Path_SR_mod_anova[, 2:4] <- round(Path_SR_mod_anova[, 2:4], digits = 2)
Path_SR_mod_anova$`Pr(>F)` <- round(Path_SR_mod_anova$`Pr(>F)`, 3)
Path_SR_mod_anova <- Path_SR_mod_anova[c(1:7,9,10,8),] # Reorder
print(Path_SR_mod_anova) 


############################# Fungal composition ###############################
Path_Bray_dist_field_no <- vegdist(Path_field_hel_no, method = 'bray')

set.seed(1234)
Path_permanova_mod_field <- GUniFrac::adonis3(Path_field_hel_no ~ Year + Site + Origin/Species + Year:Site + 
                                    Origin:Year + Origin:Site + Origin:Year:Site + 
                                    (Origin/Species):Year + (Origin/Species):Site, method = "bray", by = "margin",
                                  data = Field_group_scale, permutations = 9999)

Path_permanova_mod_sum_field <- as.data.frame(Path_permanova_mod_field$aov.tab)
Path_permanova_mod_sum_field$R2 <- round(Path_permanova_mod_sum_field$R2,3)
Path_permanova_mod_sum_field$df <- paste0(Path_permanova_mod_sum_field$Df, ",", Path_permanova_mod_sum_field$Df[length(Path_permanova_mod_sum_field$Df)-1])
Path_permanova_mod_sum_field$`F` <- round(Path_permanova_mod_sum_field$F.Model,2)
Path_permanova_mod_sum_field$`p.adj` <- round(p.adjust(Path_permanova_mod_sum_field$`Pr(>F)`, method = "BH"), 3)
print(Path_permanova_mod_sum_field)



############################ Arbuscular Mycorrhizal ############################
AMF_field_hel_no <- as.data.frame(field_hel_no)[ ,AMF$ASV_name]

dim(subset(Field_group, AMF_SR == "0"))
dim(subset(Field_group, AMF_SR != "0"))

## Richness
field_AMF_SR <- as.data.frame(specnumber(AMF_field_hel_no))
colnames(field_AMF_SR) <- "AMF_SR"
field_AMF_SR$Sample_ID  <- rownames(field_AMF_SR)
Field_group <- Field_group %>% left_join(field_AMF_SR)
rownames(Field_group) <- Field_group$Sample_ID

AMF_field_SR_mod <- lm(AMF_SR ~ Year + Site + Origin/Species + Year:Site + 
                         Origin:Year + Origin:Site + Origin:Year:Site + 
                         (Origin/Species):Year + (Origin/Species):Site, data = Field_group)
anova(AMF_field_SR_mod)
shapiro.test(residuals(AMF_field_SR_mod))
AMF_SR_mod_anova <- as.data.frame(anova(AMF_field_SR_mod))
AMF_SR_mod_anova$p_adj <- round(p.adjust(AMF_SR_mod_anova$`Pr(>F)`, method = "BH"), 3) 
AMF_SR_mod_anova[, 2:4] <- round(AMF_SR_mod_anova[, 2:4], digits = 2)
AMF_SR_mod_anova$`Pr(>F)` <- round(AMF_SR_mod_anova$`Pr(>F)`, 3)
AMF_SR_mod_anova <- AMF_SR_mod_anova[c(1:7,9,10,8),] # Reorder
print(AMF_SR_mod_anova) 


############################# Fungal composition ###############################
AMF_field_hel_no_add <- AMF_field_hel_no
AMF_field_hel_no_add[rowSums(AMF_field_hel_no_add) == 0, ] <- 1e-6


AMF_Bray_dist_field_no <- vegdist(AMF_field_hel_no, method = 'bray')

set.seed(1234)
AMF_permanova_mod_field <- GUniFrac::adonis3(AMF_field_hel_no_add ~ Year + Site + Origin/Species + Year:Site + 
                                    Origin:Year + Origin:Site + Origin:Year:Site + 
                                    (Origin/Species):Year + (Origin/Species):Site, method = "bray", by = "margin",
                                  data = Field_group_scale, permutations = 9999)

AMF_permanova_mod_sum_field <- as.data.frame(AMF_permanova_mod_field$aov.tab)
AMF_permanova_mod_sum$R2 <- round(AMF_permanova_mod_sum$R2,3)
AMF_permanova_mod_sum$df <- paste0(AMF_permanova_mod_sum$Df, ",", AMF_permanova_mod_sum$Df[length(AMF_permanova_mod_sum$Df)-1])
AMF_permanova_mod_sum$`F` <- round(AMF_permanova_mod_sum$F.Model,2)
AMF_permanova_mod_sum$`p.adj` <- round(p.adjust(AMF_permanova_mod_sum$`Pr(>F)`, method = "BH"), 3)
print(AMF_permanova_mod_sum)



###########################  Greenhouse experiment #############################
# loading database
green_otu <- read.xlsx("Greenhouse_fungi_Flattening.xlsx", sheet = "green_flattening", colNames = T, rowNames = T)
#green_otu[1:6,1:6]
green_otu <- green_otu[,-c(1)]

## load group data
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)
#colnames(Green_group)
Green_group <- Green_group[colnames(green_otu), ]
#rownames(t(green_otu)) %in% rownames(Green_group)
Green_hel_no <- t(green_otu)/9690

## Tax INFORMATION
tax_default <- read.xlsx("Greenhouse_fungi_Flattening.xlsx", sheet = "ASV_tax", rowNames = F, colNames = T)[,1:3]
#head(tax_default)
rownames(tax_default)  <- tax_default$ASV_name
tax_default <- tax_default[colnames(Green_hel_no),]

OTU_tax <- tax_default %>% tidyr::separate(col = taxonomy_all, 
                                          into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                          sep = ";\\s*", fill = "right")
head(OTU_tax[,1:6])
OTU_tax$GENUS <- sub("g__", "", OTU_tax$Genus)

# Loading FungalTraits database
FungalTraits <- read.xlsx("FungalTraits.xlsx", sheet = "FungalTraits", rowNames = F, colNames = T)
head(FungalTraits[,c(2:6,8)])
OTU_tax2 <- OTU_tax %>% left_join(FungalTraits[,c("GENUS","primary_lifestyle")], by = "GENUS")
head(OTU_tax2)

####  |Plant Pathogen|
Pathogen_id1 <- subset(OTU_tax2, primary_lifestyle == "plant_pathogen")$ASV_name
#### Fusarium
Pathogen_id2 <- subset(OTU_tax, Genus == "g__Fusarium")$ASV_name
### |Plant Pathogen| & Fusarium
Pathogens <- as.data.frame(unique(c(Pathogen_id1, Pathogen_id2)))
Pathogens$guild <- "Plant Pathogen"; colnames(Pathogens)[1] <- "ASV_name"

############################ Arbuscular Mycorrhizal ############################
####  Arbuscular Mycorrhizal
AMF_id1 <- subset(OTU_tax2, primary_lifestyle == "arbuscular_mycorrhizal")$ASV_name; length(AMF_id1)
#### Glomeromycota
AMF_id2 <- subset(OTU_tax, Phylum == "p__Glomeromycota")$ASV_name; length(AMF_id2)
### Arbuscular Mycorrhizal & Glomeromycota
AMF <- as.data.frame(unique(c(AMF_id1, AMF_id2)))
AMF$guild <- "Arbuscular Mycorrhizal"; colnames(AMF)[1] <- "ASV_name"

################################################################################
#### Plant Pathogen
Path_Green_hel_no <- as.data.frame(Green_hel_no)[ ,Pathogens$ASV_name]

## Richness
Green_path_SR <- as.data.frame(specnumber(Path_Green_hel_no))
colnames(Green_path_SR) <- "Path_SR"
Green_path_SR$Sample_ID <- rownames(Green_path_SR)
Green_group <- Green_group %>% left_join(Green_path_SR)
rownames(Green_group) <- Green_group$Sample_ID

Path_Green_SR_mod <- lm(Path_SR ~ Origin/Species, data = Green_group)
anova(Path_Green_SR_mod)
shapiro.test(residuals(Path_Green_SR_mod))
Path_SR_mod_anova <- as.data.frame(anova(Path_Green_SR_mod))
Path_SR_mod_anova$p_adj <- round(p.adjust(Path_SR_mod_anova$`Pr(>F)`, method = "BH"), 3) 
Path_SR_mod_anova[, 2:4] <- round(Path_SR_mod_anova[, 2:4], digits = 2)
Path_SR_mod_anova$`Pr(>F)` <- round(Path_SR_mod_anova$`Pr(>F)`, 3)
print(Path_SR_mod_anova) 


############################# Fungal composition ###############################
Path_Bray_dist_Green_no <- vegdist(Path_Green_hel_no, method = 'bray')

set.seed(1234)
Path_permanova_mod_Green <- GUniFrac::adonis3(Path_Green_hel_no ~ Origin/Species, method = "bray", by = "margin",
                                             data = Green_group, permutations = 9999)

Path_permanova_mod_sum_Green <- as.data.frame(Path_permanova_mod_Green$aov.tab)
Path_permanova_mod_sum_Green$R2 <- round(Path_permanova_mod_sum_Green$R2,3)
Path_permanova_mod_sum_Green$df <- paste0(Path_permanova_mod_sum_Green$Df, ",", Path_permanova_mod_sum_Green$Df[length(Path_permanova_mod_sum_Green$Df)-1])
Path_permanova_mod_sum_Green$`F` <- round(Path_permanova_mod_sum_Green$F.Model,2)
Path_permanova_mod_sum_Green$`p.adj` <- round(p.adjust(Path_permanova_mod_sum_Green$`Pr(>F)`, method = "BH"), 3)
print(Path_permanova_mod_sum_Green)


############################ Arbuscular Mycorrhizal ############################
AMF_Green_hel_no <- as.data.frame(Green_hel_no)[ ,AMF$ASV_name]

## Richness
Green_AMF_SR <- as.data.frame(specnumber(AMF_Green_hel_no))
colnames(Green_AMF_SR) <- "AMF_SR"
Green_AMF_SR$Sample_ID <- rownames(Green_AMF_SR)
Green_group <- Green_group %>% left_join(Green_AMF_SR)
rownames(Green_group) <- Green_group$Sample_ID

AMF_Green_SR_mod_Green <- lm(AMF_SR ~ Origin/Species, data = Green_group)
anova(AMF_Green_SR_mod_Green)
shapiro.test(residuals(AMF_Green_SR_mod_Green))
AMF_SR_mod_anova_Green <- as.data.frame(anova(AMF_Green_SR_mod_Green))
AMF_SR_mod_anova_Green$p_adj <- round(p.adjust(AMF_SR_mod_anova_Green$`Pr(>F)`, method = "BH"), 5) 
AMF_SR_mod_anova_Green[, 2:4] <- round(AMF_SR_mod_anova_Green[, 2:4], digits = 2)
AMF_SR_mod_anova_Green$`Pr(>F)` <- round(AMF_SR_mod_anova_Green$`Pr(>F)`, 3)
print(AMF_SR_mod_anova_Green) 


############################# Fungal composition ###############################
AMF_Green_hel_no_add  <- AMF_Green_hel_no
AMF_Green_hel_no_add[rowSums(AMF_Green_hel_no_add) == 0, ] <- 1e-6
set.seed(1234)
AMF_permanova_mod_Green <- GUniFrac::adonis3(AMF_Green_hel_no_add ~ Origin/Species, method = "bray", by = "margin",
                                            data = Green_group, permutations = 9999)

AMF_permanova_mod_sum_Green <- as.data.frame(AMF_permanova_mod_Green$aov.tab)
AMF_permanova_mod_sum_Green$R2 <- round(AMF_permanova_mod_sum_Green$R2,3)
AMF_permanova_mod_sum_Green$df <- paste0(AMF_permanova_mod_sum_Green$Df, ",", AMF_permanova_mod_sum_Green$Df[length(AMF_permanova_mod_sum_Green$Df)-1])
AMF_permanova_mod_sum_Green$`F` <- round(AMF_permanova_mod_sum_Green$F.Model,2)
AMF_permanova_mod_sum_Green$`p.adj` <- round(p.adjust(AMF_permanova_mod_sum_Green$`Pr(>F)`, method = "BH"), 3)
print(AMF_permanova_mod_sum_Green)


