################################################################################
################################## Fig.3 #######################################
################################################################################
library(ggplot2)
library(openxlsx)
library(adespatial)
library(vegan)
library(phytools)
library(picante)
library(geiger)
library(GUniFrac)
library(patchwork)
library(car)
library(dplyr)
library(patchwork)
library(forcats)
library(ggeffects)
library(MuMIn)

# Custom style
mytheme = theme(
  legend.position = "none",
  panel.grid = element_blank(), 
  strip.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"),
  strip.text.x = element_text(size = 9, color = "black"), 
  legend.title = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 10),
  legend.background = element_rect(fill = NA),
  axis.ticks = element_line(color= "black"),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(colour = "black", size = 14),
  axis.title.y = element_text(colour = "black", size = 14),
  axis.text = element_text(colour = "black", size = 12),
  plot.background = element_blank(), 
  plot.tag = element_text(size = 16, face = "bold")) 

# Loading the grouping metadata of soil samples
Field_group = read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

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

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
field_hel_no <- t(fungi_Flattening)/9477
Bray_dist_field_no <- vegdist(field_hel_no, method = 'bray')


field_hel <- as.data.frame(decostand(t(fungi_Flattening), method = 'hellinger'))
Bray_dist_field <- vegdist(field_hel, method = 'bray')

################################################################################
###########################  Greenhouse experiment #############################
green_otu <- read.xlsx("Greenhouse_fungi_Flattening.xlsx", sheet = "green_flattening", colNames = T, rowNames = T)
green_otu[1:6,1:6]
green_otu <- green_otu[,-c(1)]

## load group data
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)
colnames(Green_group)
Green_group <- Green_group[colnames(green_otu), ]
rownames(t(green_otu)) %in% rownames(Green_group)
green_hel_no <- t(green_otu)/9690
rowSums(green_hel_no)
Bray_dist_green_no <- vegdist(green_hel_no, method = 'bray')

################################################################################
Green_dist_data <- reshape2::melt(as.matrix(Bray_dist_green_no), varnames = c("Sample_ID_A", "Sample_ID_B"),
                                  value.name = "dist", na.rm = T)
colnames(Green_dist_data)[1] = "Sample_ID"
### 
Green_dist_data <- Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
colnames(Green_dist_data)[c(1,2)] <- c("Sample_ID2","Sample_ID")
Green_dist_data <- Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
Green_dist_re <- Rmisc::summarySE(Green_dist_data, measurevar = c("dist"), groupvars = c("Species.x", "Species.y"))
# 
Bray_dist_green_mean <- reshape2::dcast(Green_dist_re, Species.x  ~ Species.y , value.var = "dist")
rownames(Bray_dist_green_mean) <- Bray_dist_green_mean$Species.x
Bray_dist_green_mean <- Bray_dist_green_mean[,-1]
diag(Bray_dist_green_mean) <- 0

############################## dispersion difference ###########################
rownames(as.matrix(Bray_dist_field_no)) %in% Field_group$Sample_ID
bdisper_field <- vegan::betadisper(Bray_dist_field_no, Field_group$Origin, type = "centroid", bias.adjust = FALSE)
## Test significance of dispersion difference
anova(bdisper_field)
TukeyHSD(bdisper_field)

Field_aov_result_all <- NULL
Field_beita_data_all <- NULL
Green_aov_result_all <- NULL
Green_beita_data_all <- NULL

Year = unique(Field_group$Year)
Site = unique(Field_group$Site)

for(i in Year){
  for(ii in Site) {
    ## Field survey
    select_group = subset(Field_group, Year == i & Site == ii)
    field_select_dist = as.matrix(Bray_dist_field_no)[select_group$Sample_ID, select_group$Sample_ID]
    select_group = select_group[rownames(field_select_dist),]
    bdisper_field <- betadisper(as.dist(field_select_dist), select_group$Origin, type = "centroid", bias.adjust = FALSE)
    ## Test significance of dispersion difference
    field_aov_mod = anova(bdisper_field)
    # TukeyHSD(bdisper_field)
    Field_aov_result = data.frame(Type = "Field", Year = i, Site = ii, F_value = field_aov_mod$`F value`[1], p_value = field_aov_mod$`Pr(>F)`[1])
    ##
    Field_beita_data = as.data.frame(bdisper_field$distances) 
    Field_beita_data$Sample_ID = rownames(Field_beita_data)
    colnames(Field_beita_data)[1] = "Field_beita"
    Field_beita_data = Field_beita_data %>% left_join(select_group[,c("Sample_ID", "Species")])
    Field_beita_data = Field_beita_data[,c("Sample_ID","Species","Field_beita")]
    ## Merge data sets
    Field_aov_result_all = rbind(Field_aov_result_all, Field_aov_result)
    Field_beita_data_all = rbind(Field_beita_data_all, Field_beita_data)
    ############################### Greenhouse exp.
    select_sp <- select_group$Species
    select_Origin <- select_group$Origin
    green_select_dist <- as.matrix(Bray_dist_green_mean)[select_sp, select_sp]
    bdisper_green <- betadisper(as.dist(green_select_dist), select_Origin, type = "centroid", bias.adjust = FALSE)
    ## Test significance of dispersion difference
    green_aov_mod <- anova(bdisper_green)
    Green_aov_result <- data.frame(Type = "Greenhouse", Years = i, Site = ii, F_value = green_aov_mod$`F value`[1], p_value = green_aov_mod$`Pr(>F)`[1])
    ##
    Green_beita_data <- as.data.frame(bdisper_green$distances) 
    Green_beita_data$Species <- rownames(Green_beita_data)
    colnames(Green_beita_data)[1] <- "Green_beita"
    Green_beita_data <- Green_beita_data %>% left_join(select_group[,c("Sample_ID", "Species")])
    Green_beita_data <- Green_beita_data[,c("Sample_ID","Species","Green_beita")]
    ## Merge data sets
    Green_aov_result_all <- rbind(Green_aov_result_all, Green_aov_result)
    Green_beita_data_all <- rbind(Green_beita_data_all, Green_beita_data)
  }
}

############################### model selection ################################
beita_data_Centroid <- merge(Field_beita_data_all[,c(1,3)], Field_group_scale, by = 'Sample_ID', all.x = TRUE)
pd_attributes_variable2 <- attributes(scale(beita_data_Centroid[c("Field_beita")]))
colnames(beita_data_Centroid)

multi_model = lm(scale(Field_beita) ~ Site_pool + Soil_ph + Wcont + Soil_N + Tave + Prec + Funct_pPC1 + Funct_pPC2 + Phylo_vPC1 + Phylo_vPC2 +   
                   Funct_pPC1:Soil_ph + Funct_pPC1:Wcont + Funct_pPC1:Soil_N + Funct_pPC1:Tave + Funct_pPC1:Prec + 
                   Funct_pPC2:Soil_ph + Funct_pPC2:Wcont + Funct_pPC2:Soil_N + Funct_pPC2:Tave + Funct_pPC2:Prec + 
                   Phylo_vPC1:Soil_ph + Phylo_vPC1:Wcont + Phylo_vPC1:Soil_N + Phylo_vPC1:Tave + Phylo_vPC1:Prec +
                   Phylo_vPC2:Soil_ph + Phylo_vPC2:Wcont + Phylo_vPC2:Soil_N + Phylo_vPC2:Tave + Phylo_vPC2:Prec
                 ,data = beita_data_Centroid)

shapiro.test(residuals(multi_model))
car::Anova(multi_model)
summary(multi_model)

## Model selection procedure
simplified_model <- stepAICc(multi_model, direction = 'backward')
car::Anova(simplified_model)
performance::r2(simplified_model)

## add plant origin in the full model
multi_model = lm(scale(Field_beita) ~ Site_pool + Origin + Soil_ph + Wcont + Soil_N + Tave + Prec + Funct_pPC1 + Funct_pPC2 + Phylo_vPC1 + Phylo_vPC2 +   
                   Funct_pPC1:Soil_ph + Funct_pPC1:Wcont + Funct_pPC1:Soil_N + Funct_pPC1:Tave + Funct_pPC1:Prec + 
                   Funct_pPC2:Soil_ph + Funct_pPC2:Wcont + Funct_pPC2:Soil_N + Funct_pPC2:Tave + Funct_pPC2:Prec + 
                   Phylo_vPC1:Soil_ph + Phylo_vPC1:Wcont + Phylo_vPC1:Soil_N + Phylo_vPC1:Tave + Phylo_vPC1:Prec +
                   Phylo_vPC2:Soil_ph + Phylo_vPC2:Wcont + Phylo_vPC2:Soil_N + Phylo_vPC2:Tave + Phylo_vPC2:Prec
                 ,data = beita_data_Centroid)

shapiro.test(residuals(multi_model))
car::Anova(multi_model)
summary(multi_model)

## Model selection procedure
simplified_model <- stepAICc(multi_model, direction = 'backward')
car::Anova(simplified_model)
performance::r2(simplified_model) ## Produce the same results as the previous model.

################################################################################

options(na.action = "na.fail")
dd12 <- dredge(simplified_model, subset = ~ Site_pool &
                 dc(Funct_pPC1, Soil_N, Funct_pPC1:Soil_N) & 
                 dc(Funct_pPC1, Tave, Funct_pPC1:Tave) & 
                 ##
                 dc(Phylo_vPC1, Tave, Phylo_vPC1:Tave) & 
                 dc(Phylo_vPC1, Prec, Phylo_vPC1:Prec), trace = 2)

#  ΔAICc < 2 
avg_model <- model.avg(object = dd12, subset = delta < 2)
r_squared <- summary(avg_model)$r.squared

# get the coefficient and the standard error
resultModel <- summary(object = MuMIn::model.avg(object = dd12, 
                                                 subset = delta < 2))$coefmat.subset[c(-1), ]

################################### Table S8 ###################################
Table_S8_results <- as.data.frame(resultModel)
Table_S8_results$`q-vaules` <- p.adjust(Table_S8_results$`Pr(>|z|)`, method = "BH")
Table_S8_results$`q-vaules` <- sprintf("%.3f", Table_S8_results$`q-vaules`)
Table_S8_results$Predictors <- c("FunctPC1", "PhyloPC1", "Prec", "Site pool", "Soil N", "Soil pH", "Tave", 
                                 "FunctPC1 x Soil N", "FunctPC1 x Tave", "PhyloPC1 x Prec", "PhyloPC1 x Tave")
Table_S8_results <- Table_S8_results[c(4,6,5,7,3,1,2,8,9,10,11), c(7,1,2,4,5,6)] # reorder
rownames(Table_S8_results) <- NULL
print(Table_S8_results)

## PLOT
MegaModelSummary <- as.data.frame(resultModel)
MegaModelSummary$Parameter <- rownames(resultModel)
MegaModelSummary$Parameter2 <- c("FunctPC1", "PhyloPC1", "Prec", "Site pool", "Soil N", "Soil pH", "Tave", 
                                 "FunctPC1 × Soil N", "FunctPC1 × Tave", "PhyloPC1 × Prec", "PhyloPC1 × Tave")

MegaModelSummary$Parameter2 <- factor(MegaModelSummary$Parameter2, levels = rev(c("Site pool","Soil pH", "Soil N", "Tave", "Prec", "FunctPC1", "PhyloPC1",
                                                                                 "FunctPC1 × Soil N", "FunctPC1 × Tave", "PhyloPC1 × Prec", "PhyloPC1 × Tave")))
MegaModelSummary$`q-vaules` <- p.adjust(MegaModelSummary$`Pr(>|z|)`, method = "BH")
MegaModelSummary$`q-vaules` <- sprintf("%.3f", MegaModelSummary$`q-vaules`)
#MegaModelSummary$`q-vaules` <- ifelse(MegaModelSummary$`q-vaules` < 0.001, "<0.001", MegaModelSummary$`q-vaules`)

MegaModelSummary$Group <- c("Plant attributes", "Plant attributes", "Climate", "Site pool", "Soil properties", "Soil properties",
                            "Climate", "Interaction", "Interaction", "Interaction", "Interaction")
MegaModelSummary$Group <- factor(MegaModelSummary$Group, levels = c("Site pool", "Soil properties", "Climate", "Plant attributes", "Interaction"))


ggplot(MegaModelSummary, aes(x = Parameter2, y = Estimate, fill = Group))+
  geom_errorbar(aes(ymin=Estimate-1.96*`Std. Error`, ymax=Estimate+1.96*`Std. Error`), width=0, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 11.2), color = "black", linetype = "dashed") + 
  geom_text(aes(y = Estimate+1.96*`Std. Error`, 
                label = paste("italic(p)==", `q-vaules`)),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 3.5) + 
  annotate("text", x = 11.3, y = -0.25,
           label = expression(italic(R)^2 * " = 0.161"),
           colour = "black", size = 4) + 
  labs(x = '', y = 'Standard regression coefficients', color = '', tag = "a") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#45A0A3","#88D8BF", "#F3DBC1","#9796C8","#473C8B")) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x =  element_text(color = "black", size = 14),
        legend.text = element_text(size = 9, color = "black"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',legend.title = element_blank(),legend.key = element_blank(),
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_shape_manual(values = c(16,21)) -> mian_Fig_3a; mian_Fig_3a

MegaModelSummary_var = MegaModelSummary %>%
  mutate(Group=fct_relevel(Group,
                           (c("Site pool","Soil properties",
                              "Climate","Plant attributes",
                              "Interaction")))) %>% 
  dplyr::group_by(Group) %>%
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  mutate(new_col=sum_value/sum(sum_value),
         label = paste0(round(new_col * 100, 1), "%")) 

MegaModelSummary_var$Group = factor(MegaModelSummary_var$Group, levels = c("Site pool", "Soil properties", "Climate", "Plant attributes", "Interaction"))

ggplot()+
  geom_bar(data = MegaModelSummary_var, aes(x = "", y = new_col*100, fill = Group), 
           stat = "identity", width = 0.5, color = "black")+
  geom_text(data = MegaModelSummary_var, 
            aes(x = "", y = new_col*100, label = label, group = Group), 
            position = position_stack(vjust = 0.5), color = "black", size = 4) +
  scale_y_continuous(expand = c(0, 0), position = "right")+
  scale_fill_manual(values = c("#45A0A3","#88D8BF", "#F3DBC1","#9796C8","#473C8B")) +
  theme_classic()+
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linewidth = NA))+
  labs(y = "Relative effect of estimates (%)") -> mian_Fig_3b; mian_Fig_3b


(mian_Fig_3a+mian_Fig_3b) + plot_layout(widths = c(0.7,0.3)) -> mian_Fig_3A; mian_Fig_3A


# interaction: Tave x Funct_pPC1
predmw.m1.rsc2 = ggeffect(simplified_model, terms = c("Funct_pPC1","Tave"))
plot(predmw.m1.rsc2)

eff_mod_data <- data.frame(predmw.m1.rsc2)
colnames(eff_mod_data)[1] = "Funct_pPC1"
colnames(eff_mod_data)[2] = "Field_beita"
colnames(eff_mod_data)[6] = "Tave"

# backtransform pd_attributes_variable2

eff_mod_data["Funct_pPC1"] <- pd_attributes_variable$`scaled:center`["Funct_pPC1"] + pd_attributes_variable$`scaled:scale`["Funct_pPC1"]*eff_mod_data["Funct_pPC1"]
eff_mod_data["Field_beita"] <- pd_attributes_variable2$`scaled:center`["Field_beita"] + pd_attributes_variable2$`scaled:scale`["Field_beita"]*eff_mod_data["Field_beita"]

eff_mod_data$Tave <- ifelse(eff_mod_data$Tave == "-1", "- 1 SD", 
                            ifelse(eff_mod_data$Tave == "0", "Mean", "+ 1 SD"))
eff_mod_data$Tave <- factor(eff_mod_data$Tave, levels = c("- 1 SD", "Mean", "+ 1 SD"))

ggplot()+
  geom_line(data = eff_mod_data, mapping = aes(Funct_pPC1,Field_beita,color=factor(Tave)), size=1.5)+
  theme_bw() + mytheme + 
  guides(color = guide_legend(override.aes = list(fill = NA)))+
  scale_fill_manual(values=rev(c("#324554", "#B89088", "#E9D1CD")),guide=F)+
  scale_color_manual(values=rev(c("#324554", "#B89088", "#E9D1CD")),breaks=waiver(), name="Mean annual temperature")+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  labs(x = "The PC1 of plant traits", 
       y = "Compositional dispersions\n(distances from the centroid)", tag = "b", color = "Mean annual temperature") +
  annotate("segment", x = -0.25, xend = -1.25, y = 0.49, 
           yend = 0.49, colour = "black", size = 0.5, arrow = arrow(angle = 15, 
                                                                    length = unit(0.5,  "cm"))) + 
  annotate("segment", x = 0.5, xend = 1.5, y = 0.49, 
           yend = 0.49, colour = "black", size = 0.5, arrow = arrow(angle = 15, 
                                                                    length = unit(0.5,  "cm"))) +
  theme(legend.position = c(0.65,0.18),
        legend.key = element_blank(),
        legend.title = element_text(size=11),
        legend.text= element_text(size=10),
        legend.background = element_blank())-> mian_Fig_3B; mian_Fig_3B


### 13.20 x 5.39
(mian_Fig_3A|mian_Fig_3B) + plot_layout(widths = c(0.62,0.38))


######################### Structural equation models ###########################
beita_SEM_Centroid = beita_data_Centroid

beita_SEM_Centroid$Phylo_vPC1_soil_N = beita_SEM_Centroid$Phylo_vPC1 * beita_SEM_Centroid$Soil_N
beita_SEM_Centroid$Phylo_vPC1_Soil_ph = beita_SEM_Centroid$Phylo_vPC1 * beita_SEM_Centroid$Soil_ph
beita_SEM_Centroid$Phylo_vPC1_Wcont = beita_SEM_Centroid$Phylo_vPC1 * beita_SEM_Centroid$Wcont
beita_SEM_Centroid$Phylo_vPC1_Tave = beita_SEM_Centroid$Phylo_vPC1 * beita_SEM_Centroid$Tave
beita_SEM_Centroid$Phylo_vPC1_Prec = beita_SEM_Centroid$Phylo_vPC1 * beita_SEM_Centroid$Prec

beita_SEM_Centroid$Funct_pPC1_soil_N = beita_SEM_Centroid$Funct_pPC1 * beita_SEM_Centroid$Soil_N
beita_SEM_Centroid$Funct_pPC1_Soil_ph = beita_SEM_Centroid$Funct_pPC1 * beita_SEM_Centroid$Soil_ph
beita_SEM_Centroid$Funct_pPC1_Wcont = beita_SEM_Centroid$Funct_pPC1 * beita_SEM_Centroid$Wcont
beita_SEM_Centroid$Funct_pPC1_Tave = beita_SEM_Centroid$Funct_pPC1 * beita_SEM_Centroid$Tave
beita_SEM_Centroid$Funct_pPC1_Prec = beita_SEM_Centroid$Funct_pPC1 * beita_SEM_Centroid$Prec

### SEM
library(lavaan)
library(AICcmodavg) 

model_0 =  "
Field_beita ~ Site_pool + Prec + Tave + Soil_N + Soil_ph + Funct_pPC1 + Phylo_vPC1 + 
              Funct_pPC1_soil_N + Funct_pPC1_Tave + 
              Phylo_vPC1_Tave + Phylo_vPC1_Prec
Site_pool ~ Prec + Tave + Soil_N + Soil_ph
Soil_ph ~ Tave + Prec
Soil_N ~ Tave + Prec
"
fit_model_0 <- sem(model_0, data = beita_SEM_Centroid)
show(fit_model_0)
summary(fit_model_0, standardized = TRUE, rsq = TRUE)
# Finding the missing path
mi0<-modindices(fit_model_0);print(mi0[mi0$mi>3.0,])

## detect "Phylo_vPC1 x Prec"
model_1 =  "
Field_beita ~ Site_pool + Prec + Tave + Soil_N + Soil_ph + Funct_pPC1 + Phylo_vPC1 + 
              Funct_pPC1_soil_N + Funct_pPC1_Tave + 
              Phylo_vPC1_Tave
Site_pool ~ Prec + Tave + Soil_N + Soil_ph
Soil_ph ~ Tave + Prec
Soil_N ~ Tave + Prec
"
fit_model_1 <- sem(model_1, data = beita_SEM_Centroid)
show(fit_model_1)
summary(fit_model_1, standardized = TRUE, rsq = TRUE)


## detect "Phylo_vPC1 x Tave"
model_2 =  "
Field_beita ~ Site_pool + Prec + Tave + Soil_N + Soil_ph + Funct_pPC1 + Phylo_vPC1 + 
              Funct_pPC1_soil_N + Funct_pPC1_Tave
Site_pool ~ Prec + Tave + Soil_N + Soil_ph
Soil_ph ~ Tave + Prec
Soil_N ~ Tave + Prec
"
fit_model_2 <- sem(model_2, data = beita_SEM_Centroid)
show(fit_model_2)
summary(fit_model_2, standardized = TRUE, rsq = TRUE)


## detect "Phylo_vPC1"
model_3 =  "
Field_beita ~ Site_pool + Prec + Tave + Soil_N + Soil_ph + Funct_pPC1 + 
              Funct_pPC1_soil_N + Funct_pPC1_Tave
Site_pool ~ Prec + Tave + Soil_N + Soil_ph
Soil_ph ~ Tave + Prec
Soil_N ~ Tave + Prec
"
fit_model_3 <- sem(model_3, data = beita_SEM_Centroid)
show(fit_model_3)
summary(fit_model_3, standardized = TRUE, rsq = TRUE)


##
fitMeasures(fit_model_0, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic", "srmr"))
fitMeasures(fit_model_1, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic", "srmr"))
fitMeasures(fit_model_2, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic", "srmr"))
fitMeasures(fit_model_3, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic", "srmr"))

## fit_model_3 is the best model

## Another model
################### add plant origin -> FunctPC1 & PhyloPC1 ####################
model_0 =  "
Field_beita ~ Site_pool + Prec + Tave + Soil_N + Soil_ph + Funct_pPC1 + Phylo_vPC1 + 
              Funct_pPC1_soil_N + Funct_pPC1_Tave + 
              Phylo_vPC1_Tave + Phylo_vPC1_Prec
Site_pool ~ Prec + Tave + Soil_N + Soil_ph
Soil_ph ~ Tave + Prec
Soil_N ~ Tave + Prec
Funct_pPC1 ~ Origin
Phylo_vPC1 ~ Origin
"
fit_model_0 <- sem(model_0, data = beita_SEM_Centroid)
show(fit_model_0)
summary(fit_model_0, standardized = TRUE, rsq = TRUE)
# Finding the missing path
mi0<-modindices(fit_model_0);print(mi0[mi0$mi>3.0,])

## detect "Phylo_vPC1 x Prec"
model_1 =  "
Field_beita ~ Site_pool + Prec + Tave + Soil_N + Soil_ph + Funct_pPC1 + Phylo_vPC1 + 
              Funct_pPC1_soil_N + Funct_pPC1_Tave + 
              Phylo_vPC1_Tave
Site_pool ~ Prec + Tave + Soil_N + Soil_ph
Soil_ph ~ Tave + Prec
Soil_N ~ Tave + Prec
Funct_pPC1 ~ Origin
Phylo_vPC1 ~ Origin
"
fit_model_1 <- sem(model_1, data = beita_SEM_Centroid)
show(fit_model_1)
summary(fit_model_1, standardized = TRUE, rsq = TRUE)


## detect "Phylo_vPC1 x Tave"
model_2 =  "
Field_beita ~ Site_pool + Prec + Tave + Soil_N + Soil_ph + Funct_pPC1 + Phylo_vPC1 + 
              Funct_pPC1_soil_N + Funct_pPC1_Tave
Site_pool ~ Prec + Tave + Soil_N + Soil_ph
Soil_ph ~ Tave + Prec
Soil_N ~ Tave + Prec
Funct_pPC1 ~ Origin
Phylo_vPC1 ~ Origin
"
fit_model_2 <- sem(model_2, data = beita_SEM_Centroid)
show(fit_model_2)
summary(fit_model_2, standardized = TRUE, rsq = TRUE)

##
fitMeasures(fit_model_0, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic", "srmr"))
fitMeasures(fit_model_1, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic", "srmr"))
fitMeasures(fit_model_2, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic", "srmr"))


## Funct_pPC1 ~ Origin
summary(aov(Funct_pPC1 ~ Origin, data = beita_SEM_Centroid))
boxplot(Funct_pPC1 ~ Origin, data = beita_SEM_Centroid, 
        main = "Funct_pPC1 by Origin", ylab = "Funct_pPC1", xlab = "Origin")

stripchart(Funct_pPC1 ~ Origin, data = beita_SEM_Centroid,
           vertical = TRUE, method = "jitter", pch = 21, col = "blue", bg = "lightblue", add = TRUE)

group_means <- tapply(beita_SEM_Centroid$Funct_pPC1, beita_SEM_Centroid$Origin, mean)
points(1:length(group_means), group_means, col = "red", pch = 19, cex = 1.2)

## 均值填充相同数据后导致植物来源对PhyloPC1产生显著影响！！！
## Phylo_vPC1 ~ Origin
summary(aov(Phylo_vPC1 ~ Origin, data = beita_SEM_Centroid))
boxplot(Phylo_vPC1 ~ Origin, data = beita_SEM_Centroid, 
        main = "Phylo_vPC1 by Origin", ylab = "Phylo_vPC1", xlab = "Origin")

stripchart(Phylo_vPC1 ~ Origin, data = beita_SEM_Centroid,
           vertical = TRUE, method = "jitter", pch = 21, col = "blue", bg = "lightblue", add = TRUE)

group_means <- tapply(beita_SEM_Centroid$Phylo_vPC1, beita_SEM_Centroid$Origin, mean)
points(1:length(group_means), group_means, col = "red", pch = 19, cex = 1.2)
