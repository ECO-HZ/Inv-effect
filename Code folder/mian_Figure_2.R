################################################################################
################################## Fig 2 #######################################
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
library(ggplotify)
library(aplot)
library(deming)
library(ggtext)

# Custom style
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
  axis.title.x = element_text(colour='black', size=14),
  axis.title.y = element_text(colour='black', size=14),
  axis.text = element_text(colour='black',size=12),
  plot.title = element_textbox(
    size = 14, color = "black", fill = "grey90",
    box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
    halign = 0.5, width = grid::unit(1, "npc")), #r = unit(3, "pt")     
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
############################### Fig.2a NMDS ####################################
set.seed(1234)
field_nmds1 <- metaMDS(Bray_dist_field_no, k = 3, trymax = 100, autotransform = TRUE, wascores = TRUE)
#field_nmds1 <- metaMDS(Bray_dist_field, k = 2, trymax = 100, autotransform = TRUE, wascores = TRUE)
field_nmds1.stress <- field_nmds1$stress
field_plot_data <- data.frame(field_nmds1$point)
field_plot_data$Sample_ID <- rownames(field_plot_data)
names(field_plot_data)[1:3] <- c('NMDS1', 'NMDS2', 'NMDS3')
field_plot_data <- merge(field_plot_data, Field_group, by = 'Sample_ID', all.x = TRUE)
field_plot_data$Origin <- factor(field_plot_data$Origin, levels = c("Native", "Exotic"))
field_plot_data$Site <- factor(field_plot_data$Site, levels = c("Guangzhou", "Guilin", "Changsha", "Wuhan", "Zhengzhou", "Tai'an"))

ggplot(field_plot_data, aes(NMDS1, NMDS2, shape = Origin, fill = Year))+
  geom_point(size = 2.2, color = "black", alpha = 1, show.legend = T) +
  scale_shape_manual(values = c(21,22)) +
  scale_color_manual(values = c("#898EA1","#CF9742","#3A7C72")) + 
  scale_fill_manual(values = c("#898EA1","#CF9742","#3A7C72")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size=12))+
  theme_bw() + mytheme +
  theme(legend.position = "right") + 
  labs(x = "NMDS1", y = "NMDS2", tag = "a")+
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.01)) -> mian_Fig_2_top; mian_Fig_2_top

ggplot(field_plot_data, aes(NMDS1, NMDS2, shape = Origin, fill = Site))+
  geom_point(size = 2.2, color = "black",alpha = 1, show.legend = T) +
  scale_shape_manual(values = c(21,22)) +
  scale_color_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  scale_fill_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(text = element_text(size = 12))+
  theme_bw() + mytheme +
  theme(legend.position = "right") + 
  labs(x = "NMDS1", y = "NMDS2")+
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.01)) -> mian_Fig_2_bottom; mian_Fig_2_bottom

# Combined
(mian_Fig_2_top/mian_Fig_2_bottom) -> mian_Fig_2A; mian_Fig_2A

################################################################################
############################# Fig.2b PERMAONVA #################################
## Bray-Curtis distance matrix, 9999 Freedman-Lane permutations
set.seed(1234)
permanova_mod <- GUniFrac::adonis3(field_hel_no ~ Year + Site + Origin/Species + Year:Site + 
                                    Origin:Year + Origin:Site + Origin:Year:Site + 
                                    (Origin/Species):Year + (Origin/Species):Site, method = "bray", by = "margin",
                                  data = Field_group_scale, permutations = 9999)

permanova_mod_sum <- as.data.frame(permanova_mod$aov.tab)
permanova_mod_sum$R2 <- round(permanova_mod_sum$R2,3)
permanova_mod_sum$df <- paste0(permanova_mod_sum$Df, ",", permanova_mod_sum$Df[length(permanova_mod_sum$Df)-1])
permanova_mod_sum$`F` <- round(permanova_mod_sum$F.Model,2)
permanova_mod_sum$`p.adj` <- round(p.adjust(permanova_mod_sum$`Pr(>F)`, method = "BH"), 3)
permanova_mod_sum$Predictors <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                  "Year × Site × Origin", "Year × Species", "Site × Species","Residuals","Total")
permanova_mod_sum <- permanova_mod_sum[c(1:7,9,10,8), c("Predictors","R2","F","Pr(>F)","p.adj")] # Reorder
rownames(permanova_mod_sum) <- NULL
print(permanova_mod_sum)

# plot
permanova_data1 <- as.data.frame(permanova_mod$aov.tab)[1:10, ]
permanova_data1$Label <- c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                          "Year × Site × Origin", "Year × Species", "Site × Species")
permanova_data1$Label <- factor(permanova_data1$Label, levels = rev(c("Year", "Site", "Origin", "Species", "Year × Site", "Year × Origin", "Site × Origin", 
                                                                     "Year × Species", "Site × Species", "Year × Site × Origin")))
permanova_data1$R2_per <- sprintf("%.1f", permanova_data1$R2*100)
permanova_data1$R2 <- round(permanova_data1$R2,2)
permanova_data1$`q-vaules` <- p.adjust(permanova_data1$`Pr(>F)`, method = "BH")
permanova_data1$p_label <- ifelse(
  permanova_data1$`q-vaules` < 0.001,
  "italic(p) < 0.001",
  paste0("italic(p) == ", formatC(permanova_data1$`q-vaules`, format = "f", digits = 3, flag = "0")))
permanova_data1$p_shape <- as.factor(ifelse(permanova_data1$`q-vaules` >= 0.05, 1, 0))

ggplot(permanova_data1, aes(x = Label, y = as.numeric(R2_per))) + 
  geom_segment(aes(x = Label, xend = Label, y = 0, yend = as.numeric(R2_per)-0.015),                 
               linetype = "solid", size = 3, color = "#D3D3D3") + 
  geom_point(aes(fill = p_shape), color = "black", size = 10, pch = 22) + 
  geom_text(aes(label = R2_per), color = "black", size = 3) + 
  geom_text(aes(label = p_label),            
            hjust = -0.5,
            vjust = 0.5, angle = 0, color = "black", size = 3.5, parse = TRUE) + 
  labs(x = NULL, y = "Explained variance (%)", tag = "b") + 
  scale_fill_manual(values = c("#08ADB5", "white")) + 
  theme_bw() + mytheme +
  theme(axis.text.x = element_text(size = 11, angle = 0),
        plot.margin = margin(l = 60, r = 20, t = 10, b = 10)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), limits = c(0, 30)) +
  coord_flip() -> mian_Fig_2B; mian_Fig_2B


################################################################################
##############Fig.2c pair-wise Bray-Curtis dissimilarities #####################

# The pair-wise Bray-Curtis dissimilarities of fungal composition between exotic 
# plants and native plants was compared per site per year.
Field_aov_result_BC = NULL
Field_BC_data_all = NULL

Year = unique(Field_group$Year)
Site = unique(Field_group$Site)

for(i in Year){
  for(ii in Site) {
    ## Field survey
    select_group <- subset(Field_group, Year == i & Site == ii)
    Native_sample <- subset(select_group, Origin == "Native")$Sample_ID
    Exotic_sample <- subset(select_group, Origin == "Exotic")$Sample_ID
    ## Native
    Native_matrix <- as.matrix(Bray_dist_field_no)[Native_sample, Native_sample]
    dim(Native_matrix)
    lower_tri_mask <- lower.tri(Native_matrix, diag = FALSE)
    Native_matrix_lower <- Native_matrix
    Native_matrix_lower[!lower_tri_mask] <- NA
    Native_BC_long <- reshape2::melt(Native_matrix_lower, na.rm = TRUE, 
                                     varnames = c("Sample1", "Sample2"), value.name = "Bray_Distance")
    Native_BC_long$Origin <- "Native"
    ## Exotic
    Exotic_matrix <- as.matrix(Bray_dist_field_no)[Exotic_sample, Exotic_sample]
    dim(Exotic_matrix)
    lower_tri_mask <- lower.tri(Exotic_matrix, diag = FALSE)
    Exotic_matrix_lower <- Exotic_matrix
    Exotic_matrix_lower[!lower_tri_mask] <- NA
    Exotic_BC_long <- reshape2::melt(Exotic_matrix_lower, na.rm = TRUE, 
                                     varnames = c("Sample1", "Sample2"), value.name = "Bray_Distance")
    Exotic_BC_long$Origin <- "Exotic"
    ##
    Pairwise_BC_data <- rbind(Native_BC_long, Exotic_BC_long)
    Pairwise_BC_data$Year <- i; Pairwise_BC_data$Site <- ii
    ANOVAs_mod = aov(Bray_Distance ~ Origin, data = Pairwise_BC_data)
    aov_summary <- summary(ANOVAs_mod)[[1]]
    Field_aov_result <- data.frame(Type = "Field", Year = i, Site = ii,
                                   F_value = round(aov_summary$'F value'[1], 2),
                                   p_value = round(aov_summary$'Pr(>F)'[1], 3))
    ## Merge data sets
    Field_aov_result_BC = rbind(Field_aov_result_BC, Field_aov_result)
    Field_BC_data_all = rbind(Field_BC_data_all, Pairwise_BC_data)
  }
}

head(Field_BC_data_all)
Field_BC_data_all$Year <- as.factor(Field_BC_data_all$Year)

Field_pair_BC_sum <- Field_BC_data_all %>% group_by(Year,Site,Origin) %>%
  summarise(mean_BC = mean(Bray_Distance), sd_BC = sd(Bray_Distance),      
            n = n(), se_BC = sd_BC / sqrt(n), .groups = "drop")

Field_pair_BC_sum$Site = factor(Field_pair_BC_sum$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Field_pair_BC_sum$Origin = factor(Field_pair_BC_sum$Origin, levels = c("Native","Exotic"))

# Create a data frame with a gray background
background_data <- data.frame(
  Site = c("Guangzhou", "Changsha", "Zhengzhou"), 
  xmin = c(0.5, 2.5, 4.5),         
  xmax = c(1.5, 3.5, 5.5))

print(subset(Field_aov_result_BC, p_value <= 0.05))

ggplot(Field_pair_BC_sum, aes(x = Site, y = mean_BC, fill = Year, color = Year, shape = Origin)) + 
  geom_rect(data = background_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), 
            fill = "gray95", alpha = 0.5, inherit.aes = FALSE) + 
  geom_point(aes(fill = Year),position = position_dodge(0.9), color = "black", size = 3) +
  scale_color_manual(values = c("#898EA1","#CF9742","#3A7C72")) + 
  scale_fill_manual(values = c("#898EA1","#CF9742","#3A7C72")) +
  geom_errorbar(aes(ymin = mean_BC - se_BC, ymax = mean_BC + se_BC, color = Year), 
                width = 0, linewidth = 0.4, position = position_dodge(width = 0.9)) +
  geom_segment(aes(x = 0.9, xend = 1.1, y = 0.87, yend = 0.87), color = "black") + 
  annotate("text", x = 1, y = 0.88, label = expression(italic("p") ~ "= 0.044"), size = 3.5) +
  geom_segment(aes(x = 1.9, xend = 2.1, y = 0.87, yend = 0.87), color = "black") + 
  annotate("text", x = 2, y = 0.88, label = expression(italic("p") ~ "< 0.001"), size = 3.5) +
  geom_segment(aes(x = 2.6, xend = 2.8, y = 0.87, yend = 0.87), color = "black") + 
  annotate("text", x = 2.7, y = 0.88, label = expression(italic("p") ~ "= 0.001"), size = 3.5) +
  geom_segment(aes(x = 2.9, xend = 3.1, y = 0.89, yend = 0.89), color = "black") + 
  annotate("text", x = 3.0, y = 0.90, label = expression(italic("p") ~ "< 0.001"), size = 3.5) +
  geom_segment(aes(x = 3.2, xend = 3.4, y = 0.87, yend = 0.87), color = "black") + 
  annotate("text", x = 3.3, y = 0.88, label = expression(italic("p") ~ "< 0.001"), size = 3.5) +
  geom_segment(aes(x = 3.6, xend = 3.8, y = 0.89, yend = 0.89), color = "black") + 
  annotate("text", x = 3.7, y = 0.90, label = expression(italic("p") ~ "< 0.001"), size = 3.5) +
  geom_segment(aes(x = 4.9, xend = 5.1, y = 0.89, yend = 0.89), color = "black") + 
  annotate("text", x = 5.0, y = 0.90, label = expression(italic("p") ~ "= 0.023"), size = 3.5) +
  geom_segment(aes(x = 6.2, xend = 6.4, y = 0.87, yend = 0.87), color = "black") + 
  annotate("text", x = 6.3, y = 0.88, label = expression(italic("p") ~ "< 0.001"), size = 3.5) +
  scale_shape_manual(values = c(21, 22)) +
  theme_bw() + mytheme + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1)) +   
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.15))) + 
  labs(x = '', y = '') -> mian_Fig_2c; mian_Fig_2c


### Total database
Native_sample <- subset(Field_group, Origin == "Native")$Sample_ID
Exotic_sample <- subset(Field_group, Origin == "Exotic")$Sample_ID

## Native
Native_matrix <- as.matrix(Bray_dist_field_no)[Native_sample, Native_sample]
dim(Native_matrix)
lower_tri_mask <- lower.tri(Native_matrix, diag = FALSE)
Native_matrix_lower <- Native_matrix
Native_matrix_lower[!lower_tri_mask] <- NA
dim(Native_matrix_lower)
Native_BC_long <- reshape2::melt(Native_matrix_lower, na.rm = TRUE, 
                                 varnames = c("Sample1", "Sample2"), value.name = "Bray_Distance")
head(Native_BC_long)
Native_BC_long$Origin = "Native"

## Exotic
Exotic_matrix = as.matrix(Bray_dist_field_no)[Exotic_sample, Exotic_sample]
dim(Exotic_matrix)
lower_tri_mask <- lower.tri(Exotic_matrix, diag = FALSE)
Exotic_matrix_lower <- Exotic_matrix
Exotic_matrix_lower[!lower_tri_mask] <- NA
Exotic_BC_long <- reshape2::melt(Exotic_matrix_lower, na.rm = TRUE, 
                                 varnames = c("Sample1", "Sample2"), value.name = "Bray_Distance")
Exotic_BC_long$Origin = "Exotic"
head(Exotic_BC_long)

## One-way ANOVAs
Pairwise_BC_data_total <- rbind(Native_BC_long, Exotic_BC_long)
mod <- aov(Bray_Distance ~ Origin, data = Pairwise_BC_data_total)
summary(mod)

data_all <- Rmisc::summarySE(Pairwise_BC_data_total, groupvars = "Origin", measurevar = "Bray_Distance")
data_all$Origin <- factor(data_all$Origin, levels = c("Native","Exotic"))

ggplot(data_all, aes(x = Origin, y = Bray_Distance, shape = Origin)) + 
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = Bray_Distance - se, ymax = Bray_Distance + se, color = Origin), 
                width = 0, linewidth = 0.4, color = "black") +
  geom_segment(aes(x = 1, xend = 2, y = 0.89, yend = 0.89)) + 
  annotate("text", x = 1.5, y = 0.90, label = expression(italic("p") ~ "< 0.001"), size = 3.5) +
  scale_shape_manual(values = c(21, 22)) +
  theme_bw() + mytheme + 
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),
        legend.position = "none") + 
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.01), 
                     limits = c(ggplot_build(mian_Fig_2c)$layout$panel_scales_y[[1]]$range$range)) + 
  labs(x = '', y = 'Mean pair-wise\nBray–Curtis dissimilarities', 
       tag = "c") -> mian_Fig_2cc; mian_Fig_2cc


mian_Fig_2C <- mian_Fig_2c %>% insert_left(mian_Fig_2cc,width=0.2) %>% as.ggplot()
mian_Fig_2C # 11.25 x 4.20


################################################################################
################################## Fig.2d ######################################
## Relationships between the mean pairwise fungal compositional dissimilarity 
## estimated in the greenhouse experiment and those estimated in the field 
## survey for the same groups of plant species cooccurring at the same site 
## and time in the field

###########################  Greenhouse experiment #############################
# Loading greenhouse experimental data
green_otu <- read.xlsx("Greenhouse_fungi_Flattening.xlsx", sheet = "green_flattening", colNames = T, rowNames = T)
green_otu[1:6,1:6]
green_otu <- green_otu[,-c(1)]

## load group data
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "sample_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)
colnames(Green_group)
Green_group <- Green_group[colnames(green_otu), ]
rownames(t(green_otu)) %in% rownames(Green_group)

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) 
# in Greenhouse experiment
green_hel_no = t(green_otu)/9690
rowSums(green_hel_no)
Bray_dist_green_no <- vegdist(green_hel_no, method = 'bray')


# note:
# To match the rhizosphere fungal data with the corresponding species from the 
# field survey, we calculated the mean pairwise dissimilarity values among 
# species (averaged across three replicate samples).

Green_dist_data <- reshape2::melt(as.matrix(Bray_dist_green_no), varnames = c("Sample_ID_A", "Sample_ID_B"),
                                  value.name = "dist", na.rm = T)
colnames(Green_dist_data)[1] <- "Sample_ID"
Green_dist_data <- Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
colnames(Green_dist_data)[c(1,2)] <- c("Sample_ID2","Sample_ID")
Green_dist_data <- Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
Green_dist_re <- Rmisc::summarySE(Green_dist_data, measurevar = c("dist"), groupvars = c("Species.x", "Species.y"))
Bray_dist_green_mean <- reshape2::dcast(Green_dist_re, Species.x  ~ Species.y , value.var = "dist")
rownames(Bray_dist_green_mean) <- Bray_dist_green_mean$Species.x
Bray_dist_green_mean <- Bray_dist_green_mean[,-1]
diag(Bray_dist_green_mean) <- 0

#################### Fit a generalized Deming regression #######################
# pisewise - Native & Native
Year <- unique(Field_group_scale$Year)
Site <- unique(Field_group_scale$Site)
BC_aov_data_all <- NULL
β_Bray_Native <- NULL

for(i in Year){
  for(ii in Site) {
    ## Field
    select_group <- subset(Field_group, Year == i & Site == ii)
    ##
    native_sample <- subset(select_group, Origin == "Native")$Sample_ID
    native_latin <- subset(select_group, Origin == "Native")$Species
    ## Field
    select_dist_Field <- as.matrix(Bray_dist_field_no)[native_sample, native_sample]
    lower_tri_mask_Field <- lower.tri(select_dist_Field, diag = FALSE)
    Native_matrix_lower_Field <- select_dist_Field
    Native_matrix_lower_Field[!lower_tri_mask_Field] <- NA
    Native_Field_BC_data <- reshape2::melt(Native_matrix_lower_Field, na.rm = TRUE, 
                                           varnames = c("Sample1", "Sample2"), value.name = "BC_value")
    Native_Field_BC_data$Group <- "Field"
    ##
    Field_mean <- mean(Native_Field_BC_data$BC_value)
    Field_sd <- sd(Native_Field_BC_data$BC_value)
    Field_se <- Field_sd / sqrt(length(Native_Field_BC_data$BC_value))  

    ## Greenhouse 
    select_dist_Green <- as.matrix(Bray_dist_green_mean)[native_latin, native_latin]
    lower_tri_mask_Green <- lower.tri(select_dist_Green, diag = FALSE)
    Native_matrix_lower_Green <- select_dist_Green
    Native_matrix_lower_Green[!lower_tri_mask_Green] <- NA
    Native_Green_BC_data <- reshape2::melt(Native_matrix_lower_Green, na.rm = TRUE, 
                                           varnames = c("Sample1", "Sample2"), value.name = "BC_value")
    Native_Green_BC_data$Group <- "Greenhouse"

    Green_mean <- mean(Native_Green_BC_data$BC_value)
    Green_sd <- sd(Native_Green_BC_data$BC_value)
    Green_se <- Green_sd / sqrt(length(Native_Green_BC_data$BC_value))
    
    ## One-way ANOVAs 
    all_BC_aov_data <- rbind(Native_Field_BC_data, Native_Green_BC_data)
    aov_summary <- summary(aov(BC_value ~ Group, data = all_BC_aov_data))
    BC_aov_data_all <- rbind(BC_aov_data_all, all_BC_aov_data)
    
    # Merged data
    dist_data_all <- data.frame(Year = i , Site = ii, 
                                Field_mean = Field_mean, Field_sd = Field_sd, Field_se = Field_se,
                                Green_mean = Green_mean, Green_sd = Green_sd, Green_se = Green_se,
                                F_value = aov_summary[[1]]$"F value"[1],
                                p_value = round(aov_summary[[1]]$"Pr(>F)"[1], 3))
    β_Bray_Native <- rbind(β_Bray_Native, dist_data_all)
  }
}

head(BC_aov_data_all)
summary(aov(BC_value ~ Group, data = BC_aov_data_all))
t.test(β_Bray_Native$Field_mean, β_Bray_Native$Green_mean)
# View(β_Bray_Native)
summary(lm(β_Bray_Native$Field_mean~β_Bray_Native$Green_mean))


fit <- deming::deming(Field_mean ~ Green_mean, ystd = Field_sd, xstd = Green_sd, data = β_Bray_Native)
print(fit)
β_Bray_Native$Site <- factor(β_Bray_Native$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot()+
  geom_abline(intercept = 0, slope = 1, color = "#8B0000", linetype = 1, size = 1) +
  geom_abline(intercept = -0.1476869, slope = 1.3554311, linetype = 2, size = 1) +
  geom_point(β_Bray_Native, mapping = aes(x = Green_mean,y = Field_mean, shape = Year, fill = Site), color = "black", size = 2.5) +
  geom_errorbar(data = β_Bray_Native, mapping = aes(x = Green_mean, ymax = Field_mean + Field_se, ymin = Field_mean - Field_se, color = Site), width = 0, size = 0.5) +
  geom_errorbarh(data = β_Bray_Native, mapping = aes(y = Field_mean, xmax = Green_mean + Green_se, xmin = Green_mean - Green_se, color = Site), height = 0, size = 0.5) +
  scale_color_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  scale_fill_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  guides(col = guide_legend(ncol = 1))+
  scale_shape_manual(values = c(24,23,25)) +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.55,0.92)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.55,0.92)) + 
  theme_bw() + mytheme + theme(legend.position = "none") + 
  labs(x = NULL, 
       y = "Mean pair-wise Bray–Curtis dissimilarities\nestimated in the field",
       tag = "d", title = c("Between native plants")) -> mian_Fig_4d1; mian_Fig_4d1

## pisewise - Exotic & Exotic
Year <- unique(Field_group_scale$Year)
Site <- unique(Field_group_scale$Site)
β_Bray_Exotic <- NULL
BC_aov_data_all <- NULL

for(i in Year){
  for(ii in Site) {
    ## Field
    select_group <- subset(Field_group, Year == i & Site == ii)
    ##
    exotic_sample <- subset(select_group, Origin == "Exotic")$Sample_ID
    exotic_latin <- subset(select_group, Origin == "Exotic")$Species
    ## Field
    select_dist_Field <- as.matrix(Bray_dist_field_no)[exotic_sample, exotic_sample]

    lower_tri_mask_Field <- lower.tri(select_dist_Field, diag = FALSE)
    Exotic_matrix_lower_Field <- select_dist_Field
    Exotic_matrix_lower_Field[!lower_tri_mask_Field] <- NA
    Exotic_Field_BC_data <- reshape2::melt(Exotic_matrix_lower_Field, na.rm = TRUE, 
                                           varnames = c("Sample1", "Sample2"), value.name = "BC_value")
    Exotic_Field_BC_data$Group = "Field"
    ##
    Field_mean = mean(Exotic_Field_BC_data$BC_value)
    Field_sd = sd(Exotic_Field_BC_data$BC_value)
    Field_se = Field_sd / sqrt(length(Exotic_Field_BC_data$BC_value))

    ## Greenhouse 
    select_dist_Green <- as.matrix(Bray_dist_green_mean)[exotic_latin, exotic_latin]
    lower_tri_mask_Green <- lower.tri(select_dist_Green, diag = FALSE)
    Exotic_matrix_lower_Green <- select_dist_Green
    Exotic_matrix_lower_Green[!lower_tri_mask_Green] <- NA
    Exotic_Green_BC_data <- reshape2::melt(Exotic_matrix_lower_Green, na.rm = TRUE, 
                                           varnames = c("Sample1", "Sample2"), value.name = "BC_value")
    Exotic_Green_BC_data$Group = "Greenhouse"

    Green_mean = mean(Exotic_Green_BC_data$BC_value)
    Green_sd = sd(Exotic_Green_BC_data$BC_value)
    Green_se = Green_sd / sqrt(length(Exotic_Green_BC_data$BC_value)) 
    
    ## anova test
    all_BC_aov_data = rbind(Exotic_Field_BC_data, Exotic_Green_BC_data)
    aov_summary = summary(aov(BC_value ~ Group, data = all_BC_aov_data))
    BC_aov_data_all = rbind(BC_aov_data_all, all_BC_aov_data)
    ## 
    dist_data_all = data.frame(Year = i , Site = ii, 
                               Field_mean = Field_mean, Field_sd = Field_sd, Field_se = Field_se,
                               Green_mean = Green_mean, Green_sd = Green_sd, Green_se = Green_se,
                               F_value = aov_summary[[1]]$"F value"[1],
                               p_value = round(aov_summary[[1]]$"Pr(>F)"[1], 3))
    β_Bray_Exotic = rbind(β_Bray_Exotic, dist_data_all)
  }
}
# View(β_Bray_Exotic)
head(BC_aov_data_all)
summary(aov(BC_value ~ Group, data = BC_aov_data_all))
t.test(β_Bray_Exotic$Field_mean, β_Bray_Exotic$Green_mean)

fit <- deming::deming(Field_mean ~ Green_mean, ystd = Field_sd, xstd = Green_sd, data = β_Bray_Exotic)
print(fit)
β_Bray_Exotic$Site <- factor(β_Bray_Exotic$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot()+
  geom_abline(intercept = 0, slope = 1, color = "#96393D", linetype = 1, size = 1)+
  geom_abline(intercept = -3.337836, slope = 5.926437, linetype = 2, size = 1)+
  geom_point(β_Bray_Exotic, mapping = aes(x=Green_mean,y=Field_mean, shape = Year, fill = Site), color = "black", size=2.5)+
  geom_errorbar(data = β_Bray_Exotic,mapping = aes(x = Green_mean,ymax = Field_mean+Field_se, ymin=Field_mean-Field_se, color = Site),width=0,size=0.5,alpha = 1)+#
  geom_errorbarh(data = β_Bray_Exotic,mapping = aes(y = Field_mean,xmax=Green_mean+Green_se,xmin=Green_mean-Green_se, color = Site),height=0,size=0.5,alpha = 1)+#
  scale_color_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  scale_fill_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  guides(col = guide_legend(ncol = 1))+
  scale_shape_manual(values = c(24,23,25)) +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.55,0.92)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.55,0.92)) + 
  theme_bw() + mytheme + theme(legend.position = "none") + 
  labs(x = "Mean pair-wise Bray–Curtis dissimilarities\nestimated in the greenhouse experiment", 
       y = "", title = c("Between alien plants")) -> mian_Fig_4d2; mian_Fig_4d2


######################### pisewise - Native & Exotic ###########################
Year <- unique(Field_group_scale$Year)
Site <- unique(Field_group_scale$Site)
β_Bray_data <- NULL
BC_aov_data_all <- NULL

for(i in Year){
  for(ii in Site) {
    ## Field
    select_group = subset(Field_group, Year == i & Site == ii)
    ##
    native_sample = subset(select_group, Origin == "Native")$Sample_ID
    native_latin = subset(select_group, Origin == "Native")$Species
    exotic_sample = subset(select_group, Origin == "Exotic")$Sample_ID
    exotic_latin = subset(select_group, Origin == "Exotic")$Species
    ## Field
    select_dist_field = as.matrix(Bray_dist_field_no)[native_sample, exotic_sample]

    Field_mean = mean(select_dist_field)
    Field_sd = sd(select_dist_field)
    Field_se = Field_sd / sqrt(length(select_dist_field)) 
    Field_BC_data = data.frame(BC_value = as.vector(select_dist_field), Group = "Field")
    ## Greenhouse 
    select_dist_green = as.matrix(Bray_dist_green_mean)[native_latin, exotic_latin]

    Green_mean = mean(select_dist_green)
    Green_sd = sd(select_dist_green)
    Green_se = Green_sd / sqrt(length(select_dist_green)) 
    Green_BC_data = data.frame(BC_value = as.vector(select_dist_green), Group = "Greenhouse")
    
    ## anova test
    all_BC_aov_data = rbind(Field_BC_data, Green_BC_data)
    aov_summary = summary(aov(BC_value ~ Group, data = all_BC_aov_data))
    BC_aov_data_all = rbind(BC_aov_data_all, all_BC_aov_data)
    ## 
    dist_data_all = data.frame(Year = i , Site = ii, 
                               Field_mean = Field_mean, Field_sd = Field_sd, Field_se = Field_se,
                               Green_mean = Green_mean, Green_sd = Green_sd, Green_se = Green_se,
                               F_value = aov_summary[[1]]$"F value"[1],
                               p_value = round(aov_summary[[1]]$"Pr(>F)"[1], 3))
    β_Bray_data = rbind(β_Bray_data, dist_data_all)
  }
}

# View(β_Bray_data)
head(BC_aov_data_all)
summary(aov(BC_value ~ Group, data = BC_aov_data_all))

fit <- deming::deming(Field_mean ~ Green_mean, ystd = Field_sd, xstd = Green_sd, data = β_Bray_data)
print(fit)
β_Bray_data$Site <- factor(β_Bray_data$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))

ggplot()+
  geom_abline(intercept = 0, slope = 1, color = "#96393D", linetype = 1, size = 1)+
  geom_abline(intercept = -19.06659, slope = 28.12380, linetype = 2, size = 1)+
  geom_point(β_Bray_data, mapping = aes(x=Green_mean,y=Field_mean, shape = Year, fill = Site), color = "black", size=2.5)+
  geom_errorbar(data = β_Bray_data,mapping = aes(x = Green_mean,ymax = Field_mean+Field_se, ymin=Field_mean-Field_se, color = Site),width=0,size=0.5,alpha = 1)+#
  geom_errorbarh(data = β_Bray_data,mapping = aes(y = Field_mean,xmax=Green_mean+Green_se,xmin=Green_mean-Green_se, color = Site),height=0,size=0.5,alpha = 1)+#
  scale_color_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  scale_fill_manual(values = c("#F6DD61", "#94684E","#BF5B1D", "#3E91B7", "#0E4879","#CFBD9F"))+
  guides(col = guide_legend(ncol = 1))+
  scale_shape_manual(values = c(24,23,25)) +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.64,0.92)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.64,0.92)) + 
  theme_bw() + mytheme + theme(legend.position = c(0.85,0.45)) + 
  labs(x = NULL, 
       y = "", title = c("Between native and alien plants")) -> mian_Fig_4d3; mian_Fig_4d3


# Combined # 11.73 x 4.68
(mian_Fig_4d1|mian_Fig_4d2|mian_Fig_4d3) -> mian_Fig_2D; mian_Fig_2D

# Notice that,
# For more picture details, we have further adjusted it in Adobe illustrator.

