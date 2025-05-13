################################################################################
################################## Fig 1 #######################################
################################################################################
library(openxlsx)
library(phytools)
library(ggtree)
library(dplyr)
library(ggplot2)
library(sf)
library(ggplot2)
library(ggspatial)
library(patchwork)

# Custom style
main_theme <- theme(panel.background = element_blank(),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(size = 0.5, colour = "NA"),
                    axis.line.y = element_line(size = 0.5, colour = "NA"),
                    axis.ticks = element_line(color = "black"),
                    axis.text = element_text(color = "black", size = 11),
                    legend.position = "right",
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    legend.text = element_text(size = 9),
                    legend.title = element_text(size = 10),
                    plot.tag = element_text(size = 14, face = "bold"))

######################## Fig.1a Sampling point drawing #########################
china_map <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") 
class(china_map)

# Information of sampling site 
datasel <- data.frame(Site = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"),
                     Latitude = c(23.1, 25.2, 27.9, 30.5, 34.6, 36.2),
                     Longitude = c(113.2, 110.2, 112.9, 114.3, 113.6, 117.1))

ggplot(china_map)+
  geom_sf(data = china_map,fill = "grey95", size = 1) + 
  xlim(105,122) + ylim(18,42)+ 
  ggnewscale::new_scale_fill() + 
  geom_point(data = datasel, aes(x = Longitude, y = Latitude),
             size = 4, alpha = 1, shape = 21, color = "black", fill = "grey50") + 
  main_theme +
  annotation_scale(location = "br", style = "ticks", line_width = 0.1, pad_y = unit(0.5, "cm"), text_cex = 0.9) + 
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"), style = north_arrow_fancy_orienteering) +
  guides(fill = guide_legend(title = "Site",ncol = 1, nrow = 14, override.aes = list(shape = 22, size = 5))) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill = "#FFFFFF", colour = "black"),
        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.x.top = element_line(size = 0.5, colour = "black"),
        axis.line.y.right = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        legend.position = "right",
        panel.grid.major = element_line(color = "white",size = 0.2))+
  labs(x= "", y= "") +
  ggrepel::geom_text_repel(mapping = aes(x = Longitude, y = Latitude, label = Site), data = datasel, size = 3.5, segment.color = "black", 
                           color = "black", direction = "both", box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25)) -> mian_Fig_1_left; mian_Fig_1_left


# Loading the grouping metadata of soil samples 
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

sample_num_sum <- Field_group %>% dplyr::group_by(Year, Site, Origin) %>%
  dplyr::summarise(sample_num = n())

df_wide <- sample_num_sum %>%
  tidyr::pivot_wider(names_from = Origin, values_from = sample_num)
df_wide$Total <- df_wide$Exotic + df_wide$Native
df_wide$Site <- factor(df_wide$Site, levels = rev(c("Guangzhou","Guilin","Changsha",
                                                    "Wuhan","Zhengzhou","Tai'an")))


df_long <- df_wide %>%
  tidyr::pivot_longer(cols = c(Exotic, Native, Total), 
                      names_to = "Group",       
                      values_to = "sample_num")  

df_long$Year <- factor(df_long$Year, levels = rev(c("2018", "2020", "2021")))
df_long$Group <- factor(df_long$Group, levels = c("Native", "Exotic", "Total"))

ggplot(df_long, aes(Group, Year)) +
  geom_tile(color = "black", fill = NA) +
  geom_text(aes(label = sample_num), size = 2.8, color = 'black') + 
  ggforce::facet_col(Site~., space = "free") +
  labs(x = NULL, y = NULL) + 
  theme(strip.text.x = element_blank(), 
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, angle = 35, vjust = 1, hjust = 1, color = "black"),
        panel.grid = element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_rect(fill = 'transparent')) -> mian_Fig_1_right; mian_Fig_1_right


# Combined
(mian_Fig_1_left|mian_Fig_1_right) + plot_layout(widths = c(0.9,0.18)) -> mian_Fig_1A; mian_Fig_1A


####################### Fig.1b Plant phylogenetic tree #########################
plant_tree <- read.tree("IQ_tree_plant_2025.NEWICK")

# Loading taxonomic information of studied plant species
Species_tax <- unique(Field_group[, c("Species","Origin","Family","Genus")])

tree_df <- tibble::as_tibble(plant_tree) 
tree_df$Species <- tree_df$label
tree_df <- tree_df %>% left_join(Species_tax)

Family_color = c("Acanthaceae"="#332500","Amaranthaceae"="#542D20","Asteraceae" ="#994240","Caryophyllaceae" = "#D75C4D", "Cyperaceae" = "#E68C51",
                 "Euphorbiaceae" ="#F59D52", "Fabaceae" = "#EFBA55","Lamiaceae" ="#FCD170","Malvaceae"="#FEE1B0", "Onagraceae" = "#C5E0D0",
                 "Phytolaccaceae" = "#ABDCE0","Poaceae" ="#7FC5D8","Polygonaceae"="#73BCD5", "Solanaceae" = "#528FAC", "Urticaceae" = "#376694", "Verbenaceae" = "#1F466F", "NA" = "grey40")

tree_df$Origin <- factor(tree_df$Origin, levels = c("Native", "Exotic"))

ggtree(plant_tree,layout = "fan", branch.length = "none", ladderize = F,size = 0.8,             
       open.angle = 180) %<+% tree_df +
  geom_tiplab(aes(label = sub("_", " ", label)), size = 3, offset = 0.3, fontface = "italic") +
  geom_tree(aes(color = Family), size = 0.5) + 
  scale_color_manual(values = Family_color) +
  ggnewscale::new_scale_color() + 
  geom_tippoint(aes(color = Origin, shape = Origin), size = 1.5) +
  scale_color_manual(values = c("#70A7C3","#A67C2A")) +
  scale_shape_manual(values = c(16,15)) +
  theme(legend.position = "none") -> mian_Fig_1B

# The warning message is caused by the missing annotation information of the outer group. 
# This will not cause any errors.

# Notice that,
# For more picture details, we have further adjusted it in Adobe illustrator.


